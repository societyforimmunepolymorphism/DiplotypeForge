#!/usr/bin/env python3
"""Large-population performance benchmark (Milestone 3C).

Runs a single-population simulation at scale (default 1,000,000 individuals),
optionally with the frequency-preserving homozygosity transform enabled, and
reports runtime, peak memory, throughput, and output file sizes.

Usage:
    uv run python examples/benchmark_large.py --n 1000000 [--transform 0.3]
                                              [--compress] [--outdir DIR]

Measurement notes:
- Peak memory is measured with ``tracemalloc`` (Python allocations only; the
  numpy RNG buffers are small relative to the individual objects).
- Runtime uses ``time.perf_counter`` and is split into sampling, (optional)
  transform, and write phases.
- Throughput is individuals per second for the end-to-end run.

This is an honest measurement of the current implementation; it does not
introduce premature optimization or distributed computing.
"""

from __future__ import annotations

import argparse
import json
import time
import tracemalloc
from pathlib import Path

import numpy as np

from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from diplotypeforge.sampling.origin import OriginModel
from diplotypeforge.sampling.sampler import Sampler, make_rng
from diplotypeforge.transform import TransformationTarget, transform_homozygosity

_LOCI = (
    ("HLA-A", "g"),
    ("HLA-C", "g"),
    ("HLA-B", "g"),
    ("HLA-DRB1", "g"),
    ("HLA-DQB1", "g"),
)


def _synthetic_population(n_haplotypes: int, seed: int) -> SourcePopulation:
    """Build a synthetic 5-locus population with ``n_haplotypes`` haplotypes."""
    rng = np.random.default_rng(seed)
    loci = tuple(Locus(name, res) for name, res in _LOCI)
    weights = rng.random(n_haplotypes)
    weights /= weights.sum()
    entries = []
    for i in range(n_haplotypes):
        alleles = tuple(
            Allele(name, f"{name}*{(i % 90) + 1:02d}:{(i % 50) + 1:02d}")
            for name, _ in _LOCI
        )
        entries.append(FrequencyEntry(Haplotype(alleles), float(weights[i])))
    dist = HaplotypeFrequencyDistribution(loci=loci, entries=tuple(entries))
    return SourcePopulation(population_id="BENCH", distribution=dist)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=1_000_000)
    parser.add_argument("--haplotypes", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20120815)
    parser.add_argument(
        "--transform", type=float, default=None,
        help="if set, target homozygosity fraction for the transform",
    )
    parser.add_argument("--compress", action="store_true")
    parser.add_argument("--outdir", default="examples/output_benchmark")
    args = parser.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    pop = _synthetic_population(args.haplotypes, args.seed)
    origin = OriginModel.single_population()
    sampler = Sampler([pop], origin)

    tracemalloc.start()
    t0 = time.perf_counter()

    rng, _ = make_rng(args.seed)
    individuals = list(sampler.sample(args.n, rng, id_prefix="BENCH", replicate=0))
    t_sampled = time.perf_counter()

    transform_phase = 0.0
    tresult = None
    if args.transform is not None:
        t_rng, _ = make_rng(args.seed + 1)
        ts = time.perf_counter()
        individuals, tresult = transform_homozygosity(
            individuals, TransformationTarget(fraction=args.transform), rng=t_rng
        )
        transform_phase = time.perf_counter() - ts

    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Write via the normal result machinery.
    from diplotypeforge.output import write_diplotypes

    tw = time.perf_counter()
    dip_path = write_diplotypes(individuals, out / "diplotypes.tsv", args.compress)
    write_phase = time.perf_counter() - tw

    total = time.perf_counter() - t0
    file_size = Path(dip_path).stat().st_size

    report = {
        "n_individuals": args.n,
        "n_haplotypes": args.haplotypes,
        "transform_enabled": args.transform is not None,
        "transform_target_fraction": args.transform,
        "transform_reached": (tresult.reached if tresult else None),
        "runtime_seconds": {
            "sampling": round(t_sampled - t0, 3),
            "transform": round(transform_phase, 3),
            "write": round(write_phase, 3),
            "total": round(total, 3),
        },
        "throughput_individuals_per_second": round(args.n / total, 1) if total else None,
        "peak_memory_mb": round(peak / 1e6, 1),
        "output_file": str(dip_path),
        "output_file_size_mb": round(file_size / 1e6, 2),
        "compressed": args.compress,
    }
    (out / "benchmark.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
