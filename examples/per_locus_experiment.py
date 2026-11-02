"""Per-locus homozygosity experiment (multilocus validation).

Shows how a frequency-preserving homozygosity reassignment changes:
  * complete-haplotype homozygosity,
  * homozygosity at each individual locus,
  * the multilocus diplotype distribution,
while every haplotype count and per-locus allele count stay fixed.

Run:
    uv run python examples/per_locus_experiment.py --n 20000 --target 0.4

The printed table makes the central caveat visible: the change in multilocus
homozygosity is NOT proportional to the change at any single locus.
"""

from __future__ import annotations

import argparse

import numpy as np

from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from diplotypeforge.sampling.origin import OriginModel
from diplotypeforge.sampling.sampler import Sampler, make_rng
from diplotypeforge.transform import TransformationTarget, compare_per_locus, transform_homozygosity

LOCI = (Locus("A", "g"), Locus("B", "g"), Locus("C", "g"))


def _hap(a, b, c):
    return Haplotype((Allele("A", a), Allele("B", b), Allele("C", c)))


def population() -> SourcePopulation:
    entries = (
        FrequencyEntry(_hap("A*01:01", "B*08:01", "C*07:01"), 0.35),
        FrequencyEntry(_hap("A*01:01", "B*44:02", "C*05:01"), 0.25),
        FrequencyEntry(_hap("A*02:01", "B*08:01", "C*07:02"), 0.20),
        FrequencyEntry(_hap("A*03:01", "B*07:02", "C*07:01"), 0.20),
    )
    return SourcePopulation("MULTI", HaplotypeFrequencyDistribution(loci=LOCI, entries=entries))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--target", type=float, default=0.40)
    ap.add_argument("--seed", type=int, default=20120815)
    args = ap.parse_args()

    sampler = Sampler([population()], OriginModel.single_population())
    rng, _ = make_rng(args.seed)
    before = list(sampler.sample(args.n, rng, id_prefix="S", replicate=0))

    after, result = transform_homozygosity(
        before, TransformationTarget(fraction=args.target), rng=np.random.default_rng(args.seed)
    )
    cmp = compare_per_locus(before, after)

    print(f"N = {args.n}   requested complete-haplotype homozygosity = {args.target}")
    print(f"counts preserved: {cmp.counts_preserved}   "
          f"allele counts preserved: {cmp.allele_counts_preserved}")
    print()
    cd = cmp.complete_delta()
    print(f"{'measurement':<28}{'before':>10}{'after':>10}{'delta':>10}{'delta/complete':>16}")
    print("-" * 74)
    print(f"{'complete-haplotype':<28}"
          f"{cmp.before.complete_homozygosity:>10.4f}"
          f"{cmp.after.complete_homozygosity:>10.4f}{cd:>10.4f}{1.0:>16.3f}")
    pld = cmp.per_locus_delta()
    for locus in cmp.before.locus_names:
        b = cmp.before.per_locus_homozygosity[locus]
        a = cmp.after.per_locus_homozygosity[locus]
        d = pld[locus]
        ratio = d / cd if abs(cd) > 1e-12 else float("nan")
        print(f"{'locus ' + locus:<28}{b:>10.4f}{a:>10.4f}{d:>10.4f}{ratio:>16.3f}")
    print()
    print("Note:", cmp.note)


if __name__ == "__main__":
    main()
