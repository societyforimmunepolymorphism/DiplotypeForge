"""Command-line interface: ``diplotypeforge validate|simulate|summarize``."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

from .config import ConfigError, SimulationConfig
from .simulator import Simulator


def _cmd_validate(args: argparse.Namespace) -> int:
    try:
        config = SimulationConfig.from_file(args.config)
        sim = Simulator(config)
        sim.validate()
    except (ConfigError, ValueError) as exc:
        print(f"VALIDATION FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"OK: configuration and {len(config.populations)} population(s) are valid")
    for pop in sim.source_populations:
        print(
            f"  - {pop.population_id}: {pop.distribution.n_entries} haplotypes, "
            f"loci={list(pop.distribution.locus_names)}, "
            f"sum={pop.distribution.total_frequency():.6g}"
        )
    return 0


def _cmd_simulate(args: argparse.Namespace) -> int:
    try:
        config = SimulationConfig.from_file(args.config)
        if args.output_dir:
            config = replace(config, output_dir=args.output_dir)
        sim = Simulator(config)
        result = sim.run()
        files = result.write()
    except (ConfigError, ValueError, NotImplementedError) as exc:
        print(f"SIMULATION FAILED: {exc}", file=sys.stderr)
        return 1
    print(
        f"Simulated {result.statistics.n_individuals} individuals "
        f"(mode={config.mode}, seed={config.seed})"
    )
    for name, path in files.items():
        print(f"  {name}: {path}")
    return 0


def _cmd_transform(args: argparse.Namespace) -> int:
    import json as _json
    import tomllib
    from dataclasses import asdict

    import numpy as np

    from .transform import TransformationTarget, transform_homozygosity
    from .transform.io import load_dataset

    try:
        cfg_data = tomllib.loads(Path(args.config).read_text(encoding="utf-8"))
        tr = cfg_data.get("transformation", {})
        if not tr.get("enabled", True):
            print("transformation.enabled is false in config; nothing to do", file=sys.stderr)
            return 1
        target = TransformationTarget(
            fraction=tr.get("target_fraction"), delta=tr.get("target_delta")
        )
        dataset = load_dataset(args.input)
        before = dataset.individuals
        # transformation RNG seed: from config, else default.
        seed = int(tr.get("seed", cfg_data.get("simulation", {}).get("seed", 0)))
        rng = np.random.default_rng(seed)
        after, result = transform_homozygosity(before, target, rng=rng)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print(f"TRANSFORM FAILED: {exc}", file=sys.stderr)
        return 1

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    from .output import write_diplotypes

    compress = bool(tr.get("compress", False))
    dip_path = write_diplotypes(after, out / "diplotypes_transformed.tsv", compress)

    report = {
        "input": str(args.input),
        "transformation_result": asdict(result),
    }
    # If the sidecar manifest carried the origin matrix, add before/after diagnostics.
    if dataset.manifest and "origin_model" in dataset.manifest:
        report["origin_model"] = dataset.manifest["origin_model"]
        report["note"] = (
            "The transformed sample is NOT distributed according to the original "
            "origin matrix M; the transform is a distinct sampling mechanism."
        )
    (out / "transform_report.json").write_text(_json.dumps(report, indent=2), encoding="utf-8")

    print(
        f"Transformed {result.n_individuals} individuals: "
        f"homozygosity {result.initial_homozygosity:.4f} -> "
        f"{result.achieved_homozygosity:.4f} "
        f"(target_count={result.clamped_target_count}, reached={result.reached}, "
        f"reason={result.reason}, ops={result.n_operations})"
    )
    print(f"  feasible range: [{result.feasible_min_count}, {result.feasible_max_count}] "
          f"= [{result.feasible_min_fraction:.4f}, {result.feasible_max_fraction:.4f}]")
    print(f"  diplotypes: {dip_path}")
    print(f"  report: {out / 'transform_report.json'}")
    return 0


def _cmd_summarize(args: argparse.Namespace) -> int:
    summary_path = Path(args.output_dir) / "summary.json"
    manifest_path = Path(args.output_dir) / "manifest.json"
    if not summary_path.exists():
        print(f"no summary.json found in {args.output_dir}", file=sys.stderr)
        return 1
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    print(json.dumps(summary, indent=2))
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        print("\n--- manifest (reproducibility) ---")
        subset = {
            k: manifest[k]
            for k in ("run_hash", "software_version", "python_version", "numpy_version", "rng")
            if k in manifest
        }
        for key in ("renormalization", "duplicate_aggregation"):
            if key in manifest:
                subset[key] = manifest[key]
        print(json.dumps(subset, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="diplotypeforge",
        description="Simulate phased diploid genotypes from haplotype-frequency distributions.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_val = sub.add_parser("validate", help="validate a simulation config and its HFX inputs")
    p_val.add_argument("config", help="path to simulation.toml")
    p_val.set_defaults(func=_cmd_validate)

    p_sim = sub.add_parser("simulate", help="run a simulation and write outputs")
    p_sim.add_argument("config", help="path to simulation.toml")
    p_sim.add_argument("--output-dir", default=None, help="override output directory")
    p_sim.set_defaults(func=_cmd_simulate)

    p_sum = sub.add_parser("summarize", help="print diagnostics from an output directory")
    p_sum.add_argument("output_dir", help="path to a simulation output directory")
    p_sum.set_defaults(func=_cmd_summarize)

    p_tr = sub.add_parser(
        "transform",
        help="apply a frequency-preserving homozygosity transform to an existing dataset",
    )
    p_tr.add_argument(
        "--input", required=True, help="path to diplotypes.tsv[.gz] (with origin columns)"
    )
    p_tr.add_argument("--config", required=True, help="transformation TOML config")
    p_tr.add_argument("--output-dir", default="output_transform", help="output directory")
    p_tr.set_defaults(func=_cmd_transform)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
