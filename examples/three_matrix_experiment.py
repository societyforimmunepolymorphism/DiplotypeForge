#!/usr/bin/env python3
"""Reproducible three-matrix population-structure experiment.

Demonstrates that three origin matrices with *identical* marginal ancestry
proportions produce the **same** expected marginal haplotype frequencies but
**different** expected diplotype distributions.

Usage:
    uv run python examples/three_matrix_experiment.py [--n N] [--seed S] [--outdir DIR]

With --n > 0 the models are also simulated and observed frequencies compared to
the theoretical expectation. A PNG visualization is written if matplotlib is
installed (`uv pip install -e ".[viz]"`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from diplotypeforge.experiments import comparison_table, run_experiment


def _write_visualization(result, out_png: Path) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    names = list(result.models)
    # Union of unordered diplotype GL keys across models, stable order.
    diplos: list[str] = []
    seen: set[str] = set()
    for name in names:
        for ukey, _ in sorted(result.models[name].expected.unordered.items()):
            gl = "+".join("~".join(h) for h in ukey)
            if gl not in seen:
                seen.add(gl)
                diplos.append(gl)
    diplos.sort()

    import numpy as np

    x = np.arange(len(diplos))
    width = 0.8 / len(names)

    fig, (ax_hap, ax_dip) = plt.subplots(1, 2, figsize=(16, 6))

    # Left: marginal haplotype frequencies (identical across models).
    hap_keys = ["~".join(k) for k in result.haplotype_order()]
    xh = np.arange(len(hap_keys))
    for i, name in enumerate(names):
        vals = [result.models[name].marginal_haplotype_frequencies.get(hk, 0.0) for hk in hap_keys]
        ax_hap.bar(xh + i * width, vals, width, label=name)
    ax_hap.set_title("Expected marginal haplotype freqs f̄(h)\n(IDENTICAL across models)")
    ax_hap.set_xticks(xh + width * (len(names) - 1) / 2)
    ax_hap.set_xticklabels(hap_keys, rotation=45, ha="right")
    ax_hap.set_ylabel("frequency")
    ax_hap.legend()

    # Right: expected unordered diplotype distribution (DIFFERS across models).
    for i, name in enumerate(names):
        glmap = {
            "+".join("~".join(h) for h in k): v
            for k, v in result.models[name].expected.unordered.items()
        }
        vals = [glmap.get(gl, 0.0) for gl in diplos]
        ax_dip.bar(x + i * width, vals, width, label=name)
    ax_dip.set_title("Expected diplotype distribution under M\n(DIFFERS across models)")
    ax_dip.set_xticks(x + width * (len(names) - 1) / 2)
    ax_dip.set_xticklabels(diplos, rotation=90, fontsize=6)
    ax_dip.set_ylabel("frequency")
    ax_dip.legend()

    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n", type=int, default=50000,
        help="individuals to simulate per model (0=theory only)",
    )
    parser.add_argument("--seed", type=int, default=20120815)
    parser.add_argument("--outdir", default="examples/output_experiment")
    args = parser.parse_args()

    result = run_experiment(n_individuals=args.n, seed=args.seed)

    table = comparison_table(result)
    print(table)

    # Invariant checks (also covered by tests).
    same_marginal = result.marginal_haplotype_frequencies_identical()
    different_diplo = result.diplotype_distributions_differ()
    print()
    print(f"Marginal haplotype frequencies identical across models? {same_marginal}")
    print(f"Expected diplotype distributions differ across models?   {different_diplo}")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "comparison_table.txt").write_text(table + "\n", encoding="utf-8")

    payload = {
        "seed": args.seed,
        "n_individuals": args.n,
        "marginal_haplotype_frequencies_identical": same_marginal,
        "diplotype_distributions_differ": different_diplo,
        "models": {
            name: {
                "matrix": mr.matrix,
                "marginal_ancestry": mr.marginal_ancestry,
                "marginal_haplotype_frequencies": mr.marginal_haplotype_frequencies,
                "expected_homozygosity_m_based": mr.expected.expected_homozygosity(),
                "observed_homozygosity": mr.observed_homozygosity,
                "tv_distance_observed_vs_m": mr.tv_distance_observed_vs_m,
                "expected_unordered_diplotype_frequencies": {
                    "+".join("~".join(h) for h in k): v
                    for k, v in mr.expected.unordered.items()
                },
            }
            for name, mr in result.models.items()
        },
    }
    (out / "experiment.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    png = out / "three_matrix_experiment.png"
    if _write_visualization(result, png):
        print(f"\nVisualization written: {png}")
    else:
        print("\n(matplotlib not installed; skipped visualization. `uv pip install -e \".[viz]\"`)")
    print(f"Table + JSON written under: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
