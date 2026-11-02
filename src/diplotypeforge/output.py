"""Output writers: TSV / gzip-TSV and the machine-readable manifest.

The diplotype table and manifest constitute the output contract shared with
TissueTypeForge. See docs/DESIGN.md §6.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path
from typing import TextIO

from .glstring import diplotype_to_glstring, haplotype_to_glstring
from .model.core import SimulatedIndividual
from .stats import SummaryStatistics


def _open_tsv(path: Path, compress: bool) -> tuple[TextIO, Path]:
    if compress:
        path = path.with_suffix(path.suffix + ".gz")
        return gzip.open(path, "wt", newline="", encoding="utf-8"), path
    return path.open("w", newline="", encoding="utf-8"), path


def write_diplotypes(
    individuals: Sequence[SimulatedIndividual], path: Path, compress: bool = False
) -> Path:
    """Write individual-level table: sample_id, phased GL String, origins."""
    fh, out_path = _open_tsv(path, compress)
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["sample_id", "replicate", "diplotype_glstring", "origin_1", "origin_2"]
        )
        for ind in individuals:
            writer.writerow(
                [
                    ind.sample_id,
                    ind.replicate,
                    diplotype_to_glstring(ind.diplotype),
                    ind.origin_1,
                    ind.origin_2,
                ]
            )
    finally:
        fh.close()
    return out_path


def write_haplotype_frequencies(
    stats: SummaryStatistics, path: Path, compress: bool = False
) -> Path:
    fh, out_path = _open_tsv(path, compress)
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["haplotype", "observed_frequency"])
        for key, freq in sorted(stats.empirical_haplotype_frequencies.items()):
            writer.writerow([haplotype_key_to_glstring(key), f"{freq:.12g}"])
    finally:
        fh.close()
    return out_path


def write_diplotype_frequencies(
    stats: SummaryStatistics, path: Path, compress: bool = False
) -> Path:
    fh, out_path = _open_tsv(path, compress)
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["diplotype", "observed_frequency"])
        for key, freq in sorted(stats.empirical_diplotype_frequencies.items()):
            writer.writerow([key, f"{freq:.12g}"])
    finally:
        fh.close()
    return out_path


def haplotype_key_to_glstring(key: tuple[str, ...]) -> str:
    return "~".join(key)


def write_summary(stats: SummaryStatistics, path: Path) -> Path:
    payload = {
        "level": stats.level,
        "n_individuals": stats.n_individuals,
        "n_distinct_haplotypes_observed": stats.n_distinct_haplotypes_observed,
        "observed_homozygosity": stats.observed_homozygosity,
        "observed_heterozygosity": stats.observed_heterozygosity,
        "expected_under_M": {
            "null_model": "M-based diplotype distribution P(h,k)=sum_ij M[i,j] f_i(h) f_j(k)",
            "homozygosity": stats.expected_homozygosity_m_based,
            "heterozygosity": stats.expected_heterozygosity_m_based,
            "tv_distance_observed_vs_expected": stats.tv_distance_vs_m_based,
        },
        "expected_pooled_hwe": {
            "null_model": "pooled Hardy-Weinberg from marginal frequency f_bar(h)",
            "homozygosity": stats.expected_homozygosity_pooled_hwe,
            "heterozygosity": stats.expected_heterozygosity_pooled_hwe,
            "tv_distance_observed_vs_expected": stats.tv_distance_vs_pooled_hwe,
        },
        "origin_proportions": {
            "observed": stats.origin_proportions,
            "theoretical": stats.theoretical_origin_proportions,
        },
        "origin_pair_proportions": {
            "observed": stats.origin_pair_proportions,
            "theoretical": stats.theoretical_origin_pair_proportions,
        },
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def write_manifest(manifest: dict, path: Path) -> Path:
    path.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    return path


# re-export for convenience
__all__ = [
    "asdict",
    "diplotype_to_glstring",
    "haplotype_to_glstring",
    "write_diplotype_frequencies",
    "write_diplotypes",
    "write_haplotype_frequencies",
    "write_manifest",
    "write_summary",
]
