"""Frequency-preserving homozygosity transformation (Milestone 3A)."""

from __future__ import annotations

from .diagnostics import TransformDiagnostics, transform_diagnostics
from .feasibility import (
    HomozygosityFeasibility,
    analyze_counts,
    analyze_feasibility,
    count_homozygous,
    haplotype_counts,
)
from .per_locus import (
    PerLocusComparison,
    PerLocusDiagnostics,
    compare_per_locus,
    per_locus_diagnostics,
)
from .reassignment import (
    ReassignmentStrategy,
    TransformationResult,
    TransformationStrategy,
    TransformationTarget,
    resolve_target_count,
    transform_homozygosity,
)

__all__ = [
    "HomozygosityFeasibility",
    "PerLocusComparison",
    "PerLocusDiagnostics",
    "ReassignmentStrategy",
    "TransformDiagnostics",
    "TransformationResult",
    "TransformationStrategy",
    "TransformationTarget",
    "analyze_counts",
    "analyze_feasibility",
    "compare_per_locus",
    "count_homozygous",
    "haplotype_counts",
    "per_locus_diagnostics",
    "resolve_target_count",
    "transform_diagnostics",
    "transform_homozygosity",
]
