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
    "ReassignmentStrategy",
    "TransformDiagnostics",
    "TransformationResult",
    "TransformationStrategy",
    "TransformationTarget",
    "analyze_counts",
    "analyze_feasibility",
    "count_homozygous",
    "haplotype_counts",
    "resolve_target_count",
    "transform_diagnostics",
    "transform_homozygosity",
]
