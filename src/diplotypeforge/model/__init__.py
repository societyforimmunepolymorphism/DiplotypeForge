"""Typed domain model for DiplotypeForge."""

from __future__ import annotations

from .core import (
    HLA_RESOLUTIONS,
    KIR_RESOLUTIONS,
    Allele,
    Diplotype,
    Haplotype,
    Locus,
    SimulatedIndividual,
)
from .population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)

__all__ = [
    "HLA_RESOLUTIONS",
    "KIR_RESOLUTIONS",
    "Allele",
    "Diplotype",
    "FrequencyEntry",
    "Haplotype",
    "HaplotypeFrequencyDistribution",
    "Locus",
    "SimulatedIndividual",
    "SourcePopulation",
]
