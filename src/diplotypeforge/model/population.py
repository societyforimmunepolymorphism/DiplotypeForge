"""Frequency-table abstraction and source populations.

The :class:`HaplotypeFrequencyDistribution` is the internal, format-independent
representation of a population's haplotype frequencies. The HFX adapter maps
into it; the sampling engine reads only from it. GL String is never used as the
primary representation here.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from .core import Haplotype, Locus


@dataclass(frozen=True, slots=True)
class FrequencyEntry:
    """A single haplotype with its frequency and optional observed count."""

    haplotype: Haplotype
    frequency: float
    count: float | None = None


@dataclass(frozen=True, slots=True)
class HaplotypeFrequencyDistribution:
    """Ordered loci plus a list of haplotype/frequency entries.

    Locus order is authoritative and preserved from the source (for HFX, from
    ``outputResolution``). Entry order is preserved as given, but sampling sorts
    a derived key array to guarantee order-independence (see the sampler).
    """

    loci: tuple[Locus, ...]
    entries: tuple[FrequencyEntry, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.loci:
            raise ValueError("A distribution must declare at least one locus")

    @property
    def locus_names(self) -> tuple[str, ...]:
        return tuple(locus.name for locus in self.loci)

    @property
    def n_entries(self) -> int:
        return len(self.entries)

    def total_frequency(self) -> float:
        return sum(entry.frequency for entry in self.entries)

    def haplotypes(self) -> Iterable[Haplotype]:
        return (entry.haplotype for entry in self.entries)


@dataclass(frozen=True, slots=True)
class SourcePopulation:
    """A named source population with its frequency distribution and provenance.

    ``population_id`` is the stable internal identifier used in configuration,
    origin metadata, and output. ``label`` is a human-friendly name (often the
    HFX ``cohortDescription.population[].name``). ``provenance`` preserves input
    metadata (producer, nomenclature, cohort, original source path, etc.).
    """

    population_id: str
    distribution: HaplotypeFrequencyDistribution
    label: str = ""
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.population_id:
            raise ValueError("SourcePopulation.population_id must be non-empty")
