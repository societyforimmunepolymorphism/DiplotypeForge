"""Core typed domain objects for DiplotypeForge.

None of these types *is* a GL String. GL String is a serialization produced on
demand by :mod:`diplotypeforge.glstring`. Internally we keep structured,
phase-preserving, provenance-aware objects.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# HLA resolution values permitted by the authoritative HFX v0.1.1 schema.
HLA_RESOLUTIONS = frozenset(
    {
        "G",
        "P",
        "g",
        "allele-family",
        "protein-sequence",
        "coding-sequence",
        "genomic-sequence",
        "Serology",
    }
)
KIR_RESOLUTIONS = frozenset({"protein-sequence", "coding-sequence", "genomic-sequence"})


@dataclass(frozen=True, slots=True)
class Locus:
    """A named locus together with its declared typing resolution.

    ``name`` is the locus label as it appears in the haplotype token, e.g.
    ``"A"`` or ``"DRB1"`` (optionally prefixed with the genetic system, e.g.
    ``"HLA-A"``). ``resolution`` is an HFX resolution string.
    """

    name: str
    resolution: str

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Locus.name must be non-empty")


@dataclass(frozen=True, slots=True)
class Allele:
    """A single allele call at a single locus.

    ``locus`` is the locus name. ``name`` is the full allele token *including*
    the locus prefix as it appears in a GL String, e.g. ``"A*01:01g"``. Keeping
    the full token verbatim guarantees lossless round-tripping; structured
    accessors derive fields on demand.
    """

    locus: str
    name: str

    def __post_init__(self) -> None:
        if not self.locus:
            raise ValueError("Allele.locus must be non-empty")
        if not self.name:
            raise ValueError("Allele.name must be non-empty")

    @property
    def is_absent(self) -> bool:
        """Whether this allele represents gene absence (future CNV support)."""
        return self.name == ""


@dataclass(frozen=True, slots=True)
class Haplotype:
    """An ordered tuple of alleles across the declared loci of a distribution.

    The i-th allele corresponds to the i-th locus in the owning distribution's
    ``loci`` tuple. Ordering is therefore significant and is preserved.
    """

    alleles: tuple[Allele, ...]

    def __post_init__(self) -> None:
        if not self.alleles:
            raise ValueError("Haplotype must contain at least one allele")

    @property
    def n_loci(self) -> int:
        return len(self.alleles)

    def key(self) -> tuple[str, ...]:
        """Order-sensitive, hashable identity key (the allele tokens)."""
        return tuple(a.name for a in self.alleles)


@dataclass(frozen=True, slots=True)
class Diplotype:
    """Two haplotypes whose cross-locus correspondence (phase) is known.

    Stored as an *ordered* pair. Downstream comparisons may treat diplotypes as
    unordered genotypes via :meth:`unordered_key`, but the ordered pair is always
    retained.
    """

    haplotype_1: Haplotype
    haplotype_2: Haplotype

    def __post_init__(self) -> None:
        if self.haplotype_1.n_loci != self.haplotype_2.n_loci:
            raise ValueError("Both haplotypes of a diplotype must span the same loci")

    @property
    def is_homozygous(self) -> bool:
        """True when both haplotypes are identical across all loci."""
        return self.haplotype_1.key() == self.haplotype_2.key()

    def ordered_key(self) -> tuple[tuple[str, ...], tuple[str, ...]]:
        return (self.haplotype_1.key(), self.haplotype_2.key())

    def unordered_key(self) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """Phase-preserving but order-independent genotype identity."""
        k1, k2 = self.haplotype_1.key(), self.haplotype_2.key()
        return (k1, k2) if k1 <= k2 else (k2, k1)


@dataclass(frozen=True, slots=True)
class SimulatedIndividual:
    """One sampled individual with per-haplotype population provenance."""

    sample_id: str
    diplotype: Diplotype
    origin_1: str = field(default="")
    origin_2: str = field(default="")
    replicate: int = 0
