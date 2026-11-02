"""Per-locus homozygosity diagnostics for multilocus diplotypes.

A frequency-preserving reassignment changes *complete-haplotype* homozygosity
(both chromosomes carry the identical multilocus haplotype) while preserving
every haplotype count exactly. Because a multilocus haplotype is a tuple of
single-locus alleles, the per-locus allele counts are also preserved, but the
per-locus *homozygosity* (both chromosomes share the same allele AT THAT LOCUS)
can move differently from the complete-haplotype homozygosity.

This module quantifies, for a population of multilocus diplotypes:

  * ``complete_homozygosity`` — fraction with haplotype_1.key() == haplotype_2.key();
  * ``per_locus_homozygosity[locus]`` — fraction whose two chromosomes share the
    same allele at that locus (a marginal, single-locus view);
  * ``haplotype_counts`` and ``allele_counts`` — conserved invariants.

Key scientific caveat (documented, tested): a change in complete-haplotype
homozygosity does **not** imply a proportional change at any individual locus.
Per-locus homozygosity is bounded below by the complete-haplotype homozygosity
and depends on cross-locus allele sharing between distinct haplotypes, which the
reassignment rearranges.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

from ..model.core import SimulatedIndividual

AlleleToken = str
HaplotypeKey = tuple[str, ...]


@dataclass(slots=True)
class PerLocusDiagnostics:
    n_individuals: int
    locus_names: tuple[str, ...]
    complete_homozygosity: float
    per_locus_homozygosity: dict[str, float] = field(default_factory=dict)
    haplotype_counts: dict[HaplotypeKey, int] = field(default_factory=dict)
    per_locus_allele_counts: dict[str, dict[AlleleToken, int]] = field(default_factory=dict)


def _locus_names(individuals: Sequence[SimulatedIndividual]) -> tuple[str, ...]:
    if not individuals:
        return ()
    hap = individuals[0].diplotype.haplotype_1
    return tuple(a.locus for a in hap.alleles)


def per_locus_diagnostics(individuals: Sequence[SimulatedIndividual]) -> PerLocusDiagnostics:
    n = len(individuals)
    loci = _locus_names(individuals)
    if n == 0:
        return PerLocusDiagnostics(0, loci, 0.0)

    n_loci = len(loci)
    complete_hom = 0
    locus_hom = [0] * n_loci
    hap_counts: Counter[HaplotypeKey] = Counter()
    allele_counts: list[Counter[str]] = [Counter() for _ in range(n_loci)]

    for ind in individuals:
        h1 = ind.diplotype.haplotype_1
        h2 = ind.diplotype.haplotype_2
        if h1.key() == h2.key():
            complete_hom += 1
        hap_counts[h1.key()] += 1
        hap_counts[h2.key()] += 1
        for li in range(n_loci):
            a1 = h1.alleles[li].name
            a2 = h2.alleles[li].name
            allele_counts[li][a1] += 1
            allele_counts[li][a2] += 1
            if a1 == a2:
                locus_hom[li] += 1

    return PerLocusDiagnostics(
        n_individuals=n,
        locus_names=loci,
        complete_homozygosity=complete_hom / n,
        per_locus_homozygosity={loci[li]: locus_hom[li] / n for li in range(n_loci)},
        haplotype_counts=dict(hap_counts),
        per_locus_allele_counts={loci[li]: dict(allele_counts[li]) for li in range(n_loci)},
    )


@dataclass(slots=True)
class PerLocusComparison:
    """Before/after per-locus diagnostics for a frequency-preserving transform."""

    before: PerLocusDiagnostics
    after: PerLocusDiagnostics
    counts_preserved: bool
    allele_counts_preserved: bool
    note: str = (
        "A change in complete-haplotype homozygosity does NOT imply a proportional "
        "change at any single locus. Per-locus homozygosity is always >= the "
        "complete-haplotype homozygosity and reflects cross-locus allele sharing "
        "among distinct haplotypes, which the reassignment rearranges while keeping "
        "every haplotype count (and hence every per-locus allele count) fixed."
    )

    def complete_delta(self) -> float:
        return self.after.complete_homozygosity - self.before.complete_homozygosity

    def per_locus_delta(self) -> dict[str, float]:
        return {
            locus: self.after.per_locus_homozygosity[locus]
            - self.before.per_locus_homozygosity[locus]
            for locus in self.before.locus_names
        }


def compare_per_locus(
    before: Sequence[SimulatedIndividual], after: Sequence[SimulatedIndividual]
) -> PerLocusComparison:
    b = per_locus_diagnostics(before)
    a = per_locus_diagnostics(after)
    return PerLocusComparison(
        before=b,
        after=a,
        counts_preserved=(b.haplotype_counts == a.haplotype_counts),
        allele_counts_preserved=(b.per_locus_allele_counts == a.per_locus_allele_counts),
    )
