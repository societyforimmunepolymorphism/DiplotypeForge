"""Expected diplotype distribution: M-based vs pooled-HWE (reviewer item 5)."""

from __future__ import annotations

import numpy as np

from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from diplotypeforge.stats import (
    expected_homozygosity_pooled,
    expected_homozygosity_under_m,
)

LOCI = (Locus("A", "g"),)


def _pop(pop_id, freqs):
    entries = tuple(
        FrequencyEntry(Haplotype((Allele("A", h),)), f) for h, f in freqs.items()
    )
    return SourcePopulation(pop_id, HaplotypeFrequencyDistribution(loci=LOCI, entries=entries))


def test_single_population_m_based_equals_pooled_hwe():
    pop = _pop("P", {"A*01:01": 0.6, "A*02:01": 0.4})
    m = np.array([[1.0]])
    marginal = (m.sum(axis=1) + m.sum(axis=0)) / 2.0
    hom_m = expected_homozygosity_under_m([pop], m)
    hom_pool = expected_homozygosity_pooled([pop], marginal)
    assert abs(hom_m - hom_pool) < 1e-12
    assert abs(hom_m - (0.6**2 + 0.4**2)) < 1e-12


def test_wahlund_m_based_exceeds_pooled_hwe_under_structure():
    # Two disjoint populations, pure stratification (diagonal M): Wahlund effect
    # makes observed/M-based homozygosity exceed the pooled-HWE expectation.
    p1 = _pop("P1", {"A*01:01": 1.0})
    p2 = _pop("P2", {"A*02:01": 1.0})
    pops = [p1, p2]
    m = np.array([[0.5, 0.0], [0.0, 0.5]])  # Wahlund: both haplotypes same source
    marginal = (m.sum(axis=1) + m.sum(axis=0)) / 2.0  # [0.5, 0.5]

    hom_m = expected_homozygosity_under_m(pops, m)
    hom_pool = expected_homozygosity_pooled(pops, marginal)

    # Under full stratification every individual is homozygous at the haplotype
    # level: M-based homozygosity = 1.0.
    assert abs(hom_m - 1.0) < 1e-12
    # Pooled-HWE: f_bar = 0.5/0.5, homozygosity = 0.5^2 + 0.5^2 = 0.5.
    assert abs(hom_pool - 0.5) < 1e-12
    # The two expectations are distinct (that gap is the structure signal).
    assert hom_m > hom_pool
