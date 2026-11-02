"""Phase B: per-locus diagnostics experiment on MULTILOCUS haplotypes.

Demonstrates that a frequency-preserving reassignment changes:
  * complete-haplotype homozygosity,
  * homozygosity at each individual locus (by differing amounts),
  * the multilocus diplotype distribution,
while every haplotype count and every per-locus allele count stay fixed.

Also asserts the documented caveat: the change in multilocus homozygosity is
NOT proportional to the per-locus changes.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from diplotypeforge.sampling.origin import OriginModel
from diplotypeforge.sampling.sampler import Sampler, make_rng
from diplotypeforge.transform import (
    TransformationTarget,
    compare_per_locus,
    per_locus_diagnostics,
    transform_homozygosity,
)

LOCI = (Locus("A", "g"), Locus("B", "g"), Locus("C", "g"))


def _hap(a, b, c):
    return Haplotype((Allele("A", a), Allele("B", b), Allele("C", c)))


def _multilocus_population() -> SourcePopulation:
    # Deliberate cross-locus allele sharing so per-locus homozygosity can differ
    # from complete-haplotype homozygosity and move unevenly across loci.
    entries = (
        FrequencyEntry(_hap("A*01:01", "B*08:01", "C*07:01"), 0.35),
        FrequencyEntry(_hap("A*01:01", "B*44:02", "C*05:01"), 0.25),  # shares A*01:01
        FrequencyEntry(_hap("A*02:01", "B*08:01", "C*07:02"), 0.20),  # shares B*08:01
        FrequencyEntry(_hap("A*03:01", "B*07:02", "C*07:01"), 0.20),  # shares C*07:01
    )
    dist = HaplotypeFrequencyDistribution(loci=LOCI, entries=entries)
    return SourcePopulation("MULTI", dist)


def _sample(n: int, seed: int):
    pop = _multilocus_population()
    sampler = Sampler([pop], OriginModel.single_population())
    rng, _ = make_rng(seed)
    return list(sampler.sample(n, rng, id_prefix="S", replicate=0))


def _allele_counts(individuals):
    """Per-locus allele count multiset across both chromosomes."""
    n_loci = len(individuals[0].diplotype.haplotype_1.alleles)
    counts = [Counter() for _ in range(n_loci)]
    for ind in individuals:
        for li in range(n_loci):
            counts[li][ind.diplotype.haplotype_1.alleles[li].name] += 1
            counts[li][ind.diplotype.haplotype_2.alleles[li].name] += 1
    return counts


def _diplotype_distribution(individuals):
    c = Counter()
    for ind in individuals:
        c[ind.diplotype.unordered_key()] += 1
    return c


def test_per_locus_experiment_increase_homozygosity():
    before = _sample(5000, seed=20120815)
    after, result = transform_homozygosity(
        before, TransformationTarget(fraction=0.40), rng=np.random.default_rng(1)
    )

    cmp = compare_per_locus(before, after)

    # 1. haplotype counts and per-locus allele counts are preserved EXACTLY.
    assert cmp.counts_preserved
    assert cmp.allele_counts_preserved
    assert _allele_counts(before) == _allele_counts(after)

    # 2. complete-haplotype homozygosity actually changed toward the target.
    assert cmp.before.complete_homozygosity != cmp.after.complete_homozygosity
    assert abs(cmp.after.complete_homozygosity - result.achieved_homozygosity) < 1e-9

    # 3. per-locus homozygosity changed, and NOT uniformly across loci.
    deltas = cmp.per_locus_delta()
    assert any(abs(d) > 1e-6 for d in deltas.values())
    distinct = {round(d, 6) for d in deltas.values()}
    assert len(distinct) > 1, deltas  # loci move by different amounts

    # 4. per-locus homozygosity is always >= complete-haplotype homozygosity.
    for locus, h in cmp.after.per_locus_homozygosity.items():
        assert h >= cmp.after.complete_homozygosity - 1e-9, (locus, h)

    # 5. the change in multilocus homozygosity is NOT proportional at each locus.
    complete_delta = cmp.complete_delta()
    assert abs(complete_delta) > 1e-6
    ratios = [d / complete_delta for d in deltas.values()]
    assert not all(abs(r - ratios[0]) < 1e-6 for r in ratios), ratios

    # 6. the multilocus diplotype distribution changed.
    assert _diplotype_distribution(before) != _diplotype_distribution(after)


def test_per_locus_homozygosity_bounds():
    inds = _sample(2000, seed=7)
    diag = per_locus_diagnostics(inds)
    # each per-locus homozygosity in [complete, 1]
    for h in diag.per_locus_homozygosity.values():
        assert diag.complete_homozygosity - 1e-9 <= h <= 1.0
    assert set(diag.per_locus_homozygosity) == {"A", "B", "C"}


def test_decrease_homozygosity_also_preserves_counts():
    before = _sample(3000, seed=99)
    diag_before = per_locus_diagnostics(before)
    after, _ = transform_homozygosity(
        before, TransformationTarget(fraction=0.0), rng=np.random.default_rng(3)
    )
    cmp = compare_per_locus(before, after)
    assert cmp.counts_preserved and cmp.allele_counts_preserved
    # complete homozygosity should drop to (or toward) its feasible minimum
    assert cmp.after.complete_homozygosity <= diag_before.complete_homozygosity
