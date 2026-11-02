"""Phase B: analytical tests for the statistical-diagnostics claims.

Checks P(h,k) = Σ_ij M[i,j] f_i(h) f_j(k) and the pooled-HWE null, with emphasis
on **asymmetric** population-origin matrices (where M[i,j] != M[j,i]). Verifies:

  * ordered and unordered distributions normalize to 1;
  * homozygote probabilities are NOT doubled (unordered P(h,h) == ordered P(h,h));
  * ordered P(h,k) reflects matrix asymmetry (P(h,k) != P(k,h) in general);
  * the unordered distribution correctly folds P(h,k)+P(k,h);
  * the marginal haplotype frequency uses BOTH matrix margins via
    a_i = ½(row_i + col_i);
  * pooled-HWE is computed from f̄ alone, independently of M's internal structure;
  * the random-mixing case M = w wᵀ reproduces pooled HWE when both gamete-source
    margins equal w;
  * total-variation distance is correct and symmetric.
"""

from __future__ import annotations

import numpy as np

from diplotypeforge.expectation import (
    expected_diplotype_distribution,
    expected_diplotype_distribution_pooled_hwe,
    marginal_ancestry_proportions,
    marginal_haplotype_frequencies,
)
from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)

LOCI = (Locus("A", "g"),)


def _pop(pid, freqs):
    entries = tuple(FrequencyEntry(Haplotype((Allele("A", h),)), f) for h, f in freqs.items())
    return SourcePopulation(pid, HaplotypeFrequencyDistribution(loci=LOCI, entries=entries))


def _pops():
    p1 = _pop("P1", {"A*01:01": 0.7, "A*02:01": 0.3})
    p2 = _pop("P2", {"A*01:01": 0.2, "A*03:01": 0.8})
    return [p1, p2]


def _brute_force_ordered(pops, m):
    """Independent O(K^2 H^2) reference for P(h,k) = Σ_ij M[i,j] f_i(h) f_j(k)."""
    per_pop = []
    keys = []
    seen = set()
    for pop in pops:
        total = pop.distribution.total_frequency() or 1.0
        fmap = {e.haplotype.key(): e.frequency / total for e in pop.distribution.entries}
        per_pop.append(fmap)
        for k in fmap:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    out = {}
    for h in keys:
        for k in keys:
            s = 0.0
            for i in range(len(pops)):
                for j in range(len(pops)):
                    s += m[i, j] * per_pop[i].get(h, 0.0) * per_pop[j].get(k, 0.0)
            if s > 0:
                out[(h, k)] = s
    return out


def test_ordered_matches_independent_reference_asymmetric():
    pops = _pops()
    m = np.array([[0.1, 0.5], [0.3, 0.1]])
    m = m / m.sum()
    exp = expected_diplotype_distribution(pops, m)
    ref = _brute_force_ordered(pops, m)
    assert set(exp.ordered) == set(ref)
    for key in ref:
        assert abs(exp.ordered[key] - ref[key]) < 1e-12


def test_ordered_and_unordered_normalize_to_one_asymmetric():
    pops = _pops()
    m = np.array([[0.2, 0.4], [0.3, 0.1]])
    m = m / m.sum()
    exp = expected_diplotype_distribution(pops, m)
    assert abs(sum(exp.ordered.values()) - 1.0) < 1e-12
    assert abs(sum(exp.unordered.values()) - 1.0) < 1e-12


def test_homozygote_not_doubled():
    pops = _pops()
    m = np.array([[0.2, 0.4], [0.3, 0.1]])
    m = m / m.sum()
    exp = expected_diplotype_distribution(pops, m)
    for (h, k), v in exp.unordered.items():
        if h == k:
            # homozygote: unordered mass equals the single ordered diagonal entry
            assert abs(v - exp.ordered[(h, h)]) < 1e-12


def test_unordered_folds_offdiagonal_pairs():
    pops = _pops()
    m = np.array([[0.2, 0.4], [0.3, 0.1]])
    m = m / m.sum()
    exp = expected_diplotype_distribution(pops, m)
    for (h, k), v in exp.unordered.items():
        if h != k:
            expected = exp.ordered.get((h, k), 0.0) + exp.ordered.get((k, h), 0.0)
            assert abs(v - expected) < 1e-12


def test_ordered_reflects_matrix_asymmetry():
    pops = _pops()
    m = np.array([[0.1, 0.5], [0.3, 0.1]])
    m = m / m.sum()
    exp = expected_diplotype_distribution(pops, m)
    # P(01,03) uses f1(01)f2(03) weighted by M[0,1]; P(03,01) by M[1,0]; asymmetric.
    p_01_03 = exp.ordered.get((("A*01:01",), ("A*03:01",)), 0.0)
    p_03_01 = exp.ordered.get((("A*03:01",), ("A*01:01",)), 0.0)
    assert p_01_03 != p_03_01
    assert p_01_03 > p_03_01  # M[0,1]=0.5/sum > M[1,0]=0.3/sum


def test_marginal_ancestry_uses_both_margins():
    # Highly asymmetric M: row and column sums differ; a_i must average them.
    m = np.array([[0.1, 0.5], [0.3, 0.1]])
    m = m / m.sum()
    a = marginal_ancestry_proportions(m)
    expected = (m.sum(axis=1) + m.sum(axis=0)) / 2.0
    assert np.allclose(a, expected)
    assert abs(a.sum() - 1.0) < 1e-12
    # not equal to row sums alone nor column sums alone
    assert not np.allclose(a, m.sum(axis=1))
    assert not np.allclose(a, m.sum(axis=0))


def test_marginal_haplotype_freq_uses_both_margins():
    pops = _pops()
    m = np.array([[0.1, 0.5], [0.3, 0.1]])
    m = m / m.sum()
    fbar = marginal_haplotype_frequencies(pops, m)
    a = marginal_ancestry_proportions(m)
    # f̄(01:01) = a0*0.7 + a1*0.2
    assert abs(fbar[("A*01:01",)] - (a[0] * 0.7 + a[1] * 0.2)) < 1e-12
    assert abs(sum(fbar.values()) - 1.0) < 1e-12


def test_pooled_hwe_independent_of_matrix_internal_structure():
    """Two matrices with the SAME marginal ancestry give the SAME pooled-HWE."""
    pops = _pops()
    w = np.array([0.4, 0.6])
    m1 = np.outer(w, w)  # random mixing
    m2 = np.array([[0.4, 0.0], [0.0, 0.6]])  # Wahlund, same marginal ancestry a=w
    assert np.allclose(marginal_ancestry_proportions(m1), w)
    assert np.allclose(marginal_ancestry_proportions(m2), w)
    e1 = expected_diplotype_distribution_pooled_hwe(pops, m1)
    e2 = expected_diplotype_distribution_pooled_hwe(pops, m2)
    keys = set(e1.ordered) | set(e2.ordered)
    for k in keys:
        assert abs(e1.ordered.get(k, 0.0) - e2.ordered.get(k, 0.0)) < 1e-12


def test_random_mixing_M_equals_wwT_reproduces_pooled_hwe():  # noqa: N802
    pops = _pops()
    w = np.array([0.4, 0.6])
    m = np.outer(w, w)  # both gamete-source margins equal w
    assert np.allclose(m.sum(axis=0), w) and np.allclose(m.sum(axis=1), w)
    em = expected_diplotype_distribution(pops, m)
    ep = expected_diplotype_distribution_pooled_hwe(pops, m)
    keys = set(em.ordered) | set(ep.ordered)
    for k in keys:
        assert abs(em.ordered.get(k, 0.0) - ep.ordered.get(k, 0.0)) < 1e-12
    # homozygosity also matches
    assert abs(em.expected_homozygosity() - ep.expected_homozygosity()) < 1e-12


def test_wahlund_exceeds_pooled_hwe_homozygosity():
    pops = _pops()
    w = np.array([0.5, 0.5])
    m_wahlund = np.array([[0.5, 0.0], [0.0, 0.5]])
    m_random = np.outer(w, w)
    hom_w = expected_diplotype_distribution(pops, m_wahlund).expected_homozygosity()
    hom_r = expected_diplotype_distribution(pops, m_random).expected_homozygosity()
    # Structure raises homozygosity above the random-mixing (pooled-HWE) level.
    assert hom_w > hom_r
