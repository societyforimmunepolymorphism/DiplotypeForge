"""Milestone 2 tests: generalized K×K origin matrix, modes, expectations.

Deterministic tests use analytically tractable distributions (small, exact).
Stochastic acceptance tests use tolerances scaled by sample size.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from diplotypeforge.expectation import (
    expected_diplotype_distribution,
    expected_diplotype_distribution_pooled_hwe,
    marginal_ancestry_proportions,
    marginal_haplotype_frequencies,
    ordered_population_pair_proportions,
    unordered_population_pair_proportions,
)
from diplotypeforge.experiments import run_experiment
from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from diplotypeforge.sampling.origin import OriginModel, OriginModelError
from diplotypeforge.sampling.sampler import Sampler, make_rng

LOCUS = Locus("HLA-A", "g")


def _pop(pop_id, freqs):
    entries = tuple(
        FrequencyEntry(Haplotype((Allele("HLA-A", h),)), f) for h, f in freqs.items()
    )
    return SourcePopulation(pop_id, HaplotypeFrequencyDistribution(loci=(LOCUS,), entries=entries))


# --- OriginModel validation (dims / nonneg / normalization) -----------------


def test_matrix_must_be_square():
    with pytest.raises(OriginModelError, match="square"):
        OriginModel([[0.5, 0.5]])


def test_matrix_entries_nonnegative():
    with pytest.raises(OriginModelError, match="nonnegative"):
        OriginModel([[1.2, -0.2], [0.0, 0.0]])


def test_matrix_must_sum_to_one():
    with pytest.raises(OriginModelError, match="sum to 1"):
        OriginModel([[0.5, 0.0], [0.0, 0.4]])


def test_named_constructors_build_valid_matrices():
    assert OriginModel.single_population().k == 1
    w = OriginModel.wahlund([1, 1, 2])
    np.testing.assert_allclose(np.diag(w.matrix), [0.25, 0.25, 0.5])
    adm = OriginModel.admixture(0.8)
    np.testing.assert_allclose(adm.matrix, [[0.1, 0.4], [0.4, 0.1]])


def test_symmetric_option_symmetrizes():
    m = OriginModel.from_matrix([[0.4, 0.4], [0.0, 0.2]], symmetric=True)
    np.testing.assert_allclose(m.matrix, m.matrix.T)


# --- Theoretical marginal ancestry & pair proportions (deterministic) -------


def test_marginal_ancestry_identical_for_three_experiment_matrices():
    mats = {
        "random": np.array([[0.25, 0.25], [0.25, 0.25]]),
        "wahlund": np.array([[0.5, 0.0], [0.0, 0.5]]),
        "admixture": np.array([[0.1, 0.4], [0.4, 0.1]]),
    }
    ancestries = [marginal_ancestry_proportions(m) for m in mats.values()]
    for a in ancestries:
        np.testing.assert_allclose(a, [0.5, 0.5])


def test_asymmetric_matrix_marginal_ancestry():
    # Asymmetric: row/col sums differ; a_i = 1/2 (row_i + col_i).
    m = np.array([[0.1, 0.5], [0.2, 0.2]])  # sums to 1.0
    a = marginal_ancestry_proportions(m)
    # a_0 = 1/2((0.1+0.5)+(0.1+0.2)) = 1/2(0.6+0.3)=0.45; a_1 = 1/2(0.4+0.7)=0.55
    np.testing.assert_allclose(a, [0.45, 0.55])
    assert abs(a.sum() - 1.0) < 1e-12


def test_ordered_and_unordered_pair_proportions():
    m = np.array([[0.1, 0.4], [0.3, 0.2]])
    np.testing.assert_allclose(ordered_population_pair_proportions(m), m)
    unordered = unordered_population_pair_proportions(m)
    assert abs(unordered[(0, 0)] - 0.1) < 1e-12
    assert abs(unordered[(1, 1)] - 0.2) < 1e-12
    assert abs(unordered[(0, 1)] - 0.7) < 1e-12  # 0.4 + 0.3


# --- Expected diplotype distribution (deterministic, analytically tractable) -


def test_single_population_expectation_is_hwe():
    pop = _pop("P", {"HLA-A*01:01": 0.6, "HLA-A*02:01": 0.4})
    exp = expected_diplotype_distribution([pop], np.array([[1.0]]))
    # ordered P(h,h)=f^2
    k1 = ("HLA-A*01:01",)
    k2 = ("HLA-A*02:01",)
    assert abs(exp.ordered[(k1, k1)] - 0.36) < 1e-12
    assert abs(exp.ordered[(k2, k2)] - 0.16) < 1e-12
    # unordered heterozygote = 2 f1 f2
    assert abs(exp.unordered[(k1, k2)] - 2 * 0.6 * 0.4) < 1e-12
    # distribution sums to 1 (ordered)
    assert abs(sum(exp.ordered.values()) - 1.0) < 1e-12
    assert abs(sum(exp.unordered.values()) - 1.0) < 1e-12


def test_wahlund_full_stratification_expectation():
    # Disjoint pops, pure stratification: every individual homozygous-by-source.
    p1 = _pop("P1", {"HLA-A*01:01": 1.0})
    p2 = _pop("P2", {"HLA-A*02:01": 1.0})
    m = np.array([[0.5, 0.0], [0.0, 0.5]])
    exp = expected_diplotype_distribution([p1, p2], m)
    k1, k2 = ("HLA-A*01:01",), ("HLA-A*02:01",)
    # Only two homozygous diplotypes, each 0.5; no heterozygotes.
    assert abs(exp.ordered[(k1, k1)] - 0.5) < 1e-12
    assert abs(exp.ordered[(k2, k2)] - 0.5) < 1e-12
    assert (k1, k2) not in exp.unordered
    assert abs(exp.expected_homozygosity() - 1.0) < 1e-12


def test_m_based_differs_from_pooled_hwe_under_structure():
    p1 = _pop("P1", {"HLA-A*01:01": 1.0})
    p2 = _pop("P2", {"HLA-A*02:01": 1.0})
    m = np.array([[0.5, 0.0], [0.0, 0.5]])
    exp_m = expected_diplotype_distribution([p1, p2], m)
    exp_pool = expected_diplotype_distribution_pooled_hwe([p1, p2], m)
    assert abs(exp_m.expected_homozygosity() - 1.0) < 1e-12
    assert abs(exp_pool.expected_homozygosity() - 0.5) < 1e-12  # f_bar=0.5/0.5


def test_marginal_haplotype_frequencies_match_pooled_mix():
    p1 = _pop("P1", {"HLA-A*01:01": 0.8, "HLA-A*02:01": 0.2})
    p2 = _pop("P2", {"HLA-A*01:01": 0.2, "HLA-A*02:01": 0.8})
    m = np.array([[0.25, 0.25], [0.25, 0.25]])  # a = [0.5, 0.5]
    fbar = marginal_haplotype_frequencies([p1, p2], m)
    # 0.5*0.8 + 0.5*0.2 = 0.5 each
    assert abs(fbar[("HLA-A*01:01",)] - 0.5) < 1e-12
    assert abs(fbar[("HLA-A*02:01",)] - 0.5) < 1e-12


# --- K >= 3 matrix ----------------------------------------------------------


def test_k3_matrix_expectation_sums_to_one_and_samples():
    p1 = _pop("P1", {"HLA-A*01:01": 1.0})
    p2 = _pop("P2", {"HLA-A*02:01": 1.0})
    p3 = _pop("P3", {"HLA-A*03:01": 1.0})
    m = np.array(
        [
            [0.2, 0.1, 0.05],
            [0.1, 0.2, 0.05],
            [0.05, 0.05, 0.2],
        ]
    )
    m = m / m.sum()
    exp = expected_diplotype_distribution([p1, p2, p3], m)
    assert abs(sum(exp.ordered.values()) - 1.0) < 1e-12
    assert abs(sum(exp.unordered.values()) - 1.0) < 1e-12

    sampler = Sampler([p1, p2, p3], OriginModel(m))
    rng, _ = make_rng(123)
    inds = list(sampler.sample(2000, rng, id_prefix="K3", replicate=0))
    # all three origins appear
    origins = {o for ind in inds for o in (ind.origin_1, ind.origin_2)}
    assert origins == {"P1", "P2", "P3"}


# --- provenance preserved per haplotype -------------------------------------


def test_per_haplotype_provenance_matches_matrix():
    p1 = _pop("P1", {"HLA-A*01:01": 1.0})
    p2 = _pop("P2", {"HLA-A*02:01": 1.0})
    # Asymmetric M so origin_1/origin_2 distributions differ.
    m = np.array([[0.1, 0.6], [0.1, 0.2]])
    sampler = Sampler([p1, p2], OriginModel(m))
    rng, _ = make_rng(7)
    n = 40000
    inds = list(sampler.sample(n, rng, id_prefix="prov", replicate=0))
    # P(origin_1 = P1) = row0 sum = 0.7 ; P(origin_2 = P1) = col0 sum = 0.2
    o1_p1 = sum(1 for i in inds if i.origin_1 == "P1") / n
    o2_p1 = sum(1 for i in inds if i.origin_2 == "P1") / n
    se = math.sqrt(0.7 * 0.3 / n)
    assert abs(o1_p1 - 0.7) < 5 * se
    assert abs(o2_p1 - 0.2) < 5 * math.sqrt(0.2 * 0.8 / n)


# --- stochastic acceptance: observed matches M-based expectation ------------


@pytest.mark.parametrize("name", ["random", "wahlund", "admixture"])
def test_observed_matches_m_based_expectation_within_tolerance(name):
    n = 60000
    result = run_experiment(n_individuals=n, seed=20120815)
    mr = result.models[name]
    # TV distance between observed and M-based expectation should be small; the
    # bound is set from sample size (multinomial fluctuation ~ 1/sqrt(n)).
    tol = 8.0 / math.sqrt(n)
    assert mr.tv_distance_observed_vs_m < tol, (name, mr.tv_distance_observed_vs_m, tol)


@pytest.mark.parametrize(
    "name,expected_hom",
    [
        ("wahlund", None),   # highest homozygosity
        ("random", None),
        ("admixture", None),
    ],
)
def test_observed_homozygosity_matches_expected(name, expected_hom):
    n = 60000
    result = run_experiment(n_individuals=n, seed=99)
    mr = result.models[name]
    exp_hom = mr.expected.expected_homozygosity()
    se = math.sqrt(max(exp_hom * (1 - exp_hom), 1e-9) / n)
    assert abs(mr.observed_homozygosity - exp_hom) < 6 * se


# --- the flagship experiment invariants (deterministic) ---------------------


def test_three_matrix_experiment_same_marginal_different_diplotypes():
    result = run_experiment(n_individuals=0)  # theory only, fully deterministic
    # All three matrices have identical marginal ancestry (0.5/0.5)...
    for mr in result.models.values():
        assert abs(mr.marginal_ancestry["POP1"] - 0.5) < 1e-12
        assert abs(mr.marginal_ancestry["POP2"] - 0.5) < 1e-12
    # ...hence identical expected marginal haplotype frequencies...
    assert result.marginal_haplotype_frequencies_identical()
    # ...but DIFFERENT expected diplotype distributions.
    assert result.diplotype_distributions_differ()


def test_three_matrix_wahlund_has_highest_homozygosity():
    result = run_experiment(n_individuals=0)
    hom = {name: mr.expected.expected_homozygosity() for name, mr in result.models.items()}
    # Wahlund keeps both haplotypes within one source -> most homozygosity.
    # Admixture here has large off-diagonal mass (0.8 cross-source pairing) and
    # the two populations have largely disjoint haplotypes, so cross pairs are
    # heterozygous -> lowest homozygosity. Verified ordering:
    #   wahlund > random > admixture
    assert hom["wahlund"] > hom["random"] > hom["admixture"]


def test_random_mixing_equals_pooled_hwe():
    # When M is the outer product of the marginal ancestry (random union of
    # gametes across the mixture), the M-based expectation EQUALS pooled-HWE.
    p1 = _pop("P1", {"HLA-A*01:01": 0.8, "HLA-A*02:01": 0.2})
    p2 = _pop("P2", {"HLA-A*01:01": 0.2, "HLA-A*03:01": 0.8})
    m = np.array([[0.25, 0.25], [0.25, 0.25]])  # = a a^T with a=[0.5,0.5]
    exp_m = expected_diplotype_distribution([p1, p2], m)
    exp_pool = expected_diplotype_distribution_pooled_hwe([p1, p2], m)
    keys = set(exp_m.ordered) | set(exp_pool.ordered)
    for k in keys:
        assert abs(exp_m.ordered.get(k, 0.0) - exp_pool.ordered.get(k, 0.0)) < 1e-12
