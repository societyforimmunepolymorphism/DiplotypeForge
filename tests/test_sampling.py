"""Sampling acceptance tests: 1 (freq recovery), 2 (HWE), 7 (valid pairs), 8 (seed)."""

from __future__ import annotations

import numpy as np

from diplotypeforge.sampling.origin import OriginModel
from diplotypeforge.sampling.sampler import Sampler, make_rng
from diplotypeforge.stats import compute_summary
from tests.conftest import TOY_FREQS, build_toy_population


def _sample(n, seed, pop=None):
    pop = pop or build_toy_population()
    sampler = Sampler([pop], OriginModel.single_population())
    rng, _ = make_rng(seed)
    return list(sampler.sample(n, rng, id_prefix="DF-test", replicate=0))


def test_acceptance1_recovers_input_haplotype_frequencies():
    n = 200_000
    individuals = _sample(n, seed=42)
    stats = compute_summary(individuals, [build_toy_population()])
    emp = stats.empirical_haplotype_frequencies
    for tokens, f_in in TOY_FREQS.items():
        observed = emp.get(tokens, 0.0)
        # 2n haplotype draws; binomial SE ~ sqrt(p(1-p)/2n). Allow ~5 SE.
        se = (f_in * (1 - f_in) / (2 * n)) ** 0.5
        assert abs(observed - f_in) < 5 * se, (tokens, observed, f_in)


def test_acceptance2_diplotype_frequencies_match_hwe():
    n = 200_000
    pop = build_toy_population()
    individuals = _sample(n, seed=7, pop=pop)
    stats = compute_summary(individuals, [pop])
    # expected homozygosity = sum f_i^2; for K=1 the M-based and pooled-HWE
    # expectations coincide.
    exp_hom = sum(f * f for f in TOY_FREQS.values())
    assert abs(stats.observed_homozygosity - exp_hom) < 0.01
    assert abs(stats.expected_homozygosity_m_based - exp_hom) < 1e-9
    assert abs(stats.expected_homozygosity_pooled_hwe - exp_hom) < 1e-9
    # a specific heterozygous unordered genotype: 2 f_i f_j
    keys = list(TOY_FREQS)
    fi, fj = TOY_FREQS[keys[0]], TOY_FREQS[keys[1]]
    exp_het_pair = 2 * fi * fj
    uk = tuple(sorted((keys[0], keys[1])))
    gl = "+".join("~".join(h) for h in uk)
    observed = stats.empirical_diplotype_frequencies.get(gl, 0.0)
    assert abs(observed - exp_het_pair) < 0.01, (observed, exp_het_pair)


def test_acceptance7_all_diplotypes_are_valid_complete_pairs():
    individuals = _sample(500, seed=1)
    for ind in individuals:
        d = ind.diplotype
        assert d.haplotype_1.n_loci == 2
        assert d.haplotype_2.n_loci == 2
        # each allele token belongs to the right locus
        for hap in (d.haplotype_1, d.haplotype_2):
            assert [a.locus for a in hap.alleles] == ["A", "B"]


def test_acceptance8_fixed_seed_is_reproducible():
    a = _sample(1000, seed=12345)
    b = _sample(1000, seed=12345)
    keys_a = [(i.sample_id, i.diplotype.ordered_key()) for i in a]
    keys_b = [(i.sample_id, i.diplotype.ordered_key()) for i in b]
    assert keys_a == keys_b


def test_different_seeds_differ():
    a = _sample(1000, seed=1)
    b = _sample(1000, seed=2)
    keys_a = [i.diplotype.ordered_key() for i in a]
    keys_b = [i.diplotype.ordered_key() for i in b]
    assert keys_a != keys_b


def test_order_independence_of_input_entries():
    """Shuffling input entry order must not change sampled output (sorted internally)."""
    import random

    pop = build_toy_population()
    entries = list(pop.distribution.entries)
    random.Random(0).shuffle(entries)
    from dataclasses import replace

    shuffled_dist = replace(pop.distribution, entries=tuple(entries))
    from diplotypeforge.model.population import SourcePopulation

    shuffled_pop = SourcePopulation("TOY", shuffled_dist)

    a = _sample(2000, seed=99, pop=pop)
    b = _sample(2000, seed=99, pop=shuffled_pop)
    assert [i.diplotype.ordered_key() for i in a] == [i.diplotype.ordered_key() for i in b]


def test_origin_model_single_population_marginals():
    m = OriginModel.single_population()
    assert m.k == 1
    np.testing.assert_allclose(m.marginal_origin_proportions, [1.0])
