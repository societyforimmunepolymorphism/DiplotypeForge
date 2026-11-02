"""Milestone 3A acceptance tests: frequency-preserving homozygosity transform.

Covers brief §10 items 1-10 and 14 (items 11-13 are Milestone 3B). Includes
analytically tractable examples where the feasible min/max are known exactly.
"""

from __future__ import annotations

import copy
from collections import Counter

import numpy as np
import pytest

from diplotypeforge.glstring import parse_diplotype
from diplotypeforge.model.core import Allele, Diplotype, Haplotype, SimulatedIndividual
from diplotypeforge.transform import (
    TransformationTarget,
    analyze_feasibility,
    count_homozygous,
    transform_homozygosity,
)
from diplotypeforge.transform.reassignment import ReassignmentStrategy


def hap(*tokens: str) -> Haplotype:
    return Haplotype(tuple(Allele(t.split("*")[0], t) for t in tokens))


def ind(i: int, a: str, b: str, o1: str = "P1", o2: str = "P2") -> SimulatedIndividual:
    return SimulatedIndividual(f"s{i}", Diplotype(hap(a), hap(b)), o1, o2, 0)


def hap_counts(pop) -> dict:
    c: Counter = Counter()
    for x in pop:
        c[x.diplotype.haplotype_1.key()] += 1
        c[x.diplotype.haplotype_2.key()] += 1
    return dict(c)


def allele_counts(pop) -> dict:
    c: Counter = Counter()
    for x in pop:
        for h in (x.diplotype.haplotype_1, x.diplotype.haplotype_2):
            for al in h.alleles:
                c[al.name] += 1
    return dict(c)


def origin_label_counts(pop) -> dict:
    c: Counter = Counter()
    for x in pop:
        c[x.origin_1] += 1
        c[x.origin_2] += 1
    return dict(c)


# --- analytically tractable feasibility ------------------------------------


def test_analytic_feasibility_classic_increase_example():
    # (A,B),(A,B): counts A=2,B=2, N=2 -> H_min=0, H_max=2. Only two haplotypes
    # so H=1 is NOT reachable (parity, step 2): the (A,B),(A,B) -> (A,A),(B,B)
    # move flips BOTH individuals at once.
    pop = [ind(0, "HLA-A*01", "HLA-A*02"), ind(1, "HLA-A*01", "HLA-A*02")]
    feas = analyze_feasibility(pop)
    assert feas.h_min == 0
    assert feas.h_max == 2
    assert feas.step == 2
    assert feas.is_reachable(0) and feas.is_reachable(2)
    assert not feas.is_reachable(1)


def test_analytic_feasibility_dominant_forces_min():
    # 3 x (A,A), 1 x (B,C): A=6, B=1, C=1; N=4; c_max=6>4 -> H_min=6-4=2.
    pop = [ind(i, "HLA-A*01", "HLA-A*01") for i in range(3)] + [
        ind(3, "HLA-A*02", "HLA-A*03")
    ]
    feas = analyze_feasibility(pop)
    assert feas.max_haplotype_count == 6
    assert feas.h_min == 2  # max(0, 6-4)


def test_analytic_feasibility_two_haplotype_parity():
    # Only two haplotypes -> step 2 (parity): only every other H reachable.
    pop = [ind(i, "HLA-A*01", "HLA-A*02") for i in range(5)] + [
        ind(5, "HLA-A*01", "HLA-A*01")
    ]
    # counts: A = 2*5/... recompute: 5 hets contribute A,B each; 1 hom AA -> A=7,B=5
    feas = analyze_feasibility(pop)
    assert feas.distinct_haplotypes == 2
    assert feas.step == 2


# --- acceptance 1-3, 8, 9 (invariants) --------------------------------------


@pytest.fixture
def random_pop():
    rng = np.random.default_rng(42)
    haps = ["HLA-A*01", "HLA-A*02", "HLA-A*03", "HLA-A*04"]
    return [
        ind(i, str(rng.choice(haps)), str(rng.choice(haps)), o1="P1", o2="P2")
        for i in range(500)
    ]


def test_acceptance1_haplotype_counts_identical(random_pop):
    before = hap_counts(random_pop)
    after, _ = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.6), rng=np.random.default_rng(1)
    )
    assert hap_counts(after) == before


def test_acceptance2_allele_counts_identical(random_pop):
    before = allele_counts(random_pop)
    after, _ = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.2), rng=np.random.default_rng(1)
    )
    assert allele_counts(after) == before


def test_acceptance3_individual_count_unchanged(random_pop):
    after, _ = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.5), rng=np.random.default_rng(1)
    )
    assert len(after) == len(random_pop)


def test_acceptance4_increase_reaches_feasible_target(random_pop):
    feas = analyze_feasibility(random_pop)
    initial = count_homozygous(random_pop)
    target = (initial + feas.h_max) // 2  # a feasible count above initial
    after, result = transform_homozygosity(
        random_pop, TransformationTarget(fraction=target / len(random_pop)),
        rng=np.random.default_rng(2),
    )
    assert result.achieved_count == result.clamped_target_count
    assert count_homozygous(after) > initial


def test_acceptance5_decrease_reaches_feasible_target(random_pop):
    feas = analyze_feasibility(random_pop)
    initial = count_homozygous(random_pop)
    target = (initial + feas.h_min) // 2
    after, result = transform_homozygosity(
        random_pop, TransformationTarget(fraction=target / len(random_pop)),
        rng=np.random.default_rng(2),
    )
    assert result.achieved_count == result.clamped_target_count
    assert count_homozygous(after) < initial


def test_acceptance6_infeasible_target_reported_not_silent(random_pop):
    # Request fraction 1.0 (all homozygous) which is infeasible here.
    after, result = transform_homozygosity(
        random_pop, TransformationTarget(fraction=1.0), rng=np.random.default_rng(3)
    )
    feas = analyze_feasibility(random_pop)
    assert result.rounded_target_count == len(random_pop)        # requested
    assert result.clamped_target_count == feas.h_max             # feasible max
    assert result.achieved_count == feas.h_max
    assert result.reached is False                               # not silently "reached"
    assert result.reason == "clamped_to_feasible_max"


def test_acceptance7_terminates_at_boundaries(random_pop):
    # Both extreme targets must terminate and land on the feasible boundary.
    feas = analyze_feasibility(random_pop)
    lo, rl = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.0), rng=np.random.default_rng(4)
    )
    hi, rh = transform_homozygosity(
        random_pop, TransformationTarget(fraction=1.0), rng=np.random.default_rng(4)
    )
    assert count_homozygous(lo) == feas.h_min
    assert count_homozygous(hi) == feas.h_max


def test_acceptance8_input_unchanged(random_pop):
    snapshot = copy.deepcopy(random_pop)
    before_keys = [i.diplotype.ordered_key() for i in random_pop]
    _ = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.9), rng=np.random.default_rng(5)
    )
    after_keys = [i.diplotype.ordered_key() for i in random_pop]
    assert before_keys == after_keys
    assert [i.sample_id for i in random_pop] == [i.sample_id for i in snapshot]


def test_acceptance9_provenance_follows_haplotypes():
    # Build a population where origin is tied to haplotype so we can verify the
    # origin travels with the copy. P1 always carries A; P2 always carries B.
    pop = [ind(i, "HLA-A*01", "HLA-A*02", o1="P1", o2="P2") for i in range(10)]
    before_origin = origin_label_counts(pop)  # P1:10, P2:10
    after, _ = transform_homozygosity(
        pop, TransformationTarget(fraction=1.0), rng=np.random.default_rng(6)
    )
    # Origin-label totals are conserved (each copy keeps its origin).
    assert origin_label_counts(after) == before_origin
    # Every A-copy still carries origin P1; every B-copy still carries P2.
    for x in after:
        for h, o in ((x.diplotype.haplotype_1, x.origin_1), (x.diplotype.haplotype_2, x.origin_2)):
            if h.key() == ("HLA-A*01",):
                assert o == "P1"
            elif h.key() == ("HLA-A*02",):
                assert o == "P2"


def test_acceptance10_reproducible_same_seed(random_pop):
    a, ra = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.5), rng=np.random.default_rng(777)
    )
    b, rb = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.5), rng=np.random.default_rng(777)
    )
    assert [i.diplotype.ordered_key() for i in a] == [i.diplotype.ordered_key() for i in b]
    assert ra.achieved_count == rb.achieved_count


def test_acceptance14_output_glstrings_valid(random_pop):
    from diplotypeforge.glstring import diplotype_to_glstring

    after, _ = transform_homozygosity(
        random_pop, TransformationTarget(fraction=0.7), rng=np.random.default_rng(8)
    )
    for x in after:
        gl = diplotype_to_glstring(x.diplotype)
        reparsed = parse_diplotype(gl)
        assert reparsed.ordered_key() == x.diplotype.ordered_key()


# --- strategy robustness ----------------------------------------------------


def test_strategy_hits_every_reachable_target_small():
    # Exhaustive on a small analytically tractable population.
    pop = [ind(i, "HLA-A*01", "HLA-A*02") for i in range(6)] + [
        ind(6 + i, "HLA-A*03", "HLA-A*04") for i in range(4)
    ]
    feas = analyze_feasibility(pop)
    before = hap_counts(pop)
    for h in range(feas.h_min, feas.h_max + 1):
        if not feas.is_reachable(h):
            continue
        out, ops = ReassignmentStrategy().apply(pop, h, np.random.default_rng(h + 1))
        assert count_homozygous(out) == h
        assert hap_counts(out) == before


def test_no_target_both_or_neither_errors():
    with pytest.raises(ValueError):
        TransformationTarget()  # neither
    with pytest.raises(ValueError):
        TransformationTarget(fraction=0.5, delta=0.1)  # both
