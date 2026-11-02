"""Independent verification of the homozygosity-transform feasibility math.

Strategy: for small inventories we *brute-force enumerate every perfect matching*
of the 2N haplotype copies into N unordered pairs and record, for each matching,
how many pairs are homozygous (both copies share a key). The set of homozygous
counts observed across ALL matchings is the ground-truth achievable set. This
enumeration is completely independent of the analytical predicate in
``feasibility.py`` (it never calls ``lo_sum`` or the ceil-bound formula).

We then assert that the analytical module and the constructive reassignment
algorithm agree with this ground truth:

  * ``H_min`` / ``H_max`` equal min / max of the achievable set;
  * ``is_reachable(H)`` matches membership in the achievable set for every H;
  * the ``step`` flag correctly describes the K == 2 parity restriction;
  * the constructive algorithm reaches *exactly* every analytically-feasible H
    and preserves every haplotype count;
  * no feasible target is misclassified infeasible, and vice versa.

Hypothesis drives randomized count vectors; explicit cases cover the edge
scenarios called out in the task (single key, two keys, highly imbalanced,
all-unique, odd/even totals, boundary targets, unreachable intermediates).
"""

from __future__ import annotations

from collections import Counter
from functools import cache

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from diplotypeforge.model.core import Allele, Diplotype, Haplotype, SimulatedIndividual
from diplotypeforge.transform.feasibility import (
    analyze_counts,
    analyze_feasibility,
    count_homozygous,
    haplotype_counts,
)
from diplotypeforge.transform.reassignment import (
    ReassignmentStrategy,
    TransformationTarget,
    transform_homozygosity,
)

# --------------------------------------------------------------------------- #
# Independent brute-force ground truth
# --------------------------------------------------------------------------- #


@cache
def _achievable_set(copies: tuple[int, ...]) -> frozenset[int]:
    """Every homozygous-pair count reachable by SOME perfect matching.

    ``copies`` is a tuple of key labels (ints), length 2N. Enumerates all
    perfect matchings via recursion: pair the first remaining copy with each
    other remaining copy, recurse on the rest. Returns the set of homozygous
    counts. Independent of feasibility.py.
    """
    if not copies:
        return frozenset({0})
    first = copies[0]
    rest = copies[1:]
    out: set[int] = set()
    for i in range(len(rest)):
        partner = rest[i]
        remaining = rest[:i] + rest[i + 1 :]
        add = 1 if partner == first else 0
        for h in _achievable_set(remaining):
            out.add(h + add)
    return frozenset(out)


def _copies_from_counts(counts: dict[int, int]) -> tuple[int, ...]:
    out: list[int] = []
    for key, c in counts.items():
        out.extend([key] * c)
    return tuple(out)


def brute_force_achievable(counts: dict[int, int]) -> frozenset[int]:
    return _achievable_set(_copies_from_counts(counts))


# --------------------------------------------------------------------------- #
# Builders: turn an integer count map into SimulatedIndividuals
# --------------------------------------------------------------------------- #


def _hap(key: int) -> Haplotype:
    # single-locus haplotype whose token encodes the key
    return Haplotype(alleles=(Allele(locus="A", name=f"A*{key:02d}:01"),))


def individuals_from_counts(counts: dict[int, int]) -> list[SimulatedIndividual]:
    """Build N individuals whose 2N haplotype copies realize ``counts`` exactly.

    The initial pairing is arbitrary (we just need *some* valid pairing); the
    transform re-pairs from the copy inventory anyway.
    """
    total = sum(counts.values())
    assert total % 2 == 0, "need an even number of copies"
    flat = _copies_from_counts(counts)
    inds: list[SimulatedIndividual] = []
    for i in range(0, len(flat), 2):
        k1, k2 = flat[i], flat[i + 1]
        inds.append(
            SimulatedIndividual(
                sample_id=f"s{i // 2}",
                diplotype=Diplotype(haplotype_1=_hap(k1), haplotype_2=_hap(k2)),
                origin_1="P0",
                origin_2="P0",
            )
        )
    return inds


# --------------------------------------------------------------------------- #
# Core equivalence: analytic predicate == brute-force ground truth
# --------------------------------------------------------------------------- #


def _check_against_bruteforce(counts: dict[int, int]) -> None:
    n = sum(counts.values()) // 2
    truth = brute_force_achievable(counts)
    feas = analyze_counts(Counter(counts), n)

    assert feas.h_min == min(truth), (counts, feas.h_min, sorted(truth))
    assert feas.h_max == max(truth), (counts, feas.h_max, sorted(truth))

    # reachability matches membership for every H in [0, N]
    for h in range(n + 1):
        assert feas.is_reachable(h) == (h in truth), (counts, h, sorted(truth))

    # step flag describes the achievable stride between h_min and h_max
    reachable_sorted = sorted(h for h in truth if feas.h_min <= h <= feas.h_max)
    if len(reachable_sorted) >= 2:
        diffs = {b - a for a, b in zip(reachable_sorted, reachable_sorted[1:], strict=False)}
        expected_step = 2 if diffs == {2} else 1
        assert feas.step == expected_step, (counts, feas.step, reachable_sorted)


def _check_constructive(counts: dict[int, int]) -> None:
    """Every analytically-feasible H is reached exactly, counts preserved."""
    n = sum(counts.values()) // 2
    inds = individuals_from_counts(counts)
    original_counts = Counter(dict(haplotype_counts(inds)))
    feas = analyze_feasibility(inds)
    strat = ReassignmentStrategy()

    for h in range(feas.h_min, feas.h_max + 1):
        if not feas.is_reachable(h):
            continue
        transformed, ops = strat.apply(inds, h, rng=None)
        # exactly h homozygous
        assert count_homozygous(transformed) == h, (counts, h)
        # every haplotype count preserved exactly
        assert Counter(dict(haplotype_counts(transformed))) == original_counts, (counts, h)
        # N preserved
        assert len(transformed) == n
        # input not mutated
        assert Counter(dict(haplotype_counts(inds))) == original_counts


# --------------------------------------------------------------------------- #
# Explicit edge cases
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "counts",
    [
        {0: 2},  # single key, N=1  -> only homozygous
        {0: 6},  # single key, N=3  -> all homozygous, H in {3}
        {0: 4, 1: 2},  # two keys
        {0: 2, 1: 2},  # two keys, balanced -> parity (H in {0,2})
        {0: 8, 1: 2},  # highly imbalanced
        {0: 10, 1: 1, 2: 1},  # dominant key forces homozygotes
        {0: 1, 1: 1, 2: 1, 3: 1},  # all unique, N=2
        {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1},  # all unique, N=3
        {0: 3, 1: 3},  # two keys odd counts
        {0: 2, 1: 2, 2: 2},  # three keys balanced
        {0: 5, 1: 3, 2: 2},  # mixed
        {0: 4, 1: 4},  # two keys even -> parity stride 2
    ],
)
def test_edge_cases_match_bruteforce(counts):
    _check_against_bruteforce(counts)
    _check_constructive(counts)


def test_two_key_balanced_has_parity_stride():
    # {A:2, B:2}: matchings give H in {0, 2}; H=1 is NOT reachable.
    counts = {0: 2, 1: 2}
    truth = brute_force_achievable(counts)
    assert truth == frozenset({0, 2})
    feas = analyze_counts(Counter(counts), 2)
    assert feas.step == 2
    assert not feas.is_reachable(1)
    assert feas.is_reachable(0) and feas.is_reachable(2)


def test_single_key_forces_full_homozygosity():
    counts = {0: 8}  # N=4, all copies identical
    feas = analyze_counts(Counter(counts), 4)
    assert feas.h_min == feas.h_max == 4
    assert brute_force_achievable(counts) == frozenset({4})


def test_all_unique_spans_zero_only():
    counts = dict.fromkeys(range(6), 1)  # N=3, every copy distinct
    truth = brute_force_achievable(counts)
    assert truth == frozenset({0})  # no homozygote possible
    feas = analyze_counts(Counter(counts), 3)
    assert feas.h_min == 0 and feas.h_max == 0


def test_dominant_key_forces_minimum_homozygosity():
    # 6 individuals, 12 copies, one key has 10 copies.
    counts = {0: 10, 1: 1, 2: 1}  # N=6
    truth = brute_force_achievable(counts)
    feas = analyze_counts(Counter(counts), 6)
    # pigeonhole: 10 copies of key0 in 6 pairs -> at least 4 homozygous of key0
    assert feas.h_min == min(truth)
    assert feas.h_min >= 4


# --------------------------------------------------------------------------- #
# Property-based: random small inventories
# --------------------------------------------------------------------------- #


@st.composite
def small_count_maps(draw) -> dict[int, int]:
    """A count map with an even total of at most 10 copies (N <= 5)."""
    n_keys = draw(st.integers(min_value=1, max_value=5))
    counts = {k: draw(st.integers(min_value=1, max_value=6)) for k in range(n_keys)}
    total = sum(counts.values())
    if total % 2 == 1:  # make the total even by bumping one key
        counts[0] += 1
        total += 1
    # keep enumeration cheap
    if total > 10:
        counts = {0: 2}
    return counts


@given(small_count_maps())
@settings(max_examples=250, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_property_analytic_matches_bruteforce(counts):
    _check_against_bruteforce(counts)


@given(small_count_maps())
@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_property_constructive_reaches_every_feasible(counts):
    _check_constructive(counts)


@given(small_count_maps())
@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_property_no_feasible_target_misclassified(counts):
    """is_reachable is exactly brute-force membership (both directions)."""
    n = sum(counts.values()) // 2
    truth = brute_force_achievable(counts)
    feas = analyze_counts(Counter(counts), n)
    for h in range(n + 1):
        analytic = feas.is_reachable(h)
        ground = h in truth
        assert analytic == ground, (counts, h, analytic, ground, sorted(truth))


# --------------------------------------------------------------------------- #
# Termination: constructive algorithm never loops unboundedly
# --------------------------------------------------------------------------- #


@given(small_count_maps())
@settings(max_examples=100, deadline=2000, suppress_health_check=[HealthCheck.too_slow])
def test_transform_terminates_and_reports(counts):
    inds = individuals_from_counts(counts)
    feas = analyze_feasibility(inds)
    mid = (feas.h_min + feas.h_max) / 2
    n = len(inds)
    transformed, result = transform_homozygosity(
        inds, TransformationTarget(fraction=mid / n if n else 0.0)
    )
    assert feas.h_min <= result.achieved_count <= feas.h_max
    assert feas.is_reachable(result.achieved_count)
    assert Counter(dict(haplotype_counts(transformed))) == Counter(dict(haplotype_counts(inds)))
