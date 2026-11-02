"""Feasibility analysis for the frequency-preserving homozygosity transform.

Given ``N`` individuals (``2N`` haplotype copies) with per-haplotype counts
``c_h`` (``Σ c_h = 2N``), we rearrange the copies into ``N`` unordered pairs,
preserving every ``c_h``. A target of ``H`` homozygous pairs is **achievable**
iff there exist per-key homozygote counts ``hom_k ∈ [0, ⌊c_k/2⌋]`` with
``Σ hom_k = H`` such that the leftover copies ``leftover_k = c_k − 2·hom_k`` can
be paired heterozygously, i.e. ``max_k leftover_k ≤ N − H`` (the standard
"pair distinct items" condition: the largest leftover must not exceed half the
leftover total ``2(N−H)``).

That het-pairability condition is equivalent to a per-key **lower bound** on the
homozygotes supplied by each key:

    hom_k ≥ lo_k(H) = max(0, ⌈(c_k − (N − H)) / 2⌉)

So ``H`` is achievable **iff**

    Σ_k lo_k(H)  ≤  H  ≤  Σ_k ⌊c_k/2⌋   ( = upper bound )

``H_min`` and ``H_max`` are the smallest/largest achievable ``H`` in ``[0, N]``.
Because ``lo_k`` is nondecreasing in ``H`` and the upper bound is constant, the
achievable set is a single contiguous **stride** of ``H`` values — contiguous
(step 1) when ≥ 3 haplotypes contribute, but possibly step 2 (a parity
constraint) when only two distinct haplotypes carry all copies.

This module is **independent of any transformation algorithm**, and its
achievability predicate is exactly the one the reassignment strategy satisfies.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..model.core import SimulatedIndividual

HaplotypeKey = tuple[str, ...]


def _lo_sum(counts: Sequence[int], n: int, h: int) -> int:
    rem = n - h  # number of heterozygous pairs if H homozygotes are formed
    return sum(max(0, math.ceil((c - rem) / 2)) for c in counts)


def _achievable(counts: Sequence[int], n: int, h: int, upper: int) -> bool:
    if h < 0 or h > n or h > upper:
        return False
    return _lo_sum(counts, n, h) <= h


@dataclass(frozen=True, slots=True)
class HomozygosityFeasibility:
    n_individuals: int
    total_copies: int
    distinct_haplotypes: int
    max_haplotype_count: int
    h_min: int
    h_max: int
    step: int  # 1 if consecutive integers achievable; 2 for the K=2 parity case
    _counts: tuple[int, ...] = ()

    @property
    def min_fraction(self) -> float:
        return self.h_min / self.n_individuals if self.n_individuals else 0.0

    @property
    def max_fraction(self) -> float:
        return self.h_max / self.n_individuals if self.n_individuals else 0.0

    def is_reachable(self, h: int) -> bool:
        """Whether an exact homozygous count ``h`` is achievable."""
        if not self.h_min <= h <= self.h_max:
            return False
        upper = sum(c // 2 for c in self._counts)
        return _achievable(self._counts, self.n_individuals, h, upper)

    def nearest_reachable(self, h: int) -> int:
        """The reachable count closest to ``h`` (ties -> lower)."""
        if self.is_reachable(h):
            return h
        best = self.h_min
        best_d = abs(h - best)
        cur = self.h_min
        while cur <= self.h_max:
            if self.is_reachable(cur):
                d = abs(h - cur)
                if d < best_d:
                    best, best_d = cur, d
            cur += 1
        return best


def haplotype_counts(individuals: Sequence[SimulatedIndividual]) -> Counter[HaplotypeKey]:
    """Count every haplotype copy across both chromosomes of all individuals."""
    counts: Counter[HaplotypeKey] = Counter()
    for ind in individuals:
        counts[ind.diplotype.haplotype_1.key()] += 1
        counts[ind.diplotype.haplotype_2.key()] += 1
    return counts


def analyze_counts(count_map: Mapping[HaplotypeKey, int], n: int) -> HomozygosityFeasibility:
    """Feasibility from a haplotype-count map and the individual count ``N``."""
    counts = tuple(count_map.values())
    total = sum(counts)
    upper = sum(c // 2 for c in counts)
    c_max = max(counts) if counts else 0

    # H_min: smallest achievable H. lo_sum(H) <= H. lo_sum is nonincreasing-ish;
    # scan up from the classic lower bound max(0, c_max - N).
    classic_min = max(0, c_max - n)
    h_min = None
    for h in range(classic_min, n + 1):
        if _achievable(counts, n, h, upper):
            h_min = h
            break
    if h_min is None:
        h_min = 0  # degenerate (n == 0)

    # H_max: largest achievable H, scanning down from the upper bound.
    h_max = h_min
    for h in range(upper, h_min - 1, -1):
        if _achievable(counts, n, h, upper):
            h_max = h
            break

    # Determine the step by checking whether h_min + 1 is reachable.
    step = 1
    if h_max > h_min and not _achievable(counts, n, h_min + 1, upper):
        step = 2

    return HomozygosityFeasibility(
        n_individuals=n,
        total_copies=total,
        distinct_haplotypes=len(counts),
        max_haplotype_count=c_max,
        h_min=h_min,
        h_max=h_max,
        step=step,
        _counts=counts,
    )


def analyze_feasibility(individuals: Sequence[SimulatedIndividual]) -> HomozygosityFeasibility:
    """Compute the exact feasible homozygosity range from the haplotype counts."""
    n = len(individuals)
    counts = haplotype_counts(individuals)
    return analyze_counts(counts, n)


def count_homozygous(individuals: Sequence[SimulatedIndividual]) -> int:
    return sum(1 for ind in individuals if ind.diplotype.is_homozygous)
