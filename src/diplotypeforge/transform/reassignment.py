"""Frequency-preserving homozygosity transformation (reassignment strategy).

The transform rearranges the ``2N`` haplotype copies of a sampled population into
``N`` unordered pairs to reach a requested number of homozygous diplotypes, while
preserving every haplotype count exactly. Each haplotype copy carries its source
population origin, which travels with the copy when it moves between individuals.

Correctness guarantees (see docs/DESIGN.md §4.4 and OPEN_QUESTIONS Q13/Q15):

- Every haplotype count (hence every allele count) and ``N`` are preserved: all
  operations only *move* existing copies between individuals.
- Any feasible target in ``[H_min, H_max]`` is reached with a **bounded** number
  of elementary operations — no unbounded rejection loops. Each operation
  strictly moves the homozygous count toward the target.
- The strategy is **not** a uniform sampler over valid matchings and does not
  claim to be; an optional RNG shuffle only avoids deterministic artifacts.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol

import numpy as np

from ..model.core import Diplotype, Haplotype, SimulatedIndividual
from .feasibility import HomozygosityFeasibility, analyze_feasibility, count_homozygous

HaplotypeKey = tuple[str, ...]

REASON_EXACT = "exact"
REASON_ROUNDED = "rounded_to_1/N_grid"
REASON_CLAMPED_MIN = "clamped_to_feasible_min"
REASON_CLAMPED_MAX = "clamped_to_feasible_max"
REASON_SNAPPED_REACHABLE = "snapped_to_nearest_reachable_count"


@dataclass(frozen=True, slots=True)
class TransformationTarget:
    """A homozygosity target: exactly one of ``fraction`` or ``delta``."""

    fraction: float | None = None
    delta: float | None = None

    def __post_init__(self) -> None:
        if (self.fraction is None) == (self.delta is None):
            raise ValueError("specify exactly one of 'fraction' or 'delta'")
        if self.fraction is not None and not 0.0 <= self.fraction <= 1.0:
            raise ValueError(f"fraction must be in [0, 1]; got {self.fraction}")


@dataclass(slots=True)
class TransformationResult:
    """Full report of a transformation (never silently substitutes a target)."""

    n_individuals: int
    initial_homozygosity: float
    requested_fraction: float | None
    rounded_target_count: int
    clamped_target_count: int
    achieved_homozygosity: float
    achieved_count: int
    feasible_min_count: int
    feasible_max_count: int
    feasible_min_fraction: float
    feasible_max_fraction: float
    n_operations: int
    reached: bool
    reason: str
    strategy: str = "reassignment"


class TransformationStrategy(Protocol):
    """Interface so alternative strategies (greedy, swaps, exact) can be added."""

    name: str

    def apply(
        self,
        individuals: Sequence[SimulatedIndividual],
        target_count: int,
        rng: np.random.Generator | None,
    ) -> tuple[list[SimulatedIndividual], int]: ...


class _Copy:
    """A haplotype copy that carries its object + origin so provenance travels."""

    __slots__ = ("key", "hap", "origin")

    def __init__(self, key: HaplotypeKey, hap: Haplotype, origin: str) -> None:
        self.key = key
        self.hap = hap
        self.origin = origin


class ReassignmentStrategy:
    """One robust strategy: constructive exact pairing with bounded work.

    Rather than nudging the current pairing incrementally (which is fiddly to
    keep exact), we *construct* a valid pairing with exactly ``target_count``
    homozygous diplotypes directly from the copy inventory:

    1. Flatten the 2N copies, grouped by haplotype key.
    2. Form exactly ``H`` homozygous pairs, taken greedily from the **most
       frequent** keys first (each key contributes at most ⌊c/2⌋). Taking from
       the dominant key first guarantees the remainder is heterozygous-pairable
       for every feasible ``H ≥ H_min`` (proof in docs/DESIGN.md §5).
    3. Pair the remaining copies heterozygously by repeatedly joining the two
       most frequent distinct keys (a standard construction that succeeds iff no
       key exceeds half the remaining copies — guaranteed here).

    This is O(2N log K), deterministic, and never loops indefinitely. Every
    haplotype count is preserved because copies are only regrouped. Provenance
    rides on each ``_Copy``.
    """

    name = "reassignment"

    def apply(
        self,
        individuals: Sequence[SimulatedIndividual],
        target_count: int,
        rng: np.random.Generator | None,
    ) -> tuple[list[SimulatedIndividual], int]:
        n = len(individuals)
        # Group copies by key, preserving each copy's object + origin.
        by_key: dict[HaplotypeKey, list[_Copy]] = defaultdict(list)
        for ind in individuals:
            d = ind.diplotype
            by_key[d.haplotype_1.key()].append(
                _Copy(d.haplotype_1.key(), d.haplotype_1, ind.origin_1)
            )
            by_key[d.haplotype_2.key()].append(
                _Copy(d.haplotype_2.key(), d.haplotype_2, ind.origin_2)
            )

        if rng is not None:
            for lst in by_key.values():
                rng.shuffle(lst)

        pairs = _construct_pairs(by_key, target_count, rng)

        # Assemble individuals, reusing original sample_ids/replicate in order.
        out: list[SimulatedIndividual] = []
        ops = 0
        for ind, (c1, c2) in zip(individuals, pairs, strict=True):
            new_ind = SimulatedIndividual(
                sample_id=ind.sample_id,
                diplotype=Diplotype(haplotype_1=c1.hap, haplotype_2=c2.hap),
                origin_1=c1.origin,
                origin_2=c2.origin,
                replicate=ind.replicate,
            )
            if new_ind.diplotype.ordered_key() != ind.diplotype.ordered_key():
                ops += 1
            out.append(new_ind)
        assert len(out) == n
        return out, ops


def _construct_pairs(
    by_key: dict[HaplotypeKey, list[_Copy]],
    target_hom: int,
    rng: np.random.Generator | None,
) -> list[tuple[_Copy, _Copy]]:
    """Build the full list of N pairs with exactly ``target_hom`` homozygous."""
    pairs: list[tuple[_Copy, _Copy]] = []

    # Total copies and N.
    total = sum(len(lst) for lst in by_key.values())
    n = total // 2
    rem_after = n - target_hom  # number of heterozygous pairs we must end with
    # copies remaining after homozygous removal = 2 * rem_after

    # Step 2: choose how many homozygous pairs each key supplies, so that the
    # leftover (non-homozygous) copies are heterozygous-pairable. Het-pairability
    # requires leftover_k <= (remaining copies)/2 = rem_after for every key.
    # leftover_k = c_k - 2*hom_k  =>  hom_k >= ceil((c_k - rem_after)/2).
    counts = {k: len(lst) for k, lst in by_key.items()}
    lo = {k: max(0, -(-(c - rem_after) // 2)) for k, c in counts.items()}  # ceil div
    hi = {k: c // 2 for k, c in counts.items()}

    hom_alloc = dict(lo)
    assigned = sum(hom_alloc.values())
    remaining = target_hom - assigned
    # distribute the remaining homozygous pairs up to each key's cap, most
    # frequent first (keeps the dominant count low / remainder balanced).
    for key in sorted(by_key, key=lambda k: counts[k], reverse=True):
        if remaining <= 0:
            break
        room = hi[key] - hom_alloc[key]
        take = min(room, remaining)
        hom_alloc[key] += take
        remaining -= take

    for key, nhom in hom_alloc.items():
        lst = by_key[key]
        for _ in range(nhom):
            pairs.append((lst.pop(), lst.pop()))

    # Step 3: heterozygous pairs for the remainder, via a max-heap on counts.
    import heapq

    heap = [(-len(lst), key) for key, lst in by_key.items() if lst]
    heapq.heapify(heap)
    while heap:
        neg_a, ka = heapq.heappop(heap)
        if not heap:
            # Only one key left: by construction (lo bounds) this happens only
            # for an infeasible target; pair leftovers same-key as a fallback.
            lst = by_key[ka]
            while len(lst) >= 2:
                pairs.append((lst.pop(), lst.pop()))
            break
        neg_b, kb = heapq.heappop(heap)
        la, lb = by_key[ka], by_key[kb]
        pairs.append((la.pop(), lb.pop()))
        if la:
            heapq.heappush(heap, (-len(la), ka))
        if lb:
            heapq.heappush(heap, (-len(lb), kb))

    if rng is not None:
        rng.shuffle(pairs)
    return pairs


def _round_half_even(x: float) -> int:
    return int(np.rint(x))


def resolve_target_count(
    initial_count: int,
    n: int,
    feasibility: HomozygosityFeasibility,
    target: TransformationTarget,
) -> tuple[int, int, str]:
    """Map a target to (rounded_count, clamped_count, reason)."""
    if target.fraction is not None:
        raw = target.fraction * n
    else:
        # Keep initial_count exact; only delta*n carries float error.
        raw = initial_count + target.delta * n if n else 0.0
    # Snap away floating-point dust so an intended integer target is not nudged
    # across the half-way point by representation error.
    snapped = round(raw)
    if abs(raw - snapped) < 1e-6:
        raw = float(snapped)
    rounded = _round_half_even(raw)
    reason = REASON_EXACT if abs(raw - rounded) < 1e-9 else REASON_ROUNDED

    clamped = rounded
    if clamped < feasibility.h_min:
        clamped = feasibility.h_min
        reason = REASON_CLAMPED_MIN
    elif clamped > feasibility.h_max:
        clamped = feasibility.h_max
        reason = REASON_CLAMPED_MAX

    # Snap to the nearest *reachable* count (captures the K=2 parity constraint).
    if not feasibility.is_reachable(clamped):
        snapped_count = feasibility.nearest_reachable(clamped)
        if snapped_count != clamped:
            clamped = snapped_count
            reason = REASON_SNAPPED_REACHABLE
    return rounded, clamped, reason


def transform_homozygosity(
    individuals: Sequence[SimulatedIndividual],
    target: TransformationTarget,
    rng: np.random.Generator | None = None,
    strategy: TransformationStrategy | None = None,
) -> tuple[list[SimulatedIndividual], TransformationResult]:
    """Apply the frequency-preserving homozygosity transform.

    Returns a NEW transformed list (the input is left unchanged) and a
    :class:`TransformationResult`.
    """
    strat = strategy or ReassignmentStrategy()
    n = len(individuals)
    feas = analyze_feasibility(individuals)
    initial_count = count_homozygous(individuals)

    rounded, clamped, reason = resolve_target_count(initial_count, n, feas, target)

    transformed, ops = strat.apply(individuals, clamped, rng)
    achieved = count_homozygous(transformed)
    reached = achieved == rounded

    result = TransformationResult(
        n_individuals=n,
        initial_homozygosity=initial_count / n if n else 0.0,
        requested_fraction=target.fraction,
        rounded_target_count=rounded,
        clamped_target_count=clamped,
        achieved_homozygosity=achieved / n if n else 0.0,
        achieved_count=achieved,
        feasible_min_count=feas.h_min,
        feasible_max_count=feas.h_max,
        feasible_min_fraction=feas.min_fraction,
        feasible_max_fraction=feas.max_fraction,
        n_operations=ops,
        reached=reached,
        reason=reason,
        strategy=strat.name,
    )
    return transformed, result
