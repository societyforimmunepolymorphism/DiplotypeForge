"""Population-genetic diagnostics.

Observed quantities are computed from the sampled individuals. Expected
quantities are delegated to :mod:`diplotypeforge.expectation`, which is
**independent of the sampling implementation**. Each statistic is labeled by its
level of analysis (complete multilocus haplotype) and null model (M-based vs
pooled-HWE). See docs/DESIGN.md §4.3 / §7.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .expectation import (
    DiplotypeExpectation,
    expected_diplotype_distribution,
    expected_diplotype_distribution_pooled_hwe,
    marginal_ancestry_proportions,
    marginal_haplotype_frequencies,
    ordered_population_pair_proportions,
    unordered_population_pair_proportions,
)
from .model.core import SimulatedIndividual
from .model.population import SourcePopulation

HaplotypeKey = tuple[str, ...]


@dataclass(frozen=True, slots=True)
class HaplotypeFreqComparison:
    haplotype: HaplotypeKey
    input_frequency: float  # pooled/marginal input frequency f̄(h)
    observed_frequency: float
    observed_count: int
    expected_count: float
    standard_error: float


@dataclass(slots=True)
class SummaryStatistics:
    level: str = "multilocus-haplotype"
    n_individuals: int = 0
    n_distinct_haplotypes_observed: int = 0

    observed_homozygosity: float = 0.0
    expected_homozygosity_m_based: float = 0.0
    expected_homozygosity_pooled_hwe: float = 0.0
    observed_heterozygosity: float = 0.0
    expected_heterozygosity_m_based: float = 0.0
    expected_heterozygosity_pooled_hwe: float = 0.0

    # Origin provenance — observed vs theoretical (derived from M).
    origin_proportions: dict[str, float] = field(default_factory=dict)
    origin_pair_proportions: dict[str, float] = field(default_factory=dict)
    theoretical_origin_proportions: dict[str, float] = field(default_factory=dict)
    theoretical_origin_pair_proportions: dict[str, float] = field(default_factory=dict)

    empirical_haplotype_frequencies: dict[HaplotypeKey, float] = field(default_factory=dict)
    empirical_diplotype_frequencies: dict[str, float] = field(default_factory=dict)
    haplotype_comparisons: list[HaplotypeFreqComparison] = field(default_factory=list)

    # Total-variation distance between observed unordered diplotype freqs and
    # each expectation (0 = identical; smaller is closer).
    tv_distance_vs_m_based: float = 0.0
    tv_distance_vs_pooled_hwe: float = 0.0


# --- backward-compatible scalar helpers (delegate to expectation module) ----


def expected_homozygosity_under_m(
    populations: Sequence[SourcePopulation], matrix: np.ndarray
) -> float:
    return expected_diplotype_distribution(populations, matrix).expected_homozygosity()


def expected_homozygosity_pooled(
    populations: Sequence[SourcePopulation], marginal: np.ndarray
) -> float:
    # ``marginal`` is accepted for signature stability but recomputed from the
    # matrix inside the expectation module; callers pass marginal ancestry which
    # equals the implied matrix only for diagonal/structured cases. To keep the
    # pooled-HWE meaning exact we reconstruct f̄ directly here.
    keys: list[HaplotypeKey] = []
    seen: set[HaplotypeKey] = set()
    per_pop: list[dict[HaplotypeKey, float]] = []
    for pop in populations:
        total = pop.distribution.total_frequency() or 1.0
        fmap = {e.haplotype.key(): e.frequency / total for e in pop.distribution.entries}
        per_pop.append(fmap)
        for key in fmap:
            if key not in seen:
                seen.add(key)
                keys.append(key)
    hom = 0.0
    for key in keys:
        fbar = sum(marginal[i] * per_pop[i].get(key, 0.0) for i in range(len(per_pop)))
        hom += fbar * fbar
    return hom


def _diplotype_gl_key(unordered_key: tuple[HaplotypeKey, HaplotypeKey]) -> str:
    return "+".join("~".join(h) for h in unordered_key)


def _tv_distance(
    observed: dict[str, float], expected: DiplotypeExpectation
) -> float:
    """Total-variation distance between observed and expected unordered freqs."""
    exp_map = {_diplotype_gl_key(k): v for k, v in expected.unordered.items()}
    all_keys = set(observed) | set(exp_map)
    return 0.5 * sum(abs(observed.get(k, 0.0) - exp_map.get(k, 0.0)) for k in all_keys)


def compute_summary(
    individuals: Sequence[SimulatedIndividual],
    source_populations: Sequence[SourcePopulation] | None = None,
    origin_matrix: np.ndarray | None = None,
) -> SummaryStatistics:
    """Compute empirical frequencies and expected diagnostics.

    ``origin_matrix`` is the K×K matrix M used by the run. When omitted but a
    single population is given, M defaults to ``[[1]]``.
    """
    n = len(individuals)
    stats = SummaryStatistics(n_individuals=n)
    if n == 0:
        return stats

    hap_counter: Counter[HaplotypeKey] = Counter()
    diplo_counter: Counter[str] = Counter()
    origin_counter: Counter[str] = Counter()
    origin_pair_counter: Counter[str] = Counter()
    n_homozygous = 0

    for ind in individuals:
        d = ind.diplotype
        hap_counter[d.haplotype_1.key()] += 1
        hap_counter[d.haplotype_2.key()] += 1
        uk = d.unordered_key()
        diplo_counter[_diplotype_gl_key(uk)] += 1
        if d.is_homozygous:
            n_homozygous += 1
        origin_counter[ind.origin_1] += 1
        origin_counter[ind.origin_2] += 1
        pair = tuple(sorted((ind.origin_1, ind.origin_2)))
        origin_pair_counter[f"{pair[0]}|{pair[1]}"] += 1

    total_haps = 2 * n
    stats.empirical_haplotype_frequencies = {k: v / total_haps for k, v in hap_counter.items()}
    stats.n_distinct_haplotypes_observed = len(hap_counter)
    stats.empirical_diplotype_frequencies = {k: v / n for k, v in diplo_counter.items()}
    stats.observed_homozygosity = n_homozygous / n
    stats.observed_heterozygosity = 1.0 - stats.observed_homozygosity
    stats.origin_proportions = {k: v / total_haps for k, v in origin_counter.items()}
    stats.origin_pair_proportions = {k: v / n for k, v in origin_pair_counter.items()}

    if source_populations:
        pops = list(source_populations)
        k = len(pops)
        if origin_matrix is not None:
            matrix = np.asarray(origin_matrix, dtype=np.float64)
        elif k > 1:
            matrix = np.eye(k) / k
        else:
            matrix = np.array([[1.0]])

        # Expected diplotype distributions (independent of sampling).
        exp_m = expected_diplotype_distribution(pops, matrix)
        exp_pooled = expected_diplotype_distribution_pooled_hwe(pops, matrix)
        stats.expected_homozygosity_m_based = exp_m.expected_homozygosity()
        stats.expected_heterozygosity_m_based = exp_m.expected_heterozygosity()
        stats.expected_homozygosity_pooled_hwe = exp_pooled.expected_homozygosity()
        stats.expected_heterozygosity_pooled_hwe = exp_pooled.expected_heterozygosity()
        stats.tv_distance_vs_m_based = _tv_distance(stats.empirical_diplotype_frequencies, exp_m)
        stats.tv_distance_vs_pooled_hwe = _tv_distance(
            stats.empirical_diplotype_frequencies, exp_pooled
        )

        # Theoretical origin + origin-pair proportions (derived directly from M).
        ids = [p.population_id for p in pops]
        a = marginal_ancestry_proportions(matrix)
        stats.theoretical_origin_proportions = dict(zip(ids, a.tolist(), strict=True))
        ordered = ordered_population_pair_proportions(matrix)
        unordered_pairs = unordered_population_pair_proportions(matrix)
        stats.theoretical_origin_pair_proportions = {
            f"{ids[i]}|{ids[j]}": prop for (i, j), prop in unordered_pairs.items()
        }
        _ = ordered  # ordered proportions available via manifest/origin_model

        # Per-haplotype comparison vs the pooled/marginal input frequency.
        fbar = marginal_haplotype_frequencies(pops, matrix)
        comparisons = []
        for key, f_in in sorted(fbar.items()):
            obs_count = hap_counter.get(key, 0)
            obs_freq = obs_count / total_haps
            se = (f_in * (1.0 - f_in) / total_haps) ** 0.5
            comparisons.append(
                HaplotypeFreqComparison(
                    haplotype=key,
                    input_frequency=f_in,
                    observed_frequency=obs_freq,
                    observed_count=obs_count,
                    expected_count=f_in * total_haps,
                    standard_error=se,
                )
            )
        stats.haplotype_comparisons = comparisons

    return stats
