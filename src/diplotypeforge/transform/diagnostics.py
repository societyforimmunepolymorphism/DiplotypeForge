"""Diagnostics for a frequency-preserving transformed sample.

Reports, for a transform, the four diplotype distributions requested in the
brief, plus the total-variation distances between the relevant pairs:

1. Expected under the ORIGINAL origin matrix M (analytical, sampling-independent).
2. Observed BEFORE transformation.
3. Observed AFTER transformation.
4. Expected pooled-HWE.

The transformed sample is **explicitly not** distributed according to the
original M: the transform is a new sampling mechanism with no claimed analytical
diplotype distribution (OPEN_QUESTIONS Q16). We therefore never compare the
after-sample to M as if M generated it — we only report the distances so the
imposed deviation is visible and attributable.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from ..expectation import (
    expected_diplotype_distribution,
    expected_diplotype_distribution_pooled_hwe,
)
from ..model.core import SimulatedIndividual
from ..model.population import SourcePopulation

HaplotypeKey = tuple[str, ...]


def _gl(unordered_key: tuple[HaplotypeKey, HaplotypeKey]) -> str:
    return "+".join("~".join(h) for h in unordered_key)


def _observed_unordered(individuals: Sequence[SimulatedIndividual]) -> dict[str, float]:
    n = len(individuals)
    if n == 0:
        return {}
    counts: dict[str, int] = {}
    for ind in individuals:
        k = _gl(ind.diplotype.unordered_key())
        counts[k] = counts.get(k, 0) + 1
    return {k: v / n for k, v in counts.items()}


def _expected_unordered_gl(expectation) -> dict[str, float]:
    return {_gl(k): v for k, v in expectation.unordered.items()}


def _tv(a: dict[str, float], b: dict[str, float]) -> float:
    keys = set(a) | set(b)
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)


def _origin_pair_proportions(individuals: Sequence[SimulatedIndividual]) -> dict[str, float]:
    n = len(individuals)
    if n == 0:
        return {}
    counts: dict[str, int] = {}
    for ind in individuals:
        pair = tuple(sorted((ind.origin_1, ind.origin_2)))
        key = f"{pair[0]}|{pair[1]}"
        counts[key] = counts.get(key, 0) + 1
    return {k: v / n for k, v in counts.items()}


@dataclass(slots=True)
class TransformDiagnostics:
    n_individuals: int
    observed_homozygosity_before: float
    observed_homozygosity_after: float
    # Diplotype distributions (unordered, keyed by phased-GL string).
    expected_under_m: dict[str, float] = field(default_factory=dict)
    expected_pooled_hwe: dict[str, float] = field(default_factory=dict)
    observed_before: dict[str, float] = field(default_factory=dict)
    observed_after: dict[str, float] = field(default_factory=dict)
    # Total-variation distances.
    tv_before_vs_m: float = 0.0
    tv_after_vs_m: float = 0.0
    tv_before_vs_pooled: float = 0.0
    tv_after_vs_pooled: float = 0.0
    tv_before_vs_after: float = 0.0
    # Provenance.
    origin_pair_before: dict[str, float] = field(default_factory=dict)
    origin_pair_after: dict[str, float] = field(default_factory=dict)
    note: str = (
        "The transformed sample is NOT distributed according to the original "
        "origin matrix M. The transform is a distinct sampling mechanism with no "
        "claimed analytical diplotype distribution; distances to M and pooled-HWE "
        "are reported only to make the imposed deviation visible."
    )


def transform_diagnostics(
    before: Sequence[SimulatedIndividual],
    after: Sequence[SimulatedIndividual],
    source_populations: Sequence[SourcePopulation],
    origin_matrix: np.ndarray,
) -> TransformDiagnostics:
    pops = list(source_populations)
    m = np.asarray(origin_matrix, dtype=np.float64)

    exp_m = _expected_unordered_gl(expected_diplotype_distribution(pops, m))
    exp_pooled = _expected_unordered_gl(expected_diplotype_distribution_pooled_hwe(pops, m))
    obs_before = _observed_unordered(before)
    obs_after = _observed_unordered(after)

    n = len(before)
    hom_before = sum(1 for i in before if i.diplotype.is_homozygous) / n if n else 0.0
    hom_after = sum(1 for i in after if i.diplotype.is_homozygous) / len(after) if after else 0.0

    return TransformDiagnostics(
        n_individuals=n,
        observed_homozygosity_before=hom_before,
        observed_homozygosity_after=hom_after,
        expected_under_m=exp_m,
        expected_pooled_hwe=exp_pooled,
        observed_before=obs_before,
        observed_after=obs_after,
        tv_before_vs_m=_tv(obs_before, exp_m),
        tv_after_vs_m=_tv(obs_after, exp_m),
        tv_before_vs_pooled=_tv(obs_before, exp_pooled),
        tv_after_vs_pooled=_tv(obs_after, exp_pooled),
        tv_before_vs_after=_tv(obs_before, obs_after),
        origin_pair_before=_origin_pair_proportions(before),
        origin_pair_after=_origin_pair_proportions(after),
    )
