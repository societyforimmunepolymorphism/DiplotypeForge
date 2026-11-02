"""Reproducible replicate management and cross-replicate aggregate statistics.

The seed hierarchy (defined in :mod:`diplotypeforge.simulator`) is::

    root = SeedSequence(seed)
      children = root.spawn(R)                 # R replicates, index-stable
        child_k.spawn(2) -> [sampling_k, transformation_k]

Because ``SeedSequence.spawn`` children are deterministic *by index* and do not
depend on how many are spawned, replicate ``k`` is byte-identical whether it is
generated alone or inside a batch of any size (acceptance test 12).

This module computes a per-replicate statistics record and aggregates across
replicates (mean / sd / median / quantiles / success-fail counts), with a
machine-readable ``asdict`` representation. Finite-sampling variation and
transformation-algorithm variation are reported in clearly separated fields.
"""

from __future__ import annotations

import statistics
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field

HaplotypeKey = tuple[str, ...]


@dataclass(slots=True)
class ReplicateStatistics:
    """Statistics for a single replicate."""

    replicate: int
    n_individuals: int
    # Haplotype-frequency sampling error: RMS of (observed - input) over
    # haplotypes, a scalar summary of finite-sampling noise.
    haplotype_freq_rmse: float
    observed_homozygosity: float
    expected_homozygosity_m_based: float
    expected_homozygosity_pooled_hwe: float
    tv_distance_vs_m_based: float
    tv_distance_vs_pooled_hwe: float
    origin_pair_proportions: dict[str, float] = field(default_factory=dict)
    # Transformation (if enabled).
    transform_enabled: bool = False
    transform_reached: bool | None = None
    transform_achieved_homozygosity: float | None = None
    transform_reason: str | None = None


@dataclass(slots=True)
class AggregateField:
    mean: float
    sd: float
    median: float
    q05: float
    q25: float
    q75: float
    q95: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


def _quantile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 < len(sorted_vals):
        return sorted_vals[lo] * (1 - frac) + sorted_vals[lo + 1] * frac
    return sorted_vals[lo]


def _aggregate(values: Sequence[float]) -> AggregateField:
    vals = [float(v) for v in values if v is not None]
    if not vals:
        return AggregateField(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    s = sorted(vals)
    return AggregateField(
        mean=statistics.fmean(vals),
        sd=statistics.pstdev(vals) if len(vals) > 1 else 0.0,
        median=statistics.median(vals),
        q05=_quantile(s, 0.05),
        q25=_quantile(s, 0.25),
        q75=_quantile(s, 0.75),
        q95=_quantile(s, 0.95),
    )


@dataclass(slots=True)
class ReplicateBatchSummary:
    n_replicates: int
    n_success: int
    n_fail: int
    per_replicate: list[ReplicateStatistics]
    aggregates: dict[str, dict[str, float]]

    def to_dict(self) -> dict:
        return {
            "n_replicates": self.n_replicates,
            "n_success": self.n_success,
            "n_fail": self.n_fail,
            "per_replicate": [asdict(r) for r in self.per_replicate],
            "aggregates": self.aggregates,
            "notes": {
                "finite_sampling_variation": (
                    "haplotype_freq_rmse, tv_distance_vs_* and observed_homozygosity "
                    "reflect finite-sample multinomial noise around the expectation."
                ),
                "transformation_variation": (
                    "transform_achieved_homozygosity / transform_reached reflect the "
                    "deterministic reassignment algorithm (given a replicate's "
                    "sampled inventory), not sampling noise."
                ),
            },
        }


def aggregate_replicates(records: Sequence[ReplicateStatistics]) -> ReplicateBatchSummary:
    """Aggregate per-replicate records into mean/sd/median/quantiles + counts.

    A replicate is counted as a transformation *success* when a transform was
    enabled and ``transform_reached`` is True; when transforms are disabled every
    replicate is a (sampling) success.
    """
    numeric_fields = [
        "haplotype_freq_rmse",
        "observed_homozygosity",
        "expected_homozygosity_m_based",
        "expected_homozygosity_pooled_hwe",
        "tv_distance_vs_m_based",
        "tv_distance_vs_pooled_hwe",
    ]
    aggregates: dict[str, dict[str, float]] = {}
    for fname in numeric_fields:
        aggregates[fname] = _aggregate([getattr(r, fname) for r in records]).as_dict()

    transform_enabled = any(r.transform_enabled for r in records)
    if transform_enabled:
        achieved = [
            r.transform_achieved_homozygosity
            for r in records
            if r.transform_achieved_homozygosity is not None
        ]
        if achieved:
            aggregates["transform_achieved_homozygosity"] = _aggregate(achieved).as_dict()
        n_success = sum(1 for r in records if r.transform_reached)
        n_fail = sum(1 for r in records if r.transform_enabled and not r.transform_reached)
    else:
        n_success = len(records)
        n_fail = 0

    return ReplicateBatchSummary(
        n_replicates=len(records),
        n_success=n_success,
        n_fail=n_fail,
        per_replicate=list(records),
        aggregates=aggregates,
    )
