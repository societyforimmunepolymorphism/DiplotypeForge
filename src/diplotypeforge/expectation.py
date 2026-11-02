"""Theoretical expectations under a population-origin matrix M.

**This module is independent of the sampling implementation.** Everything here
is computed analytically from the origin matrix ``M`` and the per-population
haplotype-frequency distributions `f_i(h)`; nothing draws random numbers. The
sampler is validated *against* these expectations, not the other way round.

Definitions (see docs/DESIGN.md §4.3):

- Marginal ancestry proportion of population ``i``:
      a_i = ½ (Σ_j M[i,j] + Σ_j M[j,i])
  i.e. the expected fraction of haplotypes (over both chromosomes) originating
  from population ``i``.

- Ordered population-pair proportion: ``M[i,j]`` itself.
  Unordered population-pair proportion for ``i ≠ j``: ``M[i,j] + M[j,i]``;
  for ``i = j``: ``M[i,i]``.

- Expected ordered diplotype distribution:
      P(h, k) = Σ_i Σ_j M[i,j] · f_i(h) · f_j(k)
  Expected unordered (genotype) distribution:
      P_unordered({h,k}) = P(h,k) + P(k,h)  for h ≠ k;  P(h,h) for h = k.

- Marginal (pooled) haplotype frequency:
      f̄(h) = Σ_i a_i · f_i(h)
  Pooled-HWE expected diplotype distribution:
      P_HWE(h, k) = f̄(h) · f̄(k).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .model.population import SourcePopulation

HaplotypeKey = tuple[str, ...]


def normalized_matrix(matrix: np.ndarray | list[list[float]]) -> np.ndarray:
    """Return M as a float64 array; caller guarantees validation via OriginModel."""
    return np.asarray(matrix, dtype=np.float64)


def marginal_ancestry_proportions(matrix: np.ndarray) -> np.ndarray:
    """a_i = ½ (row_i + col_i)."""
    m = np.asarray(matrix, dtype=np.float64)
    return (m.sum(axis=1) + m.sum(axis=0)) / 2.0


def ordered_population_pair_proportions(matrix: np.ndarray) -> np.ndarray:
    """Ordered pair proportions — exactly M."""
    return np.asarray(matrix, dtype=np.float64)


def unordered_population_pair_proportions(matrix: np.ndarray) -> dict[tuple[int, int], float]:
    """Unordered pair proportions keyed by ``(i, j)`` with ``i <= j``."""
    m = np.asarray(matrix, dtype=np.float64)
    k = m.shape[0]
    out: dict[tuple[int, int], float] = {}
    for i in range(k):
        for j in range(i, k):
            out[(i, j)] = float(m[i, j]) if i == j else float(m[i, j] + m[j, i])
    return out


def _haplotype_union_and_freqs(
    populations: Sequence[SourcePopulation],
) -> tuple[list[HaplotypeKey], list[dict[HaplotypeKey, float]]]:
    """Return the ordered union of haplotype keys and per-population freq maps
    (each normalized to sum to 1)."""
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
    return keys, per_pop


def marginal_haplotype_frequencies(
    populations: Sequence[SourcePopulation], matrix: np.ndarray
) -> dict[HaplotypeKey, float]:
    """f̄(h) = Σ_i a_i f_i(h) — the pooled/marginal haplotype frequency."""
    keys, per_pop = _haplotype_union_and_freqs(populations)
    a = marginal_ancestry_proportions(matrix)
    out: dict[HaplotypeKey, float] = {}
    for key in keys:
        out[key] = float(sum(a[i] * per_pop[i].get(key, 0.0) for i in range(len(per_pop))))
    return out


@dataclass(slots=True)
class DiplotypeExpectation:
    """Expected diplotype distributions, independent of sampling."""

    haplotype_keys: list[HaplotypeKey]
    ordered: dict[tuple[HaplotypeKey, HaplotypeKey], float] = field(default_factory=dict)
    unordered: dict[tuple[HaplotypeKey, HaplotypeKey], float] = field(default_factory=dict)

    def expected_homozygosity(self) -> float:
        return float(sum(p for (h, k), p in self.ordered.items() if h == k))

    def expected_heterozygosity(self) -> float:
        return 1.0 - self.expected_homozygosity()


def expected_diplotype_distribution(
    populations: Sequence[SourcePopulation], matrix: np.ndarray
) -> DiplotypeExpectation:
    """Compute the expected ordered and unordered diplotype distributions under M.

    Returns full distributions over the union of haplotypes. Complexity is
    O(K^2 · H^2); intended for the analytically tractable cases and experiments,
    not million-haplotype panels.
    """
    keys, per_pop = _haplotype_union_and_freqs(populations)
    m = np.asarray(matrix, dtype=np.float64)
    k_pops = m.shape[0]

    # Per-population frequency vectors aligned to the haplotype union.
    idx = {key: n for n, key in enumerate(keys)}
    f = np.zeros((k_pops, len(keys)), dtype=np.float64)
    for i in range(k_pops):
        for key, val in per_pop[i].items():
            f[i, idx[key]] = val

    # Ordered: P(h,k) = Σ_ij M[i,j] f_i(h) f_j(k) = (fᵀ M f)[h,k] where f is KxH.
    # p_ordered[h,k] = Σ_i Σ_j f[i,h] M[i,j] f[j,k]
    p_ordered = f.T @ m @ f  # (H x H)

    ordered: dict[tuple[HaplotypeKey, HaplotypeKey], float] = {}
    unordered: dict[tuple[HaplotypeKey, HaplotypeKey], float] = {}
    h = len(keys)
    for a in range(h):
        for b in range(h):
            val = float(p_ordered[a, b])
            if val > 0:
                ordered[(keys[a], keys[b])] = val
    for a in range(h):
        for b in range(a, h):
            val = float(p_ordered[a, a]) if a == b else float(p_ordered[a, b] + p_ordered[b, a])
            if val > 0:
                # Canonicalize the unordered key by sorting the two haplotype
                # keys, so it matches the observed keying (Diplotype.unordered_key
                # sorts k1 <= k2) regardless of the haplotype *union* order, which
                # itself depends on input row order. Without this, the TV distance
                # to observed data would spuriously depend on input ordering.
                ka, kb = keys[a], keys[b]
                ukey = (ka, kb) if ka <= kb else (kb, ka)
                unordered[ukey] = val

    return DiplotypeExpectation(haplotype_keys=keys, ordered=ordered, unordered=unordered)


def expected_diplotype_distribution_pooled_hwe(
    populations: Sequence[SourcePopulation], matrix: np.ndarray
) -> DiplotypeExpectation:
    """Pooled-HWE expectation P_HWE(h,k) = f̄(h) f̄(k) (random union of gametes).

    This is the null that *ignores* structure: it uses only the marginal
    haplotype frequencies. Compare against :func:`expected_diplotype_distribution`
    to quantify the Wahlund/structure signal.
    """
    fbar = marginal_haplotype_frequencies(populations, matrix)
    keys = list(fbar)
    ordered: dict[tuple[HaplotypeKey, HaplotypeKey], float] = {}
    unordered: dict[tuple[HaplotypeKey, HaplotypeKey], float] = {}
    for a in keys:
        for b in keys:
            val = fbar[a] * fbar[b]
            if val > 0:
                ordered[(a, b)] = val
    for i, a in enumerate(keys):
        for b in keys[i:]:
            val = fbar[a] * fbar[b] if a == b else 2.0 * fbar[a] * fbar[b]
            if val > 0:
                unordered[(a, b) if a <= b else (b, a)] = val
    return DiplotypeExpectation(haplotype_keys=keys, ordered=ordered, unordered=unordered)
