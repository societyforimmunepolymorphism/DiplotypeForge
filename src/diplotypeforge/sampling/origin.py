"""Generalized population-origin model (the K x K matrix).

The origin model governs, per individual, which source population each of the
two haplotypes is drawn from. All simulation modes are special cases:

  * single population : M = [[1]]
  * Wahlund           : M diagonal, M[k,k] = w_k
  * 1st-gen admixture : M = [[(1-a)/2, a/2], [a/2, (1-a)/2]]
  * generalized       : arbitrary M (K x K), nonnegative, sums to 1

Only the single-population case is exercised in Milestone 1; the others are
provided as constructors and validated here so Milestone 2 can wire them in.
"""

from __future__ import annotations

import numpy as np


class OriginModelError(ValueError):
    """Raised for an invalid origin matrix specification."""


class OriginModel:
    """A validated K x K ordered population-origin probability matrix."""

    def __init__(self, matrix: np.ndarray | list[list[float]], tolerance: float = 1e-9) -> None:
        m = np.asarray(matrix, dtype=np.float64)
        if m.ndim != 2 or m.shape[0] != m.shape[1]:
            raise OriginModelError(f"origin matrix must be square 2-D; got shape {m.shape}")
        if np.any(m < 0):
            raise OriginModelError("origin matrix entries must be nonnegative")
        total = float(m.sum())
        if abs(total - 1.0) > tolerance:
            raise OriginModelError(f"origin matrix must sum to 1.0; sums to {total:.12g}")
        self.matrix = m
        self.k = m.shape[0]
        # Flattened categorical distribution over ordered (i, j) pairs.
        self._flat = m.reshape(-1)

    def sample_origin_pairs(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Return an ``(n, 2)`` int array of ordered ``(i, j)`` origin indices."""
        flat_idx = rng.choice(self._flat.size, size=n, p=self._flat)
        i = flat_idx // self.k
        j = flat_idx % self.k
        return np.stack([i, j], axis=1)

    @property
    def marginal_origin_proportions(self) -> np.ndarray:
        """Per-population marginal ancestry proportions over both haplotypes."""
        row = self.matrix.sum(axis=1)
        col = self.matrix.sum(axis=0)
        return (row + col) / 2.0

    # --- constructors for the named modes -----------------------------------

    @classmethod
    def single_population(cls) -> OriginModel:
        return cls(np.array([[1.0]]))

    @classmethod
    def wahlund(cls, weights: list[float]) -> OriginModel:
        w = np.asarray(weights, dtype=np.float64)
        if w.ndim != 1 or w.size < 1:
            raise OriginModelError("wahlund weights must be a 1-D non-empty vector")
        s = float(w.sum())
        if s <= 0:
            raise OriginModelError("wahlund weights must sum to a positive value")
        w = w / s
        return cls(np.diag(w))

    @classmethod
    def admixture(cls, a: float) -> OriginModel:
        if not 0.0 <= a <= 1.0:
            raise OriginModelError(f"admixture fraction a must be in [0, 1]; got {a}")
        off = a / 2.0
        diag = (1.0 - a) / 2.0
        return cls(np.array([[diag, off], [off, diag]]))

    @classmethod
    def from_matrix(cls, matrix: list[list[float]], *, symmetric: bool = False) -> OriginModel:
        m = np.asarray(matrix, dtype=np.float64)
        if symmetric:
            m = (m + m.T) / 2.0
        return cls(m)
