"""Seeded haplotype/diplotype sampling.

Given one or more :class:`SourcePopulation` objects and an :class:`OriginModel`,
draw individuals as ordered haplotype pairs with per-haplotype provenance.

Reproducibility guarantees (see docs/DESIGN.md §5):
  * numpy ``Generator`` backed by ``PCG64``, explicitly seeded;
  * haplotypes placed in a sorted, frozen array so results never depend on
    input/dict iteration order;
  * replicates derive independent child streams via ``SeedSequence.spawn``.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np

from ..model.core import Diplotype, Haplotype, SimulatedIndividual
from ..model.population import SourcePopulation
from .origin import OriginModel


@dataclass(frozen=True, slots=True)
class _PreparedPopulation:
    """A population reduced to a sorted haplotype array + probability vector."""

    population_id: str
    haplotypes: tuple[Haplotype, ...]
    probabilities: np.ndarray  # normalized, aligned with haplotypes

    def sample_indices(self, n: int, rng: np.random.Generator) -> np.ndarray:
        return rng.choice(len(self.haplotypes), size=n, p=self.probabilities)


def _prepare_population(pop: SourcePopulation) -> _PreparedPopulation:
    """Freeze a population into a deterministically ordered sampling table.

    Entries are sorted by their GL-String key so iteration/dict order in the
    input cannot affect sampling outcomes. Zero/negative frequencies are dropped
    and the remainder renormalized to sum to 1.
    """
    entries = [e for e in pop.distribution.entries if e.frequency > 0]
    if not entries:
        raise ValueError(f"population {pop.population_id!r} has no positive-frequency haplotypes")
    entries.sort(key=lambda e: e.haplotype.key())
    haplotypes = tuple(e.haplotype for e in entries)
    freqs = np.array([e.frequency for e in entries], dtype=np.float64)
    freqs /= freqs.sum()
    return _PreparedPopulation(pop.population_id, haplotypes, freqs)


def make_rng(seed: int) -> tuple[np.random.Generator, str]:
    """Create a seeded generator and return it with its generator name."""
    bitgen = np.random.PCG64(seed)
    return np.random.Generator(bitgen), "PCG64"


class Sampler:
    """Draws :class:`SimulatedIndividual` objects from prepared populations."""

    def __init__(self, populations: list[SourcePopulation], origin: OriginModel) -> None:
        if origin.k != len(populations):
            raise ValueError(
                f"origin model is {origin.k}x{origin.k} but {len(populations)} "
                "populations were provided"
            )
        self.populations = populations
        self.prepared = [_prepare_population(p) for p in populations]
        self.origin = origin

    def sample(
        self,
        n_individuals: int,
        rng: np.random.Generator,
        *,
        id_prefix: str,
        replicate: int = 0,
    ) -> Iterator[SimulatedIndividual]:
        """Yield ``n_individuals`` sampled individuals (streaming)."""
        if n_individuals < 0:
            raise ValueError("n_individuals must be non-negative")

        origin_pairs = self.origin.sample_origin_pairs(n_individuals, rng)

        # Pre-draw haplotype indices per population in bulk for efficiency.
        # Count how many first/second draws each population needs.
        n_pops = len(self.prepared)
        draws_1: list[np.ndarray] = [np.empty(0, dtype=int)] * n_pops
        draws_2: list[np.ndarray] = [np.empty(0, dtype=int)] * n_pops
        mask_1 = [origin_pairs[:, 0] == p for p in range(n_pops)]
        mask_2 = [origin_pairs[:, 1] == p for p in range(n_pops)]
        for p in range(n_pops):
            draws_1[p] = self.prepared[p].sample_indices(int(mask_1[p].sum()), rng)
            draws_2[p] = self.prepared[p].sample_indices(int(mask_2[p].sum()), rng)

        cursor_1 = [0] * n_pops
        cursor_2 = [0] * n_pops
        width = max(len(str(n_individuals - 1)), 1) if n_individuals else 1

        for idx in range(n_individuals):
            i, j = int(origin_pairs[idx, 0]), int(origin_pairs[idx, 1])
            hi = draws_1[i][cursor_1[i]]
            cursor_1[i] += 1
            hj = draws_2[j][cursor_2[j]]
            cursor_2[j] += 1
            hap1 = self.prepared[i].haplotypes[hi]
            hap2 = self.prepared[j].haplotypes[hj]
            diplotype = Diplotype(haplotype_1=hap1, haplotype_2=hap2)
            sample_id = f"{id_prefix}-{replicate}-{idx:0{width}d}"
            yield SimulatedIndividual(
                sample_id=sample_id,
                diplotype=diplotype,
                origin_1=self.prepared[i].population_id,
                origin_2=self.prepared[j].population_id,
                replicate=replicate,
            )
