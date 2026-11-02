"""Milestone 3C: streaming + moderate-scale sanity (fast regression guards).

The full ≥1M benchmark lives in examples/benchmark_large.py; here we only assert
the streaming generator contract and that a moderate-scale transform preserves
all counts (so the scaling path keeps working without running for minutes).
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from diplotypeforge.model.core import Allele, Haplotype, Locus, SimulatedIndividual
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from diplotypeforge.sampling.origin import OriginModel
from diplotypeforge.sampling.sampler import Sampler, make_rng
from diplotypeforge.transform import TransformationTarget, transform_homozygosity


def _population(n_haps: int, seed: int) -> SourcePopulation:
    rng = np.random.default_rng(seed)
    w = rng.random(n_haps)
    w /= w.sum()
    entries = tuple(
        FrequencyEntry(Haplotype((Allele("HLA-A", f"HLA-A*{i:03d}:01"),)), float(w[i]))
        for i in range(n_haps)
    )
    dist = HaplotypeFrequencyDistribution(loci=(Locus("HLA-A", "g"),), entries=entries)
    return SourcePopulation("BENCH", dist)


def test_sampler_is_a_streaming_generator():
    pop = _population(50, 1)
    sampler = Sampler([pop], OriginModel.single_population())
    rng, _ = make_rng(1)
    gen = sampler.sample(1000, rng, id_prefix="S", replicate=0)
    # It is lazy: consuming one item at a time works (streaming contract).
    first = next(gen)
    assert isinstance(first, SimulatedIndividual)
    rest = list(gen)
    assert len(rest) == 999


def test_moderate_scale_transform_preserves_counts():
    pop = _population(200, 7)
    sampler = Sampler([pop], OriginModel.single_population())
    rng, _ = make_rng(7)
    individuals = list(sampler.sample(20000, rng, id_prefix="S", replicate=0))

    def counts(p):
        c: Counter = Counter()
        for x in p:
            c[x.diplotype.haplotype_1.key()] += 1
            c[x.diplotype.haplotype_2.key()] += 1
        return dict(c)

    before = counts(individuals)
    transformed, result = transform_homozygosity(
        individuals, TransformationTarget(fraction=0.25), rng=np.random.default_rng(7)
    )
    assert counts(transformed) == before
    assert result.achieved_count == result.clamped_target_count
    assert len(transformed) == 20000
