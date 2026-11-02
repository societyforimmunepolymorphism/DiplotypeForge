"""The three-matrix population-structure experiment (Milestone 2 flagship).

Two synthetic source populations with *distinct* haplotype distributions are
combined under three different origin matrices that share the **same marginal
ancestry proportions** (0.5 / 0.5):

    random mixing : [[0.25, 0.25], [0.25, 0.25]]
    Wahlund       : [[0.5, 0.0 ], [0.0 , 0.5 ]]
    admixture     : [[0.1, 0.4 ], [0.4 , 0.1 ]]   (balanced, a = 0.8)

The experiment demonstrates that all three yield the **same expected marginal
haplotype frequencies** f̄(h) (because marginal ancestry is identical), yet
**different expected diplotype distributions** (because the joint origin
structure differs). Expectations come from :mod:`diplotypeforge.expectation`
and are independent of sampling; optional simulation confirms them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .expectation import (
    DiplotypeExpectation,
    expected_diplotype_distribution,
    marginal_ancestry_proportions,
    marginal_haplotype_frequencies,
)
from .model.core import Allele, Haplotype, Locus
from .model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from .sampling.origin import OriginModel
from .sampling.sampler import Sampler, make_rng
from .stats import compute_summary

# Single-locus synthetic system keeps the experiment analytically tractable and
# the diplotype table small, while still exhibiting the structure effect.
_LOCUS = Locus("HLA-A", "g")

# Two populations with DISTINCT haplotype distributions (disjoint-ish supports).
POP1_FREQS = {"HLA-A*01:01": 0.7, "HLA-A*02:01": 0.2, "HLA-A*03:01": 0.1}
POP2_FREQS = {"HLA-A*03:01": 0.2, "HLA-A*11:01": 0.5, "HLA-A*24:02": 0.3}

NAMED_MATRICES: dict[str, list[list[float]]] = {
    "random": [[0.25, 0.25], [0.25, 0.25]],
    "wahlund": [[0.5, 0.0], [0.0, 0.5]],
    "admixture": [[0.1, 0.4], [0.4, 0.1]],
}


def _make_population(pop_id: str, freqs: dict[str, float]) -> SourcePopulation:
    entries = tuple(
        FrequencyEntry(Haplotype((Allele(_LOCUS.name, hap),)), freq)
        for hap, freq in freqs.items()
    )
    dist = HaplotypeFrequencyDistribution(loci=(_LOCUS,), entries=entries)
    return SourcePopulation(population_id=pop_id, distribution=dist)


def build_experiment_populations() -> list[SourcePopulation]:
    return [_make_population("POP1", POP1_FREQS), _make_population("POP2", POP2_FREQS)]


@dataclass(slots=True)
class ModelResult:
    name: str
    matrix: list[list[float]]
    marginal_ancestry: dict[str, float]
    marginal_haplotype_frequencies: dict[str, float]
    expected: DiplotypeExpectation
    observed_diplotype_frequencies: dict[str, float] = field(default_factory=dict)
    observed_homozygosity: float | None = None
    tv_distance_observed_vs_m: float | None = None


@dataclass(slots=True)
class ExperimentResult:
    populations: list[SourcePopulation]
    models: dict[str, ModelResult]

    def haplotype_order(self) -> list[tuple[str, ...]]:
        keys: list[tuple[str, ...]] = []
        seen: set[tuple[str, ...]] = set()
        for pop in self.populations:
            for e in pop.distribution.entries:
                k = e.haplotype.key()
                if k not in seen:
                    seen.add(k)
                    keys.append(k)
        return keys

    def marginal_haplotype_frequencies_identical(self, tol: float = 1e-12) -> bool:
        """True iff all models share the same marginal haplotype frequencies."""
        names = list(self.models)
        ref = self.models[names[0]].marginal_haplotype_frequencies
        for name in names[1:]:
            other = self.models[name].marginal_haplotype_frequencies
            if set(other) != set(ref):
                return False
            if any(abs(other[k] - ref[k]) > tol for k in ref):
                return False
        return True

    def diplotype_distributions_differ(self, tol: float = 1e-9) -> bool:
        """True iff at least one pair of models has distinct expected diplotype dists."""
        names = list(self.models)
        for a in range(len(names)):
            for b in range(a + 1, len(names)):
                da = self.models[names[a]].expected.unordered
                db = self.models[names[b]].expected.unordered
                keys = set(da) | set(db)
                if any(abs(da.get(k, 0.0) - db.get(k, 0.0)) > tol for k in keys):
                    return True
        return False


def run_experiment(
    n_individuals: int = 0,
    seed: int = 20120815,
    matrices: dict[str, list[list[float]]] | None = None,
) -> ExperimentResult:
    """Run the three-matrix experiment.

    With ``n_individuals == 0`` only the theoretical expectations are computed
    (fully deterministic, no RNG). With ``n_individuals > 0`` each model is also
    simulated to produce observed diplotype frequencies for comparison.
    """
    mats = matrices or NAMED_MATRICES
    populations = build_experiment_populations()
    ids = [p.population_id for p in populations]

    models: dict[str, ModelResult] = {}
    for name, matrix in mats.items():
        m = np.asarray(matrix, dtype=np.float64)
        anc = marginal_ancestry_proportions(m)
        fbar = marginal_haplotype_frequencies(populations, m)
        expected = expected_diplotype_distribution(populations, m)

        result = ModelResult(
            name=name,
            matrix=[list(row) for row in matrix],
            marginal_ancestry=dict(zip(ids, anc.tolist(), strict=True)),
            marginal_haplotype_frequencies={"~".join(k): v for k, v in fbar.items()},
            expected=expected,
        )

        if n_individuals > 0:
            origin = OriginModel.from_matrix(matrix)
            sampler = Sampler(populations, origin)
            # Independent reproducible stream per model via SeedSequence.spawn.
            name_seed = int(
                np.random.SeedSequence([seed, abs(hash(name)) % (2**32)]).generate_state(1)[0]
            )
            rng, _ = make_rng(name_seed)
            individuals = list(
                sampler.sample(n_individuals, rng, id_prefix=f"EXP-{name}", replicate=0)
            )
            stats = compute_summary(individuals, populations, m)
            result.observed_diplotype_frequencies = stats.empirical_diplotype_frequencies
            result.observed_homozygosity = stats.observed_homozygosity
            result.tv_distance_observed_vs_m = stats.tv_distance_vs_m_based

        models[name] = result

    return ExperimentResult(populations=populations, models=models)


def comparison_table(result: ExperimentResult) -> str:
    """Render a plain-text comparison table of the three models."""
    names = list(result.models)
    lines: list[str] = []

    lines.append("=== Three-matrix population-structure experiment ===")
    lines.append("")
    lines.append("Marginal ancestry proportions (should be identical across models):")
    for name in names:
        anc = result.models[name].marginal_ancestry
        lines.append(f"  {name:10s}: " + ", ".join(f"{p}={v:.3f}" for p, v in anc.items()))
    lines.append("")

    # Marginal haplotype frequencies table (identical across models).
    hap_keys = ["~".join(k) for k in result.haplotype_order()]
    lines.append("Expected marginal haplotype frequencies f_bar(h) (identical across models):")
    header = "  haplotype".ljust(24) + "".join(f"{name:>12s}" for name in names)
    lines.append(header)
    for hk in hap_keys:
        row = f"  {hk}".ljust(24)
        for name in names:
            row += f"{result.models[name].marginal_haplotype_frequencies.get(hk, 0.0):12.5f}"
        lines.append(row)
    lines.append("")

    # Expected unordered diplotype distribution (differs across models).
    all_diplos: list[str] = []
    seen: set[str] = set()
    for name in names:
        for ukey in result.models[name].expected.unordered:
            gl = "+".join("~".join(h) for h in ukey)
            if gl not in seen:
                seen.add(gl)
                all_diplos.append(gl)
    all_diplos.sort()

    lines.append("Expected UNORDERED diplotype frequencies under M (differ across models):")
    lines.append("  diplotype".ljust(40) + "".join(f"{name:>12s}" for name in names))
    for gl in all_diplos:
        row = f"  {gl}".ljust(40)
        for name in names:
            exp = result.models[name].expected
            glmap = {
                "+".join("~".join(h) for h in k): v for k, v in exp.unordered.items()
            }
            row += f"{glmap.get(gl, 0.0):12.5f}"
        lines.append(row)
    lines.append("")

    lines.append("Expected homozygosity (M-based):")
    for name in names:
        lines.append(f"  {name:10s}: {result.models[name].expected.expected_homozygosity():.5f}")

    if result.models[names[0]].observed_homozygosity is not None:
        lines.append("")
        lines.append("Observed homozygosity (simulated) and TV distance to M-based expectation:")
        for name in names:
            mr = result.models[name]
            lines.append(
                f"  {name:10s}: obs_hom={mr.observed_homozygosity:.5f}  "
                f"TV(obs, M)={mr.tv_distance_observed_vs_m:.5f}"
            )

    return "\n".join(lines)
