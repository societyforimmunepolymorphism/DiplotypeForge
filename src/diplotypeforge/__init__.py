"""DiplotypeForge: simulate phased diploid genotypes from haplotype frequencies.

Part of the Society for Immune Polymorphism (SIP) HLA Forge initiative.
"""

from __future__ import annotations

from .config import PopulationSpec, SimulationConfig
from .expectation import (
    DiplotypeExpectation,
    expected_diplotype_distribution,
    expected_diplotype_distribution_pooled_hwe,
    marginal_ancestry_proportions,
    marginal_haplotype_frequencies,
    ordered_population_pair_proportions,
    unordered_population_pair_proportions,
)
from .glstring import (
    diplotype_to_glstring,
    haplotype_to_glstring,
    parse_diplotype,
    parse_haplotype,
)
from .model.core import (
    Allele,
    Diplotype,
    Haplotype,
    Locus,
    SimulatedIndividual,
)
from .model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)
from .sampling.origin import OriginModel, OriginModelError
from .simulator import SimulationResult, Simulator

__version__ = "0.2.0"

__all__ = [
    "Allele",
    "Diplotype",
    "DiplotypeExpectation",
    "FrequencyEntry",
    "Haplotype",
    "HaplotypeFrequencyDistribution",
    "Locus",
    "OriginModel",
    "OriginModelError",
    "PopulationSpec",
    "SimulatedIndividual",
    "SimulationConfig",
    "SimulationResult",
    "Simulator",
    "SourcePopulation",
    "__version__",
    "diplotype_to_glstring",
    "expected_diplotype_distribution",
    "expected_diplotype_distribution_pooled_hwe",
    "haplotype_to_glstring",
    "marginal_ancestry_proportions",
    "marginal_haplotype_frequencies",
    "ordered_population_pair_proportions",
    "parse_diplotype",
    "parse_haplotype",
    "unordered_population_pair_proportions",
]
