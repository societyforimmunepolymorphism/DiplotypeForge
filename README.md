# DiplotypeForge

Simulate synthetic populations of individuals with **known, phased diploid
genotypes (diplotypes)** from population **haplotype-frequency distributions**.

DiplotypeForge is a component of the Society for Immune Polymorphism (SIP)
**HLA Forge** initiative. Its companion project, **TissueTypeForge**,
independently simulates laboratory typing observations from known diplotypes.
The **phased GL String** is the output contract shared between the two projects.

> **Design philosophy: HLA-first, but not HLA-only.** The simulation engine makes
> no assumption about the number, names, or ordering of loci, and is applicable
> to other polymorphic systems (KIR, LILR) where appropriate.

## Scientific motivation

Synthetic populations with known genetic truth support:

- validation of haplotype-frequency-estimation algorithms;
- evaluation of HLA imputation algorithms;
- simulation of donor registries;
- investigation of population structure and Hardy–Weinberg disequilibrium;
- evaluation of methods for detecting mixture and admixture;
- reproducible benchmark datasets with known truth.

DiplotypeForge keeps three concepts distinct: the **population model** (input
truth), the **sampled diplotypes**, and the **derived summary statistics**.

## Status

**Milestones 1 (Foundation), 2 (Population Structure), and 3 (Controlled
Disequilibrium, Replicates & Scaling)** are implemented:

- internal typed data model;
- frequency-distribution validation (strict by default; explicit renormalization
  and duplicate-aggregation preprocessing);
- authoritative **HFX** input adapter (`societyforimmunepolymorphism/hfx` v0.1.1,
  schema-validated by default, lenient mode for historical files);
- **generalized K×K population-origin matrix** as the central sampling
  mechanism, with named configurations: single-population, Wahlund, first-gen
  admixture, and arbitrary K≥3 matrices (asymmetric supported);
- per-haplotype source-population provenance;
- **theoretical expectations independent of sampling**: marginal ancestry,
  population-pair proportions, and expected ordered/unordered diplotype
  distributions under M, plus pooled-HWE;
- diagnostics comparing observed diplotypes against both the M-based and
  pooled-HWE expectations;
- the **three-matrix experiment** (identical marginal ancestry, same marginal
  haplotype frequencies, different diplotype distributions) with a reproducible
  script, comparison table, and visualization;
- a **frequency-preserving homozygosity transform** that changes homozygosity
  while preserving every haplotype count exactly, with exact feasibility
  analysis and full reporting;
- **reproducible replicate batches** with per-replicate and aggregate statistics
  and a separate sampling/transformation RNG hierarchy;
- validated **≥ 1,000,000-individual** scaling with a documented benchmark;
- phased GL String serialization/parsing (round-trip);
- Python API + CLI (`validate` / `simulate` / `transform` / `summarize`);
  deterministic + stochastic tests; working examples.

See `docs/DEVELOPMENT_PLAN.md` for the full milestone history.

## Installation

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

```bash
uv venv --python 3.12
uv pip install -e ".[dev]"        # add ,viz for the experiment's matplotlib plot
```

## Quick start

```bash
# validate a configuration and its HFX inputs
uv run diplotypeforge validate examples/simulation.toml

# run the simulation; writes diplotypes + diagnostics + manifest
uv run diplotypeforge simulate examples/simulation.toml

# print diagnostics from an output directory
uv run diplotypeforge summarize examples/output
```

Python API:

```python
from diplotypeforge import Simulator, SimulationConfig

config = SimulationConfig.from_file("examples/simulation.toml")
result = Simulator(config).run()
result.write("examples/output/")
```

## Population models

All modes are special cases of a single **generalized K×K population-origin
matrix** `M`, where `M[i,j]` is the probability an individual's first haplotype
comes from population `i` and its second from population `j`.

- **Single-population random mating**: `M = [[1]]`; yields Hardy–Weinberg
  proportions `P(i,i)=f_i²`, `P(i≠j)=2·f_i·f_j` at the multilocus-haplotype level.
- **Wahlund stratification**: diagonal `M` (`M[k,k]=w_k`); both haplotypes from
  the same source. Produces population-level HW disequilibrium from structure.
- **First-generation admixture** (balanced two-population): `M = [[(1-a)/2, a/2],
  [a/2, (1-a)/2]]`. The admixed fraction `a` is distinct from marginal ancestry
  proportions.
- **Generalized K×K matrix**: arbitrary ordered origin matrix, K≥3, asymmetric
  supported.

Theoretical expectations (marginal ancestry, population-pair proportions,
expected ordered/unordered diplotype distributions under M, and pooled-HWE) are
computed **independently of the sampler** in `diplotypeforge.expectation`.

Full mathematical descriptions are in `docs/DESIGN.md`; the configuration
reference is in `docs/CONFIG.md`.

### Three-matrix experiment

A flagship experiment shows three matrices with identical marginal ancestry
producing the **same** marginal haplotype frequencies but **different** diplotype
distributions:

```bash
uv pip install -e ".[dev,viz]"
uv run python examples/three_matrix_experiment.py --n 50000
```

It writes a comparison table, `experiment.json`, and a PNG visualization. See
`docs/EXPERIMENT_THREE_MATRIX.md`.

Full mathematical descriptions are in `docs/DESIGN.md`.

## Controlled disequilibrium (frequency-preserving transform)

A post-sampling transform changes the number of homozygous diplotypes **while
preserving every haplotype count exactly** — a mechanism for Hardy–Weinberg
disequilibrium that is scientifically distinct from population structure. It
reports the initial / requested / achieved homozygosity, the exact feasible
`[min, max]` range, the number of operations, whether the target was reached,
and a reason for any rounding/clamping. Enable it with a `[transformation]`
block or run it standalone on an existing dataset:

```bash
diplotypeforge simulate  examples/transformation/simulate_transform.toml
diplotypeforge transform --input output/diplotypes.tsv \
                         --config examples/transformation/transform_only.toml
```

See `docs/TRANSFORM_WORKED_EXAMPLE.md` for an increase/decrease worked example
and the structure-vs-reassignment distinction.

## Replicates and scaling

Replicate batches use a `SeedSequence` hierarchy
(`root.spawn(replicates)` → each replicate `.spawn(2)` → `[sampling,
transformation]`), so replicate *k* is identical whether generated alone or in a
larger batch, and enabling the transform never changes the sampled diplotypes.
Per-replicate and aggregate statistics (mean/sd/median/quantiles, success/fail)
are written to `replicates.json`. Scaling is validated at **≥ 1,000,000
individuals**; see `docs/PERFORMANCE.md` and `examples/benchmark_large.py`.

## Output contract (shared with TissueTypeForge)

A phased diplotype is serialized as two `~`-joined haplotypes separated by `+`,
loci in declared order:

```
A*01:01g~C*07:01g~B*08:01g~DRB1*15:01~DQB1*06:02+A*02:01g~C*03:04g~B*40:01g~DRB1*04:04~DQB1*03:02g
```

Runs produce: `diplotypes.tsv[.gz]`, `haplotype_frequencies.tsv`,
`diplotype_frequencies.tsv`, `summary.json`, and a machine-readable
`manifest.json` (inputs, parameters, seed, generator, software version, files).

## Reproducibility

Identical inputs + seed + software version + configuration produce identical
output. The RNG is numpy's `Generator` (PCG64); replicates use independent
streams via `SeedSequence.spawn`; haplotypes are sorted into a frozen array so
results never depend on dictionary/JSON iteration order. The seed and generator
are recorded in every manifest.

## Limitations and assumptions

- Allele tokens are validated **syntactically**, not against an IPD database.
- One HFX file is treated as one source population (see
  `docs/OPEN_QUESTIONS.md`).
- Copy number is fixed at 2 per locus; gene absence / CNV is a future
  data-model extension (not a scalar annotation), not yet implemented.
- The frequency-preserving homozygosity transform reaches any *feasible* target
  but is **not** a uniform sampler over valid pairings (and does not claim to
  be); the transformed sample has no closed-form diplotype distribution.
- KIR copy-number modeling, recombination, and additional population models are
  out of scope.

## Citation

If you use DiplotypeForge, please cite the Wahlund and first-generation
admixture simulation designs:

> Maiers M, et al. *16th IHIW: Global analysis of registry HLA haplotypes from
> 20 million individuals.* Int J Immunogenet. 2012.

And the GL String grammar:

> Milius RP, et al. *Genotype List String: a grammar for describing HLA and KIR
> genotyping results in a text string.* Tissue Antigens. 2013.

## License

Apache-2.0.
