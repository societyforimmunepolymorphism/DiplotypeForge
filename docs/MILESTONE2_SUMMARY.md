# Milestone 2 — Implementation Summary

Milestone 2 makes the **generalized K×K population-origin matrix** the central
sampling mechanism and adds the theoretical machinery to reason about it
independently of the sampler.

## What was built

**Origin model (`sampling/origin.py`).** `OriginModel` validates a K×K matrix
(square, nonnegative, sums to 1) and samples ordered origin pairs `(i,j)`.
Named constructors: `single_population` (`[[1]]`), `wahlund` (diagonal from
weights), `admixture(a)` (balanced `[[(1-a)/2,a/2],[a/2,(1-a)/2]]`),
`from_matrix(M, symmetric=…)`. Asymmetric matrices are supported; `symmetric`
symmetrizes `(M+Mᵀ)/2`.

**Theoretical expectations (`expectation.py`, independent of sampling).** No RNG
anywhere in this module.
- `marginal_ancestry_proportions(M)` = `½(row+col)`.
- `ordered_population_pair_proportions(M)` = `M`;
  `unordered_population_pair_proportions(M)` sums both orders.
- `marginal_haplotype_frequencies(pops, M)` = `f̄(h)=Σ_i a_i f_i(h)`.
- `expected_diplotype_distribution(pops, M)` = ordered `P(h,k)=Σ_ij M[i,j]
  f_i(h) f_j(k)` (computed as `fᵀMf`) + unordered (both orders summed).
- `expected_diplotype_distribution_pooled_hwe(pops, M)` = `f̄(h)·f̄(k)`.

**Sampling.** The existing `Sampler` already drew ordered origin pairs from `M`
and recorded `origin_1`/`origin_2` per haplotype; Milestone 2 wires all four
modes through it. Per-haplotype provenance is preserved; asymmetric `M` gives
`origin_1` the row law and `origin_2` the column law.

**Diagnostics (`stats.py`).** Delegates all expectation math to
`expectation.py`. Reports observed vs **M-based** and vs **pooled-HWE**
(homozygosity/heterozygosity + total-variation distance over the diplotype
distribution), and observed vs **theoretical** origin and origin-pair
proportions.

**Config + Simulator + CLI.** `config.py` enables all four modes with per-mode
validation (population counts; `weights` length; `admixture_a∈[0,1]`; matrix
K×K dims match populations). `simulator.py` builds the right `OriginModel` and
adds an `origin_model` block (matrix, symmetric flag, population order, marginal
ancestry) to the manifest. `summary.json` now carries both expectations with TV
distances and observed-vs-theoretical origin proportions.

**Flagship experiment (`experiments.py` + `examples/three_matrix_experiment.py`).**
Two synthetic populations, three matrices with identical marginal ancestry.
Produces a comparison table, `experiment.json`, and a matplotlib PNG (optional
`[viz]` dependency). `ExperimentResult` exposes
`marginal_haplotype_frequencies_identical()` and
`diplotype_distributions_differ()`.

## Key result (verified)

All three matrices share marginal ancestry 0.5/0.5 and hence **identical**
expected marginal haplotype frequencies, but **different** diplotype
distributions. Expected homozygosity: **wahlund 0.46 > random 0.24 > admixture
0.108**. Simulated (n=50 000) observed homozygosity 0.457/0.240/0.106 with TV
distance to the M-based expectation of 0.004–0.008. A corrected modeling note:
the balanced admixture matrix `[[0.1,0.4],[0.4,0.1]]` has the *lowest*
homozygosity because its 0.8 cross-source mass pairs haplotypes from the two
(largely disjoint) populations into heterozygotes — the admixed fraction `a` is
not the ancestry proportion. `random` mixing (`M = a·aᵀ`) equals pooled-HWE
exactly.

## Verification

- `uv run ruff check .` — clean.
- `uv run pytest` — 75 passed (25 new Milestone-2 tests: matrix validation,
  named constructors, symmetric, marginal ancestry identical across the three
  matrices + an asymmetric case, ordered/unordered pair proportions, single-pop
  expectation = HWE, Wahlund full stratification, M-based ≠ pooled-HWE under
  structure, K=3 matrix sampling all three origins, asymmetric per-haplotype
  provenance, stochastic observed-vs-M TV < 8/√n, observed-vs-expected
  homozygosity within 6 SE, the experiment invariants, and random mixing =
  pooled-HWE).
- Three-matrix experiment runs end-to-end producing table, JSON, and PNG.
- Multi-population CLI path verified (a Wahlund two-population run: observed
  homozygosity 0.458 vs M-based 0.46, pooled-HWE 0.24; origin-pair proportions
  show only same-source pairs).

## Scope

Milestone 3 (frequency-preserving homozygosity transform, replicates/batch,
large-population performance, CNV groundwork) is **not** started. Milestone 2's
tests pass and the experiment runs, satisfying the gate.
