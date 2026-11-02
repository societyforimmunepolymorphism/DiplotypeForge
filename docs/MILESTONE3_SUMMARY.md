# Milestone 3 — Implementation Summary

Milestone 3 adds a complementary, **structure-independent** mechanism for
Hardy–Weinberg disequilibrium (a frequency-preserving homozygosity transform),
reproducible replicate batches, and validated large-population scaling. It was
delivered in the mandated order **M3A → M3B → M3C**, with M3A verified before
the later stages.

## M3A — Frequency-preserving homozygosity transform

**Feasibility (`transform/feasibility.py`), independent of the algorithm.** For
haplotype counts `c_h` (`Σ=2N`), a target of `H` homozygous pairs is achievable
iff `Σ_k lo_k(H) ≤ H ≤ Σ⌊c_k/2⌋` with `lo_k(H)=max(0,⌈(c_k−(N−H))/2⌉)`.
`analyze_feasibility` returns `H_min`, `H_max`, and the reachable `step`
(1 generally; **2** — a parity constraint — when only two distinct haplotypes
carry all copies). `is_reachable(H)` exposes the exact predicate.

> A correction made during implementation: the naive `H_max = Σ⌊c_k/2⌋` is only
> an *upper bound* — near the maximum the leftover copies can fail to be
> heterozygous-pairable, so the true `H_max` may be 1–2 lower. The predicate was
> validated against brute force over all count distributions for small N.

**Transform (`transform/reassignment.py`).** A `TransformationStrategy` interface
with one robust **`reassignment`** strategy. It is *constructive*: it groups the
`2N` copies by haplotype (each copy carrying its origin), forms exactly `H`
homozygous pairs chosen with per-key lower bounds that keep the remainder
heterozygous-pairable, then pairs the remainder by joining the two most frequent
distinct keys. This hits any reachable `H` with no rejection loops. Targets are a
`fraction` or a `delta`; rounding is half-even on a `1/N` grid, clamped to
`[H_min,H_max]`, then snapped to the nearest reachable value — and the
`TransformationResult` reports requested/rounded/clamped/achieved homozygosity,
feasible min/max, #operations, `reached`, and a `reason` (never a silent
substitution). The input is never mutated. **It is not a uniform sampler over
pairings and does not claim to be** (optional RNG shuffle only avoids
deterministic artifacts).

**Provenance.** Each copy carries its source-population origin, which travels
with it during reassignment; an individual's population-pair label may change.
Origin-label totals are preserved.

**Diagnostics (`transform/diagnostics.py`).** Reports expected-under-original-M,
pooled-HWE, observed-before, and observed-after diplotype distributions with all
the relevant total-variation distances, and **explicitly labels the transformed
sample as not M-distributed** (no invented analytical distribution).

**Config + CLI.** A `[transformation]` block (`enabled`, `type`, `strategy`,
`target_fraction`/`target_delta`) wires the transform into `simulate`; a
standalone `transform` subcommand reads an existing `diplotypes.tsv[.gz]`
(provenance from the `origin_*` columns, not GL Strings) plus the `manifest.json`
sidecar, and writes `diplotypes_transformed.tsv` + `transform_report.json`. When
a transform runs inside `simulate`, the untransformed sample is retained as
`diplotypes_original.tsv`.

## M3B — Replicates

The seed hierarchy is `root = SeedSequence(seed)`, `root.spawn(R)` per replicate,
each replicate `.spawn(2) → [sampling, transformation]`. Because `spawn`
children are deterministic by index, replicate *k* is identical whether
generated alone or inside a batch of any size; and because sampling is derived
independently of transformation, **enabling the transform never changes the
sampled diplotypes** (verified byte-for-byte). `replicates.py` computes a
per-replicate statistics record (sampling-error RMSE, observed/expected
homozygosity, TV distances, origin-pair proportions, transform outcome) and an
aggregate (mean/sd/median/quantiles, success/fail), written to `replicates.json`,
with finite-sampling vs transformation-algorithm variation clearly separated.

## M3C — Performance

Validated at **≥ 1,000,000 individuals** via `examples/benchmark_large.py`
(runtime split into sampling/transform/write, peak memory via `tracemalloc`,
throughput, output file sizes → `benchmark.json`). Measured at 1M (5 loci, 2000
haplotypes): ordinary ≈ 39 s / 25.5k ind/s / 217 MB peak; transformation-enabled
≈ 89 s / 11.2k ind/s / 649 MB peak (target reached exactly). This confirms the
documented memory complexity: ordinary simulation can stream with `O(1)` extra
memory, whereas the frequency-preserving transform is inherently `O(N)` because
it must hold the whole haplotype inventory to regroup it. The sampler already
addresses haplotypes by compact integer index; `docs/PERFORMANCE.md` records the
analysis and notes key interning as a future, low-risk optimization. No premature
optimization or distributed computing was introduced.

## Verification

- `uv run ruff check .` — clean.
- `uv run pytest` — **99 passed** (16 M3A, 6 M3B, 2 M3C added this milestone).
- M3A acceptance (brief §10): 1–3 counts/alleles/N preserved; 4/5 increase &
  decrease reach feasible targets; 6 infeasible reported (not silent); 7
  terminates at boundaries; 8 input unchanged; 9 provenance follows; 10
  reproducible; 14 valid phased GL round-trip; plus analytically tractable
  min/max and an exhaustive reachable-target check.
- M3B acceptance: 11 independent & reproducible; 12 replicate *k* identical alone
  vs in a batch of 100; 13 sampled diplotypes unchanged with transform on/off.
- End-to-end CLI: `simulate` with `[transformation]` (0.101 → 0.30 exact) and
  standalone `transform` (0.30 → 0.0 exact) both verified; `diplotypes_original`
  retained; diagnostics show the transformed sample moving away from M.

## Scope

No KIR copy-number modeling, recombination, or new population models were added.
No remote repositories were created and nothing was pushed or published.
