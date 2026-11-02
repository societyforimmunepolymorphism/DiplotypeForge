# Development Plan

## Approach
Test-driven within each milestone. Deterministic acceptance tests using small
synthetic frequency tables. `ruff` for lint/format, `pytest` for tests, `uv` for
environment and dependency management. Python 3.12+.

## Milestone 1 — Foundation (current)
Order of implementation (each step verified before the next). Decisions marked
**(rev)** reflect the post-review revisions.

1. **Scaffold** — `pyproject.toml` (uv, py3.12+, numpy + jsonschema runtime;
   ruff + pytest dev), `src/diplotypeforge/` layout, ruff/pytest config,
   README skeleton. Vendor the authoritative HFX v0.1.1 schema.
2. **Data models** (`model/`) — Allele, Locus, Haplotype,
   HaplotypeFrequencyDistribution, SourcePopulation, Diplotype,
   SimulatedIndividual, plus config/result types. Typed, phase-preserving,
   provenance-aware. **(rev)** Copy number fixed at 2; CNV is a documented future
   data-model extension, not a scalar annotation (Q10) — not implemented.
3. **GL String** (`glstring.py`) — serialize/parse **unambiguous phased
   diplotypes only** (`~`,`+`); reject `/`,`|`,`^`. **(rev)** Fully-qualified
   allele ids (`HLA-…`). Round-trip.
4. **Validation** (`validation.py`) — **(rev)** strict by default: frequency sum
   only renormalized within a tiny round-off tolerance; larger renormalization
   requires explicit authorization and is recorded; **duplicate haplotypes are
   errors by default**; `aggregate_duplicates` is an explicit preprocessing op
   that warns + records provenance. Structured issue list.
5. **HFX adapter** (`adapters/hfx.py`) — **(rev)** validate against the vendored
   authoritative schema with `jsonschema` **by default**; explicit **lenient
   mode** for nonconforming historical files. Normalize allele ids to
   fully-qualified form via `nomenclatureUsed.database`. Inline/CSV + header
   mapping; preserve metadata + checksums. `hfx-tools` evaluated, not depended
   upon (see DESIGN §2.2).
6. **Sampling** (`sampling/`) — seeded numpy Generator (PCG64); generalized
   origin model with `M=[[1]]` for single population; independent haplotype draws
   → ordered diplotypes; streaming; replicate seed spawning.
7. **Stats** (`stats.py`) — **(rev)** empirical haplotype/diplotype freq; expected
   diplotype distribution from the **M-based model** `Σ_ij M[i,j] f_i(h) f_j(k)`
   (ordered + unordered); observed vs expected homozygosity/heterozygosity;
   report M-based expectation distinctly from pooled-HWE; sampling-error bands.
8. **Output** (`output.py`) — TSV + gzip TSV writers; **(rev)** manifest with
   input SHA-256 checksums, canonicalized config, RNG spec, software + numpy
   versions, renormalization/aggregation provenance.
9. **Config + Simulator + CLI** (`config.py`, `simulator.py`, `cli.py`) — TOML
   config (incl. strict/lenient/renorm/aggregate keys); Simulator/SimulationResult;
   `validate`/`simulate`/`summarize`.
10. **Tests + example** — acceptance tests 1,2,7,8,9,10 (+revised defaults); a
    working 1-population example end-to-end.

### M1 acceptance subset (from brief §10)
- (1) single-population recovers input haplotype frequencies within sampling
  expectation.
- (2) diplotype frequencies agree with HWE expectation.
- (7) all diplotypes contain valid pairs of complete haplotypes.
- (8) fixed seed ⇒ reproducible results.
- (9) GL String output preserves phase and parses back to equivalent diplotypes.
- (10) invalid configurations fail with informative errors.
(Tests 3,4,5,6-swap are Milestone 2/3.)

### M1 gate
`uv run pytest` green and `uv run ruff check` clean, plus the example produces
diplotypes + manifest. **Status: met and accepted.**

## Milestone 2 — Population structure (DONE)
The generalized K×K origin matrix is the central sampling mechanism; the four
named configurations are special cases of it.

1. **Modes wired through config + simulator** — `single_population` (`M=[[1]]`),
   `wahlund` (diagonal `M[k,k]=w_k`), `admixture` (balanced 2×2
   `[[(1-a)/2,a/2],[a/2,(1-a)/2]]`), `matrix` (arbitrary K×K, K≥3, asymmetric
   allowed). Per-mode config validation (population counts, weights length,
   `admixture_a∈[0,1]`, matrix K×K dims match populations).
2. **Matrix validation** — square dims, nonnegative entries, sums to 1
   (`OriginModel`). `symmetric=True` convenience symmetrizes.
3. **Theoretical expectations** (`expectation.py`, independent of sampling):
   marginal ancestry `a_i=½(row_i+col_i)`; ordered pair proportions `M`;
   unordered pair proportions; marginal haplotype freq `f̄(h)=Σ a_i f_i(h)`;
   expected ordered diplotype `P(h,k)=Σ_ij M[i,j] f_i(h) f_j(k)` and unordered
   (sum both orders); pooled-HWE `f̄(h)f̄(k)`.
4. **Per-haplotype provenance** — `origin_1`/`origin_2` recorded per individual;
   asymmetric matrices yield distinct origin_1 (row) vs origin_2 (column) laws.
5. **Diagnostics** (`stats.py`) — observed vs M-based and vs pooled-HWE
   (homozygosity/heterozygosity + TV distance over the diplotype distribution);
   observed vs theoretical origin and origin-pair proportions.
6. **Flagship three-matrix experiment** (`experiments.py` +
   `examples/three_matrix_experiment.py`) — comparison table + JSON + matplotlib
   PNG; demonstrates identical marginal ancestry + identical marginal haplotype
   frequencies but different diplotype distributions across random/Wahlund/
   admixture.

### M2 acceptance mapping (brief §10)
- (3) Wahlund reproduces requested source proportions — origin proportions vs
  theoretical, stochastic test.
- (4) first-gen admixture reproduces source-pair proportions — origin-pair test.
- (5) generalized matrix handles ≥3 populations — K=3 test.
Plus: matrix validation, asymmetric matrices, M-based vs pooled-HWE distinction,
deterministic expectation tests, sample-size-scaled stochastic tolerances.

### M2 gate
`uv run pytest` green, `uv run ruff check` clean, and the three-matrix
experiment runs producing the table + visualization. **Do not start M3 until
met.**

## Milestone 3 — Controlled disequilibrium, replicates, scaling (DONE)
Staged; M3A was verified before M3B/M3C. All stages complete.

### M3A — Frequency-preserving homozygosity transform
- `transform/feasibility.py`: exact achievability predicate `is_reachable(H)`
  (`Σ_k lo_k(H) ≤ H ≤ Σ⌊c_h/2⌋` with `lo_k(H)=max(0,⌈(c_k−(N−H))/2⌉)`); `H_min`,
  `H_max`, and `step` (1, or 2 for the two-haplotype parity case) derived
  independently of the algorithm and validated against brute force.
- `transform/reassignment.py`: `TransformationStrategy` interface + one robust
  `reassignment` strategy; elementary operations change homozygous count in
  fixed steps with bounded work (no unbounded rejection loops); optional RNG
  shuffle; provenance follows haplotypes.
- `TransformationResult`: initial/requested/rounded/clamped/achieved
  homozygosity, `H_min`/`H_max`, #operations, `reached`, `reason`.
- Diagnostics: original-M expectation, pooled-HWE, observed before/after, TV
  distances; transformed sample explicitly not M-distributed.
- Config `[transformation]` block; CLI `transform` subcommand (reads
  `diplotypes.tsv[.gz]` + manifest sidecar; origins from columns, not GL).
- Invariant/acceptance tests (brief §10, 1–10 & 14) incl. analytically
  tractable min/max. Worked increase/decrease example.
- **Gate:** tests pass + worked example runs before M3B.

### M3B — Replicates
- Seed hierarchy: `root.spawn(replicates)`, each child `.spawn(2)` →
  `[sampling_k, transformation_k]` (toggling transform leaves sampling
  identical). Replicate k identical alone or within a batch.
- Per-replicate stats + aggregate (mean/sd/median/quantiles/success-fail),
  machine-readable. Tests 11/12/13.

### M3C — Performance
- ≥1M individuals benchmark (runtime, peak memory, throughput, file sizes).
- Compact integer haplotype ids internally; streaming preferred; document
  memory complexity (streaming O(1) extra vs transformation-enabled O(N)).
- No premature optimization, no distributed computing.

Out of scope this milestone: KIR copy-number, recombination, new population
models.

## Release hardening — v0.1.0 (DONE)

No new population models. Correctness, interoperability, reproducibility, and
packaging only.

- **HFX schema correction (PR #10)** integrated; vendored schema used verbatim
  (no `oneOf`→`anyOf` mutation); upstream commit SHA + canonical checksum
  recorded in docs and manifests; examples made strictly conformant.
- **Independent scientific verification**: brute-force matching enumeration +
  Hypothesis for the feasibility math; analytical diagnostics tests for
  asymmetric M; per-locus multilocus diagnostics experiment; provenance and
  reproducibility verified by byte-level output checksums. One bug fixed
  (order-dependent unordered diplotype keying in the M-based expectation).
- **Performance docs corrected**: current `O(N)` vs. theoretical `O(H)` pairing;
  identity/provenance/external-memory distinctions; `tracemalloc` caveat. No
  performance rewrite.
- **Packaging**: version 0.1.0; PEP 621 metadata (license, SIP attribution,
  URLs, classifiers); HFX schema bundled in the wheel; clean-environment install
  verified. `docs/API_STABILITY.md` defines the supported v0.1 surface.
- **GL String interop**: `docs/GLSTRING_INTEROP.md` + round-trip corpus; parser/
  serializer interfaces identified for a future shared `glstring` library.
- **CI**: GitHub Actions (ruff, pytest, build, wheel smoke, CLI simulation) on
  Python 3.12 + 3.13; ≥1M benchmark manual-only.
- **Docs**: README, DESIGN, CONFIG, PERFORMANCE, OPEN_QUESTIONS, CHANGELOG, and
  a DiplotypeForge↔TissueTypeForge relationship doc.

### Gate
Full suite green, ruff clean, wheel builds and installs cleanly in a fresh
environment, end-to-end CLI simulation on a strictly valid HFX file, and
`RELEASE_READINESS.md` produced. **Status: met.**
