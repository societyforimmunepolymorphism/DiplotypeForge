# Release readiness — DiplotypeForge v0.1.0

Release candidate assessment for the first public release under the Society for
Immune Polymorphism (SIP) HLA Forge initiative.

- **Branch:** `release-hardening-v0.1.0`
- **Target version:** `0.1.0`
- **Assessment date:** 2026-11-02
- **Recommendation:** **READY to tag as a release candidate** (no blockers).

## Verification results

All ten release-candidate checks were executed. Checks 4–10 were run against a
**wheel installed into a clean, isolated virtual environment** (not the source
tree), to confirm the packaged artifact is self-contained.

| # | Check | Result | Evidence |
|---|-------|--------|----------|
| 1 | Full test suite | **PASS** | `pytest` → 178 passed |
| 2 | Ruff lint | **PASS** | `ruff check src tests examples` → All checks passed |
| 3 | Package build | **PASS** | `uv build` → `diplotypeforge-0.1.0-py3-none-any.whl` + sdist |
| 4 | Clean-environment install | **PASS** | fresh venv, wheel install → `import diplotypeforge` reports `0.1.0`; CLI present; schema loads (checksum `7609fef7…`) from the wheel |
| 5 | End-to-end strict HFX simulation | **PASS** | `validate`+`simulate` on the CAU example (strict, no lenient flag) |
| 6 | Multi-population simulation | **PASS** | Wahlund, 2 populations (`examples/two_population/wahlund.toml`) |
| 7 | Frequency-preserving transform | **PASS** | target 0.35 → reached exactly (`reason: exact`); `diplotypes_original.tsv` also written |
| 8 | Replicate batch | **PASS** | 5 replicates → `n_success: 5, n_fail: 0` |
| 9 | Reproducibility by **checksum** | **PASS** | two independent runs → all 5 output files byte-identical (sha256) |
| 10 | Validate all generated GL Strings | **PASS** | 35,000 GL strings parsed + lossless round-trip, 0 invalid |

### Scientific-verification evidence (Phase B)

- **Homozygosity feasibility** independently verified against a brute-force
  enumeration of every perfect matching for small inventories, plus ~650
  Hypothesis cases: `H_min`/`H_max`, `is_reachable` (both directions), the K=2
  parity stride, constructive completeness + termination, and exact count
  preservation all confirmed.
- **Statistical diagnostics** verified analytically for asymmetric matrices:
  `P(h,k)=Σ_ij M[i,j] f_i(h) f_j(k)` matches an independent reference;
  ordered/unordered normalize to 1; homozygotes are not doubled; `a_i` uses both
  matrix margins; `M = w wᵀ` reproduces pooled HWE exactly.
- **Per-locus experiment** demonstrates (and tests) that a reassignment changing
  complete-haplotype homozygosity by +0.136 changed per-locus homozygosity by
  different amounts at each locus (−0.039 / −0.005 / −0.006), with every
  haplotype and per-locus allele count preserved — confirming non-proportionality.

## Test counts

178 tests total:

| file | tests |
|------|------:|
| test_milestone2.py | 23 |
| test_glstring_corpus.py | 22 |
| test_feasibility_property.py | 20 |
| test_milestone3a.py | 16 |
| test_hfx_schema_regression.py | 14 |
| test_glstring.py | 13 |
| test_simulator.py | 12 |
| test_provenance_reproducibility.py | 10 |
| test_diagnostics_analytical.py | 10 |
| test_validation.py | 9 |
| test_hfx.py | 9 |
| test_sampling.py | 7 |
| test_milestone3b.py | 6 |
| test_per_locus_experiment.py | 3 |
| test_milestone3c.py | 2 |
| test_expectation.py | 2 |

## Supported features

- HFX v0.1.1 input (schema-validated by default; explicit lenient mode for
  historical files), inline and external CSV frequency data.
- Four population modes via a generalized K×K origin matrix: single-population,
  Wahlund, first-generation admixture, arbitrary K≥3 (asymmetric supported).
- Sampling-independent theoretical expectations (M-based + pooled-HWE) and
  observed-vs-expected diagnostics with TV distances.
- Frequency-preserving homozygosity transform with exact feasibility analysis
  and honest reporting (no silent target substitution).
- Reproducible replicate batches with a documented `SeedSequence` hierarchy and
  per-replicate + aggregate statistics.
- Restricted phased GL String serialization/parsing (round-trip).
- Python API + CLI (`validate` / `simulate` / `summarize` / `transform`).
- Reproducibility manifest (seed, software + numpy versions, input checksums,
  HFX schema provenance) with byte-identical output under a fixed environment.

## Known limitations

- Allele tokens are validated **syntactically**, not against an IPD database.
- One HFX file is treated as one source population (OPEN_QUESTIONS Q3).
- Copy number fixed at 2 per locus; gene absence / CNV is a future extension.
- The homozygosity transform reaches any *feasible* target but is **not** a
  uniform sampler over valid pairings; the transformed sample has no closed-form
  diplotype distribution.
- The current transform holds the full inventory in memory (`O(N)`); streaming /
  external-memory variants are documented as future work, not implemented.
- KIR copy-number modeling, recombination, and additional population models are
  out of scope.

## Outstanding issues

**Blockers:** none.

**Important but non-blocking:**
- Reproducibility is guaranteed only for a fixed NumPy version (documented in the
  manifest and README). Output is not guaranteed byte-identical across NumPy
  versions.
- CI matrix covers Python 3.12 and 3.13 only (the declared support window).

**Future enhancements:**
- Extract the GL String subset into a shared SIP `glstring` library
  (`docs/GLSTRING_INTEROP.md` identifies the five-name interface and call sites).
- Streaming / external-memory homozygosity transform (`docs/PERFORMANCE.md`).
- Multi-population HFX semantics; copy-number / gene-absence data model;
  additional transform strategies (e.g. uniform matching).
- Per-locus diagnostics in the standard `summary.json`.

## Dependencies and versions (verified environment)

| component | version | license | compatibility |
|-----------|---------|---------|---------------|
| Python | 3.12.11 (CI also 3.13) | PSF | — |
| numpy | 2.5.4 (req `>=1.26`) | BSD-3-Clause | Apache-2.0 compatible |
| jsonschema | 4.26.0 (req `>=4.21`) | MIT | Apache-2.0 compatible |
| pytest (dev) | 9.2.0 | MIT | compatible |
| ruff (dev) | 0.17.3 | MIT | compatible |
| hypothesis (dev) | 6.170.14 | MPL-2.0 | compatible (test-only, not distributed) |
| matplotlib (viz, optional) | — | PSF/BSD-style | compatible |

Project license: **Apache-2.0** (unchanged). All runtime and dev dependency
licenses are compatible with Apache-2.0 redistribution.

## HFX schema provenance

- Upstream: `societyforimmunepolymorphism/hfx`, PR #10, merge commit
  `a434a217ff7c8f13a73d04b7635e34e76b3ef310`.
- Document version: `0.1.1` (unchanged by the correction).
- Canonical schema SHA-256: `7609fef72ebf381d36d1bfc49c951d6138717df5fb2ead7ff40a7176c4031bd3`.
- Vendored verbatim; no runtime dependency on the upstream repository; recorded
  in every run `manifest.json` under `hfx_schema`.

## Recommendation

DiplotypeForge **v0.1.0 is ready to tag as a release candidate.** All verification
checks pass with cited evidence, there are no blockers, the wheel is
self-contained and installs cleanly, outputs are reproducible by checksum, and
every generated GL String validates and round-trips losslessly.

Do **not** publish to PyPI or create a GitHub release without explicit
authorization (per the project's release policy).
