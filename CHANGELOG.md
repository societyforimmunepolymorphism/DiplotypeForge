# Changelog

All notable changes to DiplotypeForge are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — first public release (SIP)

First public release candidate under the Society for Immune Polymorphism (SIP)
HLA Forge initiative. Implements Milestones 1–3 plus release hardening.

### Added
- Internal typed, phase-preserving, provenance-aware data model (`Allele`,
  `Locus`, `Haplotype`, `Diplotype`, `SimulatedIndividual`, distributions).
- Authoritative **HFX v0.1.1** input adapter (`societyforimmunepolymorphism/hfx`),
  schema-validated by default, with an explicit `lenient_hfx` mode for
  nonconforming historical files.
- Generalized **K×K population-origin matrix** sampling with four named modes:
  single-population, Wahlund, first-generation admixture, and arbitrary K≥3
  (asymmetric supported). Per-haplotype source-population provenance.
- Sampling-independent **theoretical expectations** (marginal ancestry,
  population-pair proportions, ordered/unordered diplotype distributions under
  M, pooled-HWE) and observed-vs-expected diagnostics with TV distances.
- Flagship **three-matrix experiment** (identical marginal ancestry, same
  marginal haplotype frequencies, different diplotype distributions).
- **Frequency-preserving homozygosity transform** with an exact feasibility
  predicate, a constructive reassignment algorithm, and full reporting; counts
  preserved exactly.
- **Reproducible replicate batches** with a documented `SeedSequence.spawn`
  hierarchy and per-replicate + aggregate statistics.
- Restricted, phased **GL String** serialization/parsing (round-trip).
- Python API + CLI (`validate` / `simulate` / `summarize` / `transform`).
- Validated **≥ 1,000,000-individual** benchmark.

### Release-hardening changes (this release)
- **HFX schema correction integrated** from upstream PR #10 (merge commit
  `a434a217ff7c8f13a73d04b7635e34e76b3ef310`); document version unchanged
  (0.1.1). Removed the in-memory `oneOf`→`anyOf` workaround — the corrected
  `frequencyDataLocation` `oneOf` is mutually exclusive, so documents are now
  validated against the vendored schema **verbatim with no mutation**. Recorded
  the upstream commit SHA and canonical schema SHA-256
  (`7609fef7…`) in docs and in every run `manifest.json` (`hfx_schema`).
- The bundled CAU example now carries `version: "0.1.1"` and validates strictly;
  example configs no longer require `lenient_hfx`.
- **Bug fix:** `expected_diplotype_distribution` keyed unordered diplotypes in
  haplotype-union (input-order-dependent) order, so the M-based TV distance
  could vary with input `frequencyData` ordering. Unordered keys are now
  canonicalized (sorted), matching the observed keying and the pooled-HWE path.
- Independent verification of the homozygosity-feasibility mathematics via
  brute-force matching enumeration + Hypothesis; analytical diagnostics tests
  for asymmetric matrices; a per-locus multilocus diagnostics experiment.
- Reproducibility now verified by **byte-level output checksums**.
- Packaging finalized: version `0.1.0`, PEP 621 metadata (license, SIP
  attribution, URLs, classifiers), HFX schema bundled in the wheel, verified by
  a clean-environment install + CLI simulation.
- Added `docs/API_STABILITY.md`, `docs/GLSTRING_INTEROP.md`, and
  `docs/DIPLOTYPEFORGE_AND_TISSUETYPEFORGE.md`; corrected `docs/PERFORMANCE.md`
  memory-complexity claims (current `O(N)` vs. theoretical `O(H)` pairing,
  identity/provenance/external-memory distinctions, `tracemalloc` caveat).
- GitHub Actions CI (ruff, pytest, build, wheel smoke, CLI simulation) on
  Python 3.12 and 3.13; the ≥1M benchmark is manual-only.

### Known limitations
- Allele tokens are validated **syntactically**, not against an IPD database.
- One HFX file is treated as one source population.
- Copy number is fixed at 2 per locus; gene absence / CNV is a future data-model
  extension.
- The homozygosity transform reaches any *feasible* target but is **not** a
  uniform sampler over valid pairings; the transformed sample has no closed-form
  diplotype distribution.
- KIR copy-number modeling, recombination, and additional population models are
  out of scope.

[0.1.0]: https://github.com/societyforimmunepolymorphism/DiplotypeForge/releases/tag/v0.1.0
