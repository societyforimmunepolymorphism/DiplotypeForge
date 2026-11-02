# Public API stability — DiplotypeForge v0.1

This document defines what constitutes the **supported public API** of
DiplotypeForge v0.1, versus internal implementation detail that may change
without notice. v0.1 is an initial release: the API is usable and tested, but
only the surfaces listed under *Supported* carry a stability intent. Everything
else is internal.

DiplotypeForge follows semantic versioning. During the `0.x` series, minor
versions **may** introduce breaking changes to the supported API, but such
changes will be called out in `CHANGELOG.md`. From `1.0` onward, breaking
changes to the supported API will require a major-version bump.

## Supported public API (v0.1)

Everything importable directly from the top-level `diplotypeforge` package
(i.e. listed in `diplotypeforge.__all__`) is **supported**:

### Configuration loading
- `SimulationConfig` — `from_file(path)`, `from_dict(data, base_dir=...)`, and
  its public fields (`n_individuals`, `seed`, `mode`, `output_dir`, `compress`,
  `replicates`, `lenient_hfx`, `allow_renormalization`, `aggregate_duplicates`,
  `populations`, `weights`, `matrix`, `admixture_a`, `transformation`, …).
- `PopulationSpec` — `(id, hfx)`.

### Simulation execution
- `Simulator` — `Simulator(config)`, `.validate()`, `.run() -> SimulationResult`,
  `.load_populations()`.
- `SimulationResult` — `.individuals`, `.statistics`, `.manifest`,
  `.transform_results`, `.original_individuals`, `.replicate_summary`,
  `.write(output_dir=None) -> dict[str, str]`.

### Population-origin matrices
- `OriginModel` — constructors `single_population()`, `wahlund(weights)`,
  `admixture(a)`, `from_matrix(rows, symmetric=False)`; `.matrix` attribute.
- `OriginModelError`.
- Matrix-derived expectations (sampling-independent):
  `marginal_ancestry_proportions`, `ordered_population_pair_proportions`,
  `unordered_population_pair_proportions`, `marginal_haplotype_frequencies`,
  `expected_diplotype_distribution`, `expected_diplotype_distribution_pooled_hwe`,
  and the `DiplotypeExpectation` result type.

### Diplotype representation
- `Allele(locus, name)`, `Haplotype(alleles)`, `Diplotype(haplotype_1, haplotype_2)`,
  `Locus(name, resolution)`, `SimulatedIndividual(...)`.
- Stable methods: `Haplotype.key()`, `Diplotype.is_homozygous`,
  `Diplotype.ordered_key()`, `Diplotype.unordered_key()`.

### GL String serialization (restricted, phased subset)
- `haplotype_to_glstring`, `diplotype_to_glstring` (serialization);
- `parse_haplotype`, `parse_diplotype` (parsing of the accepted subset);
- `GLStringError` (raised on out-of-subset input).
  See `docs/GLSTRING_INTEROP.md` for the exact accepted/produced subset.

### Results and summary statistics
- `SummaryStatistics` (via `SimulationResult.statistics`) and its fields
  (observed vs. expected homozygosity/heterozygosity, TV distances, origin
  proportions, per-haplotype comparisons).
- The on-disk output contract: `diplotypes.tsv[.gz]`,
  `haplotype_frequencies.tsv`, `diplotype_frequencies.tsv`, `summary.json`,
  `manifest.json`, `replicates.json` (and `diplotypes_original.tsv` when a
  transform is applied). Column names and the manifest key set are part of the
  contract.

### Transformations
Supported via the `diplotypeforge.transform` subpackage:
- `transform_homozygosity(individuals, target, rng=None, strategy=None)`,
- `TransformationTarget(fraction=... | delta=...)`, `TransformationResult`,
- `analyze_feasibility` / `analyze_counts` / `HomozygosityFeasibility`,
- `transform_diagnostics` / `TransformDiagnostics`,
- `per_locus_diagnostics` / `compare_per_locus` / `PerLocusDiagnostics` /
  `PerLocusComparison`.

### Command-line interface
The `diplotypeforge` CLI and its subcommands are supported:
`validate`, `simulate`, `summarize`, `transform`. Flag names
(`--output-dir`, `--input`, `--config`) are part of the contract.

## Internal implementation detail (NOT part of the public API)

These may change at any time; do not import or depend on them:

- Any name prefixed with `_` (e.g. `Simulator._replicate_seed_sequences`,
  `_normalized_schema` — now removed, `_construct_pairs`, `_tv_distance`).
- Submodule internals not re-exported from the top-level package, e.g.
  `diplotypeforge.sampling.sampler.Sampler` / `make_rng`,
  `diplotypeforge.stats` internals, `diplotypeforge.output` writer functions,
  `diplotypeforge.replicates` internals, `diplotypeforge.experiments`,
  `diplotypeforge.transform.reassignment.ReassignmentStrategy` internals,
  `diplotypeforge.transform.io`.
- The HFX adapter beyond `read_hfx`, `HFXOptions`, `HFXError`, `load_schema`,
  `schema_checksum`, `schema_provenance` — in particular the private helpers.
- The exact text of warning/error messages and the ordering of JSON object keys
  (other than where byte-stability is explicitly guaranteed for a fixed
  software + numpy version).

## Likely future breaking changes (flagged for users)

The following are anticipated and will be handled through the changelog /
versioning policy:

1. **Copy-number / gene absence.** v0.1 fixes copy number at 2 per locus. A
   future data-model extension for gene absence / CNV will change `Haplotype`
   and the GL String subset (absence token), and may change TSV columns.
2. **Shared `glstring` library.** The restricted parser/serializer is expected
   to be replaced by an external SIP `glstring` package. The function names
   (`parse_*` / `*_to_glstring`) are intended to remain as thin shims, but the
   `GLStringError` type and the module location may change. See
   `docs/GLSTRING_INTEROP.md`.
3. **Multi-population HFX semantics.** One HFX file is currently treated as one
   source population (OPEN_QUESTIONS Q3). A future multi-population HFX mapping
   could change `read_hfx`'s return shape.
4. **Transform strategies.** Additional strategies (e.g. uniform matching) may
   be added; the `TransformationStrategy` protocol could gain members.
5. **Statistics expansion.** `SummaryStatistics` may gain fields (additive,
   non-breaking) and possibly per-locus diagnostics in the standard summary.

## Deprecation policy (0.x)

When a supported API element is slated for removal, it will first be documented
as deprecated in the changelog for at least one minor release before removal,
where practical.
