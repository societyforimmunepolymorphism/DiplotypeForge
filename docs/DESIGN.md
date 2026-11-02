# DiplotypeForge — Architecture & Design

DiplotypeForge is an open-source Python package of the Society for Immune
Polymorphism (SIP) HLA Forge initiative. It simulates synthetic populations of
individuals with **known, phased diploid genotypes (diplotypes)** drawn from
population **haplotype-frequency distributions**.

Its companion, **TissueTypeForge**, consumes DiplotypeForge output (phased GL
Strings) and independently simulates laboratory typing observations. The phased
GL String is the shared contract between the two projects.

**Design philosophy: HLA-first, but not HLA-only.** The simulation engine makes
no assumption about the number, names, or ordering of loci and is applicable to
other polymorphic systems (KIR, LILR) where appropriate.

---

## 1. Conceptual layering

The package keeps three concepts strictly distinct, per the scientific brief:

1. **Population model** — the input truth: a set of named source populations,
   each with a haplotype-frequency distribution, plus a sampling model
   describing how individuals' two haplotypes are drawn.
2. **Sampled diplotypes** — the realized synthetic individuals, each carrying an
   *ordered* pair of haplotypes with per-haplotype population provenance.
3. **Derived summary statistics** — empirical frequencies and population-genetic
   diagnostics computed *from* the sampled individuals and compared against the
   model.

These map onto three layers of code:

```
         ┌──────────────────────────────────────────────────┐
  INPUT  │  adapters/   HFX (authoritative) → internal model  │
         └───────────────────────┬──────────────────────────┘
                                 ▼
         ┌──────────────────────────────────────────────────┐
  MODEL  │  model/      typed domain objects (frequency-table │
         │              abstraction, populations, config)      │
         └───────────────────────┬──────────────────────────┘
                                 ▼
         ┌──────────────────────────────────────────────────┐
 ENGINE  │  sampling/   seeded RNG, K×K origin model,         │
         │              haplotype draws → ordered diplotypes   │
         └───────────────────────┬──────────────────────────┘
                                 ▼
         ┌──────────────────────────────────────────────────┐
 OUTPUT  │  serialize/  phased GL Strings, TSV, manifest      │
         │  stats/      empirical freq & HWE diagnostics       │
         └──────────────────────────────────────────────────┘
```

The **engine never imports the HFX adapter**; it operates only on the internal
frequency-table abstraction. HFX is one input format among potentially several.

---

## 2. Input: HFX (authoritative)

HFX is the required native interchange format. We implement an adapter against
the **authoritative specification**:

- Repository: `societyforimmunepolymorphism/hfx`
- Schema version: **0.1.1**, JSON Schema draft 2020-12, `$id https://phycus.org/hfx.schema.json`

### 2.1 Document shape

An HFX document (`*.hfx-metadata.json`) is a JSON object:

```jsonc
{
  "version": "0.1.1",                 // required (const in v0.1.1)
  "metadata": { ... },                // required
  "frequencyData": [ ... ]            // optional; present only for inline mode
}
```

`metadata` **required** keys: `outputResolution`, `hfeMethod`,
`cohortDescription`, `nomenclatureUsed`, `frequencyLocation`.
`metadata` optional keys: `creationDateTime`, `producer`, `frequencyFileHeader`,
`license`, `hfxID`.

Key points the adapter relies on:

- **`outputResolution`** is an *ordered* array of `{locus, resolution}`. **This
  array defines the canonical locus order** for every haplotype in the file.
  HLA resolution enum: `G, P, g, allele-family, protein-sequence,
  coding-sequence, genomic-sequence, Serology`. KIR:
  `protein-sequence, coding-sequence, genomic-sequence`.
- **`cohortDescription.population[].name`** is the population identifier used as
  the source-population label internally.
- **`frequencyLocation`** is either a URI (`file://name.csv`, `http://...`) or
  the literal `"inline"`.
- **Inline mode:** `frequencyData` = `[{haplotype: <GL String>, frequency:
  <number>}]`. Haplotype pattern: `^[gA-Z-~*:0-9\.;>\[\]]*$`.
- **CSV mode:** `frequencyLocation` points to a CSV. Default header is
  `haplotype,frequency`. Optional `frequencyFileHeader` maps the *expected*
  field names (`haplotype`, `frequency`, `count`) onto the *actual* column
  headers in the CSV (e.g. `{"haplotype":"Haplo","frequency":"Freq",
  "count":"Count"}`). A `count` column may be present and is preserved.

### 2.2 Conformance validation (schema-by-default)

The adapter **validates each HFX document against the authoritative JSON Schema
by default** (`hfx.schema.json` v0.1.1, vendored into the package), using
`jsonschema`. Documents that do not conform are **rejected** by default.

An explicit **lenient mode** (`strict=False` / `--lenient`) is provided for
nonconforming *historical* files: in lenient mode a schema violation is demoted
to a warning, unknown metadata keys are preserved as opaque provenance, and the
known typo `outputRResolution` (seen in official example files) is tolerated
with a warning. Lenient mode is opt-in precisely because real-world legacy files
(e.g. ones carrying an extra `checkSum` key, which the `additionalProperties:
false` schema forbids) must remain loadable without weakening the default
guarantee.

**Reuse investigation — SIP `hfx-tools` (item 3).** We evaluated
`societyforimmunepolymorphism/hfx-tools` (v0.2.0, Apache-2.0, zero runtime
deps). Its `ValidationFramework`/`ValidationResult` performs **hand-rolled field
checks** (required-field presence, `version == "0.1.1"`, frequency sanity,
file-reference existence) — it does **not** validate against the authoritative
`hfx.schema.json`. Since the review explicitly requires conformance *against the
authoritative schema*, `hfx-tools` cannot satisfy it; we therefore validate with
`jsonschema` against the vendored schema. `hfx-tools` remains the documented
**interop path** (its `read_hfx_bundle`/`HfxBundle` reader, `.hfx` archive and
Parquet support, and `ValidationResult` issue model informed our design). We do
not hard-depend on it in Milestone 1, since doing so would duplicate reading
without adding schema conformance.

### 2.3 Multi-population configuration

A single HFX file describes one cohort (which may list multiple `population`
entries, but carries one `frequencyData`/CSV). DiplotypeForge therefore treats
**one HFX file as one source population's frequency distribution** by default,
and composes **multiple HFX files** into a multi-population simulation via the
simulation configuration (`[[populations]]` entries each referencing a file and
assigning a stable internal id/label). This keeps provenance per population and
avoids guessing how to split a multi-population cohort — see Open Question Q3.

---

## 3. Internal data model

All domain types are explicit and typed (dataclasses / typed objects). None of
them is a GL String — GL String is a *serialization*, produced on demand.

| Type | Responsibility | Key fields |
|------|----------------|------------|
| `Allele` | One allele call at one locus | `locus`, `name` (field string, e.g. `01:01g`), optional raw token |
| `Locus` | A named locus + its resolution | `name`, `resolution` |
| `Haplotype` | An ordered tuple of alleles across the declared loci | `alleles: tuple[Allele, ...]`, references `loci` order |
| `HaplotypeFrequencyDistribution` | The frequency-table abstraction: ordered loci + haplotype→frequency map | `loci: tuple[Locus,...]`, `entries: list[(Haplotype, freq, count?)]`, `metadata` |
| `SourcePopulation` | One population: id/label + its distribution + provenance | `population_id`, `label`, `distribution`, `provenance` |
| `Diplotype` | Two haplotypes with known cross-locus correspondence (phase) | `haplotype_1`, `haplotype_2` (ordered pair) |
| `SimulatedIndividual` | One sampled person | `sample_id`, `diplotype`, `origin_1`, `origin_2` (per-haplotype population id) |
| `SimulationConfiguration` | Everything needed to reproduce a run | populations, mode, matrix, `n_individuals`, `seed`, `replicate`, output opts, software version |
| `SimulationResult` | Run output + provenance | individuals (or stream), manifest, empirical stats, config echo |

### 3.1 Phase and provenance invariants

- A `Diplotype` stores an **ordered pair** `(haplotype_1, haplotype_2)`.
  Downstream comparisons may treat diplotypes as unordered, but the ordered pair
  is always retained internally.
- Each haplotype's **source population is recorded separately** on the
  `SimulatedIndividual` (`origin_1`, `origin_2`), so Wahlund / admixture / K×K
  origin can be reconstructed and validated.

### 3.2 Extensibility (future-proofing)

Milestone 1 fixes copy number at **2 per locus**. The interfaces are designed so
future structure can be added *without disrupting the diploid base*, but we are
deliberate about what that future work actually requires:

- The sampling engine addresses haplotypes by **index into a frozen, sorted
  array**, independent of locus semantics, so novel locus structures affect only
  the model/serialization layers, not the engine.
- `Haplotype` holds alleles keyed by locus, leaving room for a future
  presence/absence sentinel at a locus.

**Copy-number variation is NOT merely a scalar annotation (reviewer item 6).**
We explicitly reject the assumption that CNV (KIR gene presence/absence,
DRB3/4/5, C4 copy number, gene duplications) can be represented by attaching a
copy-number count to a fixed two-haplotype, fixed-locus pair. CNV changes the
*cardinality and structure* of a genotype — the number of gene copies varies per
individual and per locus — so it requires a genotype representation that is not a
fixed-width ordered pair (variable-copy loci, presence/absence alleles). This is
documented as a **future data-model extension** and is **not implemented in
Milestone 1**, nor approximated via an annotation. See Open Question Q10.

---

## 4. Simulation engine — the generalized model

The brief requires four modes. We implement a **single generalized
population-origin model** of which the first three are special cases. (Mode D is
the base; A–C are Milestone 2 specializations expressed as matrices.)

### 4.1 The K×K origin matrix

For `K` source populations, a `K×K` matrix `M` gives the probability that an
individual's **first** haplotype originates in population `i` and its **second**
in population `j`:

```
P(origin_1 = i, origin_2 = j) = M[i, j]
```

Requirements:

- Entries nonnegative; `sum(M) == 1`.
- `M` encodes **ordered** population pairs. For unordered biological
  interpretation we offer a `symmetric=True` convenience that builds `M` from an
  upper-triangular specification by symmetrizing.
- The simulator records each haplotype's population of origin.

Sampling per individual:

1. Draw ordered origin pair `(i, j)` from the flattened `M` (categorical draw).
2. Draw `haplotype_1` from population `i`'s distribution; `haplotype_2` from
   population `j`'s distribution — **independently** (random mating within/between
   as specified by M).
3. Emit `SimulatedIndividual` with origins `(i, j)`.

### 4.2 Modes as special cases

- **A. Single-population random mating** (Milestone 1): `K = 1`, `M = [[1]]`.
  Each individual draws two haplotypes independently from `P`. This yields HWE:
  `P(i,i) = f_i^2`, `P(i≠j unordered) = 2 f_i f_j`.
- **B. Wahlund / stratification** (Milestone 2): both haplotypes come from the
  *same* source; `M` is **diagonal** with `M[k,k] = w_k` (population weights),
  off-diagonal zero. Produces population-level HWE *disequilibrium* from
  structure even though mating is random within each subpopulation.
- **C. First-generation admixture (two populations)** (Milestone 2): with
  admixed fraction `a`:
  `M = [[(1-a)/2, a/2], [a/2, (1-a)/2]]`.
  Fraction `a` get one haplotype from each source; `(1-a)/2` get both from pop 1;
  `(1-a)/2` get both from pop 2. At `a = 0.5` this reduces to independent source
  assignment with marginal ancestry 0.5 each. **The admixed fraction `a` is not
  the marginal ancestry proportion** — the engine reports both separately.
- **D. Generalized K×K** (Milestone 2): arbitrary `M`, `K ≥ 2` (tested with 3).

Milestone 1 ships mode A plus the generalized matrix scaffolding so A is literally
`GeneralizedOriginModel(M=[[1]])`.

### 4.3 Expected diplotype distribution under M (reviewer item 5)

The expected **ordered** diplotype distribution under origin matrix `M`
(K populations, per-population haplotype frequencies `f_i(h)`) is:

```
P(h, k) = Σ_i Σ_j  M[i, j] · f_i(h) · f_j(k)
```

For **unordered** diplotypes (genotypes): for `h ≠ k` sum both orders,
`P_unordered({h,k}) = P(h,k) + P(k,h)`; for `h = k`, `P_unordered = P(h,h)`.
This is the authoritative expectation used for expected frequencies and
diagnostics.

**This M-based expectation is distinct from a pooled-HWE expectation.** The
pooled-HWE null uses the marginal (mixture) haplotype frequency
`f̄(h) = Σ_i a_i f_i(h)`, where `a_i` are the marginal ancestry proportions
(`a_i = ½(Σ_j M[i,j] + Σ_j M[j,i])`), and assumes random union of gametes:
`P_HWE(h,k) = f̄(h)·f̄(k)`. When `M` encodes structure or admixture the two
expectations differ, and that difference **is** the Wahlund/structure signal.
Diagnostics report and label both. In the Milestone-1 single-population case
(`K=1`, `M=[[1]]`) they coincide: `P(h,k) = f(h)·f(k)`, i.e. Hardy–Weinberg at
the **complete multilocus-haplotype** level (`P(h,h)=f(h)²`,
`P_unordered({h,k})=2 f(h) f(k)` for `h≠k`).

### 4.4 Frequency-preserving homozygosity transform (Milestone 3)

A separate, **post-sampling** transform modifies the number of homozygous
diplotypes while preserving the empirical haplotype-frequency distribution
*exactly*. It is a complementary mechanism to the structure/admixture models:
both create Hardy–Weinberg disequilibrium, but by distinct routes (see §5.4).

Given `N` individuals = `2N` haplotype copies, the transform rearranges copies
into `N` unordered pairs, preserving all of:

- `N` individuals and `2N` haplotype copies;
- the count `c_h` of every distinct haplotype, and the distinct-haplotype set;
- the count of every allele at every locus (implied by preserving `c_h`).

The original sample is retained unchanged and remains recoverable.

**Elementary operation.** The transform works on pairs of pairs. To *increase*
homozygosity it finds two heterozygous individuals that between them hold two
copies each of two haplotypes and regroups them:

```
(A, B), (A, B)  ->  (A, A), (B, B)         # +2 homozygous, counts preserved
```

More generally, two heterozygous individuals sharing a common haplotype `A`,
`(A, B)` and `(A, C)`, can swap their non-`A` partners to make `(A, A)` and
`(B, C)` when a second `A` is available — the engine locates these configurations
via a haplotype→holders index. To *decrease* homozygosity the reverse is applied:
`(A, A), (B, B) -> (A, B), (A, B)`. Each elementary operation changes the
homozygous count by a fixed, known amount and preserves every `c_h`.

### 5. (continued) Feasibility — exact and algorithm-independent

A target of `H` homozygous pairs is **achievable** iff we can choose per-key
homozygote counts `hom_k ∈ [0, ⌊c_k/2⌋]` with `Σ hom_k = H` whose leftovers
`c_k − 2·hom_k` are heterozygous-pairable (`max_k leftover_k ≤ N − H`). That is
equivalent to a per-key lower bound `lo_k(H) = max(0, ⌈(c_k − (N−H))/2⌉)`, so:

```
H is achievable  iff  Σ_k lo_k(H) ≤ H ≤ Σ_k ⌊c_k/2⌋
```

- `H_min` and `H_max` are the smallest/largest achievable `H` (scanned from the
  classic bounds `max(0, c_max − N)` and `Σ⌊c_k/2⌋`). The naive `Σ⌊c_k/2⌋` is an
  **upper bound** but is not always itself achievable, because the leftover
  (odd-count) copies may fail to be heterozygous-pairable near the maximum — the
  exact `H_max` can be one or two below `Σ⌊c_k/2⌋`.
- The achievable set is a single contiguous **stride**: step 1 (every integer)
  in the general case, but step 2 (a **parity** constraint) when only two
  distinct haplotypes carry all copies — then only every other `H` is reachable.
- Homozygosity proportion is `H/N`, on a `1/N` grid.

These are computed directly from the haplotype-count distribution, independently
of the algorithm. `is_reachable(H)` exposes the exact predicate; a requested
count that is clamped to `[H_min, H_max]` is then **snapped to the nearest
reachable value** with reason `snapped_to_nearest_reachable_count` (we never
silently substitute — the reason is always reported).

### 5. (continued) Strategy & interface

The scientific objective (a target homozygous count) is separated from the
optimization via a `TransformationStrategy` interface. Milestone 3 ships **one
robust strategy, `reassignment`**, with clear guarantees:

- It reaches any feasible requested homozygous count with a **bounded** number
  of elementary operations (≤ `|H* − H_0|` net steps, each O(1) amortized via
  the holder index) — **no unbounded rejection loops**.
- It is **not** a uniform sampler over valid matchings and does not claim to be
  (OPEN_QUESTIONS Q15). An optional randomization step (shuffle candidate order
  with the transformation RNG) avoids deterministic artifacts.

**Target & rounding.** A target is a homozygosity *fraction* `p` or a *delta* `Δ`
from the initial fraction. It maps to an integer count `H* = round(p·N)` (round
half to even), clamped to `[H_min, H_max]`. We never silently substitute a
target: the result reports requested/rounded/clamped/achieved and a `reason`
when rounding or clamping altered the effective target (OPEN_QUESTIONS Q13).

### 5. (continued) Provenance under reassignment

Each haplotype carries its **source-population origin**, and that origin travels
with the haplotype when it moves between individuals. An individual's
population-pair label may therefore change after transformation. The result
retains the original origin matrix, the original sampled pair proportions, and
the final (post-transform) pair proportions.

### 5.4 Two distinct routes to HW disequilibrium

- **Structural (Wahlund / admixture):** disequilibrium arises because the
  *population model* is not panmictic — the joint origin matrix `M` concentrates
  or disperses same-source pairings. The M-based expectation (§4.3) predicts it
  analytically.
- **Imposed (reassignment transform):** disequilibrium is injected *after*
  sampling by rearranging a fixed haplotype inventory. It defines a **new
  sampling mechanism** with **no** claimed closed-form diplotype distribution;
  the transformed sample is explicitly *not* M-distributed. Diagnostics report
  the original-M expectation, pooled-HWE, and observed before/after separately
  (§7), and never attribute the imposed deviation to structure.

---

## 5. Reproducibility (reviewer item 7)

- **RNG specification:** `numpy.random.Generator` backed by `PCG64`, seeded
  explicitly. The manifest records the bit generator (`PCG64`), the integer
  seed, and the numpy version used.
- **Order independence:** before sampling, haplotypes are placed in a
  **sorted, frozen array** (sorted by fully-qualified GL String key) so results
  never depend on dict/JSON iteration order. Origin pairs are drawn from a fixed
  flattened index.
- **Replicates:** `SeedSequence.spawn()` derives independent, reproducible child
  streams per replicate.
- **Manifest provenance (recorded every run):**
  - **Input checksums** — SHA-256 of every input file (each HFX metadata JSON
    and any referenced CSV).
  - **Canonicalized configuration** — the fully resolved config serialized with
    sorted keys and normalized values (the same canonical form feeds the
    `run_hash`), so a run is reproducible from the manifest alone.
  - **RNG specification** — generator (`PCG64`), seed, per-replicate seed
    derivation method (`SeedSequence.spawn`).
  - **Software versions** — DiplotypeForge version, Python version, and the
    numpy version.
  - Any **renormalization** applied (original sum, tolerance, scale factor) and
    any **duplicate aggregation** provenance.
- **Reproducibility guarantee — stated precisely:**
  - Identical inputs + seed + **DiplotypeForge version** + canonical
    configuration + **the same numpy version** ⇒ **byte-identical** output.
  - We rely on numpy's documented stability of the `PCG64` stream and its
    `Generator` sampling routines. We therefore **do not guarantee identical
    output across different numpy versions**: numpy only guarantees bit-stream
    stability within a version series and may change `Generator` algorithms
    across major versions. The numpy version is recorded in the manifest so any
    cross-version difference is detectable and attributable. Reproducing results
    exactly requires pinning the recorded numpy version.
- **Streaming:** large simulations write individuals incrementally (generator →
  TSV writer) so the full output need not be held in memory.

---

## 6. Output contract (shared with TissueTypeForge)

Produced artifacts (initial formats: TSV and gzip-compressed TSV):

1. **`diplotypes.tsv[.gz]`** — true phased diplotypes as GL Strings, one per
   individual, with stable `sample_id`.
2. **Individual metadata** — `sample_id`, `origin_1`, `origin_2` (per-haplotype
   source population), replicate id.
3. **`haplotype_frequencies.tsv`** — empirical haplotype frequencies.
4. **`diplotype_frequencies.tsv`** — empirical diplotype frequencies.
5. **Population-level summary statistics** (see §7).
6. **`manifest.json`** — machine-readable: input checksums, canonicalized
   configuration, RNG specification, software/numpy versions, renormalization
   and aggregation provenance, and the generated file list (see §5).

### 6.1 Phased GL String contract

A diplotype is serialized as **two `~`-joined haplotypes separated by `+`**,
with **fully qualified allele identifiers** (genetic-system prefix included) and
loci in the declared `outputResolution` order:

```
<hapA_locus1>~<hapA_locus2>~...~<hapA_locusL>+<hapB_locus1>~...~<hapB_locusL>
```

Example (HLA 5-locus, g-group, fully qualified):

```
HLA-A*01:01g~HLA-C*07:01g~HLA-B*08:01g~HLA-DRB1*15:01~HLA-DQB1*06:02+HLA-A*02:01g~HLA-C*03:04g~HLA-B*40:01g~HLA-DRB1*04:04~HLA-DQB1*03:02g
```

This follows the GL String grammar delimiter precedence (Milius 2013; Mack 2023):
`/` (allele ambiguity) < `~` (in-phase haplotype) < `+` (chromosome/gene-copy) <
`|` (genotype ambiguity) < `^` (locus). A fully phased diplotype uses only `~`
and `+`. The string round-trips back to an equivalent (ordered-pair-preserving)
`Diplotype`.

**Milestone-1 scope (reviewer item 4).** The local serializer/parser is
**restricted to the unambiguous phased diplotypes** DiplotypeForge produces: it
handles `~` and `+` only and **rejects** allele ambiguity (`/`), genotype
ambiguity (`|`), and locus delimiters (`^`). We do **not** implement the full GL
String grammar in Milestone 1. External libraries were evaluated —
`nmdp-bioinformatics/pyglstring` is an archived (read-only), GPL/LGPL-licensed
*sanity checker* rather than a typed parser, and nomenclature tools
(`mhcgnomes`, `py-hla-match`) are not GL String diplotype serializers — so none
fit a clean Apache-2.0 dependency producing the typed objects we need. The parser
is isolated behind a module seam so a full-grammar library can replace it later.

---

## 7. Population-genetic validation / diagnostics

Each run reports statistics **explicitly labeled by level of analysis** (per-locus
vs complete multilocus haplotype) and by **null model**:

- Input vs observed haplotype frequencies + sampling error (expected ±
  binomial/multinomial SE).
- Observed vs **expected** homozygosity and heterozygosity, where "expected" is
  computed from the **M-based ordered/unordered diplotype distribution** of §4.3,
  `P(h,k)=Σ_ij M[i,j] f_i(h) f_j(k)`.
- Diplotype-frequency distribution (unordered genotypes; see Q6).
- Population-of-origin proportions and population-pair proportions.
- Deviation from the **M-based expectation** versus deviation from the
  **pooled-HWE expectation** `f̄(h)f̄(k)` — reported as two distinct,
  separately-labeled quantities so structure/admixture signal is not confused
  with sampling noise. In the single-population case the two expectations
  coincide (multilocus-haplotype HWE).

We **never** emit a single unlabeled "HWE statistic"; every statistic names its
level and its null model (M-based vs pooled-HWE).

---

## 8. Python interface

### 8.1 API

```python
from diplotypeforge import Simulator, SimulationConfig

config = SimulationConfig.from_file("simulation.toml")
simulator = Simulator(config)
result = simulator.run()
result.write("output/")
```

### 8.2 CLI

```bash
diplotypeforge validate   simulation.toml     # validate config + HFX inputs
diplotypeforge simulate   simulation.toml     # run; write outputs + manifest
diplotypeforge summarize  output/             # recompute/print diagnostics
```

### 8.3 Configuration (TOML)

```toml
[simulation]
n_individuals = 10000
seed = 20120815
mode = "single_population"     # single_population | wahlund | admixture | matrix
output_dir = "output/"
compress = false

# Validation is strict by default. Renormalization beyond floating-point
# round-off must be explicitly authorized and is recorded in the manifest.
# strict = true                 # default; HFX schema + distribution strict
# lenient_hfx = false           # set true only for nonconforming historical files
# allow_renormalization = false # set true to permit a larger renorm
# renorm_tolerance = 0.0        # max |sum-1| to accept when authorized
# aggregate_duplicates = false  # explicit preprocessing; warns + records provenance

[[populations]]
id = "CAU"
hfx = "examples/inline-frequency/cau.hfx-metadata.json"

# For Milestone 2 modes:
# [populations.weights] or [matrix] blocks, documented in CONFIG.md
```

### 8.4 Dependencies (modest footprint)

- `numpy` — RNG + vectorized sampling.
- `jsonschema` — **validate HFX against the bundled authoritative schema by
  default** (the SIP `hfx-tools` package does not provide schema-based
  conformance; see §2.2).
- Standard library for TOML (`tomllib`), CSV, JSON, gzip, `hashlib` (checksums).
- Dev: `pytest`, `ruff`. Type annotations on all public interfaces.

Target **Python 3.12+**, managed by **uv**.

---

## 9. Package structure

```
DiplotypeForge/
├── pyproject.toml
├── README.md
├── docs/
│   ├── DESIGN.md                 (this file)
│   ├── DEVELOPMENT_PLAN.md
│   ├── CONFIG.md                 configuration reference (all modes)
│   ├── EXPERIMENT_THREE_MATRIX.md flagship experiment write-up
│   ├── MILESTONE2_SUMMARY.md
│   └── OPEN_QUESTIONS.md
├── src/diplotypeforge/
│   ├── __init__.py               public API exports
│   ├── model/                    Allele, Locus, Haplotype, distribution, etc.
│   ├── adapters/
│   │   ├── hfx.py                authoritative HFX adapter
│   │   └── schema/hfx.schema.json (vendored v0.1.1)
│   ├── validation.py             frequency-distribution validation
│   ├── glstring.py               phased GL String serialize/parse
│   ├── sampling/                 RNG + origin model + samplers
│   ├── expectation.py            theoretical expectations under M (no sampling)
│   ├── experiments.py            three-matrix experiment helper
│   ├── stats.py                  diagnostics (observed vs M-based vs pooled-HWE)
│   ├── output.py                 TSV / manifest writers
│   ├── config.py                 SimulationConfig (TOML)
│   ├── simulator.py              Simulator / SimulationResult
│   └── cli.py                    validate / simulate / summarize
├── examples/                     1-, 2-, 3-population examples
└── tests/                        pytest acceptance + unit tests
```

---

## 10. Milestones

- **M1 (Foundation, this cycle):** structure, data models, frequency-distribution
  validation, HFX adapter, single-population sampling, phased GL String
  serialization, API+CLI, unit tests, a working example.
- **M2 (Population structure):** Wahlund, first-generation admixture, generalized
  K×K matrix, origin metadata, statistical validation.
- **M3 (Advanced):** frequency-preserving homozygosity transform, replicates &
  batch, large-population performance, more diagnostics/examples.

Milestone gate: do not advance past M1 until M1 tests pass and the example runs.

See `docs/OPEN_QUESTIONS.md` for scientific ambiguities requiring decisions.
