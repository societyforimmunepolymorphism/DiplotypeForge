# Configuration Reference

DiplotypeForge simulations are configured with a TOML file. The file has one
`[simulation]` table and one or more `[[populations]]` entries.

## `[simulation]`

| key | type | default | description |
|-----|------|---------|-------------|
| `n_individuals` | int ≥ 0 | 0 | number of individuals to sample per replicate |
| `seed` | int | 0 | RNG seed (PCG64) |
| `mode` | string | `single_population` | one of `single_population`, `wahlund`, `admixture`, `matrix` |
| `output_dir` | string | `output` | directory for written artifacts |
| `compress` | bool | false | gzip the TSV outputs |
| `replicates` | int ≥ 1 | 1 | independent replicate streams via `SeedSequence.spawn` |
| `strict` | bool | true | strict frequency-distribution validation |
| `lenient_hfx` | bool | false | accept nonconforming historical HFX (warn instead of reject) |
| `allow_renormalization` | bool | false | permit a frequency-sum renormalization beyond round-off |
| `renorm_tolerance` | float | 0.0 | max absolute `|sum−1|` to accept when renormalization is authorized |
| `aggregate_duplicates` | bool | false | explicitly merge duplicate haplotypes (warns + records provenance) |

### Mode parameters

| key | used by | type | description |
|-----|---------|------|-------------|
| `weights` | `wahlund` | list[float] | per-population weights (length = #populations); normalized to the diagonal of `M` |
| `admixture_a` | `admixture` | float ∈ [0,1] | admixed fraction `a`; builds `M=[[(1-a)/2, a/2],[a/2,(1-a)/2]]` |
| `matrix` | `matrix` | list[list[float]] | explicit K×K origin matrix (K = #populations) |
| `symmetric` | `matrix` | bool (default false) | if true, symmetrize the given matrix `(M+Mᵀ)/2` |

## `[[populations]]`

Each entry is one source population mapped from one HFX file.

| key | type | description |
|-----|------|-------------|
| `id` | string | stable internal population id (used in origin metadata + output) |
| `hfx` | string | path to the HFX metadata JSON (relative to the config file) |

Population count constraints by mode:

- `single_population`: exactly 1 population.
- `wahlund`: ≥ 2 populations (weights length must match if given).
- `admixture`: exactly 2 populations (the balanced two-population model).
- `matrix`: K populations, matching the K×K `matrix` dimensions.

## Examples

### Single population

```toml
[simulation]
n_individuals = 10000
seed = 20120815
mode = "single_population"
# lenient_hfx = true   # only for historical HFX files lacking `version`

[[populations]]
id = "CAU"
hfx = "cau.hfx-metadata.json"
```

### Wahlund (two populations)

```toml
[simulation]
n_individuals = 20000
seed = 20120815
mode = "wahlund"
weights = [0.5, 0.5]

[[populations]]
id = "POP1"
hfx = "pop1.hfx-metadata.json"

[[populations]]
id = "POP2"
hfx = "pop2.hfx-metadata.json"
```

### First-generation admixture

```toml
[simulation]
n_individuals = 20000
seed = 20120815
mode = "admixture"
admixture_a = 0.8     # => M = [[0.1, 0.4], [0.4, 0.1]]

[[populations]]
id = "POP1"
hfx = "pop1.hfx-metadata.json"

[[populations]]
id = "POP2"
hfx = "pop2.hfx-metadata.json"
```

### Generalized K×K matrix (K = 3)

```toml
[simulation]
n_individuals = 30000
seed = 20120815
mode = "matrix"
matrix = [
  [0.30, 0.10, 0.05],
  [0.10, 0.20, 0.05],
  [0.05, 0.05, 0.10],
]
# symmetric = true   # optional: symmetrize for unordered interpretation

[[populations]]
id = "POP1"
hfx = "pop1.hfx-metadata.json"

[[populations]]
id = "POP2"
hfx = "pop2.hfx-metadata.json"

[[populations]]
id = "POP3"
hfx = "pop3.hfx-metadata.json"
```

The matrix need not be symmetric: `M[i,j]` is the probability the **first**
haplotype comes from population `i` and the **second** from `j`, so an
asymmetric matrix gives different laws to the two sampled haplotypes while
per-haplotype provenance is retained.

## `[transformation]` (optional)

An optional post-sampling, frequency-preserving homozygosity transform. It
changes the number of homozygous diplotypes while preserving **every haplotype
count exactly** — a mechanism for Hardy–Weinberg disequilibrium distinct from
population structure. It never resamples; it rearranges the sampled copies.

| key | type | description |
|-----|------|-------------|
| `enabled` | bool | turn the transform on (default false) |
| `type` | string | `homozygosity` (only supported type) |
| `strategy` | string | `reassignment` (only supported strategy) |
| `target_fraction` | float ∈ [0,1] | target homozygous **fraction**; exactly one of fraction/delta |
| `target_delta` | float | signed change in homozygous fraction; exactly one of fraction/delta |

```toml
[simulation]
n_individuals = 5000
seed = 20120815
mode = "single_population"

[[populations]]
id = "CAU"
hfx = "cau.hfx-metadata.json"

[transformation]
enabled = true
type = "homozygosity"
strategy = "reassignment"
target_fraction = 0.30
```

The transform reports the initial / requested / achieved homozygosity, the exact
feasible `[min, max]` range, the number of operations, whether the target was
reached, and a reason for any rounding/clamping (the target is snapped to the
nearest feasible, reachable count rather than silently substituted). The
transformed sample has no closed-form diplotype distribution and is **not**
described as sampled from the origin matrix `M`. A standalone transform of an
existing dataset is available via the `diplotypeforge transform` CLI subcommand.
