# GL String interoperability note

For the proposed shared SIP `glstring` Python library and for TissueTypeForge
integration.

DiplotypeForge deliberately implements only a **restricted, phased subset** of
the Genotype List (GL) String grammar (Milius et al. 2013; Mack et al. 2023).
This note defines exactly what the current parser **accepts** and what the
serializer **produces**, so a future shared library can be dropped in without
changing the output contract, and so TissueTypeForge can rely on a precise input
shape.

DiplotypeForge does **not** implement a full GL String parser and does not plan
to in v0.1. The subset below is the stable contract.

## Grammar delimiters

GL String delimiter precedence (low → high binding):

| delim | meaning | DiplotypeForge |
|-------|---------|----------------|
| `/`   | allele ambiguity          | **rejected** (out of subset) |
| `~`   | in-phase join (haplotype) | **supported** |
| `+`   | gene-copy / chromosome separator | **supported** |
| `\|`  | genotype ambiguity        | **rejected** (out of subset) |
| `^`   | locus separator           | **rejected** (out of subset) |

## Produced subset (serializer output)

A DiplotypeForge diplotype is always serialized as exactly **two** `~`-joined
haplotypes separated by a single `+`, with loci in the distribution's declared
order:

```
<a11>~<a12>~...~<a1L>+<a21>~<a22>~...~<a2L>
```

Properties guaranteed on output:

- exactly one `+` (two chromosomes, fixed copy number 2);
- `L` `~`-joined allele tokens per haplotype, **positionally aligned by locus**
  in declared order across both haplotypes and across all individuals;
- **phase is preserved**: `haplotype_1` and `haplotype_2` keep their sampled
  order (the serializer never canonicalizes/sorts the two haplotypes);
- allele tokens are emitted **verbatim** as stored (lossless round-trip),
  normally **fully qualified** (system-prefixed, e.g. `HLA-A*01:01g`);
- no ambiguity (`/`, `|`) and no locus (`^`) delimiters ever appear.

## Accepted subset (parser input)

`parse_haplotype` / `parse_diplotype` accept exactly the produced subset, plus
these tolerances:

- surrounding whitespace is stripped;
- allele tokens need not be fully qualified for parsing (qualification is an HFX
  adapter concern), but TissueTypeForge interop expects fully-qualified tokens;
- when a `loci` tuple is supplied, haplotype tokens may be given in a different
  order and are **reordered into declared locus order** by locus name; a
  mismatch in the *set* of loci (missing / extra / duplicate locus) is an error;
- without `loci`, loci are inferred from token order.

Rejected (raise `GLStringError`):

- any of `/`, `|`, `^`;
- a haplotype containing `+`, or a diplotype not having exactly two `+`-separated
  haplotypes;
- empty strings or empty allele tokens;
- (with `loci`) wrong locus count or a locus-set mismatch.

## TissueTypeForge compatibility

TissueTypeForge consumes **phased diplotype GL Strings** as its known-truth
input. The produced subset above is that shared contract: two `~`-joined,
positionally locus-aligned, fully-qualified haplotypes separated by `+`, with
phase preserved. The `diplotypes.tsv` `diplotype_glstring` column is the
canonical hand-off. A round-trip corpus (`tests/test_glstring_corpus.py`)
exercises this subset including the exact multilocus example shared between the
two projects.

## Interfaces a shared `glstring` library could replace

The entire GL String surface is isolated in `src/diplotypeforge/glstring.py`.
A shared library could replace it by providing these four functions and one
exception, with identical signatures and semantics:

```python
haplotype_to_glstring(haplotype) -> str
diplotype_to_glstring(diplotype) -> str
parse_haplotype(gl: str, loci=None) -> Haplotype
parse_diplotype(gl: str, loci=None) -> Diplotype
GLStringError  # ValueError subclass
```

Internal call sites are limited and already funnel through these five names:

- **Serialization:** `diplotypeforge.output` (`write_diplotypes`) calls
  `diplotype_to_glstring`; `transform.diagnostics` builds GL keys directly from
  `haplotype.key()` tuples (a `_gl(...)` helper) — a candidate to route through
  the shared serializer.
- **Parsing:** `diplotypeforge.adapters.hfx` calls `parse_haplotype`;
  `diplotypeforge.transform.io` (dataset loader for standalone `transform`)
  calls `parse_diplotype`.

To avoid creating *new* dependencies on the current parser, downstream code
should import only the five names above (re-exported from the top-level
`diplotypeforge` package) and must not depend on module-private helpers such as
`_reject_unsupported`, `_locus_of_token`, or `UNSUPPORTED_DELIMS`.

## Recommended path for extraction (future)

1. Create `glstring` with the five-name surface and a superset grammar.
2. In DiplotypeForge, make `diplotypeforge.glstring` a thin shim that re-exports
   from the shared library, keeping the restricted-subset validation as a
   `strict=True` option so DiplotypeForge output stays within the subset.
3. Keep `GLStringError` as an alias of the shared library's error type.

This note is descriptive; no extraction is performed in v0.1.
