# Open Scientific Questions & Decisions

Ambiguities that require scientific decisions. Each lists the question, the
decision taken for Milestone 1 (so implementation can proceed), and what would
change it. Decisions marked **(reviewer)** were revised following design review.

## Q1. Frequency-sum tolerance and renormalization **(reviewer)**
Haplotype frequencies from EM estimation rarely sum to exactly 1.0. How close is
"valid", and do we renormalize?

**Decision:** validation is **strict by default**. A small configurable
*rounding tolerance* `sum_tolerance` (default `1e-6`) absorbs floating-point
round-off only; within it we renormalize silently (pure round-off). Any
deviation larger than `sum_tolerance` is an **error** unless the configuration
*explicitly authorizes* a larger renormalization via
`allow_renormalization = true` together with `renorm_tolerance` (the maximum
absolute deviation the user accepts). When an authorized renormalization is
applied, the original sum, the tolerance used, and the scale factor are
**recorded in the manifest** (`renormalization` block). We never auto-correct a
5%-scale deviation by default.

**Would change if:** the SIP community defines a canonical tolerance in HFX.

## Q2. Duplicate / zero / negative / missing frequencies **(reviewer)**
**Decision:**
- Duplicate haplotype identifiers are an **error by default**. Aggregating
  duplicates is an **explicit preprocessing operation** (`aggregate_duplicates`)
  that sums their frequencies, emits a **warning**, and records **provenance**
  (which keys were merged, their contributing frequencies) in the distribution
  metadata and manifest. It is never performed implicitly.
- Negative frequency → error always. Missing frequency field → error.
- Zero frequency → warning; the entry is retained in the table but cannot be
  sampled (dropped from the sampling array with a note).

## Q3. One HFX file = one population?
HFX `cohortDescription.population` is an array but the file carries a single
`frequencyData`/CSV, so it is ambiguous whether a multi-population cohort is
poolable or per-population.

**Decision:** treat **one HFX file as one source population's distribution**.
Multi-population simulations compose **multiple HFX files**, each mapped to one
internal population id in the TOML config. If a file lists multiple populations,
use the config-provided `id` (or first population name) and warn.

## Q4. GL String representation and allele identifiers **(reviewer)**
Should allele tokens be fully qualified, and do we reuse a GL String library?

**Decision:**
- Use **fully qualified allele identifiers** including the genetic-system
  prefix, e.g. `HLA-A*02:01`, `HLA-DRB1*15:01` (KIR: `KIR2DL1*0010101`). The HFX
  adapter normalizes bare tokens (e.g. `A*02:01`) to fully qualified form using
  the system implied by `nomenclatureUsed.database`
  (`IPD-IMGT/HLA` → `HLA-`, `IPD-KIR` → `KIR`), recording the normalization.
- **Library investigation (item 4):**
  - `nmdp-bioinformatics/pyglstring` — a GL String *sanity checker* (locus-block
    consistency), **not** a structured parser producing typed diplotype objects.
    It is **archived/read-only** (Jun 2025) and **LGPL/GPL-3.0** (copyleft,
    incompatible with a clean Apache-2.0 dependency). Not adopted.
  - `mhcgnomes`, `py-hla-match` — nomenclature parsers/matchers, not GL String
    diplotype serializers. Not applicable.
  - Decision: implement a **small local serializer/parser** restricted to the
    **unambiguous phased diplotypes** DiplotypeForge itself produces (delimiters
    `~` and `+` only). We explicitly **reject** allele ambiguity (`/`), genotype
    ambiguity (`|`), and locus delimiters (`^`) in Milestone 1 — we do **not**
    implement the full GL String grammar. The design keeps a clean seam so a
    full grammar library can be swapped in later.

## Q5. Expected diplotype distribution & HWE level of analysis **(reviewer)**
How are expected diplotype frequencies defined, and at what level is HWE?

**Decision:** the expected **ordered** diplotype distribution under an arbitrary
population-origin matrix `M` (K source populations, per-population haplotype
frequencies `f_i(h)`) is:

```
P(h, k) = Σ_i Σ_j  M[i, j] · f_i(h) · f_j(k)
```

For **unordered** diplotypes with `h ≠ k`, sum both orders:
`P_unordered({h,k}) = P(h,k) + P(k,h)`; for `h = k`, `P_unordered = P(h,h)`.
This is the model used for expected frequencies and simulation diagnostics.

This **M-based expectation is distinct from a pooled-HWE expectation**. The
pooled-HWE null takes the marginal (mixture) haplotype frequency
`f̄(h) = Σ_i a_i f_i(h)` (with `a_i` the marginal ancestry proportions) and
assumes random union: `P_HWE(h,k) = f̄(h) f̄(k)` (×2 for `h≠k` unordered).
Under structure/admixture (`M` not rank-1 in the mixing sense) the two differ,
and that difference is exactly the Wahlund/structure signal. Diagnostics report
**both** and label each. For the single-population Milestone-1 case `K=1`,
`M=[[1]]`, the two coincide: `P(h,k)=f(h)f(k)`, HWE at the **multilocus-
haplotype** level. No statistic is reported without naming its level and null
model.

## Q6. Ordered vs unordered diplotypes in output frequencies
**Decision:** report **unordered** diplotype (genotype) frequencies, while
retaining the ordered sampled pair per individual in the raw output. GL String
output preserves the sampled order; unordered equivalence is the documented
comparison rule (swap the two `+`-separated haplotypes).

## Q7. Allele resolution / nomenclature validation depth
**Decision:** **syntactic** validation only in M1 — tokens must parse as
`[SYSTEM-]LOCUS*fields[suffix]` consistent with declared loci/resolution. We do
**not** check tokens against an IPD database (keeps the dependency footprint
modest, avoids DB-version coupling). `nomenclatureUsed.version` is recorded for
provenance. Database-backed validation is a future option.

## Q8. Sample identifier scheme
**Decision:** stable, zero-padded, replicate-aware ids:
`DF-<run8>-<replicate>-<index>` where `run8` is a short hash of the
**canonicalized** configuration (+ seed + software version). Deterministic given
inputs; stable across reruns.

## Q9. `count` column semantics in CSV HFX
**Decision:** `frequency` is authoritative for sampling. `count` is preserved as
provenance and may inform observed-vs-expected sampling-error diagnostics. If
only `count` is present, derive frequency as `count / Σcount` with a warning.

## Q10. Copy-number variation (CNV) and gene absence **(reviewer)**
Can CNV be handled by adding a copy-number annotation to a fixed-locus diploid
representation?

**Decision:** **No — and we explicitly do not assume so.** Real CNV (e.g. KIR
gene presence/absence, DRB3/4/5 and C4 copy-number, gene duplications) changes
the *cardinality and structure* of a genotype: the number of gene copies varies
per individual and per locus, so a fixed two-haplotype, fixed-locus array cannot
represent it by merely attaching a scalar count. Properly modeling it requires a
**different data-model extension** — variable-copy loci, presence/absence
alleles, and a genotype representation that is not a fixed-width ordered pair.
This is documented as a **future data-model extension** and is **not implemented
in Milestone 1**. Milestone 1 fixes copy number at 2 per locus. The internal
interfaces avoid baking the fixed-pair assumption into the sampling engine
(which addresses haplotypes by index), so the extension can be added without
disrupting the diploid base, but it is deliberately deferred and not faked via
an annotation.

## Q11. Non-HLA systems (KIR/LILR) in M1
**Decision:** the engine is locus-agnostic and will sample any valid frequency
table; KIR resolution enums are accepted by the HFX adapter and ids are fully
qualified (`KIR...`). KIR gene presence/absence is **not** modeled in M1 (see
Q10). KIR haplotypes over the declared loci simulate as ordinary fixed-locus
haplotypes.


---

# Milestone 3 decisions

## Q12. Feasible minimum / maximum homozygosity (analytically tractable) **(M3)**
Given N individuals = 2N haplotype copies with per-haplotype counts `c_h`
(`Σ c_h = 2N`), which homozygosity counts are achievable by rearranging the 2N
copies into N unordered pairs while preserving every `c_h`?

**Decision (exact, algorithm-independent).** A target of `H` homozygous pairs is
achievable iff per-key homozygote counts `hom_k ∈ [0, ⌊c_k/2⌋]` exist with
`Σ hom_k = H` and heterozygous-pairable leftovers (`max_k (c_k−2hom_k) ≤ N−H`).
Equivalently, with `lo_k(H) = max(0, ⌈(c_k − (N−H))/2⌉)`:

```
H achievable  iff  Σ_k lo_k(H)  ≤  H  ≤  Σ_k ⌊c_k/2⌋
```

- `H_min` = smallest achievable `H` (≥ the classic bound `max(0, c_max − N)`).
- `H_max` = largest achievable `H`. **Note:** `Σ ⌊c_k/2⌋` is only an *upper
  bound* — it is not always itself achievable, because near the maximum the
  leftover (odd-count) copies can fail to be heterozygous-pairable; the exact
  `H_max` may be 1–2 below `Σ⌊c_k/2⌋`. (This corrected an initial naive
  assumption during implementation.)
- The achievable set is a contiguous **stride**: step 1 generally, but **step 2
  (parity)** when exactly two distinct haplotypes carry all copies — then only
  every other `H` is reachable.
- Homozygosity proportion `= H/N`, on a `1/N` grid.

`is_reachable(H)` exposes the exact predicate; it was validated against brute
force over all count distributions for small N. See `docs/DESIGN.md §5`.

## Q13. Rounding policy for discrete targets **(M3)**
A requested homozygosity *fraction* `p` maps to a target count `H* = round(p·N)`
(round half to even). A requested *delta* `Δ` maps to
`H* = round((p_0 + Δ)·N)` where `p_0` is the initial fraction. `H*` is then
**clamped** to `[H_min, H_max]`.

**Decision:** we **never silently substitute** a different target. We report the
requested fraction, the integer target before clamping, the achieved count, and
— if clamping or rounding changed the effective target — an explicit `reason`
(`rounded_to_1/N_grid`, `clamped_to_feasible_min`, `clamped_to_feasible_max`).
`reached` is true only if achieved == the rounded, unclamped target.

## Q14. Separate sampling vs transformation RNG **(M3)**
Enabling a transformation must not change the originally sampled diplotypes.

**Decision:** two independent child streams are spawned from the run
`SeedSequence`: a **sampling** stream and a **transformation** stream. The seed
hierarchy is:

```
root = SeedSequence(seed)
  ├─ children = root.spawn(replicates)          # one per replicate
  │    └─ child_k.spawn(2) -> [sampling_k, transformation_k]
```

Because the sampling stream is derived before and independently of the
transformation stream, toggling the transformation on/off leaves the sampling
stream — and therefore the original diplotypes — byte-identical.

## Q15. Does the transform sample uniformly from valid pairings? **(M3)**
**Decision: No, and we do not claim it does.** The implemented reassignment
strategy is a *targeted* algorithm that reaches a requested feasible homozygous
count with a bounded number of operations; it is **not** a uniform sampler over
the set of all valid matchings achieving that count. An optional randomization
step (shuffling candidate order with the transformation RNG) avoids deterministic
artifacts (always touching the same individuals/haplotypes first) but does not
confer uniformity. The documentation states this explicitly.

## Q16. Analytical distribution of the transformed sample **(M3)**
**Decision:** the frequency-preserving reassignment defines a **new sampling
mechanism**; it does **not** produce the original M-based diplotype
distribution. We do **not** invent an analytical distribution for the
transformed sample. Diagnostics report the original-M expectation, the
pooled-HWE expectation, and the observed distributions *before* and *after*
transformation, with total-variation distances, and explicitly label the
transformed sample as not M-distributed.

## Q17. Standalone transform input & provenance **(M3)**
**Decision:** the standalone `transform` command reads the `diplotypes.tsv[.gz]`
table, which already carries `origin_1`/`origin_2` columns — per-haplotype
provenance is therefore read directly, **not** reconstructed from GL Strings
(which do not encode source population). The run `manifest.json` (origin matrix,
population order, sampled pair proportions) is read as a sidecar when present so
the transform can retain and report the original provenance. A haplotype's
origin label travels with it during reassignment.
