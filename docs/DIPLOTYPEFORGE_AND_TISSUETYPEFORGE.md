# DiplotypeForge and TissueTypeForge

DiplotypeForge and TissueTypeForge are two complementary components of the
Society for Immune Polymorphism (SIP) **HLA Forge** initiative. They are
**separate, independently runnable projects** joined by a single, explicit data
contract: the **phased GL String** representing a known diploid genotype.

## Division of responsibilities

| | DiplotypeForge | TissueTypeForge |
|---|---|---|
| Input | population **haplotype-frequency** distributions (HFX) | **known** phased diplotypes (GL Strings) |
| Produces | synthetic individuals with **known, phased diplotypes** | simulated **laboratory typing observations** from those diplotypes |
| Models | population structure (K×K origin matrix), HW disequilibrium (frequency-preserving transform) | typing technology, ambiguity, error, resolution/dropout |
| Knows about | frequencies, ancestry, phase, provenance | assays, calls, ambiguity, noise |

DiplotypeForge establishes **genetic truth**; TissueTypeForge degrades that
truth into **observed data** as a laboratory would. Neither depends on the
other's internals at runtime.

## The shared contract: phased GL String

The hand-off is a phased diplotype GL String — two `~`-joined, locus-aligned,
fully-qualified haplotypes separated by `+`, with phase preserved:

```
A*01:01g~C*07:01g~B*08:01g~DRB1*15:01~DQB1*06:02+A*02:01g~C*03:04g~B*40:01g~DRB1*04:04~DQB1*03:02g
```

DiplotypeForge writes this as the `diplotype_glstring` column of
`diplotypes.tsv[.gz]`. TissueTypeForge reads the same column as its known-truth
input. The exact accepted/produced subset is specified in
[`GLSTRING_INTEROP.md`](GLSTRING_INTEROP.md). A shared SIP `glstring` library is
proposed to host this grammar for both projects; until then each project
implements the restricted subset independently against the same contract.

## A typical end-to-end pipeline

```
HFX frequencies
      │  DiplotypeForge: validate -> simulate [-> transform]
      ▼
diplotypes.tsv  (known, phased truth + per-haplotype origin + manifest)
      │  TissueTypeForge: simulate typing observations
      ▼
simulated typing results  (ambiguous / error-laden observations)
      │  evaluation
      ▼
HF-estimation / imputation / QC algorithms benchmarked against known truth
```

Because DiplotypeForge records a reproducibility manifest (seed, software and
numpy versions, input checksums, HFX schema provenance), a TissueTypeForge run
can be tied back to the exact genetic truth and the exact inputs that produced
it.

## Why keep them separate

- **Separation of concerns:** population genetics vs. assay simulation evolve
  independently.
- **Reusability:** the known-truth diplotypes are useful on their own (registry
  simulation, structure studies), and TissueTypeForge can consume truth from any
  source that emits the shared GL String contract.
- **Testability:** each side can be validated against the contract in isolation.
