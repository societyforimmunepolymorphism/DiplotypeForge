# Worked Example: Frequency-Preserving Homozygosity Transform

This example increases and then decreases homozygosity while leaving **every
haplotype count unchanged**, demonstrating the core invariant of the
Milestone-3A transform.

## Setup

Four heterozygous individuals over a single locus:

```python
import numpy as np
from diplotypeforge.model.core import Allele, Diplotype, Haplotype, SimulatedIndividual
from diplotypeforge.transform import analyze_feasibility, transform_homozygosity, TransformationTarget

def hap(t): return Haplotype((Allele("HLA-A", t),))
def ind(i, a, b): return SimulatedIndividual(f"s{i}", Diplotype(hap(a), hap(b)), "P1", "P2", 0)

pop = [
    ind(0, "HLA-A*01", "HLA-A*02"),
    ind(1, "HLA-A*01", "HLA-A*02"),
    ind(2, "HLA-A*03", "HLA-A*04"),
    ind(3, "HLA-A*03", "HLA-A*04"),
]
```

Initial state:

```
diplotypes:       HLA-A*01+HLA-A*02, HLA-A*01+HLA-A*02, HLA-A*03+HLA-A*04, HLA-A*03+HLA-A*04
haplotype counts: {A*01: 2, A*02: 2, A*03: 2, A*04: 2}   (N = 4)
feasibility:      H_min = 0, H_max = 4, step = 1
homozygotes:      0
```

## Increase homozygosity to the maximum

```python
up, result = transform_homozygosity(pop, TransformationTarget(fraction=1.0), rng=np.random.default_rng(0))
```

```
diplotypes:       HLA-A*04+HLA-A*04, HLA-A*03+HLA-A*03, HLA-A*02+HLA-A*02, HLA-A*01+HLA-A*01
haplotype counts: {A*01: 2, A*02: 2, A*03: 2, A*04: 2}   ← UNCHANGED
homozygotes:      4   (reached=True, reason="exact")
```

Every (X, Y) pair was regrouped into (X, X) and (Y, Y). This is exactly the
`(A,B),(A,B) -> (A,A),(B,B)` operation from the brief, applied across both
complementary pairs.

## Decrease it back to the minimum

```python
dn, result = transform_homozygosity(up, TransformationTarget(fraction=0.0), rng=np.random.default_rng(1))
```

```
diplotypes:       HLA-A*03+HLA-A*04, HLA-A*01+HLA-A*02, HLA-A*01+HLA-A*02, HLA-A*03+HLA-A*04
haplotype counts: {A*01: 2, A*02: 2, A*03: 2, A*04: 2}   ← STILL UNCHANGED
homozygotes:      0   (reached=True, reason="exact")
```

The `(A,A),(B,B) -> (A,B),(A,B)` reverse operation recovers a fully
heterozygous sample with the identical haplotype inventory.

## Two scientifically distinct routes to HW disequilibrium

Both this transform and the Milestone-2 structure models create departures from
Hardy–Weinberg equilibrium, but they are **not** the same mechanism and must not
be conflated:

| | Structural (Wahlund / admixture) | Imposed (reassignment transform) |
|---|---|---|
| Source of HWD | the *population model*: the origin matrix `M` concentrates/disperses same-source pairings | *post-sampling* rearrangement of a fixed haplotype inventory |
| Haplotype frequencies | set by the model | preserved **exactly** from the input sample |
| Analytical diplotype law | yes — the M-based `P(h,k)=Σ_ij M[i,j] f_i(h) f_j(k)` | **none claimed**; a new sampling mechanism |
| Diagnostic treatment | observed compared to M-based and pooled-HWE | observed-after reported *separately*; explicitly labeled not M-distributed |

The transform's diagnostics report the original-M expectation, pooled-HWE, and
the observed distributions before and after, with total-variation distances, so
the imposed deviation is visible and never attributed to population structure.

## CLI

```bash
# simulate then transform in one run (config has a [transformation] block)
diplotypeforge simulate examples/transformation/simulate_transform.toml

# or transform a previously generated dataset (provenance read from columns)
diplotypeforge transform \
    --input examples/transformation/output/diplotypes.tsv \
    --config examples/transformation/transform_only.toml \
    --output-dir examples/transformation/transform_out
```

The standalone `transform` reads per-haplotype provenance from the
`origin_1`/`origin_2` columns of the diplotypes table (never reconstructed from
the GL String) and the sibling `manifest.json` for the original origin matrix.
When a transform runs inside `simulate`, the untransformed sample is also written
to `diplotypes_original.tsv`, so the original population remains recoverable.
