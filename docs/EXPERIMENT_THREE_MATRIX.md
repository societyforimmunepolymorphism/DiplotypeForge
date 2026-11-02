# The Three-Matrix Population-Structure Experiment

This flagship experiment makes the central Milestone-2 point concrete: the
**marginal** picture of a population can be identical across very different
mating structures, while the **joint** (diplotype) picture differs sharply. It
is the signal that methods for detecting mixture and admixture rely on.

## Setup

Two synthetic single-locus (HLA-A) source populations with *distinct* haplotype
distributions:

| haplotype | POP1 | POP2 |
|-----------|------|------|
| HLA-A*01:01 | 0.7 | 0.0 |
| HLA-A*02:01 | 0.2 | 0.0 |
| HLA-A*03:01 | 0.1 | 0.2 |
| HLA-A*11:01 | 0.0 | 0.5 |
| HLA-A*24:02 | 0.0 | 0.3 |

Three origin matrices `M`, **all with marginal ancestry 0.5 / 0.5**:

| model | matrix M | interpretation |
|-------|----------|----------------|
| random | `[[0.25, 0.25], [0.25, 0.25]]` | random union of gametes over the mixture (= `a·aᵀ`) |
| wahlund | `[[0.5, 0.0], [0.0, 0.5]]` | pure stratification: both haplotypes from one source |
| admixture | `[[0.1, 0.4], [0.4, 0.1]]` | balanced first-gen admixture, `a = 0.8` |

## Theory (independent of sampling)

- Marginal ancestry `a_i = ½(Σ_j M[i,j] + Σ_j M[j,i])` = `0.5, 0.5` for all three.
- Marginal haplotype frequency `f̄(h) = Σ_i a_i f_i(h)` is therefore **identical**
  across the three models.
- Expected ordered diplotype distribution `P(h,k) = Σ_ij M[i,j] f_i(h) f_j(k)`
  (unordered sums both orders) **differs**, because the joint origin structure
  differs.

## Result

Running the experiment (`examples/three_matrix_experiment.py`, seed 20120815):

Identical marginal haplotype frequencies across all three models:

```
  haplotype        random    wahlund  admixture
  HLA-A*01:01      0.35000    0.35000    0.35000
  HLA-A*02:01      0.10000    0.10000    0.10000
  HLA-A*03:01      0.15000    0.15000    0.15000
  HLA-A*11:01      0.25000    0.25000    0.25000
  HLA-A*24:02      0.15000    0.15000    0.15000
```

Different expected diplotype distributions, e.g.:

```
  diplotype                    random    wahlund  admixture
  HLA-A*01:01+HLA-A*01:01      0.12250    0.24500    0.04900
  HLA-A*01:01+HLA-A*11:01      0.17500    0.00000    0.28000
  HLA-A*11:01+HLA-A*11:01      0.06250    0.12500    0.02500
```

Expected homozygosity (M-based): **wahlund 0.460 > random 0.240 > admixture
0.108**. (Admixture here has large off-diagonal mass — 0.8 cross-source pairing
— and because the two populations have largely disjoint haplotypes those cross
pairs are heterozygous, so admixture has the *lowest* homozygosity. The admixed
fraction `a` is not the ancestry proportion.)

Observed (simulated, n = 50 000) matches the independent expectation closely:

```
  random    : obs_hom=0.23974  TV(obs, M)=0.00804
  wahlund   : obs_hom=0.45658  TV(obs, M)=0.00394
  admixture : obs_hom=0.10642  TV(obs, M)=0.00568
```

The total-variation distances between the observed diplotype distribution and
the M-based expectation are a few thousandths — consistent with multinomial
sampling noise at this sample size.

## Running it

```bash
uv pip install -e ".[dev,viz]"
uv run python examples/three_matrix_experiment.py --n 50000 --seed 20120815
```

Outputs (under `examples/output_experiment/` by default):

- `comparison_table.txt` — the full table shown above;
- `experiment.json` — machine-readable matrices, marginals, expected and
  observed distributions, TV distances;
- `three_matrix_experiment.png` — two-panel figure: left, the overlapping
  marginal haplotype frequencies (identical); right, the diverging diplotype
  distributions.

With `--n 0` the script computes only the theoretical expectations (fully
deterministic, no RNG).

## Why it matters

The experiment is the reference demonstration that DiplotypeForge separates the
population **model** from the sampled **diplotypes** and the derived
**statistics**: three models that look identical at the allele/haplotype margin
are cleanly distinguished at the diplotype level, both in theory and in the
simulated data.
