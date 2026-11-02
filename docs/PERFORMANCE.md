# Performance & Scaling (Milestone 3C)

DiplotypeForge is validated at **≥ 1,000,000 individuals**. This document records
the benchmark methodology, measured numbers, and the memory complexity of the
ordinary (streaming) simulation versus the transformation-enabled path.

## Benchmark

```bash
# ordinary large simulation
uv run python examples/benchmark_large.py --n 1000000 --haplotypes 2000

# transformation-enabled
uv run python examples/benchmark_large.py --n 1000000 --haplotypes 2000 --transform 0.3
```

The script measures runtime (split into sampling / transform / write phases)
with `time.perf_counter`, peak Python memory with `tracemalloc`, throughput
(individuals/second), and the output file size, and writes a `benchmark.json`.

## Measured results (1,000,000 individuals, 2,000 haplotypes, 5 loci)

| | ordinary | transformation-enabled (target 0.30) |
|---|---|---|
| sampling | 35.4 s | 35.1 s |
| transform | — | 50.2 s |
| write (uncompressed TSV) | 3.7 s | 3.8 s |
| **total** | **39.3 s** | **89.3 s** |
| throughput | 25.5k ind/s | 11.2k ind/s |
| **peak memory** | **217 MB** | **649 MB** |
| output file size | 162 MB | 162 MB |
| transform reached target | — | yes (exact) |

(Measured on the reference development machine; absolute numbers vary by host.
The *relative* picture — transform-enabled uses materially more memory and time
— is the stable, portable conclusion.)

## Memory complexity

Let `N` be individuals, `H` distinct haplotypes, `L` loci.

- **Ordinary simulation (streaming-friendly):** sampling yields individuals from
  a generator; the prepared per-population sampling table is `O(H · L)` and is
  the only persistent structure needed to *produce* a stream. If the caller
  consumes the generator and streams each individual straight to the TSV writer,
  the extra memory beyond the frequency table is `O(1)` in `N`. The high-level
  `Simulator.run()` currently *materializes* the `N` individuals (`O(N)`) to
  compute summary statistics and to support replicate aggregation; the benchmark
  script above also materializes. A pure streaming write path exists at the
  `Sampler.sample(...)` generator level for callers that do not need the
  in-memory sample.
- **Transformation-enabled:** the frequency-preserving reassignment **must hold
  the full haplotype inventory** to regroup it — it groups the `2N` copies by
  haplotype and reconstructs `N` pairs. This is inherently `O(N)` memory and
  cannot be streamed, because a copy's partner may be any other individual's
  copy. The measured ~3× peak-memory increase (217 → 649 MB) reflects the
  inventory of `2N` copy records plus the rebuilt pair list. This is the
  documented, expected cost of a global frequency-preserving rearrangement
  (recognized in the brief).

So: ordinary simulation can be `O(1)`-extra streaming; a transformation-enabled
run is `O(N)` by necessity.

## Compact integer haplotype identifiers

Internally the sampler already addresses haplotypes by **integer index** into a
frozen, sorted array (`sampling/sampler.py`), so the hot sampling path does not
carry string keys. The transform groups copies by a haplotype *key* (the tuple
of allele tokens); for very large runs these keys could be interned to small
integers to reduce hashing and memory. We keep the explicit keys in Milestone 3
for clarity and correctness and note integer interning as a low-risk future
optimization — not pursued now, per the "no premature optimization" guidance.

## Streaming vs. in-memory — guidance

- For the largest *ordinary* simulations where only the TSV output is needed,
  consume `Sampler.sample(...)` as a generator and write incrementally.
- Summary statistics, replicate aggregation, and the transform all require the
  in-memory sample (or inventory) and are therefore `O(N)`.
- Compressed TSV (`compress = true`) trades a modest amount of CPU for roughly
  an order-of-magnitude smaller files.

## Non-goals

No distributed computing and no premature micro-optimization were introduced.
The benchmark characterizes the current implementation honestly so that future
optimization can be prioritized against real measurements.
