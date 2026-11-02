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
Peak memory is the `tracemalloc` Python-allocation peak — see the measurement
caveat below; it under-counts native/NumPy and RSS. The *relative* picture —
transform-enabled uses materially more memory and time — is the stable,
portable conclusion.)

## Memory complexity

Let `N` be individuals, `H` distinct haplotypes, `L` loci.

### Measurement caveat (tracemalloc)

The "peak memory" figures above are measured with Python's `tracemalloc`, which
tracks **Python-level allocations only**. It does **not** capture native
allocations made outside the CPython allocator — notably NumPy's own buffers
(backed by the system allocator), interpreter/C-extension overhead, or memory
held by the OS allocator after Python frees objects (RSS is typically higher
than the `tracemalloc` peak). Treat these numbers as a **relative, reproducible**
comparison between modes, not as an absolute resident-set-size (RSS) budget. For
RSS-accurate figures use an OS-level tool (`/usr/bin/time -v`, `tracemalloc` +
`resource.getrusage`, or a profiler such as `memray`).

### Current implementation

- **Ordinary simulation (streaming-capable):** sampling yields individuals from
  a generator; the prepared per-population sampling table is `O(H · L)` and is
  the only persistent structure needed to *produce* a stream. If the caller
  consumes the generator and streams each individual straight to the TSV writer,
  the extra memory beyond the frequency table is `O(1)` in `N`. The high-level
  `Simulator.run()` currently **materializes** the `N` individuals (`O(N)`) to
  compute summary statistics and support replicate aggregation; the benchmark
  script also materializes. A pure streaming write path exists at the
  `Sampler.sample(...)` generator level for callers that do not need the
  in-memory sample.
- **Transformation-enabled (current implementation):** the reassignment holds
  the full haplotype-copy inventory in memory — it groups the `2N` copies by
  haplotype key and reconstructs `N` pairs — so the **current implementation**
  is `O(N)` in memory. The measured ~3× peak increase (217 → 649 MB) reflects
  the `2N` copy records plus the rebuilt pair list.

### What is *inherently* required vs. an implementation choice

The earlier documentation described the transform as *inherently* `O(N)` memory.
That is too strong. The `O(N)` cost of the current code is an **implementation
choice**, not a lower bound of the problem. Distinguish four separate
requirements:

1. **Frequency-preserving pairing (counts only).** Producing *a* valid pairing
   with a feasible homozygosity count needs only the per-haplotype **counts**
   `c_h` and the per-key homozygote allocation `hom_h` (both `O(H)`). Given those
   integer allocations, individuals can be *emitted as a stream* — repeatedly
   emit `hom_h` homozygous pairs per key, then interleave the leftover copies
   heterozygously — without ever materializing all `2N` copies. This is `O(H)`
   working memory, independent of `N`.

2. **Retaining individual identities.** If each sampled individual must keep a
   stable identity (e.g. a specific `sample_id` bound to a specific pre-transform
   draw) across the rearrangement, then the mapping from identities to
   haplotype copies must be retained — `O(N)`. If identities may be assigned
   *after* pairing (the common case for synthetic data), this cost disappears.

3. **Retaining source-population provenance.** Provenance is a label per
   haplotype **copy**. If every copy's origin must travel individually, an
   `O(N)`-sized association is needed. If provenance is only required at the
   aggregate level (counts of copies per (haplotype, origin)), it collapses to
   `O(H · K)` for `K` populations and can also be streamed.

4. **External-memory / out-of-core algorithms.** Even when individual identities
   and per-copy provenance are both required, the inventory need not reside in
   RAM: copies can be grouped by key on disk (external sort / bucket by key),
   paired in a streaming merge, and written out, giving memory proportional to
   the number of distinct keys plus a bounded I/O buffer. This is the standard
   external-memory pairing approach and would make even the identity- and
   provenance-preserving transform sub-`O(N)` in RAM at the cost of disk I/O.

**Summary.** Ordinary simulation can be `O(1)`-extra streaming. The current
transform is `O(N)` because it materializes the inventory to preserve per-copy
identity *and* provenance in the simplest correct way. The *theoretical* minimum
for frequency-preserving pairing is `O(H)` (counts only); identity and
provenance requirements, and the choice to keep everything in RAM rather than
out-of-core, are what raise the current implementation to `O(N)`. A streaming or
external-memory transform is a plausible future optimization and is **not**
attempted in this release-hardening phase.

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
