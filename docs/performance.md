# Performance

This page records the measurements behind the release notes and explains where
parallelism helps.  Numbers were taken on Windows 11, Python 3.11, 16 physical
cores / 24 logical processors, with the mm10 genome cache.

Reproduce them with:

```bash
python scripts/bench.py --txdb mm10 --repeat 2
python scripts/bench.py --txdb mm10 --jobs 4 --tagmatrix --bigwig signal.bw --sample 20000
```

## Annotation

Random mm10 transcript-scale peaks (200 bp – 2 kb wide, both strands).

| Input | Before 0.7 | 0.7 sequential | 0.7, 4 jobs | 0.7, 8 jobs |
|---|---:|---:|---:|---:|
| 20,000 peaks | 15.6 s | 7.7 s | – | – |
| 100,000 peaks | 81.2 s | 45.9 s | 14.9 s | 16.2 s |

The sequential gain comes from removing per-query full-chromosome scans:

- `TSSIndex.follow_precede` / `within` and `SpanIndex.follow_precede` /
  `overlapping` use `searchsorted` on the pre-sorted feature arrays instead of
  boolean masks over an entire chromosome.  Tie-breaking still returns the
  smallest subject row index, so GenomicRanges semantics are unchanged.
- The alignment pass (`align_annotation_feature=True`) reads the feature row
  directly from cached NumPy arrays instead of building a pandas `Series` for
  every peak.
- `_nearest_feature` / `_category_hit` avoid `DataFrame.empty`, which is
  surprisingly expensive at hundreds of thousands of calls.

Annotation is still the most CPU-heavy path because it performs several
per-peak interval queries.  `n_jobs=-1` is capped at eight worker processes for
annotation: every worker holds a full genome model (about 300 MB for mm10 when
rebuilt from the cached `.npz`), and beyond eight processes memory bandwidth
and disk contention outweigh the extra cores.  Pass an explicit higher
`n_jobs` if the machine has memory to spare.

All parallel runs are byte-for-byte identical to the sequential path; the
parity tests compare full annotation tables across `level`, `same_strand`,
`ignore_overlap`, `align_annotation_feature` and `overlap` settings.

## bigWig TSS and tag matrices

`extract_tss_matrix` and the bigWig-backed `get_tag_matrix` read each bigWig
interval once per genomic block and integrate coverage at the bin edges with
cumulative sums.  `pyBigWig` handles are not picklable, so each worker opens
its own handle in the pool initializer; work units are chromosome segments so
reads stay sequential on disk while dense and sparse regions balance across
cores.

| Feature | Work unit | Parallel threshold |
|---|---|---|
| `extract_tss_matrix` | chromosome blocks | 4,000 regions per worker |
| bigWig tag matrix | chromosome blocks | 4,000 windows per worker |
| peak-coverage tag matrix | contiguous window chunks | 8,000 windows per worker |

The defaults stay sequential (`n_jobs=None`), so importing and small calls
never pay process start-up.  On Windows, spawning a pool costs roughly one to
four seconds depending on the environment; a small demo track will therefore be
slower in parallel.  The thresholds above require enough work per worker to
cover that cost, and parallel extraction only runs for vectorised value types
(`sum`, `mean`, `coverage`) or the per-window `bw.stats` fallback when the
input is genuinely large.

Peak-coverage matrices are pure-Python per window, so they scale well once the
input is large enough:

| Windows × bins | Sequential | 8 jobs |
|---|---:|---:|
| 40,000 × 100 | 4.8 s | 8.6 s |
| 100,000 × 100 | 40.1 s | 10.9 s |

Sequential and parallel matrices are compared with `numpy.testing` in
`tests/test_parallel.py`, including strand handling and dropped empty rows.

## Other paths

The multiprocessing behaviour introduced in 0.6 is unchanged:

| Operation | Work unit | Parallel threshold |
|---|---|---|
| `annotate_peak` | contiguous peak chunks | 2,000 peaks |
| BAM tag matrix / TSS enrichment | feature/window chunks | 500 regions |
| `call_peaks` (bigWig) | one task per chromosome | 2+ chromosomes |

BAM FRiP uses the C-level `pysam` counter in one streaming pass plus one
indexed query per merged peak.  On Windows, any script that uses `n_jobs > 1`
must guard its entry point with `if __name__ == "__main__":` because the
`spawn` start method re-imports the main module.
