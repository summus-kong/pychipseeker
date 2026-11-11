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
| peak-coverage tag matrix | contiguous window chunks | 50,000 windows per worker (400,000 with `spawn`) |

The defaults stay sequential (`n_jobs=None`), so importing and small calls
never pay process start-up.  On Windows, spawning a pool costs roughly one to
four seconds depending on the environment; a small demo track will therefore be
slower in parallel.  The thresholds above require enough work per worker to
cover that cost, and parallel extraction only runs for vectorised value types
(`sum`, `mean`, `coverage`) or the per-window `bw.stats` fallback when the
input is genuinely large.

## Peak-coverage tag matrices

`get_tag_matrix(peak=...)` reduces each chromosome to a piecewise-constant
coverage track and evaluates every window as a difference of the track's
cumulative integral.  Windows with the same chromosome and bin count are
evaluated in vectorised batches, so the per-window cost is proportional to the
number of bins rather than to the number of nearby peaks.  This also fixes a
latent lookup bug: the old path used `searchsorted` on peak end coordinates
without requiring them to be sorted, which could drop coverage from a long peak
that started before a shorter one.  Real ENCODE peaks and mm10 promoters give
bit-identical matrices; synthetic inputs with non-monotone ends now include the
previously skipped peaks.

Synthetic benchmark (100,000 peaks, 6 kb windows, 100 bins, Windows 11):

| Windows × bins | 0.7 | 0.8 sequential | 0.8, 8 jobs |
|---|---:|---:|---:|
| 40,000 × 100 | 51.0 s | 0.22 s | 0.22 s |
| 100,000 × 100 | 157.3 s | 0.43 s | 0.44 s |

The real mm10 promoter set (55k genes, 6 kb windows, 100 bins) drops from
16.5 s to 0.84 s on the example ENCODE peak file.  Because the sequential path
is now so fast, process parallelism only starts for very large window sets:
50,000 windows per worker with `fork`/`forkserver` and 400,000 per worker when
the start method is `spawn` (Windows and macOS).  Smaller requests keep
`n_jobs` accepted but build the matrix in-process, which is faster than paying
the pool start-up cost.

Reproduce with:

```bash
python scripts/bench_coverage.py --peaks 100000 --windows 100000 --jobs 1 8
```

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
