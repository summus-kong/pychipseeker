# Changelog

All notable changes to pychipseeker are documented here.  The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.8.4] - 2026-11-12

### Added

- Public analysis helpers: `annotation_summary`, `peaks_per_gene`,
  `top_genes`, `gene_type_summary`, `distance_summary`, `signal_groups` and
  `tss_enrichment_by_quartile`, all exported from the package root.
- Six analysis panels (M-R) in the ENCODE example: width ECDF by class, top
  genes, peaks per gene, TSS enrichment by signal quartile, width by
  TSS-distance class and gene biotypes, composed into a third A3 plate
  (`encode_analysis_plate.*`).
- `--panels` and `--list-panels` in `examples/make_encode_plate.py` compose an
  arbitrary subset of panels into `encode_custom_plate.*`, exported in every
  format including PDF.
- `gallery.html` is now interactive: grouped sections, filters, per-panel
  checkboxes, an auto-updating command and a copy button.

### Changed

- `signal_groups` assigns a score of zero to genes without an assigned peak,
  keeping them in the lowest group rather than filling the median; the
  quartile profiles and TSS-enrichment bars are now monotonic.
- The main plate's signal-quartile profiles share the same grouping as the
  analysis panel, so panels K and P agree by construction.
- README documents the analysis helpers, the eighteen panels and the
  selection workflow.

### Fixed

- Panel O statistics no longer collide with the bars, and panel P reserves
  space for its two-line method note.

## [0.8.3] - 2026-11-12

### Added

- A second A3 plate (`encode_supplementary_plate.*`) with six more views of
  the ENCODE peak set: chromosome distribution, width and signal by annotation
  class, non-promoter TSS distances, signal-quartile profiles, and a
  width-vs-signal hexbin with a median trend.
- `gallery.html`, a browsable index of both plates and all twelve panels with
  thumbnails and direct links to every export format.
- `_quartile_profiles()` in the example script groups promoters by the
  strongest peak score annotated to each gene.

### Changed

- `examples/make_encode_plate.py` now produces twelve standalone panels and a
  twelve-page `encode_figures.pdf` (previously six).
- Panel L carries a white annotation backdrop so its statistics remain legible
  over the density hexbins.

## [0.8.2] - 2026-11-11

### Fixed

- `tag_heatmap()` now sorts rows by decreasing mean signal, so the most
  enriched regions are drawn at the top of the map.  Previously the weakest
  rows were on top, which contradicted `plot_tss_heatmap(sort_rows=True)`
  and the usual genome-browser convention.

### Changed

- `examples/make_encode_plate.py` now exports the plate and every individual
  panel as PNG, PDF, SVG, EPS and TIFF, plus a 600 dpi PNG.  Vector files use
  Type 42 fonts and editable SVG text; rasters are flattened to RGB.
- Refined the plate: a larger donut, peak-coverage values normalised to a
  shared 0-100 % scale for the profile and heatmap, a lighter ECDF fill, and
  a fixed annotation-table footnote.

## [0.8.1] - 2026-11-11

### Fixed

- `plot_anno_pie(autopct=None)` no longer fails on Matplotlib 3.11, where
  `Axes.pie()` returns a two-value `PieContainer` instead of a three-tuple.

### Changed

- Redesigned the real-ENCODE figure plate: an Okabe-Ito colourblind-safe
  palette shared across panels, direct labels, a donut composition key,
  peak-width and cumulative-distance panels that expose the non-promoter tail,
  a single-hue tag heatmap, and a publication-style annotation table.
- Added small-sample counts and medians to the peak-width legend, and a
  method caption plus panel list to the plate footer.

## [0.8.0] - 2026-11-11

### Changed

- **Vectorised peak-coverage tag matrices.**  `get_tag_matrix(peak=...)` now
  reduces each chromosome to a piecewise-constant coverage track and evaluates
  every window as a difference of its cumulative integral.  Runtime scales with
  the number of bins instead of the number of nearby peaks: 100,000 synthetic
  6 kb windows drop from 157 s to 0.43 s, and the real mm10 promoter matrix on
  the example ENCODE peaks drops from 16.5 s to 0.84 s with bit-identical
  values.
- Process parallelism for peak-coverage matrices now only activates for very
  large window sets (50,000 windows per worker with `fork`, 400,000 with
  `spawn`), because the vectorised sequential path is faster than pool start-up
  below those thresholds.

### Fixed

- The peak-coverage overlap lookup no longer assumes peak end coordinates are
  sorted.  Inputs whose ends were not monotone could lose coverage from a long
  peak that started before a shorter one.
- Peak-free gaps no longer accumulate floating-point residues from the coverage
  track, so windows without signal are still dropped exactly.

### Added

- `scripts/bench_coverage.py`, a reproducible synthetic benchmark for the
  peak-coverage path.

## [0.7.0] - 2026-11-11

### Added

- **Parallel bigWig TSS and tag matrices.**  `extract_tss_matrix(n_jobs=...)`,
  `get_tag_matrix(n_jobs=...)` for `bigwig=...` inputs and
  `get_tag_matrix_from_bigwig(n_jobs=...)` process chromosome blocks with one
  `pyBigWig` handle per worker.  Peak-coverage tag matrices are chunked across
  workers as well; results are identical to the sequential run.
- **`docs/performance.md`** with the reproducible benchmark tables, parallel
  thresholds and the note that `n_jobs=-1` is capped at eight annotation
  workers because each process holds a full genome model.
- `scripts/bench.py --jobs` (plus `--tagmatrix`) benchmarks the parallel paths.

### Changed

- **Annotation hot paths use binary search.**  `TSSIndex.follow_precede` /
  `within` and `SpanIndex.follow_precede` / `overlapping` no longer build
  whole-chromosome masks, and the alignment pass reads cached arrays instead of
  constructing a pandas `Series` per peak.  Sequential annotation improves from
  15.6 s to 7.7 s at 20k peaks and from 81.2 s to 45.9 s at 100k peaks; 4 jobs
  annotate 100k peaks in 14.9 s (was 21.7 s).  Tie-breaking still returns the
  smallest subject row index.
- bigWig interval integration is shared by the plotting and tag-matrix code in
  `pychipseeker._bigwig`, removing duplicated cumulative-sum implementations.
- Parallel annotation prints one parent-side progress line instead of
  interleaved per-worker counters.
- CI adds Python 3.13 to the test matrix.

## [0.6.0] - 2026-11-10

### Added

- **Multiprocessing for large data sets.**  `annotate_peak(n_jobs=...)`,
  `call_peaks(n_jobs=...)`, `tss_enrichment(n_jobs=...)`,
  `qc_summary(n_jobs=...)`, `get_tag_matrix(n_jobs=...)`,
  `get_tag_matrix_from_bam(n_jobs=...)` and `build_report(n_jobs=...)` split
  work across processes; `-1` uses every core and results are identical to the
  sequential run.  Annotation on 100k mm10 peaks drops from 75.9 s to 21.7 s
  with 4 jobs.
- **CLI `--jobs`.**  The `annotate`, `callpeak`, `qc` and `report` commands
  accept `--jobs N` (`-1` = all cores).
- **`pychipseeker.parallel`.**  Public `available_cores`, `effective_jobs`,
  `chunk_bounds` and `parallel_map` helpers used by the parallel paths.

### Changed

- Annotation workers rebuild the `GeneAnnotation` from its GTF or cached
  `.npz` source instead of pickling the full model, keeping worker memory flat.
- BAM FRiP uses the C-level `pysam` counter: one streaming pass for the
  denominator and one indexed query per merged peak for the numerator, with the
  old streaming path kept for indexless files.
- The BAM tag matrix and BAM TSS matrix are chunked across workers for large
  feature sets; peak calling parallelises one task per chromosome.
- CI installs `pysam` in an always-succeeding optional step so the real-BAM
  tests run on Linux/macOS and skip cleanly on Windows.

## [0.5.0] - 2026-11-10

### Added

- **Chromosome harmonisation.**  `canonical_chrom`, `chrom_aliases`,
  `alias_map`, `harmonize_frame` and `chromosome_style` match UCSC (`chr1`)
  and Ensembl (`1`) names without rewriting unrelated contigs.  Annotation,
  signal matrices, tag matrices and QC now apply the mapping automatically.
- **Assembly detection and `TxDb="auto"`.**  `detect_assembly` uses chromosome
  sizes, a bigWig header or peak names; ambiguous human/mouse peak-only input
  raises instead of guessing.  On-demand specs gain cache listing, provenance
  (`cache info`) and partial-download cleanup (`cache prune`).
- **Quality control.**  New `peak_qc`, `frip`, `tss_enrichment` and
  `qc_summary` helpers report peak widths, mitochondrial/blacklist fractions,
  FRiP from bigWig or BAM, and ENCODE-style TSS enrichment.  The PDF report
  gains an optional QC page (`--qc`, `--bam`, `--blacklist`).
- **BAM-backed matrices.**  `get_tag_matrix(bam=...)`,
  `get_tag_matrix_from_bam` and `.bam`/`.cram` auto-detection build coverage
  matrices with `pysam`, including strand handling and binning.
- **Enrichment upgrades.**  Local ORA now reports odds ratios with 95%
  confidence intervals; `gsea`/`gsea_peak_genes` wrap preranked GSEA through
  `gseapy`; `plot_enrichment` draws bar and clusterProfiler-style dot plots;
  the CLI gains `enrich --plot`.
- **Peak workflow.**  `filter_peaks` removes blacklist, mitochondrial,
  width- or score-filtered peaks, and new `callpeak`, `filter` and `qc` CLI
  commands expose peak calling and QC without Python.
- **Asymmetric signal windows.**  `extract_tss_matrix`,
  `extract_tss_matrix_from_signal`, profile/heatmap helpers and
  `PeakAnnotation` methods accept independent `upstream`/`downstream`
  values, plus `drop_empty`, `return_mask` and `dtype`.
- **Streaming annotation input.**  `read_gtf` streams large files
  (`chunksize=`), truncates at `##FASTA`, explodes multi-parent GFF3 rows and
  filters by feature type (`features=`).

### Changed

- `get_tag_matrix(bigwig=...)` vectorises `sum`/`mean`/`coverage` over
  bigWig intervals; `max`/`min`/`std` keep the `pyBigWig.stats` fallback.
- Report and CLI signal panels accept asymmetric windows; the CLI cache
  listing gains `--json`, and `report` accepts `--upstream`/`--downstream`.
- The test suite adds real bigWig QC/BAM stubs, chromosome aliases, GFF3 and
  coverage gates (CI enforces 75% line coverage).
- Windows installs are no longer blocked by optional signal/BAM dependencies:
  the `bigwig`, `bam` and `all` extras skip `pyBigWig`/`pysam` on Windows,
  where they have no wheel for every Python version.  Install them via
  conda-forge/bioconda when needed.
- CI adds a core-only test job (no optional extras) and a Windows wheel
  install/smoke job; the optional `pyBigWig` install step always succeeds so
  skipped bigWig tests never appear as errors.

### Fixed

- Failed genome downloads remove their `.part` files; `cache prune` clears
  interrupted bundles and orphaned provenance records.
- GFF3 `##FASTA` sections no longer leak sequence rows into the annotation.
- Multi-parent GFF3 exons/UTRs are duplicated per parent instead of keeping a
  comma-separated transcript ID.

## [0.4.0] - 2026-11-10

### Added

- **Excel-free annotation, now substantially faster.**  Vectorised interval
  indexes, direct NumPy column access and array-backed strand handling reduce
  the first 11k-peak annotation from 45.7 s to about 7.6 s; `mm10` loading drops
  from 36.9 s to about 1.7 s.
- **Fast TSS signal matrices.**  `extract_tss_matrix` reads each chromosome's
  bigWig intervals once and integrates them with cumulative sums.  A full
  mm10 transcript set (81,540 TSS) drops from roughly 54 minutes to about 11 s,
  with results identical to `bw.stats`.
- **CLI `report` command** and `build_report()`: a multi-page PDF with the
  annotation pie/bar, TSS distance, coverage, tag-matrix and optional bigWig
  TSS profile/heatmap pages, `annoStat`, and a methods page.
- **CLI `cache list|download|clear`** for installed and downloaded genome
  bundles.
- **`--txdb`** as a mutually exclusive alternative to `--gtf` for `annotate`
  and `report`.
- **Five on-demand assemblies**: `dm6`, `danRer11`, `rn7`, `mm39` (GRCm39) and
  `sacCer3`, downloaded as Ensembl GTFs and cached as compact npz bundles.
  `ensembl:<species>[@release]`, `gencode:<human|mouse>[@release]` and
  `url:<gtf>` specs are accepted anywhere `TxDb=` is accepted.
- **Peak callers**: dependency-free threshold calling from a bigWig
  (`call_peaks`) and an optional MACS2/MACS3 wrapper (`call_peaks_macs`), plus
  `write_peaks`.
- **Input formats**: narrowPeak, broadPeak, gappedPeak (including `.gz`),
  BED headers, `track`/`browser`/comment lines, bigBed via `bigBedToBed` or
  optional `pybigtools`, and BAM/CRAM via optional `pysam`.
- **Local enrichment analysis**: `ora()` implements the hypergeometric test
  with BH/Bonferroni correction and explicit backgrounds; `read_gmt()` reads
  GMT gene sets.  `enrichr()` forwards custom backgrounds and can recompute
  adjusted p-values.
- **True UpSet layout** for `upsetplot()` (intersection bars, set sizes and
  dot matrix).
- **Packaging and quality**: `py.typed`, optional `bigwig`/`bigbed`/`bam`/
  `enrich`/`all` extras, GitHub Actions CI across Python 3.9-3.12 on Linux,
  Windows and macOS, ruff and mypy configuration, `scripts/bench.py`, docs and
  this changelog.

### Changed

- `pyBigWig` is now an optional `bigwig` extra instead of a hard dependency.
  Only direct `.bw`/`.bigWig` code paths require it.
- `supported_genomes()` now returns bundled and downloadable assemblies;
  `BUNDLED_GENOMES` and `GENOME_SOURCES` expose the distinction.
- Peak readers preserve ENCODE column names (`signalValue`, `pValue`,
  `qValue`, `peak`) instead of renaming them to `extra_*`.

### Fixed

- BED headers with fewer columns than data rows no longer trigger a pandas
  parser error.
- Comment lines beginning with `#` inside BED files are ignored consistently.
- `warnings.warn` calls now report the caller's stack level.

## [0.3.0] - 2026-10-26

- ChIPseeker 1.48-compatible annotation engine, tag matrices, profiles,
  heatmaps, compatibility helpers, bundled hg19/hg38/mm9/mm10 tables, GEO
  helpers and the real ENCODE H3K4me3 example.
