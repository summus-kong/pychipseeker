# Changelog

All notable changes to pychipseeker are documented here.  The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

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
  confidence intervals; `gsea`/`gsea_peak_genes` wrap pranked GSEA through
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
