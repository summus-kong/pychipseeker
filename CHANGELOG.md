# Changelog

All notable changes to pychipseeker are documented here.  The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

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
