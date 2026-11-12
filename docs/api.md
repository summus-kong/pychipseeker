# API overview

`pychipseeker` exposes both a small object API and ChIPseeker-style functional
helpers. Most workflows use `annotate_peak` once and then call methods on the
returned `PeakAnnotation`.

```python
import pychipseeker as pcs

ann = pcs.annotate_peak("peaks.bed", TxDb="mm10")
print(ann.annoStat)
df = ann.to_dataframe()
```

All plotting helpers return matplotlib `(fig, ax)` objects (or `(fig, axes)` for
multi-panel figures), so they compose with ordinary matplotlib code.

!!! tip "Every parameter, one page"

    The [function reference](reference/functions.md) lists all 94 public
    top-level functions with every parameter, its type, default and meaning.
    It is generated from the live signatures, so it always matches the code.

## Core objects

| Object | What it holds |
|---|---|
| `PeakAnnotation` | annotated peaks, `annoStat`, peak metadata and plotting methods |
| `GeneAnnotation` | genes, transcripts and exon/intron features built from a GTF/GFF or genome bundle |
| `TagMatrix` | a promoter/window matrix with `rows`, `upstream` and `downstream` metadata |

### `PeakAnnotation`

| Member | Description |
|---|---|
| `annoStat` | ChIPseeker-style `Feature` / `Frequency` / `Count` table |
| `detail` | detailed annotation table including genomic feature columns |
| `peaks` | the input peak table (including `score` when available) |
| `to_dataframe()` | one row per annotated peak, snake_case columns |
| `to_chipseeker()` | ChIPseeker-like column names (`peakChr`, `geneName`, …) |
| `as_granges()` | GRanges-like representation |
| `summary()` | major annotation classes |
| `anno_stat()` | detailed feature table |
| `drop_anno()` | drop annotation and re-annotate with new settings |
| `seq2gene()` / `seq2gene_links()` | peak-to-gene links for enrichment |
| `plot_anno_pie()` / `plot_anno_bar()` / `plot_dist_to_tss()` | annotation figures |
| `get_tag_matrix()` / `plot_peak_prof()` / `peak_heatmap()` | matrix helpers |
| `plot_tss_profile()` / `plot_tss_heatmap()` | bigWig TSS figures |
| `vennpie()` / `upsetplot()` | overlap views |

## Functional API

### Annotation and peaks

```python
ann = pcs.annotate_peak(
    "peaks.bed",
    "genes.gtf",              # or TxDb="mm10"
    tss_region=(-3000, 3000),
    level="transcript",       # or "gene"
    n_jobs=1,
)
```

Key helpers: `annotate_peak` / `annotatePeak`, `read_peak_file` /
`readPeakFile`, `read_peaks`, `annotate_db`, `seq2gene`, `get_flanking_genes`.

### Genomes

`load_genome`, `load_annotation`, `supported_genomes`, `download_genome`,
`download_gtf`, `cached_genomes`, `genome_provenance`, `detect_assembly`,
`clear_genome_cache`, `prune_genome_cache`.

See [Genomes and TxDb](genomes.md).

### Tag matrices and profiles

`get_promoters` / `getPromoters`, `get_bio_region` / `getBioRegion`,
`make_bio_region_from_granges`, `get_tag_matrix` / `getTagMatrix`,
`get_tag_matrix_from_bigwig`, `get_tag_matrix_from_bam`, `extract_tss_matrix`,
`plot_avg_prof` / `plotAvgProf`, `plot_peak_prof` / `plotPeakProf`,
`plot_multi_prof` / `plotMultiProf`, `tag_heatmap` / `tagHeatmap`,
`plot_tss_profile`, `plot_tss_heatmap`.

### Analysis summaries

| Function | Purpose |
|---|---|
| `annotation_summary` | counts, widths and signal per annotation class |
| `peaks_per_gene` | one row per gene with peak burden |
| `top_genes` | genes ranked by peak count |
| `gene_type_summary` | biotype counts and percentages |
| `distance_summary` | promoter fraction and distance quantiles |
| `signal_groups` | balanced signal groups for matrix rows |
| `tss_enrichment_by_quartile` | centre-to-flank enrichment per group |

See [Analysis recipes](analysis.md).

### Quality control

`qc_summary`, `peak_qc`, `frip`, `tss_enrichment`.

### Overlap and enrichment

`overlap`, `shuffle`, `vennplot`, `vennpie`, `upsetplot`, `covplot`,
`enrich_peak_overlap`, `enrich_anno_overlap`, `ora`, `enrichr`, `gsea`,
`read_gmt`, `plot_enrichment`.

See [Enrichment](enrichment.md).

### Reports

`build_report(...)` builds the same multi-page PDF as
`pychipseeker report`. See [Command line](cli.md#report).

## Output columns

`ann.to_dataframe()` returns one row per annotated peak:

| Column | Meaning |
|---|---|
| `peak_id` | peak name (BED column 4) or generated ID |
| `chrom`, `start`, `end` | 0-based, half-open peak coordinates |
| `width`, `center` | peak width and centre |
| `strand` | peak strand parsed from the input |
| `annotation` | ChIPseeker label, e.g. `Promoter (<=1kb)` |
| `gene_id`, `gene_name` | assigned gene |
| `gene_chrom`, `gene_start`, `gene_end`, `gene_strand`, `gene_length`, `gene_type` | assigned gene model |
| `transcript_id`, `transcript_name`, `transcript_type` | assigned transcript (transcript level) |
| `distance_to_tss` | signed distance from the closer peak end to the assigned TSS; negative = upstream |
| `nearest_gene_id`, `nearest_gene_name` | nearest gene by TSS |
| `nearest_transcript_id`, `nearest_distance` | nearest transcript and distance |

`ann.to_chipseeker()` renames these to `peakChr`, `peakStart`, `peakEnd`,
`geneId`, `geneName`, `distanceToTSS`, and so on. Coordinates remain 0-based and
half-open; subtract 1 from `start` when moving to 1-based formats.

## Compatibility aliases

The package intentionally ships both naming styles:

```python
ann = pcs.annotatePeak("peaks.bed", TxDb="mm10")
fig, ax = pcs.plotAnnoPie(ann)
fig, ax = pcs.plotDistToTSS(ann)
```

The camelCase names cover the common ChIPseeker entry points. The
[PARITY.md](https://github.com/summus-kong/pychipseeker/blob/main/PARITY.md)
file records the exact compatibility matrix and known differences.

## Source layout

| Module | Contents |
|---|---|
| [`annotation.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/annotation.py) | annotation engine and `PeakAnnotation` |
| [`analysis.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/analysis.py) | summary helpers |
| [`tagmatrix.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/tagmatrix.py) | promoters, matrices and profiles |
| [`plotting.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/plotting.py) | annotation and overlap figures |
| [`qc.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/qc.py) | peak QC, FRiP and TSS enrichment |
| [`enrichment.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/enrichment.py) | Enrichr, ORA and GSEA |
| [`report.py`](https://github.com/summus-kong/pychipseeker/blob/main/src/pychipseeker/report.py) | multi-page PDF report |

Python docstrings are available through `help(pcs.annotate_peak)` or
`pcs.annotate_peak?` in IPython. For the full parameter tables, see the
[function reference](reference/functions.md).
