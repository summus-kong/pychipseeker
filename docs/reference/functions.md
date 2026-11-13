# Function reference

Every public top-level function, with every parameter. The tables are
generated from the installed signatures, so required parameters, defaults
and type hints always match the code.

!!! tip "How to read this page"

    - **Required** parameters have no default (`—`).
    - `None` usually means *automatic* or *not supplied*.
    - `n_jobs=None`/`1` is sequential; `-1` uses all cores.
    - Peak coordinates are 0-based and half-open, as in BED.
    - `TxDb=`, `gtf=` and `annotation=` are alternative ways to supply the
      same gene model; use one of them.

## Parameter conventions

| Name | Meaning |
|---|---|
| `peaks` / `peak` | A BED/narrowPeak-style path, a peak DataFrame, or (where noted) a sequence. |
| `TxDb` / `gtf` / `annotation` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `level` | `transcript` (default) or `gene`; controls which TSS set is used. |
| `by` | `gene` or `transcript`; controls the window rows. |
| `type` | `start_site`, `end_site` or `body`. |
| `upstream` / `downstream` | Bases on each side of the anchor; negative upstream. |
| `nbin` / `n_bins` | Number of bins across a window. |
| `value_type` | bigWig summary: `sum`, `mean`, `max`, `min`, `coverage` or `std`. |
| `n_jobs` | `None`/`1` sequential, `-1` all cores, `n > 1` explicit workers. |

## Aliases

ChIPseeker-style aliases share the same signature and parameters as their
snake_case counterpart:

| Alias | Canonical function |
|---|---|
| `annotatePeak` | [`annotate_peak`](#annotate_peak) |
| `asGRanges` | [`as_granges`](#as_granges) |
| `combine_csAnno` | [`combine_cs_anno`](#combine_cs_anno) |
| `detailGenomicAnnotation` | [`detail_genomic_annotation`](#detail_genomic_annotation) |
| `downloadGEObedFiles` | [`download_geo_bed_files`](#download_geo_bed_files) |
| `downloadGSMbedFiles` | [`download_gsm_bed_files`](#download_gsm_bed_files) |
| `dropAnno` | [`drop_anno`](#drop_anno) |
| `enrichAnnoOverlap` | [`enrich_anno_overlap`](#enrich_anno_overlap) |
| `enrichPeakOverlap` | [`enrich_peak_overlap`](#enrich_peak_overlap) |
| `getAnnoStat` | [`get_anno_stat`](#get_anno_stat) |
| `getBioRegion` | [`get_bio_region`](#get_bio_region) |
| `getGEOgenomeVersion` | [`get_geo_genome_version`](#get_geo_genome_version) |
| `getGEOInfo` | [`get_geo_info`](#get_geo_info) |
| `getGEOspecies` | [`get_geo_species`](#get_geo_species) |
| `getPromoters` | [`get_promoters`](#get_promoters) |
| `getSampleFiles` | [`get_sample_files`](#get_sample_files) |
| `getTagMatrix` | [`get_tag_matrix`](#get_tag_matrix) |
| `makeBioRegionFromGranges` | [`make_bio_region_from_granges`](#make_bio_region_from_granges) |
| `peakHeatmap` | [`peak_heatmap`](#peak_heatmap) |
| `peakHeatmap_multiple_Sets` | [`peak_heatmap_multiple_sets`](#peak_heatmap_multiple_sets) |
| `peak_Profile_Heatmap` | [`peak_profile_heatmap`](#peak_profile_heatmap) |
| `plotAnnoBar` | [`plot_anno_bar`](#plot_anno_bar) |
| `plotAnnoPie` | [`plot_anno_pie`](#plot_anno_pie) |
| `plotAvgProf` | [`plot_avg_prof`](#plot_avg_prof) |
| `plotDistToTSS` | [`plot_dist_to_tss`](#plot_dist_to_tss) |
| `plotMultiProf` | [`plot_multi_prof`](#plot_multi_prof) |
| `plotPeakProf` | [`plot_peak_prof`](#plot_peak_prof) |
| `plotAvgProf2` | [`plot_peak_prof2`](#plot_peak_prof2) |
| `plotPeakProf2` | [`plot_peak_prof2`](#plot_peak_prof2) |
| `plotPeakProf_MultiWindows` | [`plot_peak_prof_multi_windows`](#plot_peak_prof_multi_windows) |
| `plotHeatmap` | [`plot_tss_heatmap`](#plot_tss_heatmap) |
| `plotProfile` | [`plot_tss_profile`](#plot_tss_profile) |
| `readPeakFile` | [`read_peak_file`](#read_peak_file) |
| `tagHeatmap` | [`tag_heatmap`](#tag_heatmap) |

## Annotation and result tables

Annotate peaks, read/write tables, filter peaks and work with ChIPseeker-style result objects.

### `annotate_peak`

Annotate peaks with ChIPSeeker-compatible genomic features.

!!! note
    `annotatePeak` is a bare alias. R-style keyword spellings such as `tssRegion`, `genomicAnnotationPriority`, `peak`, `annotation`, `annoDb` and `columns` are accepted through `**kwargs`.

```python
annotate_peak(peaks: 'Optional[Union[str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame, GeneAnnotation]]' = None, tss_region: 'Sequence[int]' = (-3000, 3000), level: 'str' = 'transcript', assign_genomic_annotation: 'bool' = True, genomic_annotation_priority: 'Sequence[str]' = ('Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstream', 'Intergenic'), add_flank_gene_info: 'bool' = False, flank_distance: 'int' = 5000, same_strand: 'bool' = False, ignore_overlap: 'bool' = False, ignore_upstream: 'bool' = False, ignore_downstream: 'bool' = False, overlap: 'str' = 'TSS', verbose: 'bool' = False, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, promoter_bins: 'Optional[Sequence[int]]' = None, downstream_distance: 'int' = 300, drop_unmapped: 'bool' = True, align_annotation_feature: 'bool' = True, n_jobs: 'Optional[int]' = None, TxDb: 'Optional[Union[str, Path, pd.DataFrame, GeneAnnotation]]' = None, **kwargs) -> 'PeakAnnotation'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `peaks` | Optional[Union[str, Path, DataFrame]] | `None` | Peak file or DataFrame to annotate. Required unless supplied through a compatibility alias. |
| `gtf` | Optional[Union[str, Path, DataFrame, GeneAnnotation]] | `None` | Gene model to annotate against; mutually exclusive with `TxDb`. |
| `tss_region` | Sequence[int] | `(-3000, 3000)` | Promoter window `(upstream, downstream)` in bp; negative upstream, positive downstream. |
| `level` | str | `transcript` | `transcript` (default) or `gene`; controls which feature set supplies TSSs. |
| `assign_genomic_annotation` | bool | `True` | Assign promoter/UTR/exon/intron labels when `True`. |
| `genomic_annotation_priority` | Sequence[str] | `('Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstrea...` | Feature priority used to resolve overlapping annotations. |
| `add_flank_gene_info` | bool | `False` | Add genes within `flank_distance` of each peak. |
| `flank_distance` | int | `5000` | Distance in bp used to collect flanking genes. |
| `same_strand` | bool | `False` | Only consider features on the same strand as the peak. |
| `ignore_overlap` | bool | `False` | Skip overlapping-feature assignment when `True`. |
| `ignore_upstream` | bool | `False` | Ignore upstream features when assigning annotations. |
| `ignore_downstream` | bool | `False` | Ignore downstream features when assigning annotations. |
| `overlap` | str | `TSS` | `TSS` assigns features only near the TSS; `all` considers every overlapping feature span. |
| `verbose` | bool | `False` | Print progress messages when `True`. |
| `upstream` | Optional[int] | `None` | Convenience alias for the upstream side of `tss_region`. |
| `downstream` | Optional[int] | `None` | Convenience alias for the downstream side of `tss_region`. |
| `promoter_bins` | Optional[Sequence[int]] | `None` | Promoter bin edges in bp; defaults to `tss_region` plus ChIPseeker's standard bins. |
| `downstream_distance` | int | `300` | Downstream window in bp used for the `Downstream` category. |
| `drop_unmapped` | bool | `True` | Drop peaks without a genomic annotation when `True`. |
| `align_annotation_feature` | bool | `True` | Report the isoform that supplied the genomic feature (ChIPseeker fix #252). |
| `n_jobs` | Optional[int] | `None` | `None`/`1` sequential, `-1` all cores, or an explicit worker count. |
| `TxDb` | Optional[Union[str, Path, DataFrame, GeneAnnotation]] | `None` | Assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`; mutually exclusive with `gtf`. |
| `kwargs` | — | keyword | R-style aliases such as `tssRegion`, `genomicAnnotationPriority`, `peak`, `annotation`, `annoDb` and `columns`. |

### `annotate_db`

Add OrgDb-like annotation columns to a `PeakAnnotation`.

```python
annotate_db(annotation: 'PeakAnnotation', anno_db, columns: 'Optional[Sequence[str]]' = None) -> 'PeakAnnotation'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | PeakAnnotation | — | Annotated result whose gene columns should be extended. |
| **anno_db** | — | — | Gene annotation table/mapping/path keyed by `gene_id`. |
| `columns` | Optional[Sequence[str]] | `None` | Gene metadata columns to add; `None` keeps all available columns. |

### `seq2gene`

Return unique gene IDs, following `ChIPseeker::seq2gene`.

```python
seq2gene(peaks: 'Union[str, Path, pd.DataFrame]', gene_annotation: 'Optional[GeneAnnotation]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, tss_region: 'Sequence[int]' = (-3000, 3000), flank_distance: 'int' = 5000, same_strand: 'bool' = False, level: 'str' = 'gene', upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None) -> 'list[str]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peak file path, peak DataFrame, or a sequence of peak sources. |
| `gene_annotation` | Optional[GeneAnnotation] | `None` | Prebuilt `GeneAnnotation`; mutually exclusive with `gtf`. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | Gene model used to link peaks to genes; mutually exclusive with `gene_annotation`. |
| `tss_region` | Sequence[int] | `(-3000, 3000)` | Promoter window used for the link. |
| `flank_distance` | int | `5000` | Flank window used when a peak is not in a promoter. |
| `same_strand` | bool | `False` | Only consider features on the same strand as the peak. |
| `level` | str | `gene` | Annotation level: `'transcript'` (default) or `'gene'`. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |

### `seq2gene_links`

Many-to-many peak-to-gene links used for enrichment analysis.

```python
seq2gene_links(peaks: 'Union[str, Path, pd.DataFrame]', gene_annotation: 'Optional[GeneAnnotation]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, tss_region: 'Sequence[int]' = (-3000, 3000), flank_distance: 'int' = 5000, same_strand: 'bool' = False, level: 'str' = 'gene', upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peak file path, peak DataFrame, or a sequence of peak sources. |
| `gene_annotation` | Optional[GeneAnnotation] | `None` | `GeneAnnotation` used to build promoters/links. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | Gene model used to build the peak-to-gene link table. |
| `tss_region` | Sequence[int] | `(-3000, 3000)` | Promoter window as `(upstream, downstream)` relative to the TSS. |
| `flank_distance` | int | `5000` | Distance in bp used to collect flanking genes. |
| `same_strand` | bool | `False` | Only consider features on the same strand as the peak. |
| `level` | str | `gene` | Annotation level: `'transcript'` (default) or `'gene'`. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |

### `get_flanking_genes`

Return ChIPseeker-style semicolon-delimited flanking feature info.

```python
get_flanking_genes(model: 'GeneAnnotation', peaks: 'Union[str, Path, pd.DataFrame]', level: 'str' = 'transcript', distance: 'int' = 5000) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **model** | GeneAnnotation | — | `GeneAnnotation` that supplies the gene models. |
| **peaks** | Union[str, Path, DataFrame] | — | Peaks whose flanking genes should be found. |
| `level` | str | `transcript` | Annotation level: `'transcript'` (default) or `'gene'`. |
| `distance` | int | `5000` | Distance in bp used for flanking-gene searches. |

### `get_genomic_anno_stat`

Replicate `ChIPseeker:::getGenomicAnnoStat`.

```python
get_genomic_anno_stat(annotations: 'Sequence[str]', downstream_distance: 'int' = 300, ignore_promoter_subcategory: 'bool' = False, ignore_1st_exon: 'bool' = False, ignore_1st_intron: 'bool' = False, ignore_downstream: 'bool' = False) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotations** | Sequence[str] | — | Annotation labels, usually `ann.detail['annotation']`. |
| `downstream_distance` | int | `300` | Downstream distance used by the `Downstream` label. |
| `ignore_promoter_subcategory` | bool | `False` | Merge promoter sub-bins into one `Promoter` row. |
| `ignore_1st_exon` | bool | `False` | Merge first-exon categories into `Exon`. |
| `ignore_1st_intron` | bool | `False` | Merge first-intron categories into `Intron`. |
| `ignore_downstream` | bool | `False` | Ignore downstream features when assigning annotations. |

### `promoter_label`

Return a ChIPseeker promoter subcategory label.

```python
promoter_label(distance: 'float', bins: 'Sequence[int]') -> 'str'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **distance** | float | — | Signed distance from the peak to the TSS. |
| **bins** | Sequence[int] | — | Promoter bin edges in bp. |

### `group_annotation`

Collapse detailed labels into major ChIPseeker classes.

```python
group_annotation(annotation: 'str') -> 'str'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | str | — | Detailed ChIPseeker label to collapse to a major class. |

### `read_gtf`

Read a GTF/GFF file and extract common gene/transcript attributes.

```python
read_gtf(path: 'Union[str, Path]', *, features: 'Optional[Union[str, Iterable[str]]]' = None, chunksize: 'Optional[int]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **path** | Union[str, Path] | — | Input or output path. |
| `features` | Optional[Union[str, Iterable[str]]] | `None` | Feature types to keep; `None` reads all records. |
| `chunksize` | Optional[int] | `None` | Read the GTF in chunks of this many rows; `None` reads at once. |

### `read_peaks`

Read peaks from BED, narrowPeak/broadPeak/gappedPeak, bigBed, BAM,

```python
read_peaks(path: 'Union[str, Path]', fmt: 'Optional[str]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **path** | Union[str, Path] | — | BED, narrowPeak/broadPeak/gappedPeak, bigBed, BAM, TSV or CSV file. |
| `fmt` | Optional[str] | `None` | Input format override; inferred from the suffix when omitted. |

### `read_peak_file`

Read a peak file; `as_type='GRanges'` returns a GRanges-like table.

```python
read_peak_file(peak: 'Union[str, Path, pd.DataFrame]', as_type: 'str' = 'dataframe', **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame] | — | Peak file path or peak DataFrame. |
| `as_type` | str | `dataframe` | `dataframe` (default) or `GRanges`-like table. |
| `kwargs` | — | keyword | Forwarded to `read_peaks` and the format reader. |

### `write_table`

Write a DataFrame to TSV/CSV based on the filename suffix.

```python
write_table(df: 'pd.DataFrame', path: 'Union[str, Path]') -> 'None'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **df** | DataFrame | — | Table to write as TSV. |
| **path** | Union[str, Path] | — | Destination TSV path. |

### `write_peaks`

Write a peak table as BED6+ (`chrom start end name score strand`).

```python
write_peaks(peaks: 'pd.DataFrame', path: 'Union[str, Path]') -> 'None'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | DataFrame | — | Peak table to write as BED6+. |
| **path** | Union[str, Path] | — | Destination BED path. |

### `filter_peaks`

Filter a peak table by blacklist, chromosome, width and score.

```python
filter_peaks(peaks: 'Union[str, Path, pd.DataFrame]', *, blacklist: 'Optional[Union[str, Path, pd.DataFrame]]' = None, drop_mito: 'bool' = False, min_width: 'Optional[int]' = None, max_width: 'Optional[int]' = None, min_score: 'Optional[float]' = None, chroms: 'Optional[Sequence[str]]' = None, exclude_chroms: 'Optional[Sequence[str]]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peak file path, peak DataFrame, or a sequence of peak sources. |
| `blacklist` | Optional[Union[str, Path, DataFrame]] | `None` | Blacklist BED path/DataFrame; any overlap removes the peak. |
| `drop_mito` | bool | `False` | Remove `chrM`/`MT` peaks. |
| `min_width` | Optional[int] | `None` | Keep peaks at least this wide. |
| `max_width` | Optional[int] | `None` | Keep peaks at most this wide. |
| `min_score` | Optional[float] | `None` | Keep peaks with a score at or above this value. |
| `chroms` | Optional[Sequence[str]] | `None` | Whitelist of chromosomes to keep. |
| `exclude_chroms` | Optional[Sequence[str]] | `None` | Chromosomes to drop after the whitelist step. |

### `call_peaks`

Call peaks as runs of signal above `threshold`.

!!! note
    Returns BED6+ columns with `signalValue` (mean signal) and `peak` (offset of the highest base).

```python
call_peaks(bigwig: 'Union[str, Path]', *, threshold: 'float' = 3.0, max_gap: 'int' = 1000, min_width: 'int' = 500, chroms: 'Optional[Sequence[str]]' = None, sample: 'Optional[int]' = None, seed: 'Optional[int]' = None, n_jobs: 'Optional[int]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **bigwig** | Union[str, Path] | — | Signal track to threshold. |
| `threshold` | float | `3.0` | Minimum signal for a base to belong to a peak. |
| `max_gap` | int | `1000` | Merge runs separated by at most this many bases. |
| `min_width` | int | `500` | Discard merged runs shorter than this width. |
| `chroms` | Optional[Sequence[str]] | `None` | Chromosomes to keep or process. |
| `sample` | Optional[int] | `None` | Optional downsampling of the signal track for quick tests. |
| `seed` | Optional[int] | `None` | Random seed for reproducible shuffling, sampling or permutations. |
| `n_jobs` | Optional[int] | `None` | Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit. |

### `call_peaks_macs`

Run MACS2/MACS3 (if installed) and return the called narrowPeak table.

```python
call_peaks_macs(treatment: 'Union[str, Path]', *, control: 'Optional[Union[str, Path]]' = None, genome: 'Optional[str]' = None, outdir: 'Union[str, Path]' = 'macs_output', backend: 'str' = 'macs2', extra_args: 'Iterable[str]' = (), name: 'str' = 'pychipseeker') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **treatment** | Union[str, Path] | — | Treatment BAM/CRAM file. |
| `control` | Optional[Union[str, Path]] | `None` | Optional control BAM/CRAM file. |
| `genome` | Optional[str] | `None` | Genome name understood by MACS. |
| `outdir` | Union[str, Path] | `macs_output` | Output directory. |
| `backend` | str | `macs2` | `macs2` or `macs3`. |
| `extra_args` | Iterable[str] | `()` | Extra MACS command-line arguments. |
| `name` | str | `pychipseeker` | Sample name passed to MACS. |

### `as_granges`

Return a GRanges-like table for a `PeakAnnotation` or DataFrame.

```python
as_granges(annotation) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | Peak annotation table or `PeakAnnotation` object. |

### `combine_cs_anno`

Combine several `PeakAnnotation` objects (ChIPseeker `combine_csAnno`).

```python
combine_cs_anno(*annotations) -> 'PeakAnnotation'
```

### `drop_anno`

Filter peaks by `abs(distanceToTSS) < cutoff`.

```python
drop_anno(annotation, distance_to_tss_cutoff: 'float' = 10000)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table to filter. |
| `distance_to_tss_cutoff` | float | `10000` | Keep peaks with `abs(distance_to_tss) < cutoff`. |

### `get_anno_stat`

ChIPseeker-style `getAnnoStat` for annotation objects/tables.

```python
get_anno_stat(annotation) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |

### `detail_genomic_annotation`

Return the detailed per-peak genomic annotation table.

```python
detail_genomic_annotation(annotation) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |

## Analysis summaries

Plot-ready summaries for annotated peak tables and tag matrices.

### `annotation_summary`

Summarise peak count, width and signal for each annotation category.

```python
annotation_summary(annotation: 'Union[PeakAnnotation, pd.DataFrame]') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | Union[PeakAnnotation, DataFrame] | — | `PeakAnnotation` or a DataFrame with `chrom`, `start`, `end` and `annotation`. |

### `peaks_per_gene`

Return one row per gene with its peak count and peak-level summaries.

```python
peaks_per_gene(annotation: 'Union[PeakAnnotation, pd.DataFrame]') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | Union[PeakAnnotation, DataFrame] | — | `PeakAnnotation` or annotation table with gene columns. |

### `top_genes`

Return the `n` genes carrying the most annotated peaks.

```python
top_genes(annotation: 'Union[PeakAnnotation, pd.DataFrame]', n: 'int' = 20) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | Union[PeakAnnotation, DataFrame] | — | Peak annotation table or `PeakAnnotation` object. |
| `n` | int | `20` | Number of top genes to return. |

### `gene_type_summary`

Summarise peaks by gene biotype, folding the tail into `Other`.

```python
gene_type_summary(annotation: 'Union[PeakAnnotation, pd.DataFrame]', top: 'int' = 8) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | Union[PeakAnnotation, DataFrame] | — | Peak annotation table or `PeakAnnotation` object. |
| `top` | int | `8` | Number of biotypes kept before folding the rest into `Other`. |

### `distance_summary`

Return distance-to-TSS metrics as a two-column `Metric`/`Value` table.

```python
distance_summary(annotation: 'Union[PeakAnnotation, pd.DataFrame]', thresholds: 'Sequence[int]' = (1000, 10000, 100000)) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | Union[PeakAnnotation, DataFrame] | — | Peak annotation table or `PeakAnnotation` object. |
| `thresholds` | Sequence[int] | `(1000, 10000, 100000)` | Distance cutoffs in bp used for the percentage rows. |

### `signal_groups`

Assign tag-matrix rows to balanced groups by peak signal.

!!! note
    Genes without an assigned peak receive score 0, so they stay in the lowest group instead of being filled with the median.

```python
signal_groups(tag, annotation: 'Union[PeakAnnotation, pd.DataFrame]', n_groups: 'int' = 4) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **tag** | — | — | Tag matrix whose rows are assigned to groups. |
| **annotation** | Union[PeakAnnotation, DataFrame] | — | Peak annotation table or `PeakAnnotation` object. |
| `n_groups` | int | `4` | Number of balanced groups; capped at the number of rows. |

### `tss_enrichment_by_quartile`

Summarise TSS enrichment for signal-ranked groups of tag-matrix rows.

!!! note
    Enrichment is mean coverage in the central bins divided by mean coverage in the outer flanks. The centre is derived from the tag-matrix metadata when available.

```python
tss_enrichment_by_quartile(tag, annotation: 'Union[PeakAnnotation, pd.DataFrame]', n_groups: 'int' = 4, *, center_fraction: 'float' = 0.05, flank_fraction: 'float' = 0.1) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **tag** | — | — | Tag matrix used to compute centre/flank coverage. |
| **annotation** | Union[PeakAnnotation, DataFrame] | — | Peak annotation table or `PeakAnnotation` object. |
| `n_groups` | int | `4` | Number of signal groups. |
| `center_fraction` | float | `0.05` | Fraction of bins in the central window (default 5 %). |
| `flank_fraction` | float | `0.1` | Fraction of bins in the outer flanks (default 10 %). |

## Tag matrices and profiles

Promoters, bioregions, tag matrices, profiles and heatmaps.

### `get_promoters`

ChIPseeker-compatible promoter regions (`start_site` bioregions).

```python
get_promoters(TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, upstream: 'int' = 1000, downstream: 'int' = 1000, by: 'str' = 'gene', *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model used to build promoters. |
| `upstream` | int | `1000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `1000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `by` | str | `gene` | `gene` (default) or `transcript`. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Explicit `GeneAnnotation`, GTF path or DataFrame. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path used when `TxDb`/`annotation` is omitted. |

### `get_bio_region`

Return ChIPseeker-style promoter/body regions from a GTF annotation.

```python
get_bio_region(TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, upstream: 'int' = 1000, downstream: 'int' = 1000, by: 'str' = 'gene', type: 'str' = 'start_site', *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `upstream` | int | `1000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `1000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `by` | str | `gene` | `gene` or `transcript`. |
| `type` | str | `start_site` | `start_site`, `end_site` or `body`. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |

### `make_bio_region_from_granges`

Make a bioregion from an arbitrary interval table.

```python
make_bio_region_from_granges(gr: 'pd.DataFrame', by: 'str', type: 'str', upstream: 'int' = 1000, downstream: 'int' = 1000) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gr** | DataFrame | — | Interval table with `chrom`, `start`, `end` and optionally `strand`/`tss`. |
| **by** | str | — | Feature level label stored in the output. |
| **type** | str | — | `start_site`, `end_site` or `body`. |
| `upstream` | int | `1000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `1000` | Bases downstream of the anchor (the 3' side on the feature strand). |

### `get_tag_matrix`

Build a ChIPseeker-style tag matrix.

!!! note
    When `windows` is a sequence, a list of `TagMatrix` objects is returned. With `bigwig=`/`bam=` the matrix contains signal coverage; otherwise it contains peak coverage.

```python
get_tag_matrix(peak: 'Optional[Union[str, Path, pd.DataFrame]]' = None, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, windows: 'Optional[Union[pd.DataFrame, str, Path]]' = None, type: 'str' = 'start_site', by: 'str' = 'gene', TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, weightCol: 'Optional[str]' = None, nbin: 'Optional[int]' = None, verbose: 'bool' = False, ignore_strand: 'bool' = False, *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, bigwig: 'Optional[Union[str, Path]]' = None, bam: 'Optional[Union[str, Path]]' = None, value_type: 'str' = 'sum', min_mapping_quality: 'int' = 0, exclude_duplicates: 'bool' = True, n_jobs: 'Optional[int]' = None) -> 'Union[TagMatrix, list[TagMatrix]]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `peak` | Optional[Union[str, Path, DataFrame]] | `None` | Peak file/DataFrame, or a bigWig path when `bigwig` is omitted. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `windows` | Optional[Union[DataFrame, str, Path]] | `None` | Pre-built windows; a sequence returns a list of matrices. |
| `type` | str | `start_site` | Region type used when building windows. |
| `by` | str | `gene` | `gene` or `transcript`. |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model used when `windows` is omitted. |
| `weightCol` | Optional[str] | `None` | Peak column used as weight instead of 1 per peak. |
| `nbin` | Optional[int] | `None` | Bins per window. |
| `verbose` | bool | `False` | Print progress messages when `True`. |
| `ignore_strand` | bool | `False` | Ignore window strand. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |
| `bigwig` | Optional[Union[str, Path]] | `None` | Use signal coverage instead of peak coverage. |
| `bam` | Optional[Union[str, Path]] | `None` | Use BAM/CRAM coverage instead of peak coverage. |
| `value_type` | str | `sum` | bigWig summary statistic. |
| `min_mapping_quality` | int | `0` | Minimum MAPQ for BAM coverage. |
| `exclude_duplicates` | bool | `True` | Ignore duplicate BAM reads. |
| `n_jobs` | Optional[int] | `None` | Parallel worker count for bigWig/BAM extraction. |

### `get_tag_matrix_from_bigwig`

Convenience wrapper for a bigWig-backed :class:`TagMatrix`.

```python
get_tag_matrix_from_bigwig(bigwig: 'Union[str, Path]', windows: 'Union[pd.DataFrame, str, Path]', nbin: 'Optional[int]' = None, ignore_strand: 'bool' = False, value_type: 'str' = 'sum', n_jobs: 'Optional[int]' = None) -> 'TagMatrix'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **bigwig** | Union[str, Path] | — | bigWig signal file; requires the optional `pyBigWig` dependency. |
| **windows** | Union[DataFrame, str, Path] | — | Window table (one row per region). |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `value_type` | str | `sum` | `sum`, `mean`, `max`, `min`, `coverage` or `std`. |
| `n_jobs` | Optional[int] | `None` | Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit. |

### `get_tag_matrix_from_bam`

Convenience wrapper for a BAM-backed :class:`TagMatrix`.

```python
get_tag_matrix_from_bam(bam: 'Union[str, Path]', windows: 'Union[pd.DataFrame, str, Path]', nbin: 'Optional[int]' = None, ignore_strand: 'bool' = False, value_type: 'str' = 'mean', min_mapping_quality: 'int' = 0, exclude_duplicates: 'bool' = True, n_jobs: 'Optional[int]' = None) -> 'TagMatrix'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **bam** | Union[str, Path] | — | BAM/CRAM alignment file; requires the optional `pysam` dependency. |
| **windows** | Union[DataFrame, str, Path] | — | Window table (one row per region). |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `value_type` | str | `mean` | `mean` (default) or another supported summary. |
| `min_mapping_quality` | int | `0` | Minimum BAM mapping quality (MAPQ) for covered bases. |
| `exclude_duplicates` | bool | `True` | Ignore duplicate reads/alignments when `True`. |
| `n_jobs` | Optional[int] | `None` | Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit. |

### `extract_tss_matrix`

Extract a features x bins TSS-centered signal matrix from a bigWig file.

!!! note
    `upstream`/`downstream` define an asymmetric window; otherwise `window` is used on both sides.

```python
extract_tss_matrix(bigwig: 'Union[str, Path]', genes: 'pd.DataFrame', window: 'int' = 2000, n_bins: 'int' = 100, min_value: 'float' = 0.0, *, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, drop_empty: 'bool' = False, return_mask: 'bool' = False, dtype: 'type' = <class 'numpy.float64'>, n_jobs: 'Optional[int]' = None) -> 'Union[np.ndarray, tuple[np.ndarray, np.ndarray]]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **bigwig** | Union[str, Path] | — | bigWig signal file; requires the optional `pyBigWig` dependency. |
| **genes** | DataFrame | — | TSS feature table; use `ann.gene_annotation.features('transcript')`. |
| `window` | int | `2000` | Symmetric window size in bp when `upstream`/`downstream` are omitted. |
| `n_bins` | int | `100` | Number of bins across the window. |
| `min_value` | float | `0.0` | Values below this are treated as empty. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `drop_empty` | bool | `False` | Drop features with no signal. |
| `return_mask` | bool | `False` | Return `(matrix, mask)` instead of only the matrix. |
| `dtype` | type | `<class 'numpy.float64'>` | Output dtype (default `float64`). |
| `n_jobs` | Optional[int] | `None` | Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit. |

### `extract_tss_matrix_from_signal`

Extract a TSS-centered matrix from a bedGraph-like signal table.

```python
extract_tss_matrix_from_signal(signal: 'Union[str, Path, pd.DataFrame]', genes: 'pd.DataFrame', window: 'int' = 2000, n_bins: 'int' = 100, min_value: 'Optional[float]' = 0.0, chrom_sizes: 'Optional[dict[str, int]]' = None, *, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, drop_empty: 'bool' = False, return_mask: 'bool' = False, dtype: 'type' = <class 'numpy.float64'>) -> 'Union[np.ndarray, tuple[np.ndarray, np.ndarray]]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **signal** | Union[str, Path, DataFrame] | — | bigWig path or precomputed signal table. |
| **genes** | DataFrame | — | Gene or transcript feature table; normally `ann.gene_annotation.features(...)`. |
| `window` | int | `2000` | Half-window size in bp, applied to both sides unless `upstream`/`downstream` are given. |
| `n_bins` | int | `100` | Number of bins across the window. |
| `min_value` | Optional[float] | `0.0` | Minimum signal value; lower values are clipped or treated as empty. |
| `chrom_sizes` | Optional[dict[str, int]] | `None` | Chromosome lengths used to clip windows. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `drop_empty` | bool | `False` | Drop rows with no signal when `True`. |
| `return_mask` | bool | `False` | Also return the boolean mask of retained rows. |
| `dtype` | type | `<class 'numpy.float64'>` | NumPy dtype of the returned matrix. |

### `plot_avg_prof`

Plot the mean tag-matrix profile with an optional confidence band.

```python
plot_avg_prof(tagMatrix, xlim: 'Optional[Sequence[int]]' = None, xlab: 'str' = "Genomic Region (5'->3')", ylab: 'str' = 'Peak Count Frequency', conf: 'Optional[float]' = 0.95, facet: 'str' = 'none', free_y: 'bool' = True, origin_label: 'str' = 'TSS', ax=None, figsize: 'tuple[float, float]' = (6.0, 4.0), color: 'Optional[str]' = None, label: 'Optional[str]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **tagMatrix** | — | — | Tag matrix to average. |
| `xlim` | Optional[Sequence[int]] | `None` | X-axis limits in bp relative to the anchor. |
| `xlab` | str | `Genomic Region (5'->3')` | X-axis label. |
| `ylab` | str | `Peak Count Frequency` | Y-axis label. |
| `conf` | Optional[float] | `0.95` | Bootstrap confidence level; `None` hides the band. |
| `facet` | str | `none` | `none`, `row` or `column`. |
| `free_y` | bool | `True` | Let each facet use its own y-axis range when `True`. |
| `origin_label` | str | `TSS` | Text drawn at the anchor (default `TSS`). |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `figsize` | tuple[float, float] | `(6.0, 4.0)` | Figure size in inches, as `(width, height)`. |
| `color` | Optional[str] | `None` | Line colour. |
| `label` | Optional[str] | `None` | Legend label. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_peak_prof2`

Build a tag matrix from peaks and plot its average profile.

```python
plot_peak_prof2(peak: 'Union[str, Path, pd.DataFrame]', upstream: 'int' = 3000, downstream: 'int' = 3000, conf: 'Optional[float]' = 0.95, by: 'str' = 'gene', type: 'str' = 'start_site', weightCol: 'Optional[str]' = None, TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, xlab: 'str' = "Genomic Region (5'->3')", ylab: 'str' = 'Peak Count Frequency', facet: 'str' = 'none', free_y: 'bool' = True, verbose: 'bool' = False, nbin: 'Optional[int]' = None, ignore_strand: 'bool' = False, *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame] | — | Peak file path or peak DataFrame. |
| `upstream` | int | `3000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `3000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `conf` | Optional[float] | `0.95` | Confidence level for the bootstrap band; `None` disables the band. |
| `by` | str | `gene` | Feature level used to build windows: `'gene'` or `'transcript'`. |
| `type` | str | `start_site` | Region type: `'start_site'`, `'end_site'` or `'body'`. |
| `weightCol` | Optional[str] | `None` | Optional peak column used as weight (for example `signalValue`). |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `xlab` | str | `Genomic Region (5'->3')` | X-axis label. |
| `ylab` | str | `Peak Count Frequency` | Y-axis label. |
| `facet` | str | `none` | Facet layout: `'none'`, `'row'` or `'column'`. |
| `free_y` | bool | `True` | Let each facet use its own y-axis range when `True`. |
| `verbose` | bool | `False` | Print progress messages when `True`. |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_peak_prof`

Plot a profile from a tag matrix, or build one from peaks.

```python
plot_peak_prof(tagMatrix=None, peak: 'Optional[Union[str, Path, pd.DataFrame]]' = None, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, conf: 'Optional[float]' = 0.95, by: 'str' = 'gene', type: 'str' = 'start_site', windows_name: 'Optional[Sequence[str]]' = None, weightCol: 'Optional[str]' = None, TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, xlab: 'str' = "Genomic Region (5'->3')", ylab: 'str' = 'Peak Count Frequency', facet: 'str' = 'row', free_y: 'bool' = True, verbose: 'bool' = False, nbin: 'Optional[int]' = None, ignore_strand: 'bool' = False, *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tagMatrix` | — | `None` | Optional precomputed matrix; otherwise built from `peak`. |
| `peak` | Optional[Union[str, Path, DataFrame]] | `None` | Peak file path or peak DataFrame. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `conf` | Optional[float] | `0.95` | Confidence level for the bootstrap band; `None` disables the band. |
| `by` | str | `gene` | Feature level used to build windows: `'gene'` or `'transcript'`. |
| `type` | str | `start_site` | Region type: `'start_site'`, `'end_site'` or `'body'`. |
| `windows_name` | Optional[Sequence[str]] | `None` | Names for multiple window sets. |
| `weightCol` | Optional[str] | `None` | Optional peak column used as weight (for example `signalValue`). |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `xlab` | str | `Genomic Region (5'->3')` | X-axis label. |
| `ylab` | str | `Peak Count Frequency` | Y-axis label. |
| `facet` | str | `row` | Facet layout: `'none'`, `'row'` or `'column'`. |
| `free_y` | bool | `True` | Let each facet use its own y-axis range when `True`. |
| `verbose` | bool | `False` | Print progress messages when `True`. |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_peak_prof_multi_windows`

Plot profiles for several bioregions/windows on the same peak set.

```python
plot_peak_prof_multi_windows(peak: 'Union[str, Path, pd.DataFrame]', windows: 'Sequence[pd.DataFrame]', windows_name: 'Optional[Sequence[str]]' = None, weightCol: 'Optional[str]' = None, nbin: 'Optional[int]' = None, conf: 'Optional[float]' = 0.95, ignore_strand: 'bool' = False, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame] | — | Peak file path or peak DataFrame. |
| **windows** | Sequence[DataFrame] | — | Sequence of window tables to compare. |
| `windows_name` | Optional[Sequence[str]] | `None` | Labels for the window sets. |
| `weightCol` | Optional[str] | `None` | Optional peak column used as weight (for example `signalValue`). |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `conf` | Optional[float] | `0.95` | Confidence level for the bootstrap band; `None` disables the band. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_multi_prof`

Plot several tag matrices, one line (or facet) per matrix.

```python
plot_multi_prof(tagMatrices: 'Sequence', xlab: 'str' = "Genomic Region (5'->3')", ylab: 'str' = 'Peak Count Frequency', conf: 'Optional[float]' = 0.95, facet: 'str' = 'row', free_y: 'bool' = True, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **tagMatrices** | Sequence | — | Sequence of matrices to compare. |
| `xlab` | str | `Genomic Region (5'->3')` | X-axis label. |
| `ylab` | str | `Peak Count Frequency` | Y-axis label. |
| `conf` | Optional[float] | `0.95` | Confidence level for the bootstrap band; `None` disables the band. |
| `facet` | str | `row` | Facet layout: `'none'`, `'row'` or `'column'`. |
| `free_y` | bool | `True` | Let each facet use its own y-axis range when `True`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `tag_heatmap`

Heatmap of a tag matrix, centered on the TSS/region anchor.

```python
tag_heatmap(tagMatrix, xlab: 'str' = "Genomic Region (5'->3')", ylab: 'str' = '', title: 'Optional[str]' = None, palette: 'Optional[str]' = None, ax=None, figsize: 'tuple[float, float]' = (6.0, 4.0), **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **tagMatrix** | — | — | Matrix to draw. |
| `xlab` | str | `Genomic Region (5'->3')` | X-axis label. |
| `ylab` | str | `` | Y-axis label. |
| `title` | Optional[str] | `None` | Figure or table title. |
| `palette` | Optional[str] | `None` | Colormap name. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `figsize` | tuple[float, float] | `(6.0, 4.0)` | Figure size in inches, as `(width, height)`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `peak_heatmap`

Build a tag matrix from peaks and draw a heatmap.

```python
peak_heatmap(peak: 'Union[str, Path, pd.DataFrame]', weightCol: 'Optional[str]' = None, TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, upstream: 'int' = 1000, downstream: 'int' = 1000, xlab: 'str' = '', ylab: 'str' = '', title: 'Optional[str]' = None, palette: 'Optional[str]' = None, verbose: 'bool' = False, by: 'str' = 'gene', type: 'str' = 'start_site', nbin: 'Optional[int]' = None, ignore_strand: 'bool' = False, windows: 'Optional[Union[pd.DataFrame, str, Path]]' = None, *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame] | — | Peaks used to build the matrix. |
| `weightCol` | Optional[str] | `None` | Optional peak column used as weight (for example `signalValue`). |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `upstream` | int | `1000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `1000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `xlab` | str | `` | X-axis label. |
| `ylab` | str | `` | Y-axis label. |
| `title` | Optional[str] | `None` | Figure or table title. |
| `palette` | Optional[str] | `None` | Matplotlib colormap or palette name. |
| `verbose` | bool | `False` | Print progress messages when `True`. |
| `by` | str | `gene` | Feature level used to build windows: `'gene'` or `'transcript'`. |
| `type` | str | `start_site` | Region type: `'start_site'`, `'end_site'` or `'body'`. |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `windows` | Optional[Union[DataFrame, str, Path]] | `None` | Optional pre-built windows. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `peak_heatmap_multiple_sets`

Heatmaps for several peak sets over the same bioregion.

```python
peak_heatmap_multiple_sets(peaks: 'Sequence[Union[str, Path, pd.DataFrame]]', labels: 'Optional[Sequence[str]]' = None, upstream: 'int' = 1000, downstream: 'int' = 1000, by: 'str' = 'gene', type: 'str' = 'start_site', TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, weightCol: 'Optional[str]' = None, nbin: 'Optional[int]' = None, ignore_strand: 'bool' = False, *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Sequence[Union[str, Path, DataFrame]] | — | Peak sets to compare; one heatmap row block per set. |
| `labels` | Optional[Sequence[str]] | `None` | Labels for the peak sets. |
| `upstream` | int | `1000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `1000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `by` | str | `gene` | Feature level used to build windows: `'gene'` or `'transcript'`. |
| `type` | str | `start_site` | Region type: `'start_site'`, `'end_site'` or `'body'`. |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `weightCol` | Optional[str] | `None` | Optional peak column used as weight (for example `signalValue`). |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `ignore_strand` | bool | `False` | Ignore feature strand when building windows or profiles. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `peak_profile_heatmap`

Combine average profile and heatmap in one figure.

```python
peak_profile_heatmap(peak: 'Union[str, Path, pd.DataFrame]', upstream: 'int' = 3000, downstream: 'int' = 3000, by: 'str' = 'gene', type: 'str' = 'start_site', TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, weightCol: 'Optional[str]' = None, nbin: 'Optional[int]' = None, conf: 'Optional[float]' = 0.95, *, annotation: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame] | — | Peaks used for the profile/heatmap pair. |
| `upstream` | int | `3000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `3000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `by` | str | `gene` | Feature level used to build windows: `'gene'` or `'transcript'`. |
| `type` | str | `start_site` | Region type: `'start_site'`, `'end_site'` or `'body'`. |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |
| `weightCol` | Optional[str] | `None` | Optional peak column used as weight (for example `signalValue`). |
| `nbin` | Optional[int] | `None` | Number of bins across each window; `None` uses the column widths when available. |
| `conf` | Optional[float] | `0.95` | Confidence level for the bootstrap band; `None` disables the band. |
| `annotation` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Peak annotation table or `PeakAnnotation` object. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | GTF/GFF path, feature DataFrame or `GeneAnnotation`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_tss_profile`

Plot the mean TSS-centered signal profile.

```python
plot_tss_profile(matrix: 'np.ndarray', window: 'int' = 2000, ax=None, title: 'str' = 'TSS profile', color: 'str' = '#E45756', figsize: 'tuple[float, float]' = (6.0, 4.0), show_sem: 'bool' = True, *, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **matrix** | ndarray | — | Features x bins signal matrix. |
| `window` | int | `2000` | Symmetric window size when `upstream`/`downstream` are omitted. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `title` | str | `TSS profile` | Figure or table title. |
| `color` | str | `#E45756` | Line or bar colour. |
| `figsize` | tuple[float, float] | `(6.0, 4.0)` | Figure size in inches, as `(width, height)`. |
| `show_sem` | bool | `True` | Draw the standard-error band. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_tss_heatmap`

Plot a heatmap of TSS-centered signal.

```python
plot_tss_heatmap(matrix: 'np.ndarray', window: 'int' = 2000, ax=None, title: 'str' = 'TSS heatmap', cmap: 'str' = 'viridis', figsize: 'tuple[float, float]' = (6.0, 5.0), sort_rows: 'bool' = False, *, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **matrix** | ndarray | — | Features x bins signal matrix. |
| `window` | int | `2000` | Half-window size in bp, applied to both sides unless `upstream`/`downstream` are given. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `title` | str | `TSS heatmap` | Figure or table title. |
| `cmap` | str | `viridis` | Matplotlib colormap name. |
| `figsize` | tuple[float, float] | `(6.0, 5.0)` | Figure size in inches, as `(width, height)`. |
| `sort_rows` | bool | `False` | Sort rows by mean signal. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

## Annotation plots

Pie, bar and distance-to-TSS figures for a PeakAnnotation.

### `plot_anno_pie`

Pie chart of detailed ChIPseeker annotation categories.

```python
plot_anno_pie(annotation, ax=None, title: 'str' = '', colors: 'Optional[Sequence[str]]' = None, legend: 'bool' = True, ndigit: 'int' = 2, figsize: 'tuple[float, float]' = (7.0, 4.5), **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `title` | str | `` | Figure or table title. |
| `colors` | Optional[Sequence[str]] | `None` | Sequence of colours used for the sets or bars. |
| `legend` | bool | `True` | Draw a legend when `True`. |
| `ndigit` | int | `2` | Number of digits used in labels. |
| `figsize` | tuple[float, float] | `(7.0, 4.5)` | Figure size in inches, as `(width, height)`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_anno_bar`

Horizontal stacked bar chart of detailed annotation categories.

```python
plot_anno_bar(annotation, groupby: 'Optional[str]' = None, ax=None, title: 'str' = 'Feature Distribution', xlab: 'str' = '', ylab: 'str' = 'Percentage(%)', figsize: 'tuple[float, float]' = (7.0, 4.5), colors: 'Optional[Sequence[str]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |
| `groupby` | Optional[str] | `None` | Optional metadata column used to split the summary. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `title` | str | `Feature Distribution` | Figure or table title. |
| `xlab` | str | `` | X-axis label. |
| `ylab` | str | `Percentage(%)` | Y-axis label. |
| `figsize` | tuple[float, float] | `(7.0, 4.5)` | Figure size in inches, as `(width, height)`. |
| `colors` | Optional[Sequence[str]] | `None` | Sequence of colours used for the sets or bars. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_dist_to_tss`

Diverging distance-to-TSS bar chart, following `plotDistToTSS`.

```python
plot_dist_to_tss(annotation, groupby: 'Optional[str]' = None, ax=None, distance_breaks: 'Optional[Sequence[float]]' = None, distance_column: 'str' = 'distance_to_tss', title: 'str' = 'Distribution of transcription factor-binding loci relative to TSS', xlab: 'str' = '', ylab: 'str' = "Binding sites (%) (5'->3')", figsize: 'tuple[float, float]' = (7.5, 4.5), colors: 'Optional[Sequence[str]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |
| `groupby` | Optional[str] | `None` | Optional metadata column used to split the summary. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `distance_breaks` | Optional[Sequence[float]] | `None` | Distance-axis break points in bp. |
| `distance_column` | str | `distance_to_tss` | Column that supplies TSS distances. |
| `title` | str | `Distribution of transcription factor-binding loci relative to TSS` | Figure or table title. |
| `xlab` | str | `` | X-axis label. |
| `ylab` | str | `Binding sites (%) (5'->3')` | Y-axis label. |
| `figsize` | tuple[float, float] | `(7.5, 4.5)` | Figure size in inches, as `(width, height)`. |
| `colors` | Optional[Sequence[str]] | `None` | Sequence of colours used for the sets or bars. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

## Quality control

Peak widths, blacklist/mitochondrial fractions, FRiP and TSS enrichment.

### `peak_qc`

Return a long-form table of peak quality metrics.

```python
peak_qc(peaks: 'Union[str, Path, pd.DataFrame]', *, blacklist: 'Optional[Union[str, Path, pd.DataFrame]]' = None, chrom_sizes: 'Optional[Dict[str, int]]' = None, genome: 'Optional[Union[str, Path]]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peak set to summarize. |
| `blacklist` | Optional[Union[str, Path, DataFrame]] | `None` | Blacklist intervals counted as flagged peaks. |
| `chrom_sizes` | Optional[Dict[str, int]] | `None` | Chromosome lengths for fraction calculations. |
| `genome` | Optional[Union[str, Path]] | `None` | Assembly name or chromosome-sizes source used for QC. |

### `frip`

Fraction of signal/reads falling inside peaks.

```python
frip(peaks: 'Union[str, Path, pd.DataFrame]', *, bigwig: 'Optional[Union[str, Path]]' = None, bam: 'Optional[Union[str, Path]]' = None, chrom_sizes: 'Optional[Dict[str, int]]' = None, exclude_duplicates: 'bool' = True, min_mapping_quality: 'int' = 0) -> 'float'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peaks used in the FRiP numerator. |
| `bigwig` | Optional[Union[str, Path]] | `None` | Signal track alternative to `bam`. |
| `bam` | Optional[Union[str, Path]] | `None` | Alignment file alternative to `bigwig`. |
| `chrom_sizes` | Optional[Dict[str, int]] | `None` | Mapping or DataFrame of chromosome lengths used for shuffling, FRiP or QC. |
| `exclude_duplicates` | bool | `True` | Ignore duplicate reads. |
| `min_mapping_quality` | int | `0` | Minimum BAM mapping quality (MAPQ) for covered bases. |

### `tss_enrichment`

Compute the ENCODE-style TSS enrichment score.

```python
tss_enrichment(genes: 'pd.DataFrame', *, bigwig: 'Optional[Union[str, Path]]' = None, bam: 'Optional[Union[str, Path]]' = None, upstream: 'int' = 2000, downstream: 'int' = 2000, n_bins: 'int' = 100, sample: 'Optional[int]' = 2000, seed: 'int' = 20261026, center_window: 'int' = 100, flank_window: 'Optional[int]' = None, n_jobs: 'Optional[int]' = None) -> 'dict[str, float]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **genes** | DataFrame | — | TSS feature table. |
| `bigwig` | Optional[Union[str, Path]] | `None` | bigWig signal file; requires the optional `pyBigWig` dependency. |
| `bam` | Optional[Union[str, Path]] | `None` | BAM/CRAM alignment file; requires the optional `pysam` dependency. |
| `upstream` | int | `2000` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | int | `2000` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `n_bins` | int | `100` | Number of bins across the window. |
| `sample` | Optional[int] | `2000` | Maximum number of TSSs sampled. |
| `seed` | int | `20261026` | Random seed for reproducible shuffling, sampling or permutations. |
| `center_window` | int | `100` | Half-width of the centre window. |
| `flank_window` | Optional[int] | `None` | Half-width of the flank window. |
| `n_jobs` | Optional[int] | `None` | Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit. |

### `qc_summary`

Combine peak metrics, FRiP and TSS enrichment in one long table.

```python
qc_summary(peaks: 'Union[str, Path, pd.DataFrame]', *, bigwig: 'Optional[Union[str, Path]]' = None, bam: 'Optional[Union[str, Path]]' = None, genes: 'Optional[pd.DataFrame]' = None, blacklist: 'Optional[Union[str, Path, pd.DataFrame]]' = None, chrom_sizes: 'Optional[Dict[str, int]]' = None, genome: 'Optional[Union[str, Path]]' = None, tss_upstream: 'int' = 2000, tss_downstream: 'int' = 2000, tss_sample: 'Optional[int]' = 2000, seed: 'int' = 20261026, exclude_duplicates: 'bool' = True, n_jobs: 'Optional[int]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peak file path, peak DataFrame, or a sequence of peak sources. |
| `bigwig` | Optional[Union[str, Path]] | `None` | bigWig signal file; requires the optional `pyBigWig` dependency. |
| `bam` | Optional[Union[str, Path]] | `None` | BAM/CRAM alignment file; requires the optional `pysam` dependency. |
| `genes` | Optional[DataFrame] | `None` | TSS table used for TSS enrichment. |
| `blacklist` | Optional[Union[str, Path, DataFrame]] | `None` | Blacklist BED path or DataFrame whose intervals are removed. |
| `chrom_sizes` | Optional[Dict[str, int]] | `None` | Mapping or DataFrame of chromosome lengths used for shuffling, FRiP or QC. |
| `genome` | Optional[Union[str, Path]] | `None` | Assembly name or chromosome-sizes source used for QC. |
| `tss_upstream` | int | `2000` | Upstream TSS window for QC. |
| `tss_downstream` | int | `2000` | Downstream TSS window for QC. |
| `tss_sample` | Optional[int] | `2000` | Maximum number of TSS features sampled. |
| `seed` | int | `20261026` | Random seed for reproducible shuffling, sampling or permutations. |
| `exclude_duplicates` | bool | `True` | Ignore duplicate reads in FRiP. |
| `n_jobs` | Optional[int] | `None` | Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit. |

## Enrichment

Local over-representation analysis, Enrichr, GSEA and enrichment plots.

### `read_gmt`

Read a GMT file into `{term: [genes]}`.

```python
read_gmt(path: 'Union[str, Path]') -> 'dict[str, list[str]]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **path** | Union[str, Path] | — | GMT file path. |

### `ora`

Local over-representation analysis with a hypergeometric test.

```python
ora(gene_list: 'Iterable[str]', gene_sets: 'Union[Mapping[str, Iterable[str]], str, Path]', background: 'Optional[Iterable[str]]' = None, p_adjust: 'str' = 'BH', min_overlap: 'int' = 1) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_list** | Iterable[str] | — | Query genes. |
| **gene_sets** | Union[Mapping[str, Iterable[str]], str, Path] | — | GMT path, library name or mapping of term to genes. |
| `background` | Optional[Iterable[str]] | `None` | Background genes; defaults to the union of all sets. |
| `p_adjust` | str | `BH` | `BH`, `Bonferroni` or `none`. |
| `min_overlap` | int | `1` | Minimum overlap between query and term. |

### `enrichr`

Run Enrichr through gseapy and return the result table.

```python
enrichr(gene_list: 'Iterable[str]', gene_sets: 'Union[str, Sequence[str]]' = 'GO_Biological_Process_2023', organism: 'str' = 'human', outdir: 'Optional[Union[str, Path]]' = None, background: 'Optional[Iterable[str]]' = None, p_adjust: 'Optional[str]' = None, **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_list** | Iterable[str] | — | Query genes. |
| `gene_sets` | Union[str, Sequence[str]] | `GO_Biological_Process_2023` | Enrichr library name or names. |
| `organism` | str | `human` | Enrichr organism. |
| `outdir` | Optional[Union[str, Path]] | `None` | Directory for gseapy output. |
| `background` | Optional[Iterable[str]] | `None` | Optional background genes. |
| `p_adjust` | Optional[str] | `None` | Optional adjustment method passed to gseapy. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `enrich_go`

GO Biological Process enrichment via Enrichr.

```python
enrich_go(gene_list: 'Iterable[str]', organism: 'str' = 'human', outdir: 'Optional[Union[str, Path]]' = None, **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_list** | Iterable[str] | — | Query genes. |
| `organism` | str | `human` | Enrichr organism. |
| `outdir` | Optional[Union[str, Path]] | `None` | Output directory. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `enrich_kegg`

KEGG pathway enrichment via Enrichr.

```python
enrich_kegg(gene_list: 'Iterable[str]', organism: 'str' = 'human', outdir: 'Optional[Union[str, Path]]' = None, **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_list** | Iterable[str] | — | Query genes. |
| `organism` | str | `human` | Enrichr organism. |
| `outdir` | Optional[Union[str, Path]] | `None` | Output directory. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `enrich_peak_genes`

Map peaks to genes with `seq2gene` and run Enrichr.

```python
enrich_peak_genes(annotation, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, include_gene_body: 'bool' = True, gene_column: 'str' = 'gene_name', **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |
| `upstream` | Optional[int] | `None` | Upstream promoter window for peak-to-gene links. |
| `downstream` | Optional[int] | `None` | Downstream promoter window for peak-to-gene links. |
| `include_gene_body` | bool | `True` | Include gene-body links. |
| `gene_column` | str | `gene_name` | Gene identifier column. |
| `kwargs` | — | keyword | Forwarded to `enrichr`. |

### `gsea`

Preranked Gene Set Enrichment Analysis through gseapy.

```python
gsea(gene_list: 'Union[Mapping[str, float], pd.DataFrame, Sequence[str]]', gene_sets: 'Union[str, Sequence[str], Mapping[str, Iterable[str]]]', *, min_size: 'int' = 15, max_size: 'int' = 500, permutation_num: 'int' = 1000, seed: 'int' = 20261026, outdir: 'Optional[Union[str, Path]]' = None, **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_list** | Union[Mapping[str, float], DataFrame, Sequence[str]] | — | Ranked mapping, two-column table or ranked gene sequence. |
| **gene_sets** | Union[str, Sequence[str], Mapping[str, Iterable[str]]] | — | Library name, GMT path or mapping of term to genes. |
| `min_size` | int | `15` | Minimum gene-set size. |
| `max_size` | int | `500` | Maximum gene-set size. |
| `permutation_num` | int | `1000` | Number of permutations. |
| `seed` | int | `20261026` | Random seed. |
| `outdir` | Optional[Union[str, Path]] | `None` | Output directory. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `gsea_peak_genes`

Restrict a ranked gene list to peak-linked genes and run GSEA.

```python
gsea_peak_genes(annotation, ranking: 'Union[Mapping[str, float], pd.DataFrame, Sequence[str]]', *, gene_column: 'str' = 'gene_name', upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | — | — | `PeakAnnotation` or annotation table. |
| **ranking** | Union[Mapping[str, float], DataFrame, Sequence[str]] | — | Gene ranking used for GSEA. |
| `gene_column` | str | `gene_name` | Gene identifier column. |
| `upstream` | Optional[int] | `None` | Bases upstream of the anchor (the 5' side on the feature strand). |
| `downstream` | Optional[int] | `None` | Bases downstream of the anchor (the 3' side on the feature strand). |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `plot_enrichment`

Plot an ORA/Enrichr/GSEA result table.

```python
plot_enrichment(result: 'pd.DataFrame', top: 'int' = 20, *, kind: 'str' = 'bar', ax=None, title: 'Optional[str]' = None, figsize: 'tuple[float, float]' = (7.0, 5.0), cmap: 'str' = 'Blues', **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **result** | DataFrame | — | ORA/Enrichr/GSEA result table. |
| `top` | int | `20` | Number of terms to show. |
| `kind` | str | `bar` | `bar` or `dot`. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |
| `title` | Optional[str] | `None` | Figure or table title. |
| `figsize` | tuple[float, float] | `(7.0, 5.0)` | Figure size in inches, as `(width, height)`. |
| `cmap` | str | `Blues` | Matplotlib colormap name. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

## Genomes and cache

Assembly detection, genome bundles, GTFs and the local cache.

### `canonical_chrom`

Return a comparison key for one chromosome name.

```python
canonical_chrom(value: 'object') -> 'str'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **value** | object | — | Chromosome name to normalize. |

### `chromosome_style`

Return `"ucsc"`, `"ensembl"` or `"mixed"` for a name collection.

```python
chromosome_style(names: 'Sequence[object]') -> 'str'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **names** | Sequence[object] | — | Chromosome names to classify. |

### `harmonize_frame`

Rename `frame[column]` values to matching target chromosome names.

```python
harmonize_frame(frame: 'pd.DataFrame', target: 'Iterable[object]', column: 'str' = 'chrom', *, copy: 'bool' = True) -> 'tuple[pd.DataFrame, dict[str, str]]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **frame** | DataFrame | — | DataFrame to harmonize. |
| **target** | Iterable[object] | — | Chromosome names to match. |
| `column` | str | `chrom` | Chromosome column name. |
| `copy` | bool | `True` | Copy the frame before modifying. |

### `detect_assembly`

Infer a bundled assembly from chromosome names or chromosome sizes.

```python
detect_assembly(peaks: 'Any' = None, chrom_sizes: 'Optional[Dict[str, int]]' = None, *, bigwig: 'Union[str, Path, None]' = None, tolerance: 'float' = 0.0005) -> 'str'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `peaks` | Any | `None` | Peak table used to infer the assembly. |
| `chrom_sizes` | Optional[Dict[str, int]] | `None` | Chromosome sizes used for a stronger match. |
| `bigwig` | Union[str, Path, None] | `None` | bigWig header used to infer chromosome names/sizes. |
| `tolerance` | float | `0.0005` | Relative size tolerance for the match. |

### `load_genome`

Load a bundled assembly table and return a :class:`GeneAnnotation`.

```python
load_genome(assembly: 'Union[str, None]', refresh: 'bool' = False) -> 'GeneAnnotation'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **assembly** | Union[str, None] | — | Assembly name; `None` loads the default supported bundle. |
| `refresh` | bool | `False` | Rebuild/re-download even when cached. |

### `load_annotation`

Load a GTF from an `ensembl:`/`gencode:`/`url:` spec.

```python
load_annotation(spec: 'str', refresh: 'bool' = False) -> 'GeneAnnotation'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **spec** | str | — | Assembly name or `ensembl:`/`gencode:`/`url:` spec. |
| `refresh` | bool | `False` | Re-download/reload the genome bundle even when it is cached. |

### `supported_genomes`

Assemblies accepted by `TxDb=`.

```python
supported_genomes() -> 'tuple[str, ...]'
```

### `is_genome_spec`

True for `ensembl:`, `gencode:` and `url:` annotation specs.

```python
is_genome_spec(value: 'Any') -> 'bool'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **value** | Any | — | Value to normalize or test. |

### `parse_genome_spec`

Resolve a download spec to provider, release, URL and description.

```python
parse_genome_spec(spec: 'str') -> 'dict[str, str]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **spec** | str | — | Genome specification: assembly name, `ensembl:`, `gencode:` or `url:` spec. |

### `download_genome`

Fetch a genome annotation into the cache and return the npz path.

```python
download_genome(assembly: 'Union[str, None]') -> 'Path'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **assembly** | Union[str, None] | — | Assembly to download. |

### `download_gtf`

Download the GTF referenced by an `ensembl:`/`gencode:`/`url:` spec.

```python
download_gtf(spec: 'str', dest: 'Optional[Union[str, Path]]' = None) -> 'Path'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **spec** | str | — | Ensembl/GENCODE/URL spec. |
| `dest` | Optional[Union[str, Path]] | `None` | Destination GTF path. |

### `cached_genomes`

Describe installed and downloaded bundles for `pychipseeker cache`.

```python
cached_genomes() -> 'list[dict[str, object]]'
```

### `clear_genome_cache`

Delete downloaded bundles from the cache directory.

```python
clear_genome_cache(assembly: 'Union[str, None]' = None) -> 'list[Path]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `assembly` | Union[str, None] | `None` | Assembly to clear; `None` clears all user-cached genomes. |

### `prune_genome_cache`

Remove incomplete cache entries and orphaned provenance JSON files.

```python
prune_genome_cache() -> 'list[Path]'
```

### `genome_provenance`

Return the provenance record of a genome annotation.

```python
genome_provenance(assembly: 'Union[str, None]') -> 'dict'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **assembly** | Union[str, None] | — | Assembly whose provenance record should be returned. |

## Overlap and coverage

Set overlap, shuffling, Venn/UpSet plots, coverage plots and overlap tests.

### `overlap`

Calculate exact membership counts for a list/dict of sets.

```python
overlap(Sets) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **Sets** | — | — | Peak sets or mapping of set names to peaks. |

### `shuffle`

Randomly relocate peaks within chromosomes, preserving width.

```python
shuffle(peak: 'Union[str, Path, pd.DataFrame]', chrom_sizes: 'Optional[Union[dict, pd.DataFrame]]' = None, TxDb: 'Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]]' = None, n: 'Optional[int]' = None, seed: 'Optional[int]' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame] | — | Peaks to shuffle. |
| `chrom_sizes` | Optional[Union[dict, DataFrame]] | `None` | Chromosome lengths used to draw new positions. |
| `TxDb` | Optional[Union[GeneAnnotation, str, Path, DataFrame]] | `None` | Gene model alternative to `chrom_sizes`. |
| `n` | Optional[int] | `None` | Number of shuffled peak sets to generate. |
| `seed` | Optional[int] | `None` | Random seed for reproducible shuffling, sampling or permutations. |

### `vennplot`

Plot a Venn diagram for 1-3 sets, or an overlap bar chart for more.

```python
vennplot(Sets, labels: 'Optional[Sequence[str]]' = None, title: 'str' = '', colors: 'Optional[Sequence[str]]' = None, figsize: 'tuple[float, float]' = (5.5, 5.0), ax=None)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **Sets** | — | — | Peak sets or mapping of names to peaks. |
| `labels` | Optional[Sequence[str]] | `None` | Set labels. |
| `title` | str | `` | Figure or table title. |
| `colors` | Optional[Sequence[str]] | `None` | Set colours. |
| `figsize` | tuple[float, float] | `(5.5, 5.0)` | Figure size in inches, as `(width, height)`. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |

### `vennplot_peakfile`

`vennplot` for a list of peak files (peaks compared by coordinates).

```python
vennplot_peakfile(peakfiles: 'Sequence[Union[str, Path, pd.DataFrame]]', labels: 'Optional[Sequence[str]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peakfiles** | Sequence[Union[str, Path, DataFrame]] | — | Peak files to compare. |
| `labels` | Optional[Sequence[str]] | `None` | Display labels for the sets or windows. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `vennpie`

Draw nested pies summarising ChIPseeker's genomic annotation detail.

```python
vennpie(annotation: 'PeakAnnotation', r: 'float' = 0.2, figsize=(6.0, 5.0), ax=None)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | PeakAnnotation | — | `PeakAnnotation` whose detailed labels are summarized. |
| `r` | float | `0.2` | Inner radius of the nested pie. |
| `figsize` | — | `(6.0, 5.0)` | Figure size in inches, as `(width, height)`. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |

### `upsetplot`

Draw a real UpSet plot for ChIPseeker's detailed annotation flags.

```python
upsetplot(annotation: 'PeakAnnotation', order_by: 'str' = 'freq', vennpie: 'bool' = False, n_intersections: 'int' = 20, figsize: 'tuple[float, float]' = (8.0, 4.8), ax=None)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **annotation** | PeakAnnotation | — | `PeakAnnotation` whose annotation flags are compared. |
| `order_by` | str | `freq` | `freq` or `degree`. |
| `vennpie` | bool | `False` | Add a nested pie alongside the UpSet plot. |
| `n_intersections` | int | `20` | Maximum intersections shown. |
| `figsize` | tuple[float, float] | `(8.0, 4.8)` | Figure size in inches, as `(width, height)`. |
| `ax` | — | `None` | Existing matplotlib axes; a new figure is created when omitted. |

### `covplot`

Plot peak coverage along chromosomes (ChIPseeker `covplot`).

```python
covplot(peak: 'Union[str, Path, pd.DataFrame, Sequence]', weightCol: 'Optional[str]' = None, xlab: 'str' = 'Genomic Coordinates', ylab: 'str' = 'Coverage', title: 'str' = '', chrs: 'Optional[Sequence[str]]' = None, xlim: 'Optional[Sequence[int]]' = None, lower: 'Optional[float]' = None, fill_color: 'Optional[str]' = None, figsize: 'Optional[tuple[float, float]]' = None)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peak** | Union[str, Path, DataFrame, Sequence] | — | Peaks or peak file to draw. |
| `weightCol` | Optional[str] | `None` | Optional weight column. |
| `xlab` | str | `Genomic Coordinates` | X-axis label. |
| `ylab` | str | `Coverage` | Y-axis label. |
| `title` | str | `` | Figure or table title. |
| `chrs` | Optional[Sequence[str]] | `None` | Chromosomes to draw. |
| `xlim` | Optional[Sequence[int]] | `None` | X-axis limits. |
| `lower` | Optional[float] | `None` | Lower y bound. |
| `fill_color` | Optional[str] | `None` | Coverage fill colour. |
| `figsize` | Optional[tuple[float, float]] | `None` | Figure size in inches, as `(width, height)`. |

### `enrich_peak_overlap`

Permutation test for overlap between query and target peak sets.

```python
enrich_peak_overlap(queryPeak: 'Union[str, Path, pd.DataFrame]', targetPeak, chrom_sizes: 'Optional[Union[dict, pd.DataFrame]]' = None, nShuffle: 'int' = 1000, pAdjustMethod: 'str' = 'BH', seed: 'Optional[int]' = None, pool: 'bool' = True, **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **queryPeak** | Union[str, Path, DataFrame] | — | Query peak set. |
| **targetPeak** | — | — | Target peak set. |
| `chrom_sizes` | Optional[Union[dict, DataFrame]] | `None` | Chromosome lengths for shuffling. |
| `nShuffle` | int | `1000` | Number of shuffles. |
| `pAdjustMethod` | str | `BH` | `BH`, `Bonferroni` or `none`. |
| `seed` | Optional[int] | `None` | Random seed for reproducible shuffling, sampling or permutations. |
| `pool` | bool | `True` | Pool the two sets before shuffling. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

### `enrich_anno_overlap`

Hypergeometric enrichment of genes annotated to two peak sets.

```python
enrich_anno_overlap(queryPeak: 'Union[str, Path, pd.DataFrame]', targetPeak, TxDb=None, gtf=None, annotation=None, level: 'str' = 'gene', tss_region: 'Sequence[int]' = (-3000, 3000), pAdjustMethod: 'str' = 'BH', **kwargs) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **queryPeak** | Union[str, Path, DataFrame] | — | Query peaks. |
| **targetPeak** | — | — | Target peaks. |
| `TxDb` | — | `None` | Gene model used for annotation. |
| `gtf` | — | `None` | GTF/GFF path alternative to `TxDb`. |
| `annotation` | — | `None` | Peak annotation table or `PeakAnnotation` object. |
| `level` | str | `gene` | Annotation level: `'transcript'` (default) or `'gene'`. |
| `tss_region` | Sequence[int] | `(-3000, 3000)` | Promoter window used for gene links. |
| `pAdjustMethod` | str | `BH` | Multiple-testing correction: `BH`, `Bonferroni` or `none`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

## GEO helpers

ChIPseeker-compatible helpers backed by the bundled GEO snapshot.

### `get_geo_species`

Species statistics collected from the bundled GEO snapshot.

```python
getGEOspecies(gsminfo=None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gsminfo` | — | `None` | Optional GEO snapshot path/DataFrame. |

### `get_geo_genome_version`

Genome-version statistics collected from GEO ChIP-seq samples.

```python
getGEOgenomeVersion(gsminfo=None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gsminfo` | — | `None` | Optional GEO snapshot path/DataFrame. |

### `get_geo_info`

Subset the GEO snapshot by genome version (ChIPseeker `getGEOInfo`).

```python
getGEOInfo(genome, simplify: 'bool' = True, gsminfo=None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **genome** | — | — | Genome version to subset. |
| `simplify` | bool | `True` | Keep only key columns. |
| `gsminfo` | — | `None` | Optional GEO snapshot path or DataFrame; defaults to the bundled table. |

### `download_geo_bed_files`

Download all BED supplementary files for `genome` from GEO.

```python
downloadGEObedFiles(genome, destDir=None, gsminfo=None) -> 'list[Path]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **genome** | — | — | Genome version to download. |
| `destDir` | — | `None` | Destination directory. |
| `gsminfo` | — | `None` | Optional GEO snapshot path or DataFrame; defaults to the bundled table. |

### `download_gsm_bed_files`

Download BED supplementary files for the given GSM accessions.

```python
downloadGSMbedFiles(GSM, destDir=None, gsminfo=None) -> 'list[Path]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **GSM** | — | — | GSM accession(s) to download. |
| `destDir` | — | `None` | Destination directory. |
| `gsminfo` | — | `None` | Optional GEO snapshot path or DataFrame; defaults to the bundled table. |

### `get_sample_files`

Paths of the bundled sample peak files (ChIPseeker `getSampleFiles`).

```python
get_sample_files(directory=None) -> 'dict[str, Path]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `directory` | — | `None` | Directory to scan for sample peaks. |

## Reports

Build the multi-page PDF report.

### `build_report`

Build a multi-page PDF report and return its path.

!!! note
    The PDF includes title/summary, annotation pie/bar, TSS distance, coverage, tag-matrix pages, optional signal pages and optional QC pages.

```python
build_report(peaks: 'Union[str, Path, pd.DataFrame]', *, gtf: 'Optional[Union[str, Path, pd.DataFrame]]' = None, TxDb: 'Optional[Union[str, Path, pd.DataFrame]]' = None, bigwig: 'Optional[Union[str, Path]]' = None, bam: 'Optional[Union[str, Path]]' = None, blacklist: 'Optional[Union[str, Path, pd.DataFrame]]' = None, output: 'Union[str, Path]' = 'pychipseeker_report.pdf', level: 'str' = 'transcript', tss_region: 'Sequence[int]' = (-3000, 3000), upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, sample: 'int' = 3000, seed: 'int' = 20261026, title: 'Optional[str]' = None, include_tag_matrix: 'bool' = True, qc: 'bool' = False, n_jobs: 'Optional[int]' = None) -> 'Path'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | Union[str, Path, DataFrame] | — | Peaks included in the report. |
| `gtf` | Optional[Union[str, Path, DataFrame]] | `None` | Gene model for annotation. |
| `TxDb` | Optional[Union[str, Path, DataFrame]] | `None` | Assembly or gene model alternative to `gtf`. |
| `bigwig` | Optional[Union[str, Path]] | `None` | Add signal TSS profile/heatmap pages. |
| `bam` | Optional[Union[str, Path]] | `None` | Add BAM-based FRiP/TSS QC pages. |
| `blacklist` | Optional[Union[str, Path, DataFrame]] | `None` | Blacklist used by the QC pages. |
| `output` | Union[str, Path] | `pychipseeker_report.pdf` | Destination PDF path. |
| `level` | str | `transcript` | `transcript` or `gene` annotation level. |
| `tss_region` | Sequence[int] | `(-3000, 3000)` | Promoter window for annotation. |
| `upstream` | Optional[int] | `None` | Asymmetric signal window upstream. |
| `downstream` | Optional[int] | `None` | Asymmetric signal window downstream. |
| `sample` | int | `3000` | Number of transcripts used for signal pages. |
| `seed` | int | `20261026` | Sampling seed for signal pages. |
| `title` | Optional[str] | `None` | Report title. |
| `include_tag_matrix` | bool | `True` | Include tag-matrix pages. |
| `qc` | bool | `False` | Include QC/FRiP pages. |
| `n_jobs` | Optional[int] | `None` | Parallel worker count for annotation/QC. |

## Classes and methods

The functional API returns three public objects. Their methods accept the
same parameters as the corresponding top-level helpers and are the usual
entry point for interactive work:

### `PeakAnnotation`

Returned by `annotate_peak()` / `annotatePeak()`. It stores the annotated peaks, the gene model and the settings used for annotation, and exposes the ChIPseeker-style plotting and summary methods.

#### Constructor {#peakannotation-constructor}

```python
PeakAnnotation(peaks: 'pd.DataFrame', annotation: 'pd.DataFrame', gene_annotation: 'GeneAnnotation', tss_region: 'tuple[int, int]' = (-3000, 3000), level: 'str' = 'transcript', promoter_bins: 'tuple[int, ...]' = (), genomic_annotation_priority: 'tuple[str, ...]' = ('Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstream', 'Intergenic'), downstream_distance: 'int' = 300, detail_genomic_annotation: 'pd.DataFrame' = <factory>, has_genomic_annotation: 'bool' = True, add_flank_gene_info: 'bool' = False, flank_distance: 'int' = 5000) -> None
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **peaks** | DataFrame | — | Peak table to annotate. |
| **annotation** | DataFrame | — | One row per annotated peak. |
| **gene_annotation** | GeneAnnotation | — | Gene model used for annotation. |
| `tss_region` | tuple[int, int] | `(-3000, 3000)` | Promoter window `(upstream, downstream)` in bp. |
| `level` | str | `transcript` | `transcript` (default) or `gene`. |
| `promoter_bins` | tuple[int, ...] | `()` | Promoter bin edges in bp. |
| `genomic_annotation_priority` | tuple[str, ...] | `('Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstrea...` | Feature priority used for overlapping annotations. |
| `downstream_distance` | int | `300` | Downstream window used by the `Downstream` label. |
| `detail_genomic_annotation` | DataFrame | `<factory>` | Detailed genomic-feature table. |
| `has_genomic_annotation` | bool | `True` | Whether genomic features were assigned. |
| `add_flank_gene_info` | bool | `False` | Whether flanking genes were added. |
| `flank_distance` | int | `5000` | Flank distance used for flanking-gene searches. |

#### Data attributes

| Attribute | Type | Description |
|---|---|---|
| `peaks` | DataFrame | Input peak table, including `score` when available. |
| `annotation` | DataFrame | One row per annotated peak with snake_case columns. |
| `gene_annotation` | GeneAnnotation | Gene model used for annotation. |
| `tss_region` | tuple[int, int] | Promoter window `(upstream, downstream)` in bp. |
| `level` | str | `transcript` (default) or `gene`. |
| `promoter_bins` | tuple[int, ...] | Promoter bin edges used for labels. |
| `genomic_annotation_priority` | tuple[str, ...] | Feature priority used for overlapping annotations. |
| `downstream_distance` | int | Downstream window used by the `Downstream` label. |
| `detail_genomic_annotation` | DataFrame | Detailed genomic-feature table. |
| `has_genomic_annotation` | bool | Whether genomic features were assigned. |
| `add_flank_gene_info` | bool | Whether flanking genes were added. |
| `flank_distance` | int | Flank distance used for flanking-gene searches. |
| `upstream` | int | Read-only property derived from `tss_region`. |
| `downstream` | int | Read-only property derived from `tss_region`. |
| `detail` | DataFrame | Alias for `detail_genomic_annotation`. |
| `annoStat` | DataFrame | Alias for `anno_stat()`. |
| `peakNum` | int | Number of annotated peaks. |
| `hasGenomicAnnotation` | bool | Alias for `has_genomic_annotation`. |

#### Methods

#### `to_dataframe` {#peakannotation-to_dataframe}

Return one row per annotated peak with snake_case columns.

```python
to_dataframe(include_peak_metadata: 'bool' = False) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `include_peak_metadata` | bool | `False` | Merge peak metadata columns (`name`, `score`, `strand`) into the result. |

#### `to_chipseeker` {#peakannotation-to_chipseeker}

Return a table with ChIPseeker-style column names.

```python
to_chipseeker() -> 'pd.DataFrame'
```

#### `as_data_frame` {#peakannotation-as_data_frame}

ChIPseeker `as.data.frame(csAnno)` analogue.

```python
as_data_frame() -> 'pd.DataFrame'
```

Aliases: `asDataFrame`

#### `as_granges` {#peakannotation-as_granges}

Return a GRanges-like table with `seqnames/start/end/strand`.

```python
as_granges() -> 'pd.DataFrame'
```

Aliases: `asGRanges`

#### `get_anno_stat` {#peakannotation-get_anno_stat}

Return the ChIPseeker `getAnnoStat` table.

```python
get_anno_stat() -> 'pd.DataFrame'
```

Aliases: `getAnnoStat`

#### `drop_anno` {#peakannotation-drop_anno}

Keep peaks within `abs(distanceToTSS) < cutoff`.

```python
drop_anno(distance_to_tss_cutoff: 'float' = 10000) -> "'PeakAnnotation'"
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `distance_to_tss_cutoff` | float | `10000` | Keep peaks with `abs(distance_to_tss) < cutoff`. |

Aliases: `dropAnno`

#### `anno_stat` {#peakannotation-anno_stat}

Summarize major annotation classes.

```python
anno_stat(ignore_promoter_subcategory: 'bool' = False, ignore_1st_exon: 'bool' = False, ignore_1st_intron: 'bool' = False, ignore_downstream: 'bool' = False) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ignore_promoter_subcategory` | bool | `False` | Merge promoter sub-bins into one `Promoter` row. |
| `ignore_1st_exon` | bool | `False` | Merge first-exon categories into `Exon`. |
| `ignore_1st_intron` | bool | `False` | Merge first-intron categories into `Intron`. |
| `ignore_downstream` | bool | `False` | Ignore downstream features when assigning annotations. |

#### `summary` {#peakannotation-summary}

Count peaks per class and optionally per metadata group.

```python
summary(groupby: 'Optional[str]' = None, normalize: 'bool' = False) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `groupby` | Optional[str] | `None` | Optional metadata column used to split the summary. |
| `normalize` | bool | `False` | Return percentages instead of raw counts. |

#### `seq2gene` {#peakannotation-seq2gene}

Return gene IDs linked to peaks through promoters or flanks.

!!! note
    Promoter and flank logic uses the annotation's `tss_region` and `flank_distance` unless overridden here.

```python
seq2gene(tss_region: 'Optional[Sequence[int]]' = None, flank_distance: 'int' = 5000, same_strand: 'bool' = False) -> 'list[str]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tss_region` | Optional[Sequence[int]] | `None` | Promoter window used for the link. |
| `flank_distance` | int | `5000` | Flank window used when a peak is not in a promoter. |
| `same_strand` | bool | `False` | Only consider features on the same strand as the peak. |

#### `seq2gene_links` {#peakannotation-seq2gene_links}

Return a peak-to-gene link table for enrichment.

```python
seq2gene_links(tss_region: 'Optional[Sequence[int]]' = None, flank_distance: 'int' = 5000, same_strand: 'bool' = False) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tss_region` | Optional[Sequence[int]] | `None` | Promoter window as `(upstream, downstream)` relative to the TSS. |
| `flank_distance` | int | `5000` | Distance in bp used to collect flanking genes. |
| `same_strand` | bool | `False` | Only consider features on the same strand as the peak. |

#### `plot_anno_pie` {#peakannotation-plot_anno_pie}

Draw the annotation composition pie chart.

```python
plot_anno_pie(**kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

Aliases: `plotAnnoPie`

#### `plot_anno_bar` {#peakannotation-plot_anno_bar}

Draw annotation composition as grouped bars.

```python
plot_anno_bar(**kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

Aliases: `plotAnnoBar`

#### `plot_dist_to_tss` {#peakannotation-plot_dist_to_tss}

Draw the TSS-distance distribution.

```python
plot_dist_to_tss(**kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

Aliases: `plotDistToTSS`

#### `vennpie` {#peakannotation-vennpie}

Draw a nested Venn-pie summary of annotation flags.

```python
vennpie(**kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

#### `upsetplot` {#peakannotation-upsetplot}

Draw an UpSet plot of annotation flags.

```python
upsetplot(**kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

#### `get_tag_matrix` {#peakannotation-get_tag_matrix}

Build a ChIPseeker-style tag matrix from this annotation.

!!! note
    When `windows` is a sequence, a list of `TagMatrix` objects is returned. With `bigwig=`/`bam=` the matrix contains signal coverage; otherwise it contains peak coverage.

```python
get_tag_matrix(peak=None, upstream: 'int' = 3000, downstream: 'int' = 3000, type: 'str' = 'start_site', by: 'Optional[str]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `peak` | — | `None` | Peaks used to build the matrix; defaults to `ann.peaks`. |
| `upstream` | int | `3000` | Upstream window in bp. |
| `downstream` | int | `3000` | Downstream window in bp. |
| `type` | str | `start_site` | Region type used when building windows. |
| `by` | Optional[str] | `None` | `gene` or `transcript`. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

Aliases: `getTagMatrix`

#### `plot_peak_prof` {#peakannotation-plot_peak_prof}

Plot the average peak-centered profile.

```python
plot_peak_prof(**kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

Aliases: `plotPeakProf`

#### `peak_heatmap` {#peakannotation-peak_heatmap}

Draw a peak-centered coverage heatmap.

```python
peak_heatmap(peak=None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `peak` | — | `None` | Peaks used to build the matrix. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

Aliases: `peakHeatmap`

#### `plot_tss_profile` {#peakannotation-plot_tss_profile}

Plot a bigWig, BAM or matrix TSS profile.

!!! note
    Pass exactly one of `bigwig`, `matrix` or `signal`; extra keyword arguments are forwarded to the top-level `plot_tss_profile`.

```python
plot_tss_profile(bigwig: 'Optional[Union[str, Path]]' = None, *, matrix: 'Optional[np.ndarray]' = None, signal: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bigwig` | Optional[Union[str, Path]] | `None` | bigWig signal file; requires the optional `pyBigWig` dependency. |
| `matrix` | Optional[ndarray] | `None` | Features x bins signal matrix. |
| `signal` | Optional[Union[str, Path, DataFrame]] | `None` | bigWig path or precomputed signal table. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

#### `plot_tss_heatmap` {#peakannotation-plot_tss_heatmap}

Draw a bigWig, BAM or matrix TSS heatmap.

!!! note
    Pass exactly one of `bigwig`, `matrix` or `signal`; extra keyword arguments are forwarded to the top-level `plot_tss_heatmap`.

```python
plot_tss_heatmap(bigwig: 'Optional[Union[str, Path]]' = None, *, matrix: 'Optional[np.ndarray]' = None, signal: 'Optional[Union[str, Path, pd.DataFrame]]' = None, **kwargs)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bigwig` | Optional[Union[str, Path]] | `None` | bigWig signal file; requires the optional `pyBigWig` dependency. |
| `matrix` | Optional[ndarray] | `None` | Features x bins signal matrix. |
| `signal` | Optional[Union[str, Path, DataFrame]] | `None` | bigWig path or precomputed signal table. |
| `kwargs` | — | keyword | Extra keyword arguments forwarded to matplotlib or the underlying helper. |

#### Method aliases

CamelCase aliases share the same signature and parameters as the
canonical method:

| Alias | Canonical method |
|---|---|
| `asDataFrame` | [`as_data_frame`](#peakannotation-as_data_frame) |
| `asGRanges` | [`as_granges`](#peakannotation-as_granges) |
| `getAnnoStat` | [`get_anno_stat`](#peakannotation-get_anno_stat) |
| `dropAnno` | [`drop_anno`](#peakannotation-drop_anno) |
| `plotAnnoPie` | [`plot_anno_pie`](#peakannotation-plot_anno_pie) |
| `plotAnnoBar` | [`plot_anno_bar`](#peakannotation-plot_anno_bar) |
| `plotDistToTSS` | [`plot_dist_to_tss`](#peakannotation-plot_dist_to_tss) |
| `getTagMatrix` | [`get_tag_matrix`](#peakannotation-get_tag_matrix) |
| `plotPeakProf` | [`plot_peak_prof`](#peakannotation-plot_peak_prof) |
| `peakHeatmap` | [`peak_heatmap`](#peakannotation-peak_heatmap) |

### `GeneAnnotation`

Gene, transcript, exon, intron and UTR tables plus interval indexes. Build one with `GeneAnnotation.from_gtf()` or let `annotate_peak()` create it from `TxDb=`, `gtf=` or `annotation=`.

#### Constructor {#geneannotation-constructor}

```python
GeneAnnotation(genes: 'pd.DataFrame', transcripts: 'pd.DataFrame', exons: 'pd.DataFrame', introns: 'pd.DataFrame', utr5: 'pd.DataFrame', utr3: 'pd.DataFrame', *, copy: 'bool' = True) -> 'None'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **genes** | DataFrame | — | Gene or transcript feature table; normally `ann.gene_annotation.features(...)`. |
| **transcripts** | DataFrame | — | Transcript feature table. |
| **exons** | DataFrame | — | Exon feature table. |
| **introns** | DataFrame | — | Intron feature table. |
| **utr5** | DataFrame | — | 5' UTR feature table. |
| **utr3** | DataFrame | — | 3' UTR feature table. |
| `copy` | bool | `True` | Return a copy instead of modifying the input frame in place. |

#### Data attributes

| Attribute | Type | Description |
|---|---|---|
| `genes` | DataFrame | Gene feature table. |
| `transcripts` | DataFrame | Transcript feature table. |
| `exons` | DataFrame | Exon feature table. |
| `introns` | DataFrame | Intron feature table. |
| `utr5` | DataFrame | 5' UTR feature table. |
| `utr3` | DataFrame | 3' UTR feature table. |
| `gene_index` | BinnedIntervalIndex | Binned interval index over genes. |
| `transcript_index` | BinnedIntervalIndex | Binned interval index over transcripts. |
| `exon_index` | BinnedIntervalIndex | Binned interval index over exons. |
| `intron_index` | BinnedIntervalIndex | Binned interval index over introns. |
| `utr5_index` | BinnedIntervalIndex | Binned interval index over 5' UTRs. |
| `utr3_index` | BinnedIntervalIndex | Binned interval index over 3' UTRs. |
| `tss_index` | TSSIndex | Backwards-compatible gene-level TSS index. |

#### Methods

#### `from_gtf` {#geneannotation-from_gtf}

Build a `GeneAnnotation` from a GTF/GFF path or feature table.

!!! note
    Coordinates are normalized to 0-based, half-open intervals.

```python
from_gtf(gtf: 'Union[str, Path, pd.DataFrame]') -> "'GeneAnnotation'"
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gtf** | Union[str, Path, DataFrame] | — | GTF/GFF path or a feature DataFrame. |

#### `features` {#geneannotation-features}

Return the gene or transcript feature table.

```python
features(level: 'str' = 'transcript') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `level` | str | `transcript` | `transcript` (default) or `gene`. |

#### `tss_index_for` {#geneannotation-tss_index_for}

Return (and cache) the TSS index for a level.

```python
tss_index_for(level: 'str') -> 'TSSIndex'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **level** | str | — | `transcript` or `gene`. |

#### `span_index_for` {#geneannotation-span_index_for}

Return (and cache) the span index for a level.

```python
span_index_for(level: 'str') -> 'SpanIndex'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **level** | str | — | `transcript` or `gene`. |

#### `feature_by_id` {#geneannotation-feature_by_id}

Look up one feature row by ID.

```python
feature_by_id(level: 'str', feature_id: 'str') -> 'Optional[pd.Series]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **level** | str | — | `transcript` or `gene`. |
| **feature_id** | str | — | Feature identifier to look up. |

#### `feature_row` {#geneannotation-feature_row}

Return a feature row by positional index.

```python
feature_row(level: 'str', row_index: 'int') -> 'pd.Series'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **level** | str | — | `transcript` or `gene`. |
| **row_index** | int | — | Zero-based row index in the feature table. |

#### `category_table` {#geneannotation-category_table}

Return the exon, intron or UTR table for a category.

```python
category_table(category: 'str') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **category** | str | — | `Exon`, `Intron`, `5UTR` or `3UTR`. |

#### `category_index` {#geneannotation-category_index}

Return the binned index for a feature category.

```python
category_index(category: 'str') -> 'BinnedIntervalIndex'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **category** | str | — | `Exon`, `Intron`, `5UTR` or `3UTR`. |

#### `gene_row` {#geneannotation-gene_row}

Look up a gene row by gene ID.

```python
gene_row(gene_id: 'str') -> 'Optional[pd.Series]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_id** | str | — | Gene identifier. |

#### `nearest_gene` {#geneannotation-nearest_gene}

Return the nearest gene and signed distance.

```python
nearest_gene(chrom: 'str', pos: 'int') -> 'tuple[Optional[pd.Series], float]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **chrom** | str | — | Chromosome name. |
| **pos** | int | — | 0-based genomic position. |

#### `choose_gene` {#geneannotation-choose_gene}

Choose the closest gene among candidate IDs.

```python
choose_gene(gene_ids: 'Iterable[str]', pos: 'int') -> 'tuple[Optional[pd.Series], float]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **gene_ids** | Iterable[str] | — | Candidate gene identifiers to compare. |
| **pos** | int | — | 0-based genomic position. |

#### `promoter_gene` {#geneannotation-promoter_gene}

Return the gene whose promoter contains `pos`.

```python
promoter_gene(chrom: 'str', pos: 'int', upstream: 'int', downstream: 'int') -> 'tuple[Optional[pd.Series], float]'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **chrom** | str | — | Chromosome name. |
| **pos** | int | — | 0-based genomic position. |
| **upstream** | int | — | Bases upstream of the anchor (the 5' side on the feature strand). |
| **downstream** | int | — | Bases downstream of the anchor (the 3' side on the feature strand). |

#### `promoter_label` {#geneannotation-promoter_label}

Map a TSS distance to a promoter-bin label.

```python
promoter_label(distance: 'float', bins: 'Sequence[int]') -> 'str'
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **distance** | float | — | Signed distance from the peak to the TSS. |
| **bins** | Sequence[int] | — | Promoter bin edges in bp. |

### `TagMatrix`

A `numpy.ndarray` subclass carrying ChIPseeker-style region metadata. It is returned by `get_tag_matrix()` and accepted by the profile and heatmap helpers.

#### Constructor {#tagmatrix-constructor}

```python
TagMatrix(matrix: 'np.ndarray', *, upstream: 'Optional[int]' = None, downstream: 'Optional[int]' = None, type: 'Optional[str]' = None, by: 'Optional[str]' = None, label: 'Optional[str]' = None, is_binning: 'bool' = False, rows: 'Optional[pd.DataFrame]' = None) -> "'TagMatrix'"
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| **matrix** | ndarray | — | Numeric features x bins matrix. |
| `upstream` | Optional[int] | `None` | Upstream window size in bp. |
| `downstream` | Optional[int] | `None` | Downstream window size in bp. |
| `type` | Optional[str] | `None` | Region type (`start_site`, `end_site` or `body`). |
| `by` | Optional[str] | `None` | Feature level (`gene` or `transcript`). |
| `label` | Optional[str] | `None` | Optional label used in plots. |
| `is_binning` | bool | `False` | Whether the matrix contains binned windows. |
| `rows` | Optional[DataFrame] | `None` | Region metadata with one row per matrix row. |

#### Data attributes

| Attribute | Type | Description |
|---|---|---|
| `upstream` | int \| None | Upstream window size in bp. |
| `downstream` | int \| None | Downstream window size in bp. |
| `type` | str \| None | Region type (`start_site`, `end_site` or `body`). |
| `by` | str \| None | Feature level (`gene` or `transcript`). |
| `label` | str \| None | Optional label used in plots. |
| `is_binning` | bool | Whether the matrix contains binned windows. |
| `rows` | DataFrame \| None | Region metadata with one row per matrix row. |
| `matrix` | ndarray | The underlying float matrix as a plain NumPy array. |

#### Methods

`TagMatrix` has no public methods. Use the attributes above or pass it
to `plot_avg_prof`, `tag_heatmap` and the other top-level helpers.
