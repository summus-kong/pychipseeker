# ChIPseeker compatibility matrix

Reference: **ChIPseeker 1.48.0** `NAMESPACE` (Bioconductor release 3.21).
Status is intentionally conservative.

Legend:

- **Core**: implemented with the same ChIPseeker semantics and covered by
  tests.
- **Approx.**: implemented in Python, but the R object model or plot details
  differ; exact numerical/visual identity is not guaranteed.
- **N/A**: an R/Bioconductor concept with no direct Python equivalent.
- **Missing**: not implemented yet.

| ChIPseeker API | pychipseeker | Notes |
|---|---|---|
| `annotatePeak` | Core | `annotatePeak` camelCase wrapper and `annotate_peak`; GTF/GFF or `GeneAnnotation` as `TxDb`; `annoDb` accepts a DataFrame/mapping |
| `as.GRanges` | Approx. | `asGRanges`/`as_granges` returns a GRanges-like pandas table |
| `as.data.frame,csAnno` | Core | `to_dataframe()` / `to_chipseeker()` |
| `show,csAnno` | Approx. | Python `repr`/printing |
| `combine_csAnno` | Core | `combine_cs_anno` / `combine_csAnno` |
| `covplot` | Approx. | matplotlib coverage plot |
| `downloadGEObedFiles` | Core | `downloadGEObedFiles`/`download_geo_bed_files`; bundled `gsminfo` snapshot, `ftp://`→`https://`, skips existing/failed downloads |
| `downloadGSMbedFiles` | Core | `downloadGSMbedFiles`/`download_gsm_bed_files`; filters by GSM accession from the bundled snapshot |
| `dropAnno` | Core | Filters by `abs(distanceToTSS) < cutoff`, does not drop columns |
| `enrichAnnoOverlap` | Approx. | Hypergeometric test on annotated gene sets |
| `enrichPeakOverlap` | Approx. | Permutation test; p-values can differ |
| `getAnnoStat` | Core | `getAnnoStat` / `get_anno_stat` |
| `getBioRegion` | Core | Gene/transcript/exon/intron/UTR bioregions from GTF |
| `getGEOInfo` | Core | `getGEOInfo`/`get_geo_info`; bundled `gsminfo.csv.gz` snapshot, `simplify=` supported, custom snapshot via `gsminfo=` |
| `getGEOgenomeVersion` | Core | `getGEOgenomeVersion`/`get_geo_genome_version`; `organism`/`genomeVersion`/`Freq` columns |
| `getGEOspecies` | Core | `getGEOspecies`/`get_geo_species`; `species`/`Freq` count table |
| `getPromoters` | Core | `getBioRegion(type="start_site")` |
| `getSampleFiles` | Core | `getSampleFiles`/`get_sample_files`; same GEO sample BED files as R extdata, same name regex |
| `getTagMatrix` | Core | Peak-coverage tag matrix, binning, weights, strand handling |
| `makeBioRegionFromGranges` | Core | DataFrame input instead of GRanges |
| `overlap` | Core | Exact membership counts for generic sets/peaks |
| `peakHeatmap` | Approx. | matplotlib heatmap |
| `peakHeatmap_multiple_Sets` | Approx. | Multi-panel matplotlib heatmaps |
| `peak_Profile_Heatmap` | Approx. | Combined matplotlib profile/heatmap |
| `plotAnnoBar` | Core | Same categories/percentages and stacked bar semantics |
| `plotAnnoPie` | Core | Same detailed annotation categories |
| `plotAnnoPie.csAnno` | Core | `PeakAnnotation.plot_anno_pie()` |
| `plotAvgProf` | Approx. | Mean profile with bootstrap-style band |
| `plotAvgProf2` | Approx. | Alias of the peak-to-matrix profile path |
| `plotDistToTSS` | Core | Diverging 5'->3' distance plot |
| `plotPeakProf` | Approx. | matplotlib profile, single or multiple matrices |
| `plotPeakProf2` | Approx. | Builds matrix from peaks and plots |
| `plotPeakProf_MultiWindows` | Approx. | One peak set, multiple bioregions |
| `plotMultiProf` | Approx. | Multiple tag matrices |
| `readPeakFile` | Core | `readPeakFile` / `read_peak_file` |
| `rel` | N/A | Re-exported from ggplot2 in ChIPseeker |
| `seq2gene` | Core | Many-to-many peak-to-gene assignment |
| `shuffle` | Approx. | Uniform relocation within chromosome bounds |
| `tagHeatmap` | Approx. | matplotlib heatmap |
| `upsetplot` | Approx. | Combination bar chart instead of ggupset |
| `vennpie` | Approx. | Nested matplotlib pies |
| `vennplot` | Approx. | 1-3 set Venn diagrams; bar chart for more sets |
| `vennplot.peakfile` | Approx. | `vennplot_peakfile`, peaks compared by coordinates |
| `GRangesList` | N/A | Re-exported Bioconductor class |

## Additional Python-only helpers

- `annotate_db`: apply a DataFrame/mapping of gene annotation columns.
- `get_tag_matrix_from_bigwig`: direct bigWig tag-matrix extraction.
- `plot_tss_profile` / `plot_tss_heatmap`: existing TSS signal helpers.
- `enrich_go`, `enrich_kegg`, `enrichr`: optional `gseapy` wrappers.

## Known semantic differences

1. **Coordinates**: pychipseeker uses 0-based, half-open intervals throughout;
   ChIPseeker/GRanges use 1-based, closed intervals.
2. **Annotation engine input**: GTF/GFF, a `GeneAnnotation` object or a
   bundled assembly name (`TxDb="hg38"`, `"hg19"`, `"mm10"`, `"mm9"` --
   GENCODE-based bundles fetched on demand into a local cache) replaces
   a Bioconductor `TxDb`.
3. **Isoform alignment**: the default follows ChIPseeker's development fix
   #252.  Use `align_annotation_feature=False` for release-1.48 behavior.
4. **`p_shared_exon`**: the devel fix intentionally differs from the 1.48
   release in shared-exon cases.
5. **Visualization**: matplotlib is used instead of ggplot2.  Labels,
   categories, percentages and axis direction follow ChIPseeker, but theme,
   spacing and exact aesthetics differ.
6. **`annoDb`**: a DataFrame or mapping is used instead of an OrgDb package.
7. **GEO functions**: the ChIPseeker `gsminfo` snapshot (originally collected
   from GEOmetadb) is bundled as `pychipseeker/data/gsminfo.csv.gz`
   (upstream data last update: Mar 03, 2015; extended by subsequent upstream
   releases).  `getGEOspecies`, `getGEOgenomeVersion`, `getGEOInfo`,
   `downloadGEObedFiles`, `downloadGSMbedFiles` and `getSampleFiles` operate
   on it exactly like the R functions; the downloaders rewrite `ftp://` to
   `https://` like upstream issue #254 and return the list of created
   destination paths (R returns `NULL`).  A custom snapshot can be passed
   through `gsminfo=...`.
   Upstream's internal `prepareGSMInfo()` maintainer batch job has an
   offline Python equivalent: `scripts/build_gsminfo.py` rebuilds a
   compatible `gsminfo` CSV from a locally downloaded `GEOmetadb.sqlite`
   (zero-dependency, stdlib only) plus the upstream
   `scripts/ucsc_release.csv.gz` mapping.  See the script docstring.
8. **bigWig reading**: pyBigWig is a declared dependency and the supported
   reader.  The libBigWig `fopen()` calls are patched to use `_wfopen()` so
   UTF-8 (non-ASCII, e.g. Chinese) local paths work on Windows on a build
   matching `pyBigWig` upstream; the tests use upstream's own
   `pyBigWigTest/test.bw` file.

## Re-running the R comparison

The core annotation comparison harness lives in `comparison/`.  It requires an
existing R installation with ChIPseeker and txdbmaker/Bioconductor packages;
it does not install or modify R libraries.  See the main README for usage.
