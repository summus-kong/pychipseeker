#!/usr/bin/env python3
# ruff: noqa: E402, I001
"""Generate the detailed function/parameter reference page.

The page is written to ``docs/reference/functions.md``.  Run it after adding
or changing a public function::

    python scripts/build_function_reference.py
"""

from __future__ import annotations

import inspect
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pychipseeker as pcs  # noqa: E402


CATEGORIES: Sequence[Tuple[str, str, Sequence[str]]] = (
    (
        "Annotation and result tables",
        "Annotate peaks, read/write tables, filter peaks and work with "
        "ChIPseeker-style result objects.",
        (
            "annotate_peak",
            "annotate_db",
            "seq2gene",
            "seq2gene_links",
            "get_flanking_genes",
            "get_genomic_anno_stat",
            "promoter_label",
            "group_annotation",
            "read_gtf",
            "read_peaks",
            "read_peak_file",
            "write_table",
            "write_peaks",
            "filter_peaks",
            "call_peaks",
            "call_peaks_macs",
            "as_granges",
            "combine_cs_anno",
            "drop_anno",
            "get_anno_stat",
            "detail_genomic_annotation",
        ),
    ),
    (
        "Analysis summaries",
        "Plot-ready summaries for annotated peak tables and tag matrices.",
        (
            "annotation_summary",
            "peaks_per_gene",
            "top_genes",
            "gene_type_summary",
            "distance_summary",
            "signal_groups",
            "tss_enrichment_by_quartile",
        ),
    ),
    (
        "Tag matrices and profiles",
        "Promoters, bioregions, tag matrices, profiles and heatmaps.",
        (
            "get_promoters",
            "get_bio_region",
            "make_bio_region_from_granges",
            "get_tag_matrix",
            "get_tag_matrix_from_bigwig",
            "get_tag_matrix_from_bam",
            "extract_tss_matrix",
            "extract_tss_matrix_from_signal",
            "plot_avg_prof",
            "plot_peak_prof2",
            "plot_peak_prof",
            "plot_peak_prof_multi_windows",
            "plot_multi_prof",
            "tag_heatmap",
            "peak_heatmap",
            "peak_heatmap_multiple_sets",
            "peak_profile_heatmap",
            "plot_tss_profile",
            "plot_tss_heatmap",
        ),
    ),
    (
        "Annotation plots",
        "Pie, bar and distance-to-TSS figures for a PeakAnnotation.",
        ("plot_anno_pie", "plot_anno_bar", "plot_dist_to_tss"),
    ),
    (
        "Quality control",
        "Peak widths, blacklist/mitochondrial fractions, FRiP and TSS enrichment.",
        ("peak_qc", "frip", "tss_enrichment", "qc_summary"),
    ),
    (
        "Enrichment",
        "Local over-representation analysis, Enrichr, GSEA and enrichment plots.",
        (
            "read_gmt",
            "ora",
            "enrichr",
            "enrich_go",
            "enrich_kegg",
            "enrich_peak_genes",
            "gsea",
            "gsea_peak_genes",
            "plot_enrichment",
        ),
    ),
    (
        "Genomes and cache",
        "Assembly detection, genome bundles, GTFs and the local cache.",
        (
            "canonical_chrom",
            "chromosome_style",
            "harmonize_frame",
            "detect_assembly",
            "load_genome",
            "load_annotation",
            "supported_genomes",
            "is_genome_spec",
            "parse_genome_spec",
            "download_genome",
            "download_gtf",
            "cached_genomes",
            "clear_genome_cache",
            "prune_genome_cache",
            "genome_provenance",
        ),
    ),
    (
        "Overlap and coverage",
        "Set overlap, shuffling, Venn/UpSet plots, coverage plots and overlap tests.",
        (
            "overlap",
            "shuffle",
            "vennplot",
            "vennplot_peakfile",
            "vennpie",
            "upsetplot",
            "covplot",
            "enrich_peak_overlap",
            "enrich_anno_overlap",
        ),
    ),
    (
        "GEO helpers",
        "ChIPseeker-compatible helpers backed by the bundled GEO snapshot.",
        (
            "get_geo_species",
            "get_geo_genome_version",
            "get_geo_info",
            "download_geo_bed_files",
            "download_gsm_bed_files",
            "get_sample_files",
        ),
    ),
    (
        "Reports",
        "Build the multi-page PDF report.",
        ("build_report",),
    ),
)


PARAM_DOCS: Dict[str, str] = {
    "annotation": "Peak annotation table or `PeakAnnotation` object.",
    "kwargs": "Extra keyword arguments forwarded to matplotlib or the underlying helper.",
    "upstream": "Bases upstream of the anchor (the 5' side on the feature strand).",
    "downstream": "Bases downstream of the anchor (the 3' side on the feature strand).",
    "gtf": "GTF/GFF path, feature DataFrame or `GeneAnnotation`.",
    "TxDb": "Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`.",
    "peaks": "Peak file path, peak DataFrame, or a sequence of peak sources.",
    "peak": "Peak file path or peak DataFrame.",
    "figsize": "Figure size in inches, as `(width, height)`.",
    "title": "Figure or table title.",
    "ax": "Existing matplotlib axes; a new figure is created when omitted.",
    "n_jobs": "Worker count: `None`/`1` sequential, `-1` all cores, `n > 1` explicit.",
    "bigwig": "bigWig signal file; requires the optional `pyBigWig` dependency.",
    "xlab": "X-axis label.",
    "ylab": "Y-axis label.",
    "by": "Feature level used to build windows: `'gene'` or `'transcript'`.",
    "nbin": "Number of bins across each window; `None` uses the column widths when available.",
    "weightCol": "Optional peak column used as weight (for example `signalValue`).",
    "type": "Region type: `'start_site'`, `'end_site'` or `'body'`.",
    "ignore_strand": "Ignore feature strand when building windows or profiles.",
    "seed": "Random seed for reproducible shuffling, sampling or permutations.",
    "chrom_sizes": "Mapping or DataFrame of chromosome lengths used for shuffling, FRiP or QC.",
    "level": "Annotation level: `'transcript'` (default) or `'gene'`.",
    "bam": "BAM/CRAM alignment file; requires the optional `pysam` dependency.",
    "conf": "Confidence level for the bootstrap band; `None` disables the band.",
    "tss_region": "Promoter window as `(upstream, downstream)` relative to the TSS.",
    "verbose": "Print progress messages when `True`.",
    "path": "Input or output path.",
    "genome": "Assembly name or chromosome-sizes source used for QC.",
    "outdir": "Output directory.",
    "gsminfo": "Optional GEO snapshot path or DataFrame; defaults to the bundled table.",
    "windows": "Pre-built windows (DataFrame or path) or a sequence of window tables.",
    "gene_list": "Iterable of gene identifiers or a ranked mapping.",
    "blacklist": "Blacklist BED path or DataFrame whose intervals are removed.",
    "exclude_duplicates": "Ignore duplicate reads/alignments when `True`.",
    "genes": "Gene or transcript feature table; normally `ann.gene_annotation.features(...)`.",
    "assembly": "Assembly name such as `mm10`, `hg38` or `dm6`; `None` means auto/cache-driven.",
    "colors": "Sequence of colours used for the sets or bars.",
    "window": "Half-window size in bp, applied to both sides unless `upstream`/`downstream` are given.",
    "facet": "Facet layout: `'none'`, `'row'` or `'column'`.",
    "free_y": "Let each facet use its own y-axis range when `True`.",
    "flank_distance": "Distance in bp used to collect flanking genes.",
    "same_strand": "Only consider features on the same strand as the peak.",
    "sample": "Maximum number of features/rows to sample; `None` or `0` keeps all.",
    "min_mapping_quality": "Minimum BAM mapping quality (MAPQ) for covered bases.",
    "n_bins": "Number of bins across the window.",
    "spec": "Genome specification: assembly name, `ensembl:`, `gencode:` or `url:` spec.",
    "labels": "Display labels for the sets or windows.",
    "value_type": "bigWig summary: `sum`, `mean`, `max`, `min`, `coverage` or `std`.",
    "tagMatrix": "A `TagMatrix` or a numeric matrix with one row per region.",
    "gene_sets": "GMT path, library name, or mapping of term to gene list.",
    "organism": "Enrichr organism name, for example `human` or `mouse`.",
    "top": "Number of rows/terms to keep.",
    "tag": "A `TagMatrix` (or numeric matrix) with one row per region.",
    "n_groups": "Number of signal groups; capped at the number of matrix rows.",
    "n": "Number of items to return or shuffle.",
    "ignore_downstream": "Ignore downstream features when assigning annotations.",
    "downstream_distance": "Downstream window in bp used for the `Downstream` category.",
    "gene_annotation": "`GeneAnnotation` used to build promoters/links.",
    "distance": "Distance in bp used for flanking-gene searches.",
    "annotations": "Sequence of annotation labels to summarize.",
    "min_width": "Minimum peak width in bp.",
    "chroms": "Chromosomes to keep or process.",
    "value": "Value to normalize or test.",
    "refresh": "Re-download/reload the genome bundle even when it is cached.",
    "Sets": "Sequence of peak sets or a mapping of set name to peaks.",
    "destDir": "Destination directory for downloaded files.",
    "xlim": "X-axis limits.",
    "queryPeak": "Query peak set.",
    "targetPeak": "Target peak set.",
    "pAdjustMethod": "Multiple-testing correction: `BH`, `Bonferroni` or `none`.",
    "groupby": "Optional metadata column used to split the summary.",
    "min_value": "Minimum signal value; lower values are clipped or treated as empty.",
    "drop_empty": "Drop rows with no signal when `True`.",
    "return_mask": "Also return the boolean mask of retained rows.",
    "dtype": "NumPy dtype of the returned matrix.",
    "matrix": "Numeric features x bins matrix.",
    "color": "Line or bar colour.",
    "cmap": "Matplotlib colormap name.",
    "windows_name": "Names for the window sets, used as facet labels.",
    "palette": "Matplotlib colormap or palette name.",
    "background": "Background gene set for enrichment; defaults to the gene-set union.",
    "p_adjust": "Multiple-testing correction: `BH`, `Bonferroni` or `none`.",
    "gene_column": "Annotation column that supplies gene identifiers.",
    "thresholds": "Distance cutoffs in bp for the percentage rows.",
    "center_fraction": "Fraction of bins treated as the centre window.",
    "flank_fraction": "Fraction of bins treated as the outer flanks.",
    "assign_genomic_annotation": "Assign promoter/UTR/exon/intron labels when `True`.",
    "genomic_annotation_priority": "Feature priority used to resolve overlapping annotations.",
    "add_flank_gene_info": "Add genes within `flank_distance` of each peak.",
    "ignore_overlap": "Skip overlapping-feature assignment when `True`.",
    "ignore_upstream": "Ignore upstream features when assigning annotations.",
    "overlap": "Overlap mode: `'TSS'` (default) or `'all'` feature spans.",
    "promoter_bins": "Promoter bin edges in bp, e.g. `(1000, 2000, 3000)`.",
    "drop_unmapped": "Drop peaks without a genomic annotation when `True`.",
    "align_annotation_feature": "Report the isoform that supplied the genomic feature (ChIPseeker fix #252).",
    "anno_db": "Gene-annotation table, mapping or path used to add gene metadata.",
    "columns": "Gene metadata columns to add; `None` keeps all available columns.",
    "model": "`GeneAnnotation` used for the flanking-gene search.",
    "ignore_promoter_subcategory": "Merge promoter sub-bins into one `Promoter` row.",
    "ignore_1st_exon": "Merge first-exon categories into `Exon`.",
    "ignore_1st_intron": "Merge first-intron categories into `Intron`.",
    "bins": "Promoter bin edges in bp.",
    "features": "Feature types to keep (for example `gene`, `transcript`, `exon`).",
    "chunksize": "Read the GTF in chunks of this many rows; `None` reads at once.",
    "fmt": "Explicit input format (for example `bed`, `narrowPeak`, `csv`); inferred when omitted.",
    "as_type": "Return type: `'dataframe'` (default) or `'GRanges'`.",
    "df": "Table to write or summarize.",
    "threshold": "Signal threshold for peak calling.",
    "max_gap": "Maximum gap in bp between consecutive signal runs that are merged.",
    "treatment": "Treatment alignment file for MACS.",
    "control": "Optional control alignment file for MACS.",
    "backend": "MACS executable: `macs2` (default) or `macs3`.",
    "extra_args": "Extra command-line arguments forwarded to MACS.",
    "name": "Sample name passed to MACS.",
    "drop_mito": "Remove mitochondrial chromosomes when `True`.",
    "max_width": "Maximum peak width in bp.",
    "min_score": "Minimum peak score.",
    "exclude_chroms": "Chromosomes to remove.",
    "tss_upstream": "Upstream TSS window used for TSS enrichment.",
    "tss_downstream": "Downstream TSS window used for TSS enrichment.",
    "tss_sample": "Maximum number of TSS features sampled for QC.",
    "center_window": "Half-width in bp of the centre window for TSS enrichment.",
    "flank_window": "Half-width in bp of the flank window; defaults to the centre width.",
    "names": "Sequence of chromosome names to classify.",
    "frame": "DataFrame whose chromosome column is harmonized.",
    "target": "Target chromosome names or values to harmonize against.",
    "column": "Column containing chromosome names.",
    "copy": "Return a copy instead of modifying the input frame in place.",
    "tolerance": "Relative chromosome-size tolerance for assembly detection.",
    "dest": "Destination path for the downloaded GTF.",
    "output": "Output PDF/table path.",
    "include_tag_matrix": "Include tag-matrix profile and heatmap pages in the report.",
    "qc": "Include the QC/FRiP pages in the report.",
    "distance_to_tss_cutoff": "Keep peaks whose absolute TSS distance is below this cutoff.",
    "simplify": "Keep only the most useful GEO columns when `True`.",
    "GSM": "GSM accession or sequence of accessions.",
    "directory": "Directory to scan; defaults to the bundled sample files.",
    "peakfiles": "Sequence of peak files to compare.",
    "r": "Inner pie radius for the nested Venn pie.",
    "order_by": "Intersection ordering: `'freq'` or `'degree'`.",
    "vennpie": "Add a nested pie summary when `True`.",
    "n_intersections": "Maximum number of intersections shown in the UpSet plot.",
    "chrs": "Chromosomes to draw.",
    "lower": "Lower y-axis bound for coverage plots.",
    "fill_color": "Coverage fill colour.",
    "nShuffle": "Number of shuffles used for the empirical p-value.",
    "pool": "Pool query and target peaks before shuffling when `True`.",
    "legend": "Draw a legend when `True`.",
    "ndigit": "Number of digits used in labels.",
    "distance_breaks": "Distance-axis break points in bp.",
    "distance_column": "Column that supplies TSS distances.",
    "signal": "Signal track or matrix used for the profile.",
    "show_sem": "Draw the standard-error band when `True`.",
    "sort_rows": "Sort heatmap rows by mean signal when `True`.",
    "gr": "Interval table with `chrom`, `start`, `end` and optional `strand`/`tss`.",
    "origin_label": "Label drawn at the anchor position (default `TSS`).",
    "label": "Legend label for the line or set.",
    "tagMatrices": "Sequence of tag matrices to compare.",
    "min_overlap": "Minimum overlap size required for an enrichment result.",
    "include_gene_body": "Include gene-body links in addition to promoter links.",
    "min_size": "Minimum gene-set size for GSEA.",
    "max_size": "Maximum gene-set size for GSEA.",
    "permutation_num": "Number of GSEA permutations.",
    "ranking": "Ranked gene mapping/table or ranked gene sequence.",
    "result": "Enrichment result table from `ora`, `enrichr` or `gsea`.",
    "kind": "Plot kind: `'bar'` (default) or `'dot'`.",
}


OVERRIDES: Dict[Tuple[str, str], str] = {
    ("annotate_peak", "peaks"): "Peak file or DataFrame to annotate. Required unless supplied through a compatibility alias.",
    ("annotate_peak", "gtf"): "Gene model to annotate against; mutually exclusive with `TxDb`.",
    ("annotate_peak", "TxDb"): "Assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`; mutually exclusive with `gtf`.",
    ("annotate_peak", "tss_region"): "Promoter window `(upstream, downstream)` in bp; negative upstream, positive downstream.",
    ("annotate_peak", "level"): "`transcript` (default) or `gene`; controls which feature set supplies TSSs.",
    ("annotate_peak", "upstream"): "Convenience alias for the upstream side of `tss_region`.",
    ("annotate_peak", "downstream"): "Convenience alias for the downstream side of `tss_region`.",
    ("annotate_peak", "promoter_bins"): "Promoter bin edges in bp; defaults to `tss_region` plus ChIPseeker's standard bins.",
    ("annotate_peak", "overlap"): "`TSS` assigns features only near the TSS; `all` considers every overlapping feature span.",
    ("annotate_peak", "n_jobs"): "`None`/`1` sequential, `-1` all cores, or an explicit worker count.",
    ("annotate_peak", "kwargs"): "R-style aliases such as `tssRegion`, `genomicAnnotationPriority`, `peak`, `annotation`, `annoDb` and `columns`.",
    ("annotate_db", "annotation"): "Annotated result whose gene columns should be extended.",
    ("annotate_db", "anno_db"): "Gene annotation table/mapping/path keyed by `gene_id`.",
    ("seq2gene", "gtf"): "Gene model used to link peaks to genes; mutually exclusive with `gene_annotation`.",
    ("seq2gene", "gene_annotation"): "Prebuilt `GeneAnnotation`; mutually exclusive with `gtf`.",
    ("seq2gene", "tss_region"): "Promoter window used for the link.",
    ("seq2gene", "flank_distance"): "Flank window used when a peak is not in a promoter.",
    ("seq2gene_links", "gtf"): "Gene model used to build the peak-to-gene link table.",
    ("get_flanking_genes", "model"): "`GeneAnnotation` that supplies the gene models.",
    ("get_flanking_genes", "peaks"): "Peaks whose flanking genes should be found.",
    ("get_genomic_anno_stat", "annotations"): "Annotation labels, usually `ann.detail['annotation']`.",
    ("get_genomic_anno_stat", "downstream_distance"): "Downstream distance used by the `Downstream` label.",
    ("promoter_label", "distance"): "Signed distance from the peak to the TSS.",
    ("promoter_label", "bins"): "Promoter bin edges in bp.",
    ("group_annotation", "annotation"): "Detailed ChIPseeker label to collapse to a major class.",
    ("read_gtf", "features"): "Feature types to keep; `None` reads all records.",
    ("read_peaks", "fmt"): "Input format override; inferred from the suffix when omitted.",
    ("read_peaks", "path"): "BED, narrowPeak/broadPeak/gappedPeak, bigBed, BAM, TSV or CSV file.",
    ("read_peak_file", "as_type"): "`dataframe` (default) or `GRanges`-like table.",
    ("read_peak_file", "kwargs"): "Forwarded to `read_peaks` and the format reader.",
    ("write_table", "df"): "Table to write as TSV.",
    ("write_table", "path"): "Destination TSV path.",
    ("write_peaks", "peaks"): "Peak table to write as BED6+.",
    ("write_peaks", "path"): "Destination BED path.",
    ("filter_peaks", "blacklist"): "Blacklist BED path/DataFrame; any overlap removes the peak.",
    ("filter_peaks", "drop_mito"): "Remove `chrM`/`MT` peaks.",
    ("filter_peaks", "min_width"): "Keep peaks at least this wide.",
    ("filter_peaks", "max_width"): "Keep peaks at most this wide.",
    ("filter_peaks", "min_score"): "Keep peaks with a score at or above this value.",
    ("filter_peaks", "chroms"): "Whitelist of chromosomes to keep.",
    ("filter_peaks", "exclude_chroms"): "Chromosomes to drop after the whitelist step.",
    ("call_peaks", "bigwig"): "Signal track to threshold.",
    ("call_peaks", "threshold"): "Minimum signal for a base to belong to a peak.",
    ("call_peaks", "max_gap"): "Merge runs separated by at most this many bases.",
    ("call_peaks", "min_width"): "Discard merged runs shorter than this width.",
    ("call_peaks", "sample"): "Optional downsampling of the signal track for quick tests.",
    ("call_peaks_macs", "treatment"): "Treatment BAM/CRAM file.",
    ("call_peaks_macs", "control"): "Optional control BAM/CRAM file.",
    ("call_peaks_macs", "genome"): "Genome name understood by MACS.",
    ("call_peaks_macs", "backend"): "`macs2` or `macs3`.",
    ("call_peaks_macs", "extra_args"): "Extra MACS command-line arguments.",
    ("combine_cs_anno", "annotations"): "Two or more `PeakAnnotation` objects to combine.",
    ("drop_anno", "annotation"): "`PeakAnnotation` or annotation table to filter.",
    ("drop_anno", "distance_to_tss_cutoff"): "Keep peaks with `abs(distance_to_tss) < cutoff`.",
    ("get_anno_stat", "annotation"): "`PeakAnnotation` or annotation table.",
    ("detail_genomic_annotation", "annotation"): "`PeakAnnotation` or annotation table.",
    ("annotation_summary", "annotation"): "`PeakAnnotation` or a DataFrame with `chrom`, `start`, `end` and `annotation`.",
    ("peaks_per_gene", "annotation"): "`PeakAnnotation` or annotation table with gene columns.",
    ("top_genes", "n"): "Number of top genes to return.",
    ("gene_type_summary", "top"): "Number of biotypes kept before folding the rest into `Other`.",
    ("distance_summary", "thresholds"): "Distance cutoffs in bp used for the percentage rows.",
    ("signal_groups", "tag"): "Tag matrix whose rows are assigned to groups.",
    ("signal_groups", "n_groups"): "Number of balanced groups; capped at the number of rows.",
    ("tss_enrichment_by_quartile", "tag"): "Tag matrix used to compute centre/flank coverage.",
    ("tss_enrichment_by_quartile", "n_groups"): "Number of signal groups.",
    ("tss_enrichment_by_quartile", "center_fraction"): "Fraction of bins in the central window (default 5 %).",
    ("tss_enrichment_by_quartile", "flank_fraction"): "Fraction of bins in the outer flanks (default 10 %).",
    ("get_promoters", "TxDb"): "Gene model used to build promoters.",
    ("get_promoters", "by"): "`gene` (default) or `transcript`.",
    ("get_promoters", "type"): "Anchor type; promoters use `start_site`.",
    ("get_promoters", "annotation"): "Explicit `GeneAnnotation`, GTF path or DataFrame.",
    ("get_promoters", "gtf"): "GTF/GFF path used when `TxDb`/`annotation` is omitted.",
    ("get_bio_region", "type"): "`start_site`, `end_site` or `body`.",
    ("get_bio_region", "by"): "`gene` or `transcript`.",
    ("make_bio_region_from_granges", "gr"): "Interval table with `chrom`, `start`, `end` and optionally `strand`/`tss`.",
    ("make_bio_region_from_granges", "by"): "Feature level label stored in the output.",
    ("make_bio_region_from_granges", "type"): "`start_site`, `end_site` or `body`.",
    ("get_tag_matrix", "peak"): "Peak file/DataFrame, or a bigWig path when `bigwig` is omitted.",
    ("get_tag_matrix", "windows"): "Pre-built windows; a sequence returns a list of matrices.",
    ("get_tag_matrix", "TxDb"): "Gene model used when `windows` is omitted.",
    ("get_tag_matrix", "type"): "Region type used when building windows.",
    ("get_tag_matrix", "by"): "`gene` or `transcript`.",
    ("get_tag_matrix", "weightCol"): "Peak column used as weight instead of 1 per peak.",
    ("get_tag_matrix", "nbin"): "Bins per window.",
    ("get_tag_matrix", "ignore_strand"): "Ignore window strand.",
    ("get_tag_matrix", "bigwig"): "Use signal coverage instead of peak coverage.",
    ("get_tag_matrix", "bam"): "Use BAM/CRAM coverage instead of peak coverage.",
    ("get_tag_matrix", "value_type"): "bigWig summary statistic.",
    ("get_tag_matrix", "min_mapping_quality"): "Minimum MAPQ for BAM coverage.",
    ("get_tag_matrix", "exclude_duplicates"): "Ignore duplicate BAM reads.",
    ("get_tag_matrix", "n_jobs"): "Parallel worker count for bigWig/BAM extraction.",
    ("get_tag_matrix_from_bigwig", "windows"): "Window table (one row per region).",
    ("get_tag_matrix_from_bigwig", "value_type"): "`sum`, `mean`, `max`, `min`, `coverage` or `std`.",
    ("get_tag_matrix_from_bam", "windows"): "Window table (one row per region).",
    ("get_tag_matrix_from_bam", "value_type"): "`mean` (default) or another supported summary.",
    ("extract_tss_matrix", "genes"): "TSS feature table; use `ann.gene_annotation.features('transcript')`.",
    ("extract_tss_matrix", "window"): "Symmetric window size in bp when `upstream`/`downstream` are omitted.",
    ("extract_tss_matrix", "min_value"): "Values below this are treated as empty.",
    ("extract_tss_matrix", "drop_empty"): "Drop features with no signal.",
    ("extract_tss_matrix", "return_mask"): "Return `(matrix, mask)` instead of only the matrix.",
    ("extract_tss_matrix", "dtype"): "Output dtype (default `float64`).",
    ("extract_tss_matrix_from_signal", "signal"): "bigWig path or precomputed signal table.",
    ("extract_tss_matrix_from_signal", "chrom_sizes"): "Chromosome lengths used to clip windows.",
    ("plot_avg_prof", "tagMatrix"): "Tag matrix to average.",
    ("plot_avg_prof", "xlim"): "X-axis limits in bp relative to the anchor.",
    ("plot_avg_prof", "conf"): "Bootstrap confidence level; `None` hides the band.",
    ("plot_avg_prof", "facet"): "`none`, `row` or `column`.",
    ("plot_avg_prof", "origin_label"): "Text drawn at the anchor (default `TSS`).",
    ("plot_avg_prof", "color"): "Line colour.",
    ("plot_avg_prof", "label"): "Legend label.",
    ("plotAvgProf2", "peak"): "Peaks used to build the matrix.",
    ("plotAvgProf2", "upstream"): "Upstream window in bp.",
    ("plotAvgProf2", "downstream"): "Downstream window in bp.",
    ("plot_peak_prof", "tagMatrix"): "Optional precomputed matrix; otherwise built from `peak`.",
    ("plot_peak_prof", "windows_name"): "Names for multiple window sets.",
    ("plot_peak_prof_multi_windows", "windows"): "Sequence of window tables to compare.",
    ("plot_peak_prof_multi_windows", "windows_name"): "Labels for the window sets.",
    ("plot_multi_prof", "tagMatrices"): "Sequence of matrices to compare.",
    ("tag_heatmap", "tagMatrix"): "Matrix to draw.",
    ("tag_heatmap", "palette"): "Colormap name.",
    ("peak_heatmap", "peak"): "Peaks used to build the matrix.",
    ("peak_heatmap", "windows"): "Optional pre-built windows.",
    ("peak_heatmap_multiple_sets", "peaks"): "Peak sets to compare; one heatmap row block per set.",
    ("peak_heatmap_multiple_sets", "labels"): "Labels for the peak sets.",
    ("peak_profile_heatmap", "peak"): "Peaks used for the profile/heatmap pair.",
    ("plot_tss_profile", "matrix"): "Features x bins signal matrix.",
    ("plot_tss_profile", "window"): "Symmetric window size when `upstream`/`downstream` are omitted.",
    ("plot_tss_profile", "show_sem"): "Draw the standard-error band.",
    ("plot_tss_heatmap", "matrix"): "Features x bins signal matrix.",
    ("plot_tss_heatmap", "sort_rows"): "Sort rows by mean signal.",
    ("plot_anno_pie", "annotation"): "`PeakAnnotation` or annotation table.",
    ("plot_anno_bar", "annotation"): "`PeakAnnotation` or annotation table.",
    ("plot_dist_to_tss", "annotation"): "`PeakAnnotation` or annotation table.",
    ("peak_qc", "peaks"): "Peak set to summarize.",
    ("peak_qc", "blacklist"): "Blacklist intervals counted as flagged peaks.",
    ("peak_qc", "chrom_sizes"): "Chromosome lengths for fraction calculations.",
    ("frip", "peaks"): "Peaks used in the FRiP numerator.",
    ("frip", "bigwig"): "Signal track alternative to `bam`.",
    ("frip", "bam"): "Alignment file alternative to `bigwig`.",
    ("frip", "exclude_duplicates"): "Ignore duplicate reads.",
    ("tss_enrichment", "genes"): "TSS feature table.",
    ("tss_enrichment", "sample"): "Maximum number of TSSs sampled.",
    ("tss_enrichment", "center_window"): "Half-width of the centre window.",
    ("tss_enrichment", "flank_window"): "Half-width of the flank window.",
    ("qc_summary", "genes"): "TSS table used for TSS enrichment.",
    ("qc_summary", "tss_upstream"): "Upstream TSS window for QC.",
    ("qc_summary", "tss_downstream"): "Downstream TSS window for QC.",
    ("qc_summary", "tss_sample"): "Maximum number of TSS features sampled.",
    ("qc_summary", "exclude_duplicates"): "Ignore duplicate reads in FRiP.",
    ("read_gmt", "path"): "GMT file path.",
    ("ora", "gene_list"): "Query genes.",
    ("ora", "gene_sets"): "GMT path, library name or mapping of term to genes.",
    ("ora", "background"): "Background genes; defaults to the union of all sets.",
    ("ora", "p_adjust"): "`BH`, `Bonferroni` or `none`.",
    ("ora", "min_overlap"): "Minimum overlap between query and term.",
    ("enrichr", "gene_list"): "Query genes.",
    ("enrichr", "gene_sets"): "Enrichr library name or names.",
    ("enrichr", "organism"): "Enrichr organism.",
    ("enrichr", "outdir"): "Directory for gseapy output.",
    ("enrichr", "background"): "Optional background genes.",
    ("enrichr", "p_adjust"): "Optional adjustment method passed to gseapy.",
    ("enrich_go", "gene_list"): "Query genes.",
    ("enrich_go", "organism"): "Enrichr organism.",
    ("enrich_go", "outdir"): "Output directory.",
    ("enrich_kegg", "gene_list"): "Query genes.",
    ("enrich_kegg", "organism"): "Enrichr organism.",
    ("enrich_kegg", "outdir"): "Output directory.",
    ("enrich_peak_genes", "annotation"): "`PeakAnnotation` or annotation table.",
    ("enrich_peak_genes", "upstream"): "Upstream promoter window for peak-to-gene links.",
    ("enrich_peak_genes", "downstream"): "Downstream promoter window for peak-to-gene links.",
    ("enrich_peak_genes", "include_gene_body"): "Include gene-body links.",
    ("enrich_peak_genes", "gene_column"): "Gene identifier column.",
    ("enrich_peak_genes", "kwargs"): "Forwarded to `enrichr`.",
    ("gsea", "gene_list"): "Ranked mapping, two-column table or ranked gene sequence.",
    ("gsea", "gene_sets"): "Library name, GMT path or mapping of term to genes.",
    ("gsea", "min_size"): "Minimum gene-set size.",
    ("gsea", "max_size"): "Maximum gene-set size.",
    ("gsea", "permutation_num"): "Number of permutations.",
    ("gsea", "seed"): "Random seed.",
    ("gsea", "outdir"): "Output directory.",
    ("gsea_peak_genes", "annotation"): "`PeakAnnotation` or annotation table.",
    ("gsea_peak_genes", "ranking"): "Gene ranking used for GSEA.",
    ("gsea_peak_genes", "gene_column"): "Gene identifier column.",
    ("plot_enrichment", "result"): "ORA/Enrichr/GSEA result table.",
    ("plot_enrichment", "top"): "Number of terms to show.",
    ("plot_enrichment", "kind"): "`bar` or `dot`.",
    ("detect_assembly", "peaks"): "Peak table used to infer the assembly.",
    ("detect_assembly", "chrom_sizes"): "Chromosome sizes used for a stronger match.",
    ("detect_assembly", "bigwig"): "bigWig header used to infer chromosome names/sizes.",
    ("detect_assembly", "tolerance"): "Relative size tolerance for the match.",
    ("load_genome", "assembly"): "Assembly name; `None` loads the default supported bundle.",
    ("load_genome", "refresh"): "Rebuild/re-download even when cached.",
    ("load_annotation", "spec"): "Assembly name or `ensembl:`/`gencode:`/`url:` spec.",
    ("download_genome", "assembly"): "Assembly to download.",
    ("download_gtf", "spec"): "Ensembl/GENCODE/URL spec.",
    ("download_gtf", "dest"): "Destination GTF path.",
    ("clear_genome_cache", "assembly"): "Assembly to clear; `None` clears all user-cached genomes.",
    ("genome_provenance", "assembly"): "Assembly whose provenance record should be returned.",
    ("canonical_chrom", "value"): "Chromosome name to normalize.",
    ("chromosome_style", "names"): "Chromosome names to classify.",
    ("harmonize_frame", "frame"): "DataFrame to harmonize.",
    ("harmonize_frame", "target"): "Chromosome names to match.",
    ("harmonize_frame", "column"): "Chromosome column name.",
    ("harmonize_frame", "copy"): "Copy the frame before modifying.",
    ("overlap", "Sets"): "Peak sets or mapping of set names to peaks.",
    ("shuffle", "peak"): "Peaks to shuffle.",
    ("shuffle", "chrom_sizes"): "Chromosome lengths used to draw new positions.",
    ("shuffle", "TxDb"): "Gene model alternative to `chrom_sizes`.",
    ("shuffle", "n"): "Number of shuffled peak sets to generate.",
    ("vennplot", "Sets"): "Peak sets or mapping of names to peaks.",
    ("vennplot", "labels"): "Set labels.",
    ("vennplot", "colors"): "Set colours.",
    ("vennplot_peakfile", "peakfiles"): "Peak files to compare.",
    ("vennpie", "annotation"): "`PeakAnnotation` whose detailed labels are summarized.",
    ("vennpie", "r"): "Inner radius of the nested pie.",
    ("upsetplot", "annotation"): "`PeakAnnotation` whose annotation flags are compared.",
    ("upsetplot", "order_by"): "`freq` or `degree`.",
    ("upsetplot", "vennpie"): "Add a nested pie alongside the UpSet plot.",
    ("upsetplot", "n_intersections"): "Maximum intersections shown.",
    ("covplot", "peak"): "Peaks or peak file to draw.",
    ("covplot", "weightCol"): "Optional weight column.",
    ("covplot", "chrs"): "Chromosomes to draw.",
    ("covplot", "lower"): "Lower y bound.",
    ("enrich_peak_overlap", "queryPeak"): "Query peak set.",
    ("enrich_peak_overlap", "targetPeak"): "Target peak set.",
    ("enrich_peak_overlap", "chrom_sizes"): "Chromosome lengths for shuffling.",
    ("enrich_peak_overlap", "nShuffle"): "Number of shuffles.",
    ("enrich_peak_overlap", "pAdjustMethod"): "`BH`, `Bonferroni` or `none`.",
    ("enrich_peak_overlap", "pool"): "Pool the two sets before shuffling.",
    ("enrich_anno_overlap", "queryPeak"): "Query peaks.",
    ("enrich_anno_overlap", "targetPeak"): "Target peaks.",
    ("enrich_anno_overlap", "TxDb"): "Gene model used for annotation.",
    ("enrich_anno_overlap", "gtf"): "GTF/GFF path alternative to `TxDb`.",
    ("enrich_anno_overlap", "tss_region"): "Promoter window used for gene links.",
    ("get_geo_species", "gsminfo"): "Optional GEO snapshot path/DataFrame.",
    ("get_geo_genome_version", "gsminfo"): "Optional GEO snapshot path/DataFrame.",
    ("get_geo_info", "genome"): "Genome version to subset.",
    ("get_geo_info", "simplify"): "Keep only key columns.",
    ("download_geo_bed_files", "genome"): "Genome version to download.",
    ("download_geo_bed_files", "destDir"): "Destination directory.",
    ("download_gsm_bed_files", "GSM"): "GSM accession(s) to download.",
    ("download_gsm_bed_files", "destDir"): "Destination directory.",
    ("get_sample_files", "directory"): "Directory to scan for sample peaks.",
    ("build_report", "peaks"): "Peaks included in the report.",
    ("build_report", "gtf"): "Gene model for annotation.",
    ("build_report", "TxDb"): "Assembly or gene model alternative to `gtf`.",
    ("build_report", "bigwig"): "Add signal TSS profile/heatmap pages.",
    ("build_report", "bam"): "Add BAM-based FRiP/TSS QC pages.",
    ("build_report", "blacklist"): "Blacklist used by the QC pages.",
    ("build_report", "output"): "Destination PDF path.",
    ("build_report", "level"): "`transcript` or `gene` annotation level.",
    ("build_report", "tss_region"): "Promoter window for annotation.",
    ("build_report", "upstream"): "Asymmetric signal window upstream.",
    ("build_report", "downstream"): "Asymmetric signal window downstream.",
    ("build_report", "sample"): "Number of transcripts used for signal pages.",
    ("build_report", "seed"): "Sampling seed for signal pages.",
    ("build_report", "title"): "Report title.",
    ("build_report", "include_tag_matrix"): "Include tag-matrix pages.",
    ("build_report", "qc"): "Include QC/FRiP pages.",
    ("build_report", "n_jobs"): "Parallel worker count for annotation/QC.",
}


FUNCTION_NOTES: Dict[str, str] = {
    "annotate_peak": (
        "`annotatePeak` is a bare alias. R-style keyword spellings such as "
        "`tssRegion`, `genomicAnnotationPriority`, `peak`, `annotation`, "
        "`annoDb` and `columns` are accepted through `**kwargs`."
    ),
    "get_tag_matrix": (
        "When `windows` is a sequence, a list of `TagMatrix` objects is returned. "
        "With `bigwig=`/`bam=` the matrix contains signal coverage; otherwise it "
        "contains peak coverage."
    ),
    "extract_tss_matrix": (
        "`upstream`/`downstream` define an asymmetric window; otherwise "
        "`window` is used on both sides."
    ),
    "signal_groups": (
        "Genes without an assigned peak receive score 0, so they stay in the "
        "lowest group instead of being filled with the median."
    ),
    "tss_enrichment_by_quartile": (
        "Enrichment is mean coverage in the central bins divided by mean coverage "
        "in the outer flanks. The centre is derived from the tag-matrix metadata "
        "when available."
    ),
    "call_peaks": (
        "Returns BED6+ columns with `signalValue` (mean signal) and `peak` "
        "(offset of the highest base)."
    ),
    "build_report": (
        "The PDF includes title/summary, annotation pie/bar, TSS distance, "
        "coverage, tag-matrix pages, optional signal pages and optional QC pages."
    ),
}


def _signature_text(func: Any) -> str:
    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):
        return "(...)"
    return f"{func.__name__}{signature}"


def _clean_type(annotation: Any) -> str:
    text = str(annotation)
    if text in {"<class 'inspect._empty'>", "inspect._empty", "empty"}:
        return "—"
    text = text.strip("'\"")
    for old, new in (
        ("pd.DataFrame", "DataFrame"),
        ("np.ndarray", "ndarray"),
        ("typing.", ""),
        ("NoneType", "None"),
    ):
        text = text.replace(old, new)
    return text


def _clean_default(value: Any) -> str:
    if value is inspect._empty:
        return "—"
    if value is None:
        return "`None`"
    if isinstance(value, str):
        return f"`{value}`"
    text = repr(value)
    if len(text) > 60:
        text = text[:57] + "..."
    return f"`{text}`"


def _escape(text: str) -> str:
    return text.replace("|", "\\|")


def _summary(func: Any) -> str:
    doc = inspect.getdoc(func) or ""
    if not doc:
        return "No description available."
    first = doc.split("\n\n", 1)[0].strip().splitlines()[0].strip()
    first = re.sub(r"``([^`]+)``", r"`\1`", first)
    return first


def _param_description(name: str, func_name: str) -> str:
    return OVERRIDES.get((func_name, name)) or PARAM_DOCS.get(name, "")


def _render_function(name: str) -> str:
    func = getattr(pcs, name)
    signature = inspect.signature(func)
    lines = [f"### `{name}`", "", _summary(func), ""]
    note = FUNCTION_NOTES.get(name)
    if note:
        lines += ["!!! note", f"    {note}", ""]
    lines += ["```python", _signature_text(func), "```", ""]
    lines += ["| Parameter | Type | Default | Description |", "|---|---|---|---|"]
    for param in signature.parameters.values():
        if param.kind is inspect.Parameter.VAR_POSITIONAL:
            continue
        description = _param_description(param.name, name)
        if not description:
            description = "*(description not yet written)*"
        ptype = _clean_type(param.annotation)
        default = "keyword" if param.kind is inspect.Parameter.VAR_KEYWORD else _clean_default(param.default)
        label = f"**{param.name}**" if param.default is inspect._empty and param.kind is not inspect.Parameter.VAR_KEYWORD else f"`{param.name}`"
        lines.append(
            f"| {label} | {_escape(ptype)} | {_escape(default)} | {_escape(description)} |"
        )
    lines.append("")
    return "\n".join(lines)


RENDERED_NAMES = {name for _title, _intro, names in CATEGORIES for name in names}


def _aliases() -> Dict[str, List[str]]:
    by_id: Dict[int, List[str]] = {}
    for name in pcs.__all__:
        obj = getattr(pcs, name, None)
        if not callable(obj) or inspect.isclass(obj):
            continue
        by_id.setdefault(id(obj), []).append(name)
    result: Dict[str, List[str]] = {}
    for names in by_id.values():
        names = sorted(
            names,
            key=lambda item: (item not in RENDERED_NAMES, item.lower() != item, item),
        )
        canonical = names[0]
        if len(names) > 1:
            result[canonical] = names[1:]
    return result


def main() -> int:
    lines = [
        "# Function reference",
        "",
        "Every public top-level function, with every parameter. The tables are",
        "generated from the installed signatures, so required parameters, defaults",
        "and type hints always match the code.",
        "",
        "!!! tip \"How to read this page\"",
        "",
        "    - **Required** parameters have no default (`—`).",
        "    - `None` usually means *automatic* or *not supplied*.",
        "    - `n_jobs=None`/`1` is sequential; `-1` uses all cores.",
        "    - Peak coordinates are 0-based and half-open, as in BED.",
        "    - `TxDb=`, `gtf=` and `annotation=` are alternative ways to supply the",
        "      same gene model; use one of them.",
        "",
        "## Parameter conventions",
        "",
        "| Name | Meaning |",
        "|---|---|",
        "| `peaks` / `peak` | A BED/narrowPeak-style path, a peak DataFrame, or (where noted) a sequence. |",
        "| `TxDb` / `gtf` / `annotation` | Gene model: assembly name, GTF/GFF path, DataFrame or `GeneAnnotation`. |",
        "| `level` | `transcript` (default) or `gene`; controls which TSS set is used. |",
        "| `by` | `gene` or `transcript`; controls the window rows. |",
        "| `type` | `start_site`, `end_site` or `body`. |",
        "| `upstream` / `downstream` | Bases on each side of the anchor; negative upstream. |",
        "| `nbin` / `n_bins` | Number of bins across a window. |",
        "| `value_type` | bigWig summary: `sum`, `mean`, `max`, `min`, `coverage` or `std`. |",
        "| `n_jobs` | `None`/`1` sequential, `-1` all cores, `n > 1` explicit workers. |",
        "",
        "## Aliases",
        "",
        "ChIPseeker-style aliases share the same signature and parameters as their",
        "snake_case counterpart:",
        "",
    ]
    aliases = _aliases()
    lines += ["| Alias | Canonical function |", "|---|---|"]
    for canonical in sorted(aliases):
        for alias in aliases[canonical]:
            lines.append(f"| `{alias}` | [`{canonical}`](#{canonical.lower()}) |")
    lines.append("")

    for title, intro, names in CATEGORIES:
        lines += [f"## {title}", "", intro, ""]
        for name in names:
            if not hasattr(pcs, name):
                raise SystemExit(f"missing public function: {name}")
            lines.append(_render_function(name))

    output = ROOT / "docs" / "reference" / "functions.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"wrote {output.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
