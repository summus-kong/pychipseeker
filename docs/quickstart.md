# Quick start

This page walks through a complete small analysis: annotate peaks, inspect the
result, draw figures, build a promoter matrix, summarise the peak set and write
a PDF report.

## 1. Annotate

```python
import pychipseeker as pcs

ann = pcs.annotate_peak("peaks.bed", TxDb="mm10")
```

`TxDb="mm10"` downloads the compact genome bundle once and caches it. You can
also pass a GTF/GFF path, a `GeneAnnotation`, or another assembly name:

```python
ann = pcs.annotate_peak("peaks.bed", "genes.gtf", tss_region=(-3000, 3000))
ann = pcs.annotate_peak("peaks.bed", TxDb="hg38")
ann = pcs.annotate_peak("peaks.bed", TxDb="auto")   # infer from the peaks
```

Large peak sets can use all cores; results are identical to the sequential run:

```python
ann = pcs.annotate_peak("peaks.bed", TxDb="mm10", n_jobs=-1)
```

## 2. Inspect the result

```python
print(ann.annoStat)          # Feature / Frequency / Count
df = ann.to_dataframe()      # one row per annotated peak
print(df.head())
```

Useful columns include `annotation`, `gene_name`, `gene_type`,
`distance_to_tss`, `nearest_gene_name` and `nearest_distance`. The full list is
in the [API overview](api.md#output-columns).

ChIPseeker-style names are one call away:

```python
cs = ann.to_chipseeker()
print(cs[["peakChr", "geneName", "distanceToTSS"]].head())
```

## 3. Draw and export figures

```python
import matplotlib.pyplot as plt

fig, ax = ann.plot_anno_pie()
fig.savefig("annotation-pie.pdf", bbox_inches="tight")   # vector
fig.savefig("annotation-pie.svg", bbox_inches="tight")   # editable text
fig.savefig("annotation-pie.png", dpi=300, bbox_inches="tight")
plt.close(fig)

fig, ax = ann.plot_anno_bar()
fig, ax = ann.plot_dist_to_tss()
```

Every plotting helper returns a matplotlib `(fig, ax)` pair, so you can keep
editing titles, colours and layout before saving. PDF, SVG and EPS outputs are
vector files.

## 4. Build a promoter matrix

Tag matrices use **peak coverage** by default and do not need a bigWig:

```python
promoters = pcs.get_promoters(
    TxDb=ann.gene_annotation, upstream=3000, downstream=3000
)
tag = pcs.get_tag_matrix(ann.peaks, windows=promoters, nbin=100)

fig, ax = pcs.plotAvgProf(tag)
fig, ax = pcs.tagHeatmap(tag)
```

With the optional `bigwig` extra, use real signal instead:

```python
fig, ax = ann.plot_tss_profile(
    "signal.bw", upstream=2000, downstream=4000, n_bins=100
)
fig, ax = ann.plot_tss_heatmap(
    "signal.bw", upstream=2000, downstream=4000, n_bins=100
)
```

!!! tip "Which matrix should I use?"

    Peak-coverage matrices answer “where do my called peaks sit relative to a
    TSS?”. bigWig matrices answer “how strong is the signal around each TSS?”.
    Both are valid; the first works without any signal file.

## 5. Summarise the peak set

```python
df = ann.to_dataframe()
peaks = df.merge(ann.peaks[["peak_id", "score"]], on="peak_id", how="left")

print(pcs.annotation_summary(df))
print(pcs.top_genes(df, n=12))
print(pcs.gene_type_summary(df, top=8))
print(pcs.distance_summary(df))
```

`signal_groups` and `tss_enrichment_by_quartile` split promoters by peak
strength and compare centre-to-flank coverage:

```python
groups = pcs.signal_groups(tag, peaks, n_groups=4)
print(pcs.tss_enrichment_by_quartile(tag, peaks, n_groups=4))
```

See [Analysis recipes](analysis.md) for the returned columns and plotting
examples.

## 6. Quality control

```python
qc = pcs.qc_summary(
    ann.peaks,
    bigwig="signal.bw",            # or bam="alignments.bam"
    genes=ann.gene_annotation.features("transcript"),
    blacklist="ENCFF356LFX.bed.gz",
    genome="mm10",
)
print(qc)
```

The table includes peak counts and widths, mitochondrial/blacklist fractions,
FRiP and TSS enrichment.

## 7. One PDF report

```bash
pychipseeker report \
    --peaks peaks.narrowPeak \
    --txdb mm10 \
    --bigwig signal.bw \
    --qc \
    --output report.pdf
```

The report contains title/summary pages, annotation pie and bar charts, TSS
distance, coverage, tag-matrix profile and heatmap, optional signal and QC
pages, and methods text. See [Command line](cli.md#report).

## Next steps

- [Figure gallery](gallery.md) — 18 panels, three A3 plates, custom selection.
- [Analysis recipes](analysis.md) — every summary helper with examples.
- [Genomes and TxDb](genomes.md) — cache, offline mirrors, assembly detection.
- [API overview](api.md) — main objects and functions.
