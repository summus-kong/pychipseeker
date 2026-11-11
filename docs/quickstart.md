# Quick start

```python
import pychipseeker as pcs

# Prebuilt genome tables: hg38, hg19, mm10, mm9 (downloaded once)
ann = pcs.annotatePeak("peaks.bed", TxDb="mm10")
print(ann.annoStat)

# A GTF/GFF file, DataFrame or GeneAnnotation object works too
ann = pcs.annotate_peak("peaks.bed", "genes.gtf", tss_region=(-3000, 3000))

# Large peak sets: annotate on all cores (results match the sequential run)
ann = pcs.annotate_peak("peaks.bed", TxDb="mm10", n_jobs=-1)

# ChIPseeker-style result table and plots
df = ann.to_chipseeker()
fig, ax = ann.plot_anno_pie()
fig, ax = ann.plot_anno_bar()
fig, ax = ann.plot_dist_to_tss()
```

## Tag matrices

```python
promoters = pcs.getPromoters(TxDb=ann.gene_annotation,
                             upstream=3000, downstream=2000)
tag = pcs.getTagMatrix(ann.peaks, windows=promoters, nbin=100)
fig, ax = pcs.plotAvgProf(tag, conf=0.95)
fig, ax = pcs.tagHeatmap(tag)
```

Pass `n_jobs=-1` to `getTagMatrix`/`getPromoters`-based pipelines to process
large window sets in parallel.  Peak-coverage matrices are vectorised per
chromosome in 0.8 and only switch to processes for very large window sets; see
[Performance](performance.md).

## bigWig signal around TSS

The `bigwig` extra is required.

```python
matrix = pcs.extract_tss_matrix(
    "signal.bw",
    ann.gene_annotation.features("transcript"),
    upstream=2000,
    downstream=4000,
    n_bins=100,
    n_jobs=-1,
)
fig, ax = pcs.plot_tss_profile(matrix, upstream=2000, downstream=4000)
fig, ax = pcs.plot_tss_heatmap(matrix, upstream=2000, downstream=4000)
```

## Quality control

```python
qc = pcs.qc_summary(
    ann.peaks,
    bam="alignments.bam",          # or bigwig="signal.bw"
    genes=ann.gene_annotation.features("transcript"),
    blacklist="blacklist.bed",
    genome="mm10",
)
print(qc)
```

## One PDF report

```bash
pychipseeker report \
    --peaks peaks.narrowPeak \
    --txdb mm10 \
    --bigwig signal.bw \
    --output report.pdf
```
