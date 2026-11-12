# Analysis recipes

The `analysis` helpers turn an annotation table into small, plot-ready
DataFrames. They need only NumPy and pandas, accept either a `PeakAnnotation`
or a DataFrame, and never mutate the input.

```python
import pychipseeker as pcs

ann = pcs.annotate_peak("peaks.bed", TxDb="mm10")
df = ann.to_dataframe()
```

For signal-aware summaries, make sure the table has a `score` column. A
`PeakAnnotation` gets one automatically; a plain DataFrame needs the merge:

```python
df = df.merge(ann.peaks[["peak_id", "score"]], on="peak_id", how="left")
```

| Helper | Returns |
|---|---|
| `annotation_summary` | counts, widths and optional median signal per annotation class |
| `peaks_per_gene` | one row per gene with peak count and peak-level medians |
| `top_genes` | the `n` genes with the most assigned peaks |
| `gene_type_summary` | counts and percentages per gene biotype, tail folded into `Other` |
| `distance_summary` | promoter fraction and distance quantiles as a Metric/Value table |
| `signal_groups` | balanced signal groups for tag-matrix rows |
| `tss_enrichment_by_quartile` | centre-to-flank enrichment per signal group |

## Annotation composition

```python
summary = pcs.annotation_summary(ann)
print(summary)
```

Columns: `annotation`, `count`, `percent`, `median_width_bp`, `mean_width_bp`,
`total_bp`, and `median_signal` when a score is available. `percent` is the
share of all peaks, so the column sums to 100.

```python
import matplotlib.pyplot as plt

top = summary.head(8).iloc[::-1]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.barh(top["annotation"], top["count"], color="#0072B2")
ax.set_xlabel("Peaks")
fig.savefig("annotation-summary.pdf", bbox_inches="tight")
```

## Peaks per gene and top genes

```python
per_gene = pcs.peaks_per_gene(ann)
top = pcs.top_genes(ann, n=12)
```

`peaks_per_gene` returns `gene_id`, `gene_name`, `n_peaks`, `total_bp`,
`median_width_bp`, and `median_signal` / `median_distance_to_tss` when those
columns exist. Genes are sorted by peak count, then signal.

`top_genes` adds a `rank` column and keeps the same columns:

| rank | gene_id | gene_name | n_peaks | median_width_bp | median_signal |
|---:|---|---|---:|---:|---:|
| 1 | ENSMUSG… | Bcor | 8 | … | … |
| 2 | ENSMUSG… | Cbx8 | 6 | … | … |

To count genes by burden:

```python
counts = per_gene["n_peaks"].value_counts().sort_index()
print(counts.head())
```

## Gene biotypes

```python
biotypes = pcs.gene_type_summary(ann, top=8)
```

Columns: `gene_type`, `count`, `percent`. The eight largest biotypes are kept
and everything else is folded into one `Other` row; percentages still refer to
the complete peak set, so the column sums to about 100.

## Distance to TSS

```python
metrics = pcs.distance_summary(ann, thresholds=(1000, 10000, 100000))
print(metrics)
```

The result is a two-column `Metric` / `Value` table:

| Metric | Value |
|---|---|
| Peaks | 11291 |
| Median distance (bp) | … |
| 90th percentile (bp) | … |
| At TSS (distance = 0) | … |
| Within 1 kb | 92.0 |
| Within 10 kb | 96.5 |
| Within 100 kb | 99.5 |

Values for the “Within …” rows and “At TSS” are percentages. Distances use the
absolute signed distance reported by the annotation.

## Signal groups

```python
groups = pcs.signal_groups(tag, ann, n_groups=4)
print(groups.head())
```

Columns: `row` (tag-matrix row index), `group` (0-based), `label`
(`Q1 lowest` … `Q4 highest`) and `score`.

How scores are assigned:

1. Each promoter inherits the strongest peak score annotated to its gene.
2. Genes without an assigned peak receive a score of **zero**, so they stay in
   the lowest group instead of being silently dropped or filled with the
   median.
3. Groups are balanced by rank, which avoids empty bins when many scores tie.
4. If row metadata or scores are unavailable, row-mean coverage is used
   instead, so synthetic matrices still work.

## TSS enrichment by quartile

```python
enrichment = pcs.tss_enrichment_by_quartile(tag, ann, n_groups=4)
print(enrichment)
```

Columns: `quartile`, `n`, `median_signal`, `tss_enrichment`, `max_coverage`.

Enrichment is the mean coverage in the central 5 % of bins divided by the mean
coverage in the outer 10 % of bins. A value of 1 means the centre looks like the
flanks; larger values mean a sharper TSS peak.

```python
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.bar(enrichment["quartile"], enrichment["tss_enrichment"], color="#1F78B4")
ax.set_ylabel("TSS enrichment (centre / flank)")
fig.savefig("tss-enrichment.pdf", bbox_inches="tight")
```

!!! tip "Asymmetric windows"

    `tss_enrichment_by_quartile` finds the centre from the tag matrix metadata
    when available, so asymmetric upstream/downstream windows are handled
    correctly. Override the window fractions with `center_fraction=` and
    `flank_fraction=` if your bins need a different definition.

## A complete two-panel example

```python
import matplotlib.pyplot as plt
import pychipseeker as pcs

ann = pcs.annotate_peak("peaks.bed", TxDb="mm10")
df = ann.to_dataframe()

promoters = pcs.get_promoters(
    TxDb=ann.gene_annotation, upstream=3000, downstream=3000
)
tag = pcs.get_tag_matrix(ann.peaks, windows=promoters, nbin=120)

summary = pcs.tss_enrichment_by_quartile(tag, df, n_groups=4)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
pcs.top_genes(df, n=10).iloc[::-1].plot.barh(
    x="gene_name", y="n_peaks", ax=axes[0], legend=False, color="#0072B2"
)
axes[1].bar(summary["quartile"], summary["tss_enrichment"], color="#1F78B4")
axes[1].set_ylabel("TSS enrichment")
fig.savefig("gene-summary.pdf", bbox_inches="tight")
```

## Validation

All helpers raise clear errors for missing columns (`chrom`, `start`, `end`,
`annotation`; `gene_id`/`gene_name`; `gene_type`; `distance_to_tss`) and for
invalid arguments such as `n <= 0` or `n_groups < 2`. Empty or invalid peak
tables are rejected rather than returning misleading zeros.

The full test suite is in `tests/test_analysis.py`.
