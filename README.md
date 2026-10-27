# pychipseeker

A Python implementation of **ChIPseeker-compatible peak annotation,
summarization and visualization**.

`pychipseeker` follows the R/Bioconductor
[ChIPseeker](https://bioconductor.org/packages/ChIPseeker/) rules for
transcript/gene-level annotation, promoter bins, annotation priority and signed
TSS distances.  The core annotation engine has been checked against
**ChIPseeker 1.48.0**.  The package also provides Python equivalents for the
common tag-matrix, profile, heatmap, overlap and summary helpers.

It is **not** an official ChIPseeker port.  Bioconductor `OrgDb` classes and
the GEO database are not bundled; instead `TxDb=` accepts either a GTF/GFF
path, a `GeneAnnotation` object, or a bundled assembly name (`"hg38"`,
`"hg19"`, `"mm10"`, `"mm9"`); see [PARITY.md](PARITY.md) for the exact
coverage and the remaining gaps.

## Features

- Annotate ChIP-seq / ATAC-seq peaks against a GTF/GFF annotation.
- Assign promoters, 5' UTR, 3' UTR, exons, introns, downstream and distal
  intergenic regions using ChIPseeker's default priority.
- Compute signed distance to the nearest TSS and assign nearest gene/transcript.
- ChIPseeker-style `annoStat`, `plotAnnoPie`, `plotAnnoBar` and
  `plotDistToTSS`.
- Promoter/bioregion construction (`getPromoters`, `getBioRegion`,
  `makeBioRegionFromGranges`).
- Peak-coverage tag matrices (`getTagMatrix`), average profiles, heatmaps and
  multi-window profiles.
- bigWig signal matrices and TSS profiles using `pyBigWig`.
- Peak-set overlap, Venn/UpSet summaries, coverage plots and overlap
  enrichment tests.
- Thin Enrichr wrapper via optional `gseapy`.

## Dependencies

Core installation requires:

- `numpy`
- `pandas`
- `matplotlib`
- `pyBigWig` (declared dependency; used for direct `.bw`/`.bigWig` access)

If a platform has no prebuilt `pyBigWig` wheel, install a C compiler or use a
conda-forge binary.  On Windows, pip ships no native wheel, but the C
extension builds from source with a MinGW toolchain (`NOCURL` disables
remote-file access; local files work, including non-ASCII paths).  R is not
required for any Python functionality.

Optional:

- `gseapy` for Enrichr integration: `pip install -e ".[enrich]"`
- `pytest` for development: `pip install -e ".[dev]"`

## Installation

```bash
pip install -e .
```

## Quick start

```python
import pychipseeker as pcs

# ChIPseeker-style API: annotatePeak is a bare alias of annotate_peak,
# so camelCase keywords and TxDb= work on the same function
ann = pcs.annotatePeak(
    "peaks.bed",
    TxDb="genes.gtf",          # a GTF/GFF path or a GeneAnnotation object
    tssRegion=(-3000, 3000),
    level="transcript",        # ChIPseeker default
)

print(ann)
print(ann.annoStat)            # ChIPseeker annoStat table
print(ann.peakNum)

# Snake_case API is also available
ann2 = pcs.annotate_peak("peaks.bed", "genes.gtf", upstream=3000, downstream=3000)

# Bundled genome lookup — TxDb="assembly", like ChIPseeker's TxDb= packages
ann = pcs.annotatePeak("peaks.bed", TxDb="mm10")   # also hg38, hg19, mm9
# aliases (case-insensitive) work too: "GRCh38", "GRCh37", "GRCm38", "GRCm37"
# bundles are fetched on demand into a local cache on first use; see
# pychipseeker.genomes docs and NOTICE.md for provenance (GENCODE releases).

# ChIPseeker-style column names
chip_df = ann.to_chipseeker()

# Plots
fig, ax = ann.plot_anno_pie()
fig, ax = ann.plot_anno_bar()
fig, ax = ann.plot_dist_to_tss()

# Tag matrix / profile / heatmap
promoters = pcs.getPromoters(TxDb=ann.gene_annotation, upstream=3000, downstream=3000)
tag = pcs.getTagMatrix(peak=ann.peaks, windows=promoters, nbin=100)
fig, ax = pcs.plotAvgProf(tag, conf=0.95)
fig, ax = pcs.tagHeatmap(tag)

# Or let the helper build the matrix and plot it
fig, ax = pcs.plotPeakProf2(
    "peaks.bed", TxDb="genes.gtf", upstream=3000, downstream=3000, nbin=100
)
fig, ax = pcs.peakHeatmap(
    "peaks.bed", TxDb="genes.gtf", upstream=3000, downstream=3000, nbin=100
)
```

`upstream=` and `downstream=` are convenience aliases for an asymmetric
`tss_region`.  `tss_region` controls the promoter window (default ±3 kb).
`downstream_distance` controls the separate **Downstream** annotation category
and defaults to ChIPseeker's 300 bp.

Extra peak metadata columns (for example a cluster label) are preserved and can
be used for grouped summaries and plots:

```python
ann = pcs.annotate_peak(peaks_with_cluster, "genes.gtf")
print(ann.summary(groupby="cluster"))
fig, ax = ann.plot_anno_bar(groupby="cluster")
```

### Real ENCODE example

`examples/call_encode_peaks.py` extracts peaks from a real ENCODE signal
track; `examples/make_encode_plate.py` and `examples/make_encode_report.py`
turn them into a single-page A3 plate and a multi-page PDF report:

```bash
export ENCODE_BIGWIG=/path/to/ENCFF893ZGR.bigWig   # H3K4me3, mm10, ENCSR000CDQ
python examples/call_encode_peaks.py --bigwig "$ENCODE_BIGWIG" \
    --out examples/encode_peaks.bed --threshold 3
python examples/make_encode_plate.py    # encode_figure_plate.pdf + 6 figures
python examples/make_encode_report.py   # 11-page encode_report.pdf
```

The committed `examples/encode_peaks.bed` holds the resulting 11,291 peaks
(92% within 1 kb of a TSS, median width 958 bp).  The track is a
fold-change/signal track, so these peaks are **threshold-called from the
signal**, not ENCODE's own narrowPeak calls.

### OrgDb-like columns

ChIPseeker's `annoDb=` expects an R OrgDb package.  Pure Python does not have
that object, so `pychipseeker` accepts a DataFrame or mapping with a `gene_id`
key instead:

```python
import pandas as pd

db = pd.DataFrame({
    "gene_id": ["G1", "G2"],
    "SYMBOL": ["GeneA", "GeneB"],
    "GENENAME": ["gene A", "gene B"],
})
ann = pcs.annotatePeak("peaks.bed", TxDb="genes.gtf", annoDb=db,
                       columns=["SYMBOL", "GENENAME"])
```

This is an explicit interoperability hook, not a replacement for
TxDb/OrgDb/AnnotationDbi.

## Tag matrices and profiles

`getTagMatrix` follows ChIPseeker's default semantics: it builds a matrix from
**peak coverage** over promoter/bioregion windows.  It does not require a
bigWig file.

```python
promoters = pcs.getPromoters(TxDb=ann.gene_annotation, upstream=1000, downstream=1000)
tag = pcs.getTagMatrix(ann.peaks, windows=promoters, nbin=100)
```

A bigWig-backed matrix is available as an extension when you want to use an
external signal track:

```python
tag = pcs.getTagMatrix(
    bigwig="signal.bw",
    windows=promoters,
    nbin=100,
    value_type="mean",   # sum/mean/max/min/coverage/std
)

# Or use the TSS convenience path
fig, ax = ann.plot_tss_profile("signal.bw", window=2000, n_bins=100)
fig, ax = ann.plot_tss_heatmap("signal.bw", window=2000, n_bins=100)
```

The older bedGraph/precomputed-matrix path remains available:

```python
fig, ax = ann.plot_tss_profile(matrix=matrix, window=2000)
fig, ax = ann.plot_tss_profile(signal="signal.bedgraph", window=2000)
```

Main profile/heatmap functions:

- `getPromoters`, `getBioRegion`, `makeBioRegionFromGranges`
- `getTagMatrix` / `getTagMatrixFromBigWig`
- `plotAvgProf`, `plotAvgProf2`
- `plotPeakProf`, `plotPeakProf2`, `plotPeakProf_MultiWindows`
- `plotMultiProf`
- `tagHeatmap`, `peakHeatmap`, `peakHeatmap_multiple_Sets`
- `peak_Profile_Heatmap`

All plotting helpers return matplotlib `(fig, ax)` (or `(fig, axes)` for
multi-panel plots) and use the same 5'->3' orientation convention as
ChIPseeker.

## Compatibility helpers

Names with dots in R cannot be Python identifiers, so `as.GRanges` and
`vennplot.peakfile` are exported as `asGRanges` / `as_granges` and
`vennplot_peakfile`.

```python
from pychipseeker import (
    readPeakFile, asGRanges, combine_csAnno, dropAnno, getAnnoStat,
    overlap, covplot, vennplot, vennplot_peakfile, vennpie, upsetplot,
    shuffle, enrichPeakOverlap, enrichAnnoOverlap,
)

peaks = readPeakFile("peaks.bed", as_type="GRanges")
gr = asGRanges(ann)
filtered = dropAnno(ann, distanceToTSS_cutoff=10000)
combined = combine_csAnno(ann1, ann2)
stat = getAnnoStat(ann)

fig, ax = covplot(ann.peaks)
fig, ax = vennplot([set_a, set_b], labels=["A", "B"])
fig, ax = vennpie(ann)
fig, ax = upsetplot(ann)

random_peaks = shuffle(ann.peaks, seed=1)
res = enrichPeakOverlap(ann.peaks, [ann.peaks], nShuffle=100, seed=1)
res = enrichAnnoOverlap(ann.peaks, [ann.peaks], TxDb=ann.gene_annotation)
```

`overlap` accepts generic iterables, DataFrames of peaks (compared by
coordinates), or `PeakAnnotation` objects (compared by peak ID).  The overlap
enrichment functions are Python approximations of ChIPseeker's permutation and
hypergeometric tests; p-values are not expected to be byte-identical.

## Command line

```bash
pychipseeker annotate \
    --peaks peaks.bed \
    --gtf genes.gtf \
    --output peaks.annotated.tsv \
    --plot-dir plots \
    --summary
```

Plot from an existing annotation table:

```bash
pychipseeker plot \
    --annotated peaks.annotated.tsv \
    --outdir plots
```

## Output columns

`ann.to_dataframe()` returns one row per annotated peak:

| Column | Meaning |
|---|---|
| `peak_id` | Peak name (BED column 4) or generated ID |
| `chrom`, `start`, `end` | Peak coordinates (0-based, half-open) |
| `width`, `center` | Peak width and center |
| `strand` | Peak strand parsed from the input |
| `annotation` | ChIPseeker label, e.g. `Promoter (<=1kb)`, `5' UTR`, `3' UTR`, `Exon (...)`, `Intron (...)`, `Downstream (<=300bp)`, `Distal Intergenic` |
| `gene_id`, `gene_name` | Gene assigned by the genomic annotation |
| `gene_chrom`, `gene_start`, `gene_end`, `gene_strand`, `gene_length`, `gene_type` | Assigned gene model |
| `transcript_id`, `transcript_name`, `transcript_type` | Assigned transcript (transcript level only) |
| `distance_to_tss` | Signed distance from the closer peak end to the assigned TSS; negative = upstream |
| `nearest_gene_id`, `nearest_gene_name` | Nearest gene by TSS |
| `nearest_transcript_id`, `nearest_distance` | Nearest transcript and its signed distance |

`ann.to_chipseeker()` renames these to ChIPseeker-like names (`peakChr`,
`peakStart`, `peakEnd`, `geneId`, `geneName`, `distanceToTSS`, ...).
Coordinates remain 0-based and half-open; convert ChIPseeker/GRanges
coordinates with `start - 1` when moving between conventions.

## Annotation rules

Feature priority follows ChIPseeker's default:

1. Promoter
2. 5' UTR
3. 3' UTR
4. Exon
5. Intron
6. Downstream
7. Intergenic

Promoter annotation uses the distance to the **nearest TSS**, not "any TSS in
the window".  Distances are strand-aware and use the closer of the peak start
and inclusive end, matching ChIPseeker's GRanges calculations.

`anno_stat()` produces the detailed `Feature`/`Frequency` table used by
ChIPseeker (`Promoter (<=1kb)`, `1st Exon`, `Other Exon`, `1st Intron`,
`Downstream (<=300)`, ...).  `summary()` is a convenience table of major
annotation classes.

### Isoform alignment

By default, when a peak overlaps an exon/intron/UTR from one isoform but the
nearest TSS belongs to another isoform, `pychipseeker` reports the
gene/transcript and distance of the isoform that supplied the genomic
annotation.  This follows the ChIPseeker development fix for issue #252 and
keeps `annotation`, `transcript_id` and `distance_to_tss` consistent.

Set `align_annotation_feature=False` (or CLI
`--no-align-annotation-feature`) to reproduce ChIPseeker 1.48 and earlier,
where the annotation text and transcript columns can refer to different
isoforms.  The nearest feature is still available in `nearest_*` columns.

## Supported annotation formats

- **GTF / GTF.gz**: fully supported, including GENCODE/Ensembl attributes.
- **GFF3**: best-effort.  Gene IDs are recovered from `gene`/`mRNA` rows when
  present.  If child rows only have `Parent=transcript:...` and the file has no
  mRNA rows, convert it to GTF first or use an Ensembl/GENCODE GTF.

## Comparison with R ChIPseeker

The `comparison/` directory contains a small GTF/BED fixture and optional R and
Python runners.  The core annotation engine was checked against
**ChIPseeker 1.48.0** on 22 fixture peaks (transcript level, gene level, small
and default promoter windows, `overlap="TSS"` and `overlap="all"`).  The
comparison matched:

- annotation labels and subcategories;
- signed `distanceToTSS`;
- gene/transcript IDs and gene coordinates/length;
- `annoStat` labels, order and frequencies;
- `detailGenomicAnnotation` booleans.

One intentional difference remains: for a peak in a shared exon, ChIPseeker
1.48 reports the nearest-TSS isoform even though the annotation label came from
a different isoform.  `pychipseeker` defaults to the ChIPseeker development fix
(#252) and reports the aligned isoform; `align_annotation_feature=False`
reproduces the release behavior.

The comparison scripts require an existing R installation with ChIPseeker and
its Bioconductor dependencies.  They are not part of the Python test suite.
Set `RSCRIPT` to your `Rscript.exe` if it is not on `PATH`:

```bash
python comparison/compare.py          # requires R + ChIPseeker + txdbmaker
python comparison/compare.py --no-r   # only re-check existing r_*/py_* files
```

## Limitations

- Not a byte-for-byte port and not an official Bioconductor package.
- No R `TxDb`/`OrgDb`/`GRanges` classes.  Use GTF/GFF, `GeneAnnotation`, pandas
  DataFrames and the `annoDb` DataFrame/mapping hook instead.
- ChIPseeker's GEO helper functions are implemented against a bundled snapshot
  of the upstream `gsminfo` table (`pychipseeker/data/gsminfo.csv.gz`):
  `getGEOspecies`, `getGEOgenomeVersion`, `getGEOInfo`, `downloadGEObedFiles`,
  `downloadGSMbedFiles` and `getSampleFiles` all work, and `gsminfo=...`
  accepts a custom snapshot.  To refresh the snapshot locally, use the
  zero-dependency offline converter
  `python scripts/build_gsminfo.py GEOmetadb.sqlite -o gsminfo.csv.gz`
  (see the script's docstring; the very large `GEOmetadb.sqlite` download is
  up to the user).
- GFF3 support is best-effort; GTF is the primary annotation input.
- Visualization is matplotlib-based.  Category labels, percentages and axis
  directions follow ChIPseeker, but ggplot2 aesthetics are not reproduced
  exactly.
- `enrichPeakOverlap` / `enrichAnnoOverlap` implement the same statistical
  ideas with Python code; exact p-values, multiple-testing details and
  permutation counts can differ.
- Peak annotation uses a bin-based interval index and Python loops.  It is
  intended for typical peak sets, not whole-genome single-base processing.
- `pyBigWig` is a declared dependency.  On platforms without a prebuilt wheel,
  a C compiler or conda-forge is required.

See [PARITY.md](PARITY.md) for the function-by-function status.

## Development

```bash
pip install -e ".[dev]"
pytest -q
```

The core test suite uses synthetic GTF/BED data.  Real bigWig tests use
upstream pyBigWig's own `pyBigWigTest/test.bw` file (bundled as
`tests/data/test.bw`) and exercise `get_tag_matrix(bigwig=...)`,
`get_tag_matrix_from_bigwig` and `extract_tss_matrix` against a real file;
a `pyBigWig` stub is still used to cover the optional-import code paths.

## License

MIT.  See `LICENSE`.  Bundled ChIPseeker-derived datasets (`gsminfo`,
`ucsc_release`, the GEO sample BED files) and the patched `windows_build/io.c`
come from differently licensed third parties; see [NOTICE.md](NOTICE.md).
