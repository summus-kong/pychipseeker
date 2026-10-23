"""Run a small end-to-end example.

The promoter window is intentionally small (10 bp) so that the demo peaks
exercise 5' UTR / exon / intron categories; with the default 3 kb window most
of these peaks would correctly be classified as promoters.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd

import pychipseeker as pcs

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
OUT.mkdir(exist_ok=True)

ann = pcs.annotate_peak(
    peaks=HERE / "example_peaks.bed",
    gtf=HERE / "example.gtf",
    upstream=10,
    downstream=10,
)

print(ann.summary().to_string(index=False))
print()
print(
    ann.to_dataframe()[
        ["peak_id", "annotation", "gene_name", "distance_to_tss"]
    ].to_string(index=False)
)

fig, _ = ann.plot_anno_pie()
fig.savefig(OUT / "annotation_pie.png", dpi=200, bbox_inches="tight")
fig.clf()

fig, _ = ann.plot_anno_bar()
fig.savefig(OUT / "annotation_bar.png", dpi=200, bbox_inches="tight")
fig.clf()

fig, _ = ann.plot_dist_to_tss()
fig.savefig(OUT / "distance_to_tss.png", dpi=200, bbox_inches="tight")
fig.clf()

# The bedGraph/precomputed-matrix path remains available alongside pyBigWig.
signal = pd.DataFrame(
    {
        "chrom": ["chr1", "chr2"],
        "start": [0, 0],
        "end": [10000, 20000],
        "value": [1.0, 2.0],
    }
)
fig, _ = ann.plot_tss_profile(signal=signal, window=100, n_bins=50)
fig.savefig(OUT / "tss_profile.png", dpi=200, bbox_inches="tight")
fig.clf()

fig, _ = ann.plot_tss_heatmap(signal=signal, window=100, n_bins=50)
fig.savefig(OUT / "tss_heatmap.png", dpi=200, bbox_inches="tight")
fig.clf()

# ChIPseeker-style peak-coverage tag matrix and profile/heatmap.
promoters = pcs.getPromoters(
    TxDb=ann.gene_annotation, upstream=1000, downstream=1000
)
tag = pcs.getTagMatrix(peak=ann.peaks, windows=promoters, nbin=50)
fig, _ = pcs.plotAvgProf(tag, conf=None)
fig.savefig(OUT / "tag_profile.png", dpi=200, bbox_inches="tight")
fig.clf()

fig, _ = pcs.tagHeatmap(tag)
fig.savefig(OUT / "tag_heatmap.png", dpi=200, bbox_inches="tight")
fig.clf()

print(f"\nSaved plots to {OUT}")
