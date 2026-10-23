"""Build a complete multi-page PDF report from the bundled example data.

Run::

    python examples/make_report_pdf.py

Uses only the virtual data shipped with pychipseeker —
``examples/example_peaks.bed`` (11 hand-wired peaks covering every genomic
category) and ``examples/example.gtf`` — and writes every visualization as
one page of ``examples/output/pychipseeker_report.pdf``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot as plt
import pandas as pd

import pychipseeker as pcs
from pychipseeker._version import __version__

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"


def _title_page(stat: pd.DataFrame, df: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.93, "pychipseeker", ha="center", fontsize=28, weight="bold")
    fig.text(0.5, 0.895, "ChIPseeker-style peak annotation report", ha="center", fontsize=13)
    fig.text(0.5, 0.86, f"pychipseeker {__version__} · virtual example data", ha="center", fontsize=10, color="0.4")

    lines = [
        "Input data (bundled, virtual):",
        f"  peaks : examples/example_peaks.bed  ({len(df)} synthetic peaks)",
    ]
    fig.text(0.12, 0.78, "\n".join(lines), fontsize=11, family="monospace")

    fig.text(0.12, 0.68, "Annotation statistics (ChIPseeker annoStat):", fontsize=12, weight="bold")
    ax = fig.add_axes((0.12, 0.36, 0.76, 0.28))
    ax.axis("off")
    table = ax.table(
        cellText=stat.round(2).astype(str).values,
        colLabels=[c.replace("_", " ").title() for c in stat.columns],
        loc="upper center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)

    genes = df.groupby("gene_name").size().sort_values(ascending=False)
    fig.text(
        0.12, 0.30,
        "Peaks assigned per gene:\n  "
        + ";  ".join(f"{g}: {n}" for g, n in genes.items()),
        fontsize=11, family="monospace",
    )

    kinds = df["annotation"].value_counts()
    fig.text(
        0.12, 0.18,
        "Categories hit by the demo peaks:\n  "
        + "\n  ".join(f"{k}: {n}" for k, n in kinds.items()),
        fontsize=11, family="monospace",
    )
    fig.text(
        0.5, 0.06,
        "promoter window = 10 bp (deliberately narrow so UTR/exon/intron categories appear)",
        ha="center", fontsize=9, style="italic", color="0.35",
    )
    return fig


def _final_page() -> plt.Figure:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.85, "Reproduce this report", ha="center", fontsize=16, weight="bold")
    fig.text(
        0.10, 0.60,
        "python examples/make_report_pdf.py\n\n"
        "from pychipseeker.report steps:\n"
        "  ann = pcs.annotate_peak(peaks=..., gtf=..., upstream=10, downstream=10)\n"
        "  ann.plot_anno_pie() / ann.plot_anno_bar() / ann.plot_dist_to_tss()\n"
        "  pcs.covplot(ann.peaks)\n"
        "  promoters = pcs.getPromoters(TxDb=ann.gene_annotation, upstream=1000, downstream=1000)\n"
        "  tag = pcs.getTagMatrix(peak=ann.peaks, windows=promoters, nbin=50)\n"
        "  pcs.plotAvgProf(tag, conf=None)  ·  pcs.tagHeatmap(tag)\n"
        "  ann.plot_tss_profile(signal=..., window=100, n_bins=50)\n"
        "  ann.plot_tss_heatmap(signal=..., window=100, n_bins=50)",
        fontsize=10.5, family="monospace", va="top",
    )
    fig.text(0.5, 0.06, f"pychipseeker {__version__}", ha="center", fontsize=9, color="0.4")
    return fig


def main() -> int:
    OUT.mkdir(exist_ok=True)
    pdf_path = OUT / "pychipseeker_report.pdf"

    ann = pcs.annotate_peak(
        peaks=HERE / "example_peaks.bed",
        gtf=HERE / "example.gtf",
        upstream=10,
        downstream=10,
    )
    stat = ann.summary()
    df = ann.to_dataframe()
    print(stat.to_string(index=False))

    # Virtual signal for the signal-based TSS profile/heatmap pages
    # (ChIPseeker accepts plain bedGraph-like DataFrames here).
    signal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [10000, 20000],
            "value": [1.0, 2.0],
        }
    )

    promoters = pcs.getPromoters(TxDb=ann.gene_annotation, upstream=1000, downstream=1000)
    tag = pcs.getTagMatrix(peak=ann.peaks, windows=promoters, nbin=50)

    figures = [
        (_title_page(stat, df), "Title / overview"),
        (ann.plot_anno_pie()[0], "Annotation distribution (pie)"),
        (ann.plot_anno_bar()[0], "Annotation distribution (bar)"),
        (ann.plot_dist_to_tss()[0], "Distance to TSS"),
        (pcs.covplot(ann.peaks, title="Peak coverage over chromosomes")[0], "Peak coverage (covplot)"),
        (ann.plot_tss_profile(signal=signal, window=100, n_bins=50)[0], "TSS avg profile (signal)"),
        (ann.plot_tss_heatmap(signal=signal, window=100, n_bins=50)[0], "TSS heatmap (signal)"),
        (pcs.plot_avg_prof(tag, conf=None)[0], "Promoter avg profile (tagMatrix)"),
        (pcs.tag_heatmap(tag)[0], "Promoter tag heatmap (tagMatrix)"),
        (_final_page(), "How to reproduce"),
    ]

    with PdfPages(str(pdf_path)) as pdf:
        for fig, name in figures:
            fig.suptitle(name, fontsize=10, color="0.3", y=0.995)
            pdf.savefig(fig, bbox_inches="tight")
            print("page:", name)
            plt.close(fig)

    print(f"PDF written: {pdf_path} ({pdf_path.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
