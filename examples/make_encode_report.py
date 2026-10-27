"""Build a multi-page PDF report from real ENCODE H3K4me3 peaks.

Run::

    python examples/make_encode_report.py

Inputs:

- ``examples/encode_peaks.bed`` -- 11,291 peaks called from the real ENCODE
  track ENCFF893ZGR (H3K4me3, mm10, ENCSR000CDQ)
- ``C:/Users/houwe/AppData/Local/Temp/ENCFF893ZGR.bigWig`` -- the source track
  (override with ``--bigwig``)

Output: ``examples/output/encode/encode_report.pdf``.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pychipseeker as pcs
from pychipseeker._version import __version__

HERE = Path(__file__).resolve().parent
PAGE = (11.0, 8.5)


def _default_bigwig() -> Path:
    """Resolve ENCFF893ZGR.bigWig from the environment or the usual places."""

    env = os.environ.get("ENCODE_BIGWIG")
    if env:
        return Path(env)
    for candidate in (
        HERE / "ENCFF893ZGR.bigWig",
        Path(tempfile.gettempdir()) / "ENCFF893ZGR.bigWig",
    ):
        if candidate.exists():
            return candidate
    return HERE / "ENCFF893ZGR.bigWig"


def _title_page(ann, bigwig: Path, n_all: int) -> plt.Figure:
    fig = plt.figure(figsize=PAGE)
    fig.text(0.06, 0.90, "pychipseeker", fontsize=30, weight="bold")
    fig.text(
        0.06,
        0.845,
        "ChIPseeker-style annotation of real ENCODE H3K4me3 peaks",
        fontsize=15,
    )
    fig.text(
        0.06,
        0.795,
        "ENCSR000CDQ · H3K4me3 · Mus musculus placenta · mm10",
        fontsize=11,
        color="0.35",
    )

    lines = [
        ("Signal track", f"{bigwig.name}  ({bigwig.stat().st_size / 1e6:.1f} MB)"),
        ("Experiment", "ENCSR000CDQ (ENCODE), target H3K4me3"),
        ("Assembly / annotation", "mm10 (UCSC), bundled pychipseeker mm10 bundle"),
        ("Peaks", f"{n_all:,} threshold-called peaks (signal ≥ 3, ≥ 500 bp)"),
        ("Peak width", "median 958 bp, mean 1,119 bp, max 8.9 kb"),
        ("Annotation", f"TxDb='mm10', promoter window −1 kb / +1 kb"),
        ("Promoter (≤1 kb)", f"{ann.annoStat.loc[0, 'Frequency']:.1f}%"),
        ("Distal intergenic", f"{ann.annoStat.iloc[-1]['Frequency']:.2f}%"),
        ("pychipseeker", __version__),
    ]
    y = 0.70
    for label, value in lines:
        fig.text(0.06, y, label, fontsize=10, color="0.35")
        fig.text(0.31, y, value, fontsize=10)
        y -= 0.045

    fig.text(
        0.06,
        0.17,
        "Peaks are extracted from the real signal track by thresholding, because the\n"
        "ENCODE portal was not reachable from this machine. They are not ENCODE's\n"
        "narrowPeak calls; the annotation, plotting and statistics code paths are the\n"
        "same ones used for any peak table.",
        fontsize=9,
        color="0.45",
        linespacing=1.5,
    )
    return fig


def _stats_page(ann) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PAGE)
    ax.axis("off")
    ax.set_title("annoStat — detailed annotation categories", fontsize=14, loc="left")
    stat = ann.annoStat
    table = ax.table(
        cellText=stat[["Feature", "Frequency", "Count"]].round(3).astype(str).values,
        colLabels=["Feature", "Frequency (%)", "Count"],
        loc="center",
        cellLoc="center",
        colWidths=[0.42, 0.25, 0.18],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.9)
    for cell in table.get_celld().values():
        cell.set_edgecolor("white")
    fig.text(
        0.12,
        0.12,
        "Promoter sub-categories follow ChIPseeker (<=1 kb, 1-2 kb, 2-3 kb).",
        fontsize=9,
        color="0.45",
    )
    return fig


def _methods_page(ann, bigwig: Path) -> plt.Figure:
    fig = plt.figure(figsize=PAGE)
    fig.text(0.06, 0.90, "Methods / reproducibility", fontsize=17, weight="bold")
    body = (
        "1. Signal source\n"
        f"   {bigwig}\n"
        "   ENCFF893ZGR, H3K4me3 fold-change-over-control / signal p-value track,\n"
        "   mm10, experiment ENCSR000CDQ (8-week mouse placenta).\n\n"
        "2. Peak calling (threshold caller, not MACS2)\n"
        "   Runs of bases with signal >= 3.0, bridged across gaps <= 1 kb and kept\n"
        "   when >= 500 bp wide. 11,291 autosomal/sex-chromosome peaks.\n\n"
        "3. Annotation\n"
        "   pychipseeker.annotate_peak(peaks=..., TxDb='mm10')\n"
        "   Promoter = TSS -1 kb / +1 kb; categories follow ChIPseeker.\n\n"
        "4. Figures\n"
        "   plot_anno_pie / plot_anno_bar / plot_dist_to_tss / plot_avg_prof /\n"
        "   tag_heatmap / covplot / plot_tss_profile / plot_tss_heatmap.\n"
        "   Signal panels use 3,000 transcripts sampled with seed 20261026\n"
        "   (all 81,540 transcripts need ~1 h of per-window bigWig reads).\n\n"
        "5. Reproduce\n"
        "   python examples/call_encode_peaks.py --bigwig <ENCFF893ZGR.bigWig>\n"
        "   python examples/make_encode_report.py\n"
        "   python examples/make_encode_plate.py"
    )
    fig.text(0.06, 0.84, body, fontsize=10, va="top", family="monospace", linespacing=1.45)
    fig.text(0.06, 0.06, f"pychipseeker {__version__}", fontsize=9, color="0.4")
    return fig


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", default=str(HERE / "encode_peaks.bed"))
    parser.add_argument("--bigwig", default=str(_default_bigwig()))
    parser.add_argument("--txdb", default="mm10")
    parser.add_argument("--out", default=str(HERE / "output" / "encode" / "encode_report.pdf"))
    args = parser.parse_args()

    bigwig = Path(args.bigwig)
    if not bigwig.exists():
        parser.error(
            f"bigWig not found: {bigwig}\n"
            "pass --bigwig PATH or set ENCODE_BIGWIG (ENCFF893ZGR.bigWig)"
        )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    ann = pcs.annotate_peak(peaks=args.peaks, TxDb=args.txdb, verbose=False)
    n_all = len(ann.peaks)
    promoters = pcs.get_promoters(
        TxDb=ann.gene_annotation, upstream=3000, downstream=3000
    )
    tag = pcs.get_tag_matrix(peak=ann.peaks, windows=promoters, nbin=60)

    # Reading 81k transcript windows one by one from the bigWig takes ~1 h, so
    # the signal panels use a fixed random subset (the peak panels use all
    # peaks).  This is the "sampled proportion" allowed by the task.
    transcripts = ann.gene_annotation.features("transcript")
    n_signal = min(3000, len(transcripts))
    tss_subset = transcripts.sample(n=n_signal, random_state=20261026)
    signal_matrix = pcs.extract_tss_matrix(
        bigwig, tss_subset, window=2000, n_bins=100
    )

    pages = [(ann.plot_anno_pie()[0], "Annotation distribution (pie)")]
    pages.append((ann.plot_anno_bar()[0], "Annotation distribution (bar)"))
    pages.append((ann.plot_dist_to_tss()[0], "Distance to TSS"))
    pages.append((pcs.covplot(ann.peaks, title="Peak coverage")[0], "Peak coverage"))
    pages.append(
        (
            pcs.plot_tss_profile(
                signal_matrix,
                window=2000,
                title="H3K4me3 signal around TSS (±2 kb)",
                show_sem=False,
            )[0],
            f"Signal profile around TSS ({n_signal:,} sampled transcripts)",
        )
    )
    pages.append(
        (
            pcs.plot_tss_heatmap(
                signal_matrix,
                window=2000,
                cmap="Reds",
                sort_rows=True,
                title="H3K4me3 signal around TSS",
            )[0],
            f"Signal heatmap around TSS ({n_signal:,} sampled transcripts)",
        )
    )
    pages.append((pcs.plot_avg_prof(tag, conf=None)[0], "Peak profile (tagMatrix)"))
    pages.append((pcs.tag_heatmap(tag)[0], "Peak heatmap (tagMatrix)"))
    pages.append((_stats_page(ann), "Annotation statistics"))
    pages.insert(0, (_title_page(ann, bigwig, n_all), "Title"))
    pages.append((_methods_page(ann, bigwig), "Methods"))

    with PdfPages(str(out)) as pdf:
        for fig, name in pages:
            fig.suptitle(name, fontsize=9, color="0.45", y=0.985)
            pdf.savefig(fig, bbox_inches="tight")
            print("page:", name)
            plt.close(fig)

    print(f"PDF written: {out} ({out.stat().st_size / 1024:.0f} KiB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
