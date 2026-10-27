"""Build the figure plate from real ENCODE H3K4me3 peaks.

Run::

    python examples/call_encode_peaks.py \
        --bigwig C:/path/ENCFF893ZGR.bigWig \
        --out examples/encode_peaks.bed
    python examples/make_encode_plate.py

``examples/encode_peaks.bed`` holds 11,291 peaks called from the real
ENCODE track ENCFF893ZGR (H3K4me3, mm10, experiment ENCSR000CDQ) at
signal >= 3, bridged over 1 kb gaps and at least 500 bp wide.

Outputs (``examples/output/encode``):

- ``encode_figure_plate.pdf`` / ``.png`` -- one A3 portrait plate, panels A-F
- ``encode_figures.pdf`` -- six pages, one figure each
- ``figures/*.pdf`` / ``figures/*.png`` -- individual vector figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure

import pychipseeker as pcs

HERE = Path(__file__).resolve().parent


def _save_one(fig: plt.Figure, stem: str, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figdir / f"{stem}.png", dpi=220, bbox_inches="tight")


def _stats_axes(ann, ax) -> None:
    stat = ann.annoStat.copy()
    ax.axis("off")
    ax.set_title("Annotation statistics (annoStat)", loc="left", fontsize=11)
    table = ax.table(
        cellText=stat[["Feature", "Frequency", "Count"]].round(2).astype(str).values,
        colLabels=["Feature", "Frequency (%)", "Count"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(0.92, 1.5)
    for cell in table.get_celld().values():
        cell.set_edgecolor("white")


def _panel_figures(ann, tag):
    yield ann.plot_anno_pie()[0], "A_annotation_pie", "A"
    yield ann.plot_anno_bar()[0], "B_annotation_bar", "B"
    yield ann.plot_dist_to_tss()[0], "C_distance_to_tss", "C"
    yield pcs.plot_avg_prof(tag, conf=None)[0], "D_tag_profile", "D"
    yield pcs.tag_heatmap(tag)[0], "E_tag_heatmap", "E"
    yield pcs.covplot(ann.peaks, title="Peak coverage")[0], "F_peak_coverage", "F"


def _composite_plate(ann, tag, caption: str) -> plt.Figure:
    original = Figure.tight_layout
    Figure.tight_layout = lambda self, *args, **kwargs: None  # type: ignore[assignment]
    try:
        fig = plt.figure(figsize=(12.0, 15.5))
        gs = fig.add_gridspec(
            3,
            2,
            left=0.05,
            right=0.93,
            top=0.90,
            bottom=0.09,
            hspace=0.42,
            wspace=0.78,
        )
        axes = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(2)]

        ann.plot_anno_pie(ax=axes[0])
        legend = axes[0].get_legend()
        if legend is not None:
            legend.set_frame_on(False)
            for text in legend.get_texts():
                text.set_fontsize(8)
        ann.plot_anno_bar(ax=axes[1])
        ann.plot_dist_to_tss(ax=axes[2])
        pcs.plot_avg_prof(tag, conf=None, ax=axes[3])
        pcs.tag_heatmap(tag, ax=axes[4])
        _stats_axes(ann, axes[5])

        for letter, ax in zip("ABCDEF", axes):
            ax.text(
                -0.06,
                1.045,
                letter,
                transform=ax.transAxes,
                fontsize=19,
                weight="bold",
                ha="right",
                va="bottom",
            )

        fig.suptitle(
            "pychipseeker — real ENCODE H3K4me3 peaks (mm10)",
            fontsize=18,
            weight="bold",
            y=0.975,
        )
        fig.text(0.5, 0.022, caption, ha="center", fontsize=8.5, color="0.4")
        return fig
    finally:
        Figure.tight_layout = original  # type: ignore[assignment]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", default=str(HERE / "encode_peaks.bed"))
    parser.add_argument("--txdb", default="mm10")
    parser.add_argument("--outdir", default=str(HERE / "output" / "encode"))
    parser.add_argument(
        "--sample", type=int, default=0, help="subsample peaks (fixed seed); 0 = all"
    )
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--upstream", type=int, default=3000)
    parser.add_argument("--downstream", type=int, default=3000)
    args = parser.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    figdir = out / "figures"

    ann = pcs.annotate_peak(peaks=args.peaks, TxDb=args.txdb, verbose=False)
    used = len(ann.peaks)
    if args.sample and used > args.sample:
        peaks = ann.peaks.sample(n=args.sample, random_state=args.seed)
        peaks = peaks.sort_values(["chrom", "start"], kind="mergesort")
        ann = pcs.annotate_peak(peaks=peaks, TxDb=args.txdb, verbose=False)
        used = len(ann.peaks)

    promoters = pcs.get_promoters(
        TxDb=ann.gene_annotation,
        upstream=args.upstream,
        downstream=args.downstream,
    )
    tag = pcs.get_tag_matrix(peak=ann.peaks, windows=promoters, nbin=60)

    caption = (
        f"ENCODE ENCSR000CDQ · H3K4me3 · mm10 · {used:,} peaks "
        f"(signal ≥ 3, 500 bp–8.9 kb) · promoter window ±{args.upstream / 1000:g} kb"
    )

    with PdfPages(str(out / "encode_figures.pdf")) as pdf:
        for fig, stem, _letter in _panel_figures(ann, tag):
            _save_one(fig, stem, figdir)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    plate = _composite_plate(ann, tag, caption)
    plate.savefig(out / "encode_figure_plate.pdf", bbox_inches="tight")
    plate.savefig(out / "encode_figure_plate.png", dpi=220, bbox_inches="tight")
    plt.close(plate)

    print(ann.annoStat.to_string(index=False))
    print(f"peaks   : {used:,}")
    print(f"promoter: {ann.annoStat.loc[0, 'Frequency']:.1f}% (<=1 kb)")
    print(f"plate   : {out / 'encode_figure_plate.pdf'}")
    print(f"figures : {out / 'encode_figures.pdf'}")
    print(f"singles : {figdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
