"""Generate a ChIPseeker-style figure plate and individual vector figures.

Run::

    python examples/make_figure_plate.py

Outputs (under ``examples/output``):

- ``pychipseeker_figure_plate.pdf`` / ``.png`` -- one A3-portrait plate with
  panels A-F (annotation pie/bar, distance to TSS, tag profile, tag heatmap,
  annoStat table)
- ``pychipseeker_figures.pdf`` -- the same figures one per page (vector)
- ``figures/*.pdf`` / ``figures/*.png`` -- individual figures for reuse

The example promoter window is deliberately narrow (10 bp) so that the
virtual peaks exercise every genomic category, exactly as in
``run_example.py``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure

import pychipseeker as pcs

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
FIGDIR = OUT / "figures"


def _save_one(fig: plt.Figure, stem: str) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{stem}.png", dpi=220, bbox_inches="tight")


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
    table.set_fontsize(9)
    table.scale(1.0, 1.55)
    for cell in table.get_celld().values():
        cell.set_edgecolor("white")


def _panel_figures(ann, tag):
    """Yield ``(figure, stem, panel letter)`` for the six standard panels."""

    yield ann.plot_anno_pie()[0], "A_annotation_pie", "A"
    yield ann.plot_anno_bar()[0], "B_annotation_bar", "B"
    yield ann.plot_dist_to_tss()[0], "C_distance_to_tss", "C"
    yield pcs.plot_avg_prof(tag, conf=None)[0], "D_tag_profile", "D"
    yield pcs.tag_heatmap(tag)[0], "E_tag_heatmap", "E"
    yield pcs.covplot(ann.peaks, title="Peak coverage")[0], "F_peak_coverage", "F"


def _composite_plate(ann, tag) -> plt.Figure:
    """One vector plate with panels A-F.

    The plotting helpers call ``fig.tight_layout()`` internally, which would
    override the GridSpec positions of a shared figure.  Disable it only while
    composing the plate; the individual figures still use the normal layout.
    """

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
            hspace=0.40,
            wspace=0.55,
        )
        axes = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(2)]

        ann.plot_anno_pie(ax=axes[0])
        ann.plot_anno_bar(ax=axes[1])
        ann.plot_dist_to_tss(ax=axes[2])
        pcs.plot_avg_prof(tag, conf=None, ax=axes[3])
        pcs.tag_heatmap(tag, ax=axes[4])
        _stats_axes(ann, axes[5])

        for letter, ax in zip("ABCDEF", axes):
            ax.text(
                -0.02,
                1.01,
                letter,
                transform=ax.transAxes,
                fontsize=19,
                weight="bold",
                ha="right",
                va="bottom",
            )

        fig.suptitle(
            "pychipseeker — ChIPseeker-style peak annotation",
            fontsize=18,
            weight="bold",
            y=0.975,
        )
        fig.text(
            0.5,
            0.022,
            "Virtual example peaks + GTF (promoter window 10 bp) · panels A-E are "
            "vector matplotlib output; panel F: annotation statistics",
            ha="center",
            fontsize=8.5,
            color="0.4",
        )
        return fig
    finally:
        Figure.tight_layout = original  # type: ignore[assignment]


def main() -> int:
    OUT.mkdir(exist_ok=True)
    ann = pcs.annotate_peak(
        peaks=HERE / "example_peaks.bed",
        gtf=HERE / "example.gtf",
        upstream=10,
        downstream=10,
    )
    promoters = pcs.get_promoters(
        TxDb=ann.gene_annotation, upstream=1000, downstream=1000
    )
    tag = pcs.get_tag_matrix(peak=ann.peaks, windows=promoters, nbin=50)

    # Individual vector figures + a one-figure-per-page PDF.
    with PdfPages(str(OUT / "pychipseeker_figures.pdf")) as pdf:
        for fig, stem, _letter in _panel_figures(ann, tag):
            _save_one(fig, stem)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    plate = _composite_plate(ann, tag)
    plate.savefig(OUT / "pychipseeker_figure_plate.pdf", bbox_inches="tight")
    plate.savefig(
        OUT / "pychipseeker_figure_plate.png", dpi=220, bbox_inches="tight"
    )
    plt.close(plate)

    print(ann.annoStat.to_string(index=False))
    print(f"plate   : {OUT / 'pychipseeker_figure_plate.pdf'}")
    print(f"figures : {OUT / 'pychipseeker_figures.pdf'}")
    print(f"singles : {FIGDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
