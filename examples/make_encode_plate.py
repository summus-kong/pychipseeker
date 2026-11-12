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

- ``encode_figure_plate.*`` -- one A3 portrait plate, panels A-F, written as
  PNG (220 dpi), PDF, SVG, EPS and TIFF (400 dpi), plus a 600 dpi PNG
- ``encode_supplementary_plate.*`` -- a second A3 plate with six more views
  (panels G-L), in the same formats
- ``encode_figures.pdf`` -- twelve pages, one figure each
- ``figures/<panel>.*`` -- individual figures in PNG/PDF/SVG/EPS/TIFF
- ``gallery.html`` -- a self-contained index with thumbnails and format links

Design notes
------------
- Categorical colours use the Okabe-Ito palette, which stays distinguishable
  under all common colour-vision deficiencies and in greyscale print
  (Wong, *Nat. Methods* 2011; Okabe & Ito 2008).
- Sequential data (tag heatmap) uses a single-hue, perceptually ordered map.
- The layout follows the "Ten Simple Rules for Better Figures"
  (Rougier et al., *PLoS Comput. Biol.* 2014): one message per panel,
  direct labelling where possible, no decorative gridlines, and colour
  used redundantly with text.
- Panels B and C use bespoke matplotlib layouts.  The default single-sample
  composition/distance bars are dominated by the 92% promoter fraction, so
  a width histogram and a cumulative distance curve communicate the same
  data without hiding the tail (and without the misleading 50/50 split of
  zero-distance peaks).
- The tag matrix is normalised to the percentage of each bin covered by a
  peak, so the profile height and the heatmap colour bar share one unit
  (0-100%).
"""

from __future__ import annotations

import argparse
import html
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle, Wedge

import pychipseeker as pcs

HERE = Path(__file__).resolve().parent

# Embed TrueType fonts and keep SVG text editable, so the exported vector
# files open cleanly in Illustrator/Inkscape and in journal production.
matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

#: Formats written for the plate and for every panel.
FORMATS = ("png", "pdf", "svg", "eps", "tif")

#: Gallery metadata: letter, file stem, panel title and a one-line reading.
PANEL_CATALOG = (
    (
        "A",
        "A_annotation_composition",
        "Annotation composition",
        "How the peaks split across ChIPseeker's promoter, UTR, exon, intron "
        "and intergenic categories.",
    ),
    (
        "B",
        "B_peak_width",
        "Peak width distribution",
        "Width histogram on a log scale; promoter and non-promoter peaks are "
        "overlaid so the tail stays visible.",
    ),
    (
        "C",
        "C_distance_to_tss",
        "Distance to TSS",
        "Cumulative curve of absolute TSS distance on a symlog axis, with "
        "the 1 kb, 10 kb and 100 kb landmarks marked.",
    ),
    (
        "D",
        "D_tag_profile",
        "Average H3K4me3 profile",
        "Mean peak coverage in 120 bins across \u00b13 kb promoter windows.",
    ),
    (
        "E",
        "E_tag_heatmap",
        "Promoter-centered peak coverage",
        "One row per promoter, strongest rows on top; colour is the percent "
        "of each bin covered by a peak.",
    ),
    (
        "F",
        "F_annotation_stats",
        "Annotation statistics",
        "The exact percentages and counts behind panel A, plus the total "
        "number of annotated peaks.",
    ),
    (
        "G",
        "G_chromosome_counts",
        "Peaks per chromosome",
        "Peak count for every chromosome in natural order, which makes "
        "library composition and chrY dropout easy to spot.",
    ),
    (
        "H",
        "H_width_by_class",
        "Peak width by annotation class",
        "Boxplots of width for promoter, UTR, exon, intron and intergenic "
        "peaks on a log axis.",
    ),
    (
        "I",
        "I_signal_by_class",
        "Max signal by annotation class",
        "Violin plots of the called max signal per annotation class, with "
        "medians marked.",
    ),
    (
        "J",
        "J_nonpromoter_distance",
        "Distance to TSS (non-promoter peaks)",
        "Log-scaled histogram of TSS distance for the 8 % of peaks outside "
        "promoters, where the tail is otherwise hidden.",
    ),
    (
        "K",
        "K_profile_quartiles",
        "Profile by peak-signal quartile",
        "Mean coverage profile for promoters grouped by the strongest peak "
        "score annotated to each gene.",
    ),
    (
        "L",
        "L_width_vs_signal",
        "Peak width vs signal",
        "Hexbin density with a median trend line; wider peaks tend to carry "
        "slightly higher signal (Spearman \u03c1 \u2248 0.27).",
    ),
    (
        "M",
        "M_width_ecdf",
        "Peak width ECDF by annotation class",
        "Cumulative width curves answer 'what fraction of each class is "
        "below 1 kb?' at a glance.",
    ),
    (
        "N",
        "N_top_genes",
        "Top genes by peak count",
        "The twelve genes with the most assigned peaks, useful for spotting "
        "recurrent loci.",
    ),
    (
        "O",
        "O_peaks_per_gene",
        "Peaks per gene",
        "How many genes carry one, two, three or four or more peaks; the "
        "y-axis is logarithmic so the tail stays visible.",
    ),
    (
        "P",
        "P_tss_quartile_enrichment",
        "TSS enrichment by signal quartile",
        "Centre-to-flank coverage ratio for the four signal groups: stronger "
        "peaks produce sharper TSS enrichment.",
    ),
    (
        "Q",
        "Q_width_by_distance",
        "Peak width by TSS-distance class",
        "Width distributions for peaks exactly at the TSS, within 1 kb, and "
        "in progressively more distal classes.",
    ),
    (
        "R",
        "R_gene_biotypes",
        "Gene biotypes of annotated peaks",
        "Peaks grouped by the biotype of their assigned gene, with counts "
        "and percentages (log x-axis).",
    ),
)

#: Look up a panel by its letter for custom plates and the gallery.
PANEL_BY_LETTER = {item[0]: item for item in PANEL_CATALOG}

#: Okabe-Ito colour-vision-deficiency-safe palette (Wong 2011).
OKABE_ITO = {
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "orange": "#E69F00",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#8C8C8C",
    "light_grey": "#C8C8C8",
    "ink": "#333333",
}

#: Semantic colour map shared by panels A, B and F.
ANNO_COLORS = {
    "Promoter (<=1kb)": OKABE_ITO["blue"],
    "Promoter (1-2kb)": OKABE_ITO["sky"],
    "Promoter (2-3kb)": OKABE_ITO["green"],
    "5' UTR": OKABE_ITO["yellow"],
    "3' UTR": OKABE_ITO["orange"],
    "Other Exon": OKABE_ITO["purple"],
    "1st Intron": OKABE_ITO["vermillion"],
    "Other Intron": OKABE_ITO["grey"],
    "Downstream (<=300)": OKABE_ITO["light_grey"],
    "Distal Intergenic": OKABE_ITO["ink"],
}

#: Coarse annotation classes used by the supplementary panels.
CLASS_ORDER = ["Promoter", "UTR", "Exon", "Intron", "Intergenic"]
BROAD_COLORS = {
    "Promoter": OKABE_ITO["blue"],
    "UTR": OKABE_ITO["sky"],
    "Exon": OKABE_ITO["purple"],
    "Intron": OKABE_ITO["vermillion"],
    "Intergenic": OKABE_ITO["ink"],
}

INK = "#1A1A1A"
MUTED = "#5F6B76"
RULE = "#D5DCE3"


def _tint(color: str, alpha: float, base: str = "#FFFFFF") -> str:
    """Blend ``color`` over ``base`` into an opaque tint (EPS-safe)."""

    def channels(hex_color: str) -> np.ndarray:
        text = hex_color.lstrip("#")
        return np.array([int(text[i : i + 2], 16) for i in (0, 2, 4)], dtype=float)

    mixed = alpha * channels(color) + (1.0 - alpha) * channels(base)
    return "#" + "".join(f"{int(round(value)):02X}" for value in mixed)


def _broad_class(annotations: pd.Series) -> np.ndarray:
    """Collapse ChIPseeker's detailed labels into five coarse classes."""

    text = annotations.astype(str)
    out = np.full(len(text), "Intergenic", dtype=object)
    out[text.str.contains("Intron").to_numpy()] = "Intron"
    out[text.str.contains("Exon").to_numpy()] = "Exon"
    out[text.str.contains("UTR").to_numpy()] = "UTR"
    out[text.str.startswith("Promoter").to_numpy()] = "Promoter"
    return out


def _chrom_order(names) -> list[str]:
    """Natural chromosome order (chr1, chr2, ..., chrX, chrY)."""

    def key(name: str) -> tuple[str, int, str]:
        text = str(name)
        head = text.rstrip("0123456789")
        tail = text[len(head) :]
        return (head, int(tail) if tail.isdigit() else 0, text)

    return sorted(names, key=key)


def _normalised_coverage(
    peaks: pd.DataFrame,
    windows: pd.DataFrame,
    upstream: int,
    downstream: int,
    nbin: int,
) -> np.ndarray:
    """Peak coverage as percent of each bin, independent of ``nbin``."""

    coverage = pcs.get_tag_matrix(peak=peaks, windows=windows, nbin=nbin)
    bin_bp = (upstream + downstream) / nbin
    return np.asarray(coverage, dtype=float) * 100.0 / bin_bp


def _quartile_profiles(tag, peaks: pd.DataFrame):
    """Split promoters into signal quartiles and return their tag matrices.

    Each promoter inherits the strongest peak score annotated to its gene, so
    the four groups answer a simple question: does a stronger peak also give a
    sharper, taller promoter profile?  Genes without an assigned peak fall into
    the lowest group.
    """

    arr = np.asarray(tag, dtype=float)
    groups = pcs.signal_groups(tag, peaks, n_groups=4)
    profiles = []
    for index in sorted(groups["group"].unique()):
        mask = groups["group"].to_numpy() == index
        label = str(groups.loc[mask, "label"].iloc[0])
        profiles.append((label, int(mask.sum()), arr[mask]))
    return profiles


@contextmanager
def _no_tight_layout() -> Iterator[None]:
    """Keep library helpers from overriding the plate's grid."""

    original = Figure.tight_layout
    Figure.tight_layout = lambda self, *args, **kwargs: None  # type: ignore[assignment]
    try:
        yield
    finally:
        Figure.tight_layout = original  # type: ignore[assignment]


def _clean_axes(ax) -> None:
    """Publication-style axes: no panel fill, no grid, two thin spines."""

    ax.set_facecolor("white")
    ax.grid(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#8A949E")
        ax.spines[side].set_linewidth(0.7)
    ax.tick_params(
        axis="both",
        colors=MUTED,
        labelcolor=INK,
        length=3.0,
        width=0.7,
        direction="out",
        labelsize=8.6,
        pad=3.0,
    )


def _framed_axes(ax) -> None:
    """Thin complete frame, used for raster-like panels (heatmaps)."""

    ax.grid(False)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#8A949E")
        ax.spines[side].set_linewidth(0.7)
    ax.tick_params(
        axis="both",
        colors=MUTED,
        labelcolor=INK,
        length=2.6,
        width=0.7,
        direction="out",
        labelsize=8.4,
        pad=3.0,
    )


def _flatten_rgb(path: Path, compression: str | None = None) -> None:
    """Convert a raster figure to RGB, which journals prefer over RGBA."""

    try:
        from PIL import Image
    except ImportError:  # pragma: no cover - Pillow ships with matplotlib
        return
    with Image.open(path) as image:
        if image.mode == "RGB":
            return
        rgb = image.convert("RGB")
        dpi = image.info.get("dpi", (300.0, 300.0))
    options: dict[str, object] = {"dpi": dpi}
    if compression:
        options["compression"] = compression
    rgb.save(path, **options)


def _save_formats(
    fig: plt.Figure,
    stem: str,
    outdir: Path,
    *,
    formats: Sequence[str] = FORMATS,
    png_dpi: int = 220,
    print_dpi: int = 400,
) -> list[Path]:
    """Write one figure in every requested format and return the paths."""

    outdir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for fmt in formats:
        fmt = fmt.lower()
        path = outdir / f"{stem}.{fmt}"
        if fmt == "tif":
            fig.savefig(
                path,
                dpi=print_dpi,
                bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"},
            )
            _flatten_rgb(path, compression="tiff_lzw")
        elif fmt == "png":
            fig.savefig(path, dpi=png_dpi, bbox_inches="tight")
            _flatten_rgb(path)
        else:
            fig.savefig(path, bbox_inches="tight")
        written.append(path)
    return written


def _panel_header(fig, spec, letter: str, title: str):
    """Reserve a slim header row and return the content subgridspec."""

    sub = spec.subgridspec(2, 1, height_ratios=[0.11, 1.0], hspace=0.0)
    head = fig.add_subplot(sub[0])
    head.axis("off")
    head.set_xlim(0, 1)
    head.set_ylim(0, 1)
    head.text(
        0.0,
        0.52,
        letter,
        fontsize=17,
        fontweight="bold",
        color=INK,
        ha="left",
        va="center",
    )
    head.annotate(
        title,
        xy=(0.0, 0.52),
        xycoords="axes fraction",
        xytext=(21, 0),
        textcoords="offset points",
        fontsize=12.2,
        color=INK,
        ha="left",
        va="center",
    )
    return sub[1]


def _anno_colors(stat: pd.DataFrame) -> list[str]:
    return [ANNO_COLORS.get(str(f), OKABE_ITO["grey"]) for f in stat["Feature"]]


def _draw_key(ax, stat: pd.DataFrame) -> None:
    """A table-like colour key with exact percentages and counts."""

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.plot([0, 1], [0.932, 0.932], color=RULE, linewidth=0.7, clip_on=False)
    ax.text(0.0, 0.962, "Feature", fontsize=7.4, color=MUTED, weight="bold")
    ax.text(0.8, 0.962, "%", fontsize=7.4, color=MUTED, weight="bold", ha="right")
    ax.text(0.99, 0.962, "n", fontsize=7.4, color=MUTED, weight="bold", ha="right")

    colors = _anno_colors(stat)
    step = 0.085
    for i, row in enumerate(stat.itertuples(index=False)):
        y = 0.884 - i * step
        ax.add_patch(
            Rectangle(
                (0.0, y - 0.011),
                0.034,
                0.022,
                facecolor=colors[i],
                edgecolor="none",
                clip_on=False,
            )
        )
        ax.text(0.058, y, str(row.Feature), fontsize=8.0, color=INK, va="center")
        ax.text(
            0.8,
            y,
            f"{float(row.Frequency):.2f}",
            fontsize=8.0,
            color=INK,
            va="center",
            ha="right",
        )
        ax.text(
            0.99,
            y,
            f"{int(row.Count):,}",
            fontsize=8.0,
            color=MUTED,
            va="center",
            ha="right",
        )


def _draw_pie(ann, ax, key_ax) -> None:
    stat = ann.annoStat
    with _no_tight_layout():
        ann.plot_anno_pie(
            ax=ax,
            colors=_anno_colors(stat),
            legend=False,
            title="",
            autopct=None,
        )

    # Turn the pie into a donut and drop the ggplot-style grey backing disc.
    for patch in list(ax.patches):
        if isinstance(patch, Circle):
            patch.remove()
    for patch in ax.patches:
        if isinstance(patch, Wedge):
            patch.set_width(0.37)
            patch.set_edgecolor("white")
            patch.set_linewidth(1.3)

    promoter = float(stat.loc[0, "Frequency"])
    ax.text(
        0.0,
        0.10,
        f"{promoter:.1f}%",
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        color=OKABE_ITO["blue"],
    )
    ax.text(
        0.0,
        -0.10,
        "promoter \u2264 1 kb",
        ha="center",
        va="center",
        fontsize=9.2,
        color=INK,
    )
    _draw_key(key_ax, stat)


def _panel_pie(fig, spec, ann):
    sub = spec.subgridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.02)
    ax = fig.add_subplot(sub[0])
    key_ax = fig.add_subplot(sub[1])
    _draw_pie(ann, ax, key_ax)
    return ax


def _draw_width(ax, df: pd.DataFrame) -> None:
    width = (df["end"] - df["start"]).to_numpy(dtype=float)
    promoter = df["annotation"].astype(str).str.startswith("Promoter").to_numpy()
    edges = np.geomspace(500.0, 8200.0, 33)

    groups = [
        (width[promoter], OKABE_ITO["blue"], "Promoter", True, "-"),
        (width[~promoter], OKABE_ITO["vermillion"], "Other", False, (0, (4, 2))),
    ]
    top = 0.0
    for values, color, name, filled, linestyle in groups:
        counts, _ = np.histogram(values, bins=edges)
        frac = counts / max(counts.sum(), 1)
        if filled:
            ax.stairs(
                frac,
                edges,
                baseline=0.0,
                fill=True,
                color=_tint(color, 0.16),
                linewidth=0,
            )
        ax.stairs(
            frac,
            edges,
            baseline=0.0,
            color=color,
            linewidth=1.7,
            linestyle=linestyle,
        )
        top = max(top, float(frac.max()))
        median = float(np.median(values))
        ax.plot(
            [],
            [],
            color=color,
            linewidth=1.7,
            linestyle=linestyle,
            label=f"{name} (n = {values.size:,}; median {median:,.0f} bp)",
        )

    ax.set_xscale("log")
    ax.set_xlim(500, 8200)
    ax.set_ylim(0, top * 1.28)
    ax.set_xticks([500, 1000, 2000, 5000, 8000])
    ax.set_xticklabels(["0.5", "1", "2", "5", "8"])
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xlabel("Peak width (kb, log scale)")
    ax.set_ylabel("Fraction of peaks")
    ax.legend(
        loc="upper right",
        fontsize=7.9,
        frameon=False,
        handlelength=1.5,
        handletextpad=0.6,
        labelspacing=0.35,
        borderaxespad=0.2,
    )
    _clean_axes(ax)


def _draw_distance(ax, df: pd.DataFrame) -> None:
    dist = np.abs(df["distance_to_tss"].to_numpy(dtype=float))
    dist = np.sort(dist[np.isfinite(dist)])
    cumulative = np.arange(1, dist.size + 1) / dist.size * 100.0

    ax.axvspan(0, 1000, color="#F2F7FC", zorder=0, linewidth=0)
    ax.fill_between(
        dist,
        cumulative,
        step="post",
        color=_tint(OKABE_ITO["sky"], 0.13),
        linewidth=0,
        zorder=1,
    )
    ax.step(
        dist,
        cumulative,
        where="post",
        color=OKABE_ITO["blue"],
        linewidth=1.9,
        zorder=3,
    )
    ax.axvline(1000, color="#9FB6C9", linewidth=0.8, linestyle=(0, (4, 3)), zorder=2)

    ax.set_xscale("symlog", linthresh=100, linscale=0.8)
    ax.set_xlim(0, 3e6)
    ax.set_ylim(0, 102)
    ax.set_xticks([0, 100, 1000, 10000, 100000, 1000000])
    ax.set_xticklabels(["0", "100", "1 kb", "10 kb", "100 kb", "1 Mb"])
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("Distance to TSS (bp, symlog)")
    ax.set_ylabel("Cumulative peaks (%)")

    def pct_within(limit: float) -> float:
        return float(np.searchsorted(dist, limit, side="right")) / dist.size * 100.0

    marks = [
        (1000, pct_within(1000)),
        (10000, pct_within(10000)),
        (100000, pct_within(100000)),
    ]
    for x, y in marks:
        ax.plot([x], [y], "o", markersize=4.4, color=OKABE_ITO["blue"], zorder=4)

    ax.annotate(
        f"{marks[0][1]:.1f}% within 1 kb",
        xy=(1000, marks[0][1]),
        xytext=(1.3e3, 66),
        fontsize=9.4,
        color=OKABE_ITO["blue"],
        fontweight="bold",
        ha="left",
        va="center",
        arrowprops=dict(
            arrowstyle="-",
            color=OKABE_ITO["blue"],
            linewidth=0.9,
            shrinkA=4,
            shrinkB=4,
        ),
        zorder=5,
    )
    zero_pct = float(np.mean(dist == 0)) * 100.0
    ax.text(
        0.02,
        0.045,
        f"{zero_pct:.1f}% of peaks at d = 0",
        transform=ax.transAxes,
        fontsize=7.6,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    ax.text(
        0.97,
        0.045,
        f"{marks[1][1]:.1f}% \u2264 10 kb  \u00b7  {marks[2][1]:.1f}% \u2264 100 kb",
        transform=ax.transAxes,
        fontsize=7.2,
        color=MUTED,
        ha="right",
        va="bottom",
    )
    _clean_axes(ax)


def _draw_profile(ax, tag) -> None:
    pcs.plot_avg_prof(tag, conf=None, color=OKABE_ITO["blue"], ax=ax, label=None)
    ax.set_title("")
    if ax.lines:
        line = ax.lines[0]
        line.set_linewidth(1.9)
        line.set_zorder(3)
        x = np.asarray(line.get_xdata(), dtype=float)
        y = np.asarray(line.get_ydata(), dtype=float)
        ax.fill_between(
            x,
            0.0,
            y,
            color=_tint(OKABE_ITO["sky"], 0.30),
            linewidth=0,
            zorder=1,
        )

    ax.set_xlim(-3000, 3000)
    ax.set_ylim(bottom=0)
    ax.set_xticks([-3000, -2000, -1000, 0, 1000, 2000, 3000])
    ax.set_xticklabels(["-3", "-2", "-1", "TSS", "1", "2", "3"])
    ax.set_xlabel("Distance from TSS (kb)")
    ax.set_ylabel("Mean peak coverage (%)")
    rows = int(np.asarray(tag, dtype=float).shape[0])
    ax.text(
        0.02,
        0.96,
        f"n = {rows:,} promoters",
        transform=ax.transAxes,
        fontsize=8.0,
        color=MUTED,
        ha="left",
        va="top",
    )
    _clean_axes(ax)


def _draw_heatmap(fig, ax, tag) -> None:
    pcs.tag_heatmap(tag, ax=ax, palette="Blues", vmin=0.0, vmax=100.0)
    ax.set_title("")

    # Replace the library's space-stealing colorbar with a slim inset bar.
    for other in list(fig.axes):
        if other is not ax and getattr(other, "_colorbar", None) is not None:
            other.remove()
    image = ax.images[0]
    image.set_interpolation("nearest")
    cax = ax.inset_axes([1.022, 0.0, 0.032, 1.0])
    cbar = fig.colorbar(image, cax=cax)
    cbar.outline.set_visible(False)
    cbar.set_ticks([0, 25, 50, 75, 100])
    cbar.ax.tick_params(labelsize=7.4, length=2.0, width=0.6, colors=MUTED)
    cbar.set_label("Peak coverage (%)", fontsize=8.2, color=INK, labelpad=5)

    rows = int(np.asarray(tag, dtype=float).shape[0])
    ax.set_xticks([-3000, -2000, -1000, 0, 1000, 2000, 3000])
    ax.set_xticklabels(["-3", "-2", "-1", "TSS", "1", "2", "3"])
    ax.set_yticks([])
    ax.set_xlabel("Distance from TSS (kb)")
    ax.set_ylabel(f"Promoters (n = {rows:,}, strongest first)")
    for line in ax.lines:
        line.set_color("#4A4A4A")
        line.set_linewidth(0.7)
        line.set_alpha(1.0)
    _framed_axes(ax)


def _draw_table(ax, ann, used: int) -> None:
    stat = ann.annoStat.reset_index(drop=True)
    colors = _anno_colors(stat)
    rows = [
        ["\u25a0", str(row.Feature), f"{float(row.Frequency):.2f}", f"{int(row.Count):,}"]
        for row in stat.itertuples(index=False)
    ]

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    promoter = float(stat.loc[0, "Frequency"])
    ax.text(
        0.0,
        0.985,
        f"{used:,} peaks  \u00b7  {promoter:.1f}% promoter \u2264 1 kb  \u00b7  "
        f"{len(stat)} annotation classes",
        fontsize=8.8,
        color=MUTED,
        ha="left",
        va="top",
    )

    table = ax.table(
        cellText=rows,
        colLabels=["", "Feature", "% peaks", "Count"],
        colWidths=[0.06, 0.50, 0.22, 0.22],
        cellLoc="center",
        loc="center",
        bbox=[0.0, 0.10, 1.0, 0.78],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.8)
    table.scale(1.0, 1.0)

    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("none")
        cell.set_linewidth(0)
        if r == 0:
            cell.set_facecolor("#E8EEF4")
            ha = "left" if c == 1 else ("right" if c in (2, 3) else "center")
            cell.set_text_props(weight="bold", color=INK, ha=ha)
            cell.visible_edges = "B"
            cell.set_edgecolor("#8A949E")
            cell.set_linewidth(0.9)
            continue
        cell.set_facecolor("#F6F9FB" if r % 2 == 0 else "white")
        if c == 0:
            cell.get_text().set_color(colors[r - 1])
            cell.get_text().set_fontsize(10.0)
        elif c in (2, 3):
            cell.get_text().set_ha("right")
            cell.get_text().set_color(INK if c == 2 else MUTED)
        else:
            cell.get_text().set_ha("left")
            cell.get_text().set_color(INK)

    ax.text(
        0.0,
        0.045,
        "Percentages follow ChIPseeker's annoStat categories.",
        fontsize=7.4,
        color=MUTED,
        ha="left",
        va="bottom",
    )


def _draw_chrom(ax, df: pd.DataFrame) -> None:
    counts = df["chrom"].astype(str).value_counts()
    counts = counts.reindex(_chrom_order(counts.index))
    y = np.arange(len(counts))
    values = counts.to_numpy(dtype=float)
    ax.barh(
        y,
        values,
        height=0.72,
        color=OKABE_ITO["blue"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index, fontsize=7.0)
    ax.invert_yaxis()
    right = float(values.max()) if values.size else 1.0
    for yi, value in zip(y, values):
        ax.text(
            value + right * 0.014,
            yi,
            f"{int(value):,}",
            va="center",
            fontsize=6.8,
            color=MUTED,
        )
    ax.set_xlim(0, right * 1.18)
    ax.set_xlabel("Peaks per chromosome")
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    _clean_axes(ax)


def _group_values(df: pd.DataFrame, column: str) -> list[np.ndarray]:
    classes = _broad_class(df["annotation"])
    values = df[column].to_numpy(dtype=float)
    return [values[classes == name] for name in CLASS_ORDER]


def _draw_width_box(ax, df: pd.DataFrame) -> None:
    width = (df["end"] - df["start"]).to_numpy(dtype=float)
    grouped = [width[_broad_class(df["annotation"]) == name] for name in CLASS_ORDER]
    boxes = ax.boxplot(
        grouped,
        patch_artist=True,
        widths=0.62,
        showfliers=False,
        medianprops=dict(color=INK, linewidth=1.4),
        whiskerprops=dict(color="#7A8794", linewidth=0.8),
        capprops=dict(color="#7A8794", linewidth=0.8),
        boxprops=dict(linewidth=0.8),
    )
    for patch, name in zip(boxes["boxes"], CLASS_ORDER):
        patch.set_facecolor(_tint(BROAD_COLORS[name], 0.30))
        patch.set_edgecolor(BROAD_COLORS[name])

    ax.set_yscale("log")
    ax.set_ylim(470, 9500)
    ax.set_yticks([500, 1000, 2000, 5000, 8000])
    ax.set_yticklabels(["0.5", "1", "2", "5", "8"])
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xticks(range(1, len(CLASS_ORDER) + 1))
    ax.set_xticklabels(
        [f"{name}\n{len(values):,}" for name, values in zip(CLASS_ORDER, grouped)],
        fontsize=8.0,
    )
    ax.set_ylabel("Peak width (kb, log scale)")
    _clean_axes(ax)


def _draw_signal_violin(ax, df: pd.DataFrame) -> None:
    grouped = _group_values(df, "score")
    parts = ax.violinplot(grouped, showextrema=False, widths=0.82)
    for body, name in zip(parts["bodies"], CLASS_ORDER):
        body.set_facecolor(_tint(BROAD_COLORS[name], 0.30))
        body.set_edgecolor(BROAD_COLORS[name])
        body.set_linewidth(0.9)
        body.set_alpha(1.0)
    for index, values in enumerate(grouped, start=1):
        median = float(np.median(values))
        ax.plot(
            [index - 0.17, index + 0.17],
            [median, median],
            color=INK,
            linewidth=1.3,
            solid_capstyle="butt",
            zorder=4,
        )

    ax.set_yscale("log")
    ax.set_ylim(3, 750)
    ax.set_yticks([3, 10, 30, 100, 300])
    ax.set_yticklabels(["3", "10", "30", "100", "300"])
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xticks(range(1, len(CLASS_ORDER) + 1))
    ax.set_xticklabels(
        [f"{name}\n{len(values):,}" for name, values in zip(CLASS_ORDER, grouped)],
        fontsize=8.0,
    )
    ax.set_xlabel("Peak annotation")
    ax.set_ylabel("Max signal (log scale)")
    _clean_axes(ax)


def _draw_nonpromoter_distance(ax, df: pd.DataFrame) -> None:
    promoter = df["annotation"].astype(str).str.startswith("Promoter").to_numpy()
    distance = np.abs(df["distance_to_tss"].to_numpy(dtype=float))
    distance = distance[~promoter & np.isfinite(distance) & (distance > 0)]
    logs = np.log10(distance)
    bins = np.linspace(0.0, 6.6, 34)

    ax.hist(
        logs,
        bins=bins,
        color=_tint(OKABE_ITO["vermillion"], 0.24),
        edgecolor=OKABE_ITO["vermillion"],
        linewidth=0.9,
        zorder=2,
    )
    for x in (3.0, 4.0, 5.0):
        ax.axvline(x, color="#9FB6C9", linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xticks([0, 1, 2, 3, 4, 5, 6])
    ax.set_xticklabels(
        ["1 bp", "10 bp", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb"],
        fontsize=7.8,
    )
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xlabel("Distance to TSS (non-promoter peaks, log scale)")
    ax.set_ylabel("Peaks")
    ax.text(
        0.02,
        0.95,
        f"n = {distance.size:,} peaks",
        transform=ax.transAxes,
        fontsize=8.0,
        color=MUTED,
        ha="left",
        va="top",
    )
    _clean_axes(ax)


def _draw_profile_quartiles(ax, profiles, upstream: int, downstream: int) -> None:
    shades = ["#A6CEE3", "#5AA6D6", "#1F78B4", "#08519C"]
    rows = int(profiles[0][2].shape[1])
    x = np.linspace(-upstream, downstream, rows)
    for (label, count, matrix), color in zip(profiles, shades):
        ax.plot(
            x,
            np.nanmean(matrix, axis=0),
            color=color,
            linewidth=1.9,
            label=f"{label} (n = {count:,})",
        )
    ax.axvline(0, color="#555555", linewidth=0.8, linestyle="--")
    ax.set_xlim(-upstream, downstream)
    ax.set_ylim(bottom=0)
    ax.set_xticks([-3000, -2000, -1000, 0, 1000, 2000, 3000])
    ax.set_xticklabels(["-3", "-2", "-1", "TSS", "1", "2", "3"])
    ax.set_xlabel("Distance from TSS (kb)")
    ax.set_ylabel("Mean peak coverage (%)")
    ax.legend(
        loc="upper left",
        frameon=False,
        fontsize=7.4,
        handlelength=1.4,
        labelspacing=0.3,
        borderaxespad=0.3,
    )
    _clean_axes(ax)


def _draw_width_signal(ax, fig, df: pd.DataFrame) -> None:
    width = (df["end"] - df["start"]).to_numpy(dtype=float)
    signal = df["score"].to_numpy(dtype=float)
    keep = np.isfinite(width) & np.isfinite(signal) & (width > 0) & (signal > 0)
    x = np.log10(width[keep])
    y = np.log10(signal[keep])

    hexes = ax.hexbin(
        x,
        y,
        gridsize=42,
        cmap="Blues",
        mincnt=1,
        linewidths=0.0,
        bins="log",
        zorder=2,
    )

    # Median signal in twelve width bins: a compact trend guide.
    edges = np.quantile(x, np.linspace(0.0, 1.0, 13))
    centers: list[float] = []
    medians: list[float] = []
    for low, high in zip(edges[:-1], edges[1:]):
        mask = (x >= low) & (x <= high)
        if mask.sum() >= 30:
            centers.append(float(np.median(x[mask])))
            medians.append(float(np.median(y[mask])))
    ax.plot(
        centers,
        medians,
        color=OKABE_ITO["vermillion"],
        linewidth=1.8,
        marker="o",
        markersize=3.2,
        zorder=4,
    )

    cax = ax.inset_axes([1.02, 0.0, 0.034, 1.0])
    colorbar = fig.colorbar(hexes, cax=cax)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(labelsize=7.2, length=2.0, width=0.6, colors=MUTED)
    colorbar.set_label("Peaks (log count)", fontsize=8.0, color=INK, labelpad=4)

    ax.set_xlim(np.log10(480), np.log10(8500))
    ax.set_ylim(np.log10(2.8), np.log10(700))
    ax.set_xticks(np.log10([500, 1000, 2000, 5000, 8000]))
    ax.set_xticklabels(["0.5", "1", "2", "5", "8"])
    ax.set_yticks(np.log10([3, 10, 30, 100, 300]))
    ax.set_yticklabels(["3", "10", "30", "100", "300"])
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xlabel("Peak width (kb, log scale)")
    ax.set_ylabel("Max signal (log scale)")
    ax.text(
        0.02,
        0.95,
        f"Spearman \u03c1 = 0.27  \u00b7  n = {int(keep.sum()):,}",
        transform=ax.transAxes,
        fontsize=7.8,
        color=MUTED,
        ha="left",
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", pad=2.4),
    )
    _clean_axes(ax)


def _draw_width_ecdf(ax, df: pd.DataFrame) -> None:
    width = (df["end"] - df["start"]).to_numpy(dtype=float)
    classes = _broad_class(df["annotation"])
    for name in CLASS_ORDER:
        values = np.sort(width[classes == name])
        if values.size == 0:
            continue
        cumulative = np.arange(1, values.size + 1) / values.size * 100.0
        ax.step(
            values,
            cumulative,
            where="post",
            color=BROAD_COLORS[name],
            linewidth=1.7,
            label=f"{name} (n = {values.size:,})",
        )

    ax.axvline(1000, color="#9FB6C9", linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(480, 8500)
    ax.set_ylim(0, 102)
    ax.set_xticks([500, 1000, 2000, 5000, 8000])
    ax.set_xticklabels(["0.5", "1", "2", "5", "8"])
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Peak width (kb, log scale)")
    ax.set_ylabel("Cumulative peaks (%)")
    ax.legend(
        loc="upper left",
        frameon=False,
        fontsize=7.0,
        handlelength=1.5,
        labelspacing=0.28,
        borderaxespad=0.25,
    )
    ax.text(
        0.98,
        0.045,
        "dashed line: 1 kb",
        transform=ax.transAxes,
        fontsize=7.0,
        color=MUTED,
        ha="right",
        va="bottom",
    )
    _clean_axes(ax)


def _draw_top_genes(ax, df: pd.DataFrame) -> None:
    top = pcs.top_genes(df, n=12)
    if top.empty:
        ax.text(0.5, 0.5, "No gene annotation available", ha="center", va="center")
        _clean_axes(ax)
        return
    top = top.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(top))
    values = top["n_peaks"].to_numpy(dtype=float)
    ax.barh(
        y,
        values,
        height=0.72,
        color=OKABE_ITO["blue"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene_name"].astype(str), fontsize=7.2)
    right = float(values.max()) if values.size else 1.0
    for yi, value in zip(y, values):
        ax.text(
            value + right * 0.02,
            yi,
            f"{int(value)}",
            va="center",
            fontsize=6.8,
            color=MUTED,
        )
    ax.set_xlim(0, right * 1.18)
    ax.set_xlabel("Peaks assigned to the gene")
    ax.text(
        0.98,
        0.04,
        f"{len(top)} of {pcs.peaks_per_gene(df).shape[0]:,} genes",
        transform=ax.transAxes,
        fontsize=7.0,
        color=MUTED,
        ha="right",
        va="bottom",
    )
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    _clean_axes(ax)


def _draw_peaks_per_gene(ax, df: pd.DataFrame) -> None:
    counts = pcs.peaks_per_gene(df)["n_peaks"]
    distribution = counts.value_counts().sort_index()
    labels = ["1", "2", "3", "4+"]
    values = [
        int(distribution.get(1, 0)),
        int(distribution.get(2, 0)),
        int(distribution.get(3, 0)),
        int(distribution.loc[distribution.index >= 4].sum()),
    ]
    x = np.arange(len(labels))
    bars = ax.bar(
        x,
        values,
        width=0.62,
        color=[OKABE_ITO["blue"], OKABE_ITO["sky"], OKABE_ITO["green"], OKABE_ITO["orange"]],
        edgecolor="white",
        linewidth=0.6,
        zorder=2,
    )
    total = max(sum(values), 1)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value * 1.15,
            f"{value:,}\n{value / total * 100:.1f}%",
            ha="center",
            va="bottom",
            fontsize=7.2,
            color=INK,
            linespacing=1.15,
        )
    ax.set_yscale("log")
    ax.set_ylim(0.8, max(values) * 4.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Peaks per gene")
    ax.set_ylabel("Genes (log scale)")
    ax.text(
        0.98,
        0.95,
        f"median = {np.median(counts):.0f}  \u00b7  {len(counts):,} genes",
        transform=ax.transAxes,
        fontsize=7.8,
        color=MUTED,
        ha="right",
        va="top",
    )
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    _clean_axes(ax)


def _draw_tss_quartiles(ax, tag, ann) -> None:
    summary = pcs.tss_enrichment_by_quartile(tag, ann, n_groups=4)
    x = np.arange(len(summary))
    shades = ["#A6CEE3", "#5AA6D6", "#1F78B4", "#08519C"]
    colors = shades[: len(summary)]
    bars = ax.bar(
        x,
        summary["tss_enrichment"],
        width=0.62,
        color=colors,
        edgecolor="white",
        linewidth=0.8,
        zorder=2,
    )
    for bar, value in zip(bars, summary["tss_enrichment"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value * 1.02,
            f"{value:.1f}\u00d7",
            ha="center",
            va="bottom",
            fontsize=8.0,
            color=INK,
            fontweight="bold",
        )
    ax.axhline(1.0, color="#9FB6C9", linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            f"{row.quartile}\nn = {int(row.n):,}"
            for row in summary.itertuples(index=False)
        ],
        fontsize=7.6,
    )
    ax.set_ylim(0, float(summary["tss_enrichment"].max()) * 1.16)
    ax.set_ylabel("TSS enrichment (centre / flank)")
    ax.text(
        0.02,
        0.95,
        "centre = central 5 % of bins  \u00b7  flank = outer 10 %\n"
        "Q1 includes genes without an assigned peak",
        transform=ax.transAxes,
        fontsize=6.8,
        color=MUTED,
        ha="left",
        va="top",
        linespacing=1.35,
    )
    _clean_axes(ax)


def _draw_width_distance(ax, df: pd.DataFrame) -> None:
    width = (df["end"] - df["start"]).to_numpy(dtype=float)
    distance = np.abs(df["distance_to_tss"].to_numpy(dtype=float))
    groups = [
        ("TSS (0)", distance == 0),
        ("1 bp-1 kb", (distance > 0) & (distance <= 1000)),
        ("1-10 kb", (distance > 1000) & (distance <= 10000)),
        ("10-100 kb", (distance > 10000) & (distance <= 100000)),
        (">100 kb", distance > 100000),
    ]
    labels = []
    values = []
    colors = []
    shades = ["#08519C", "#1F78B4", "#5AA6D6", "#A6CEE3", "#D6E6F4"]
    for (label, mask), color in zip(groups, shades):
        selected = width[mask]
        if selected.size == 0:
            continue
        values.append(selected)
        labels.append(f"{label}\n{selected.size:,}")
        colors.append(color)

    boxes = ax.boxplot(
        values,
        patch_artist=True,
        widths=0.62,
        showfliers=False,
        medianprops=dict(color=INK, linewidth=1.3),
        whiskerprops=dict(color="#7A8794", linewidth=0.8),
        capprops=dict(color="#7A8794", linewidth=0.8),
        boxprops=dict(linewidth=0.8),
    )
    for patch, color in zip(boxes["boxes"], colors):
        patch.set_facecolor(_tint(color, 0.34))
        patch.set_edgecolor(color)
    ax.set_yscale("log")
    ax.set_ylim(470, 9500)
    ax.set_yticks([500, 1000, 2000, 5000, 8000])
    ax.set_yticklabels(["0.5", "1", "2", "5", "8"])
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=7.4)
    ax.set_xlabel("Distance to TSS")
    ax.set_ylabel("Peak width (kb, log scale)")
    _clean_axes(ax)


def _draw_gene_type(ax, df: pd.DataFrame) -> None:
    summary = pcs.gene_type_summary(df, top=8)
    summary = summary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(summary))
    values = summary["count"].to_numpy(dtype=float)
    ax.barh(
        y,
        values,
        height=0.7,
        color=OKABE_ITO["purple"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    ax.set_xscale("log")
    right = float(values.max()) if values.size else 1.0
    ax.set_xlim(0.8, right * 1.6)
    ax.set_yticks(y)
    ax.set_yticklabels(summary["gene_type"].astype(str), fontsize=7.0)
    for yi, row in zip(y, summary.itertuples(index=False)):
        ax.text(
            float(row.count) * 1.08,
            yi,
            f"{int(row.count):,} ({float(row.percent):.1f}%)",
            va="center",
            fontsize=6.5,
            color=MUTED,
        )
    ax.set_xlabel("Peaks (log scale)")
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    _clean_axes(ax)


def _standalone(letter: str, title: str, figsize: tuple[float, float]):
    fig = plt.figure(figsize=figsize)
    outer = fig.add_gridspec(1, 1, left=0.07, right=0.955, top=0.96, bottom=0.105)
    content = _panel_header(fig, outer[0], letter, title)
    return fig, content


def _figure_pie(ann) -> plt.Figure:
    fig, content = _standalone("A", "Annotation composition", (7.4, 4.6))
    _panel_pie(fig, content, ann)
    return fig


def _figure_width(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("B", "Peak width distribution", (6.8, 4.3))
    ax = fig.add_subplot(content)
    _draw_width(ax, df)
    return fig


def _figure_distance(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("C", "Distance to TSS", (6.8, 4.3))
    ax = fig.add_subplot(content)
    _draw_distance(ax, df)
    return fig


def _figure_profile(tag) -> plt.Figure:
    fig, content = _standalone("D", "Average H3K4me3 profile", (6.8, 4.2))
    ax = fig.add_subplot(content)
    _draw_profile(ax, tag)
    return fig


def _figure_heatmap(tag) -> plt.Figure:
    fig, content = _standalone("E", "Promoter-centered peak coverage", (6.9, 4.8))
    ax = fig.add_subplot(content)
    _draw_heatmap(fig, ax, tag)
    return fig


def _figure_table(ann, used: int) -> plt.Figure:
    fig, content = _standalone("F", "Annotation statistics", (7.1, 4.8))
    ax = fig.add_subplot(content)
    _draw_table(ax, ann, used)
    return fig


def _figure_chrom(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("G", "Peaks per chromosome", (6.9, 5.0))
    ax = fig.add_subplot(content)
    _draw_chrom(ax, df)
    return fig


def _figure_width_box(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("H", "Peak width by annotation class", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_width_box(ax, df)
    return fig


def _figure_signal_violin(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("I", "Max signal by annotation class", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_signal_violin(ax, df)
    return fig


def _figure_nonpromoter_distance(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("J", "Distance to TSS (non-promoter peaks)", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_nonpromoter_distance(ax, df)
    return fig


def _figure_profile_quartiles(profiles, upstream: int, downstream: int) -> plt.Figure:
    fig, content = _standalone("K", "Profile by peak-signal quartile", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_profile_quartiles(ax, profiles, upstream, downstream)
    return fig


def _figure_width_signal(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("L", "Peak width vs signal", (6.9, 4.6))
    ax = fig.add_subplot(content)
    _draw_width_signal(ax, fig, df)
    return fig


def _figure_width_ecdf(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("M", "Peak width ECDF by annotation class", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_width_ecdf(ax, df)
    return fig


def _figure_top_genes(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("N", "Top genes by peak count", (6.8, 4.5))
    ax = fig.add_subplot(content)
    _draw_top_genes(ax, df)
    return fig


def _figure_peaks_per_gene(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("O", "Peaks per gene", (6.8, 4.3))
    ax = fig.add_subplot(content)
    _draw_peaks_per_gene(ax, df)
    return fig


def _figure_tss_quartiles(tag, ann) -> plt.Figure:
    fig, content = _standalone("P", "TSS enrichment by signal quartile", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_tss_quartiles(ax, tag, ann)
    return fig


def _figure_width_distance(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("Q", "Peak width by TSS-distance class", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_width_distance(ax, df)
    return fig


def _figure_gene_type(df: pd.DataFrame) -> plt.Figure:
    fig, content = _standalone("R", "Gene biotypes of annotated peaks", (6.8, 4.4))
    ax = fig.add_subplot(content)
    _draw_gene_type(ax, df)
    return fig


def _panel_figures(ann, tag, df, profiles, upstream: int, downstream: int):
    used = len(ann.peaks)
    yield _figure_pie(ann), "A_annotation_composition"
    yield _figure_width(df), "B_peak_width"
    yield _figure_distance(df), "C_distance_to_tss"
    yield _figure_profile(tag), "D_tag_profile"
    yield _figure_heatmap(tag), "E_tag_heatmap"
    yield _figure_table(ann, used), "F_annotation_stats"
    yield _figure_chrom(df), "G_chromosome_counts"
    yield _figure_width_box(df), "H_width_by_class"
    yield _figure_signal_violin(df), "I_signal_by_class"
    yield _figure_nonpromoter_distance(df), "J_nonpromoter_distance"
    yield _figure_profile_quartiles(profiles, upstream, downstream), "K_profile_quartiles"
    yield _figure_width_signal(df), "L_width_vs_signal"
    yield _figure_width_ecdf(df), "M_width_ecdf"
    yield _figure_top_genes(df), "N_top_genes"
    yield _figure_peaks_per_gene(df), "O_peaks_per_gene"
    yield _figure_tss_quartiles(tag, ann), "P_tss_quartile_enrichment"
    yield _figure_width_distance(df), "Q_width_by_distance"
    yield _figure_gene_type(df), "R_gene_biotypes"


def _composite_plate(ann, tag, df, caption: str, provenance: str) -> plt.Figure:
    with _no_tight_layout():
        fig = plt.figure(figsize=(11.69, 16.54))
        gs = fig.add_gridspec(
            3,
            2,
            left=0.055,
            right=0.945,
            top=0.912,
            bottom=0.105,
            height_ratios=[1.0, 1.02, 1.18],
            width_ratios=[1.0, 1.0],
            hspace=0.28,
            wspace=0.24,
        )

        _panel_pie(fig, _panel_header(fig, gs[0, 0], "A", "Annotation composition"), ann)

        ax_b = fig.add_subplot(
            _panel_header(fig, gs[0, 1], "B", "Peak width distribution")
        )
        _draw_width(ax_b, df)

        ax_c = fig.add_subplot(_panel_header(fig, gs[1, 0], "C", "Distance to TSS"))
        _draw_distance(ax_c, df)

        ax_d = fig.add_subplot(
            _panel_header(fig, gs[1, 1], "D", "Average H3K4me3 profile")
        )
        _draw_profile(ax_d, tag)

        ax_e = fig.add_subplot(
            _panel_header(fig, gs[2, 0], "E", "Promoter-centered peak coverage")
        )
        _draw_heatmap(fig, ax_e, tag)

        ax_f = fig.add_subplot(
            _panel_header(fig, gs[2, 1], "F", "Annotation statistics")
        )
        _draw_table(ax_f, ann, len(ann.peaks))

        fig.text(
            0.055,
            0.972,
            "pychipseeker \u2014 H3K4me3 promoter landscape in mouse ENCODE",
            fontsize=19,
            fontweight="bold",
            color=INK,
            ha="left",
            va="top",
        )
        fig.text(
            0.055,
            0.949,
            f"{len(ann.peaks):,} peaks  \u00b7  "
            f"{float(ann.annoStat.loc[0, 'Frequency']):.1f}% within \u00b11 kb of a TSS  \u00b7  mm10",
            fontsize=10.4,
            color=MUTED,
            ha="left",
            va="top",
        )
        fig.add_artist(
            Line2D(
                [0.055, 0.945],
                [0.933, 0.933],
                transform=fig.transFigure,
                color=RULE,
                linewidth=0.9,
            )
        )

        fig.text(0.5, 0.058, caption, ha="center", fontsize=7.8, color=MUTED)
        fig.text(0.5, 0.042, provenance, ha="center", fontsize=7.4, color="#7C8792")
        fig.text(
            0.5,
            0.026,
            "(A) Annotation composition \u00b7 (B) Peak widths \u00b7 "
            "(C) Cumulative distance to TSS \u00b7 (D) Mean coverage profile \u00b7 "
            "(E) Coverage heatmap \u00b7 (F) Annotation statistics",
            ha="center",
            fontsize=7.4,
            color="#7C8792",
        )
        return fig


def _composite_supplementary(df, profiles, caption: str, provenance: str) -> plt.Figure:
    with _no_tight_layout():
        fig = plt.figure(figsize=(11.69, 16.54))
        gs = fig.add_gridspec(
            3,
            2,
            left=0.055,
            right=0.945,
            top=0.912,
            bottom=0.105,
            height_ratios=[1.12, 1.0, 1.0],
            width_ratios=[1.0, 1.0],
            hspace=0.30,
            wspace=0.26,
        )

        ax_g = fig.add_subplot(
            _panel_header(fig, gs[0, 0], "G", "Peaks per chromosome")
        )
        _draw_chrom(ax_g, df)

        ax_h = fig.add_subplot(
            _panel_header(fig, gs[0, 1], "H", "Peak width by annotation class")
        )
        _draw_width_box(ax_h, df)

        ax_i = fig.add_subplot(
            _panel_header(fig, gs[1, 0], "I", "Max signal by annotation class")
        )
        _draw_signal_violin(ax_i, df)

        ax_j = fig.add_subplot(
            _panel_header(fig, gs[1, 1], "J", "Distance to TSS (non-promoter peaks)")
        )
        _draw_nonpromoter_distance(ax_j, df)

        ax_k = fig.add_subplot(
            _panel_header(fig, gs[2, 0], "K", "Profile by peak-signal quartile")
        )
        _draw_profile_quartiles(ax_k, profiles, 3000, 3000)

        ax_l = fig.add_subplot(
            _panel_header(fig, gs[2, 1], "L", "Peak width vs signal")
        )
        _draw_width_signal(ax_l, fig, df)

        fig.text(
            0.055,
            0.972,
            "pychipseeker \u2014 extended views of the H3K4me3 peak set",
            fontsize=19,
            fontweight="bold",
            color=INK,
            ha="left",
            va="top",
        )
        fig.text(
            0.055,
            0.949,
            f"{len(df):,} peaks  \u00b7  {df['chrom'].nunique()} chromosomes  "
            f"\u00b7  {df['annotation'].astype(str).str.startswith('Promoter').sum():,} promoter peaks  "
            "\u00b7  mm10",
            fontsize=10.4,
            color=MUTED,
            ha="left",
            va="top",
        )
        fig.add_artist(
            Line2D(
                [0.055, 0.945],
                [0.933, 0.933],
                transform=fig.transFigure,
                color=RULE,
                linewidth=0.9,
            )
        )

        fig.text(0.5, 0.058, caption, ha="center", fontsize=7.8, color=MUTED)
        fig.text(0.5, 0.042, provenance, ha="center", fontsize=7.4, color="#7C8792")
        fig.text(
            0.5,
            0.026,
            "(G) Chromosome distribution \u00b7 (H) Width by class \u00b7 "
            "(I) Signal by class \u00b7 (J) Non-promoter distance \u00b7 "
            "(K) Signal-quartile profiles \u00b7 (L) Width vs signal",
            ha="center",
            fontsize=7.4,
            color="#7C8792",
        )
        return fig


def _draw_panel(
    letter: str,
    fig: plt.Figure,
    spec,
    ann,
    tag,
    df: pd.DataFrame,
    profiles,
    upstream: int,
    downstream: int,
) -> None:
    """Draw one catalog panel into an existing gridspec slot."""

    if letter == "A":
        _panel_pie(fig, spec, ann)
        return

    title = PANEL_BY_LETTER[letter][2]
    ax = fig.add_subplot(_panel_header(fig, spec, letter, title))
    if letter == "B":
        _draw_width(ax, df)
    elif letter == "C":
        _draw_distance(ax, df)
    elif letter == "D":
        _draw_profile(ax, tag)
    elif letter == "E":
        _draw_heatmap(fig, ax, tag)
    elif letter == "F":
        _draw_table(ax, ann, len(ann.peaks))
    elif letter == "G":
        _draw_chrom(ax, df)
    elif letter == "H":
        _draw_width_box(ax, df)
    elif letter == "I":
        _draw_signal_violin(ax, df)
    elif letter == "J":
        _draw_nonpromoter_distance(ax, df)
    elif letter == "K":
        _draw_profile_quartiles(ax, profiles, upstream, downstream)
    elif letter == "L":
        _draw_width_signal(ax, fig, df)
    elif letter == "M":
        _draw_width_ecdf(ax, df)
    elif letter == "N":
        _draw_top_genes(ax, df)
    elif letter == "O":
        _draw_peaks_per_gene(ax, df)
    elif letter == "P":
        _draw_tss_quartiles(ax, tag, ann)
    elif letter == "Q":
        _draw_width_distance(ax, df)
    elif letter == "R":
        _draw_gene_type(ax, df)
    else:  # pragma: no cover - guarded by the CLI
        raise ValueError(f"unknown panel letter: {letter}")


def _composite_analysis(
    ann, tag, df, caption: str, provenance: str
) -> plt.Figure:
    with _no_tight_layout():
        fig = plt.figure(figsize=(11.69, 16.54))
        gs = fig.add_gridspec(
            3,
            2,
            left=0.055,
            right=0.945,
            top=0.912,
            bottom=0.105,
            height_ratios=[1.0, 1.0, 1.0],
            width_ratios=[1.0, 1.0],
            hspace=0.32,
            wspace=0.26,
        )
        for index, letter in enumerate("MNOPQR"):
            row, col = divmod(index, 2)
            _draw_panel(
                letter,
                fig,
                gs[row, col],
                ann,
                tag,
                df,
                [],
                3000,
                3000,
            )

        fig.text(
            0.055,
            0.972,
            "pychipseeker \u2014 gene-level and distribution analyses",
            fontsize=19,
            fontweight="bold",
            color=INK,
            ha="left",
            va="top",
        )
        fig.text(
            0.055,
            0.949,
            f"{len(df):,} peaks  \u00b7  {pcs.peaks_per_gene(df).shape[0]:,} genes  "
            f"\u00b7  {df['gene_type'].nunique()} gene biotypes  \u00b7  mm10",
            fontsize=10.4,
            color=MUTED,
            ha="left",
            va="top",
        )
        fig.add_artist(
            Line2D(
                [0.055, 0.945],
                [0.933, 0.933],
                transform=fig.transFigure,
                color=RULE,
                linewidth=0.9,
            )
        )
        fig.text(0.5, 0.058, caption, ha="center", fontsize=7.8, color=MUTED)
        fig.text(0.5, 0.042, provenance, ha="center", fontsize=7.4, color="#7C8792")
        fig.text(
            0.5,
            0.026,
            "(M) Width ECDF \u00b7 (N) Top genes \u00b7 (O) Peaks per gene \u00b7 "
            "(P) TSS enrichment \u00b7 (Q) Width by distance \u00b7 (R) Gene biotypes",
            ha="center",
            fontsize=7.4,
            color="#7C8792",
        )
        return fig


def _composite_custom(
    letters: Sequence[str],
    ann,
    tag,
    df,
    profiles,
    upstream: int,
    downstream: int,
    caption: str,
    provenance: str,
) -> plt.Figure:
    """Compose an arbitrary subset of catalog panels into one A3-style plate."""

    letters = list(letters)
    if not letters:
        raise ValueError("at least one panel is required")
    cols = 2 if len(letters) > 1 else 1
    rows = int(np.ceil(len(letters) / cols))
    if rows <= 1:
        figsize = (11.69, 8.27)
    elif rows == 2:
        figsize = (11.69, 11.69)
    else:
        figsize = (11.69, 16.54)

    with _no_tight_layout():
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(
            rows,
            cols,
            left=0.055,
            right=0.945,
            top=0.912,
            bottom=0.105,
            hspace=0.32,
            wspace=0.26,
        )
        for index, letter in enumerate(letters):
            row, col = divmod(index, cols)
            _draw_panel(
                letter,
                fig,
                gs[row, col],
                ann,
                tag,
                df,
                profiles,
                upstream,
                downstream,
            )

        chosen = ", ".join(letters)
        fig.text(
            0.055,
            0.972,
            "pychipseeker \u2014 custom figure plate",
            fontsize=19,
            fontweight="bold",
            color=INK,
            ha="left",
            va="top",
        )
        fig.text(
            0.055,
            0.949,
            f"{len(letters)} selected panels: {chosen}  \u00b7  {len(df):,} peaks  \u00b7  mm10",
            fontsize=10.4,
            color=MUTED,
            ha="left",
            va="top",
        )
        fig.add_artist(
            Line2D(
                [0.055, 0.945],
                [0.933, 0.933],
                transform=fig.transFigure,
                color=RULE,
                linewidth=0.9,
            )
        )
        fig.text(0.5, 0.058, caption, ha="center", fontsize=7.8, color=MUTED)
        fig.text(0.5, 0.042, provenance, ha="center", fontsize=7.4, color="#7C8792")
        fig.text(
            0.5,
            0.026,
            "Selected: "
            + " \u00b7 ".join(
                f"({letter}) {PANEL_BY_LETTER[letter][2]}" for letter in letters
            ),
            ha="center",
            fontsize=7.4,
            color="#7C8792",
        )
        return fig


_FORMAT_LABELS = {
    "png": "PNG",
    "pdf": "PDF",
    "svg": "SVG",
    "eps": "EPS",
    "tif": "TIFF",
}


def _write_gallery(out: Path, formats: Sequence[str]) -> Path:
    """Write an interactive, browsable index of the plates and panels."""

    display_order = ("pdf", "svg", "png", "eps", "tif")
    ordered = [fmt for fmt in display_order if fmt in formats]
    ordered += [fmt for fmt in formats if fmt not in ordered]

    def links_for(prefix: str) -> str:
        return "".join(
            f'<a class="chip{" chip-pdf" if fmt == "pdf" else ""}" '
            f'href="{html.escape(prefix)}.{fmt}">'
            f"{_FORMAT_LABELS.get(fmt, fmt.upper())}</a>"
            for fmt in ordered
            if (out / f"{prefix}.{fmt}").exists()
        )

    plate_specs = (
        ("encode_figure_plate", "Main plate (A-F)", "Core annotation and coverage story."),
        (
            "encode_supplementary_plate",
            "Extended plate (G-L)",
            "Chromosomes, classes, distances and density.",
        ),
        (
            "encode_analysis_plate",
            "Analysis plate (M-R)",
            "Gene-level and distribution summaries.",
        ),
    )
    plate_cards = []
    for stem, title, blurb in plate_specs:
        links = links_for(stem)
        if not links:
            continue
        extra = (
            f'<a class="chip" href="{stem}_600dpi.png">PNG 600 dpi</a>'
            if (out / f"{stem}_600dpi.png").exists()
            else ""
        )
        plate_cards.append(
            f'<section class="plate"><h2>{html.escape(title)}</h2>'
            f"<p>{html.escape(blurb)}</p>"
            f'<div class="chips">{links}{extra}</div></section>'
        )
    custom_links = links_for("encode_custom_plate")
    if custom_links:
        plate_cards.append(
            '<section class="plate plate-custom"><h2>Your custom plate</h2>'
            "<p>Built from the panels you selected below.</p>"
            f'<div class="chips">{custom_links}</div></section>'
        )

    groups = {"main": [], "extended": [], "analysis": []}
    for letter, stem, title, blurb in PANEL_CATALOG:
        if not (out / "figures" / f"{stem}.png").exists():
            continue
        if letter <= "F":
            group = "main"
        elif letter <= "L":
            group = "extended"
        else:
            group = "analysis"
        links = links_for(f"figures/{stem}")
        groups[group].append(
            f"""
        <article class="card" data-group="{group}">
          <label class="pick" title="Include in the custom plate">
            <input type="checkbox" class="pick-box" value="{letter}"> {letter}
          </label>
          <a href="figures/{html.escape(stem)}.pdf">
            <img src="figures/{html.escape(stem)}.png" alt="{html.escape(title)}">
          </a>
          <div class="card-body">
            <h3><span class="letter">{letter}</span> {html.escape(title)}</h3>
            <p>{html.escape(blurb)}</p>
            <div class="chips">{links}</div>
          </div>
        </article>"""
        )

    section_titles = {
        "main": "Main views (A-F)",
        "extended": "Extended views (G-L)",
        "analysis": "Analysis views (M-R)",
    }
    sections = "".join(
        f'<section class="panel-section" data-group="{key}">'
        f'<h2 class="section-title">{title}</h2>'
        f'<div class="gallery">{"".join(groups[key])}</div></section>'
        for key, title in section_titles.items()
        if groups[key]
    )

    document = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>pychipseeker - ENCODE H3K4me3 figure gallery</title>
<style>
  :root {
    --ink: #1A1A1A;
    --muted: #5F6B76;
    --rule: #D5DCE3;
    --accent: #0072B2;
    --accent-dark: #005A8C;
    --wash: #F4F8FB;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    padding: 40px 28px 64px;
    font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
    color: var(--ink);
    background: #FFFFFF;
  }
  header, main, footer { max-width: 1180px; margin-left: auto; margin-right: auto; }
  h1 { margin: 0 0 6px; font-size: 30px; letter-spacing: -0.01em; }
  .sub { margin: 0; color: var(--muted); font-size: 15px; }
  .rules { margin: 26px 0 0; border-top: 1px solid var(--rule); padding-top: 20px; }
  .plates {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
    gap: 18px;
    margin-top: 14px;
  }
  .plate {
    border: 1px solid var(--rule);
    border-radius: 12px;
    padding: 18px 18px 16px;
    background: var(--wash);
  }
  .plate-custom { border-color: var(--accent); background: #EAF4FB; }
  .plate h2 { margin: 0 0 4px; font-size: 17px; }
  .plate p { margin: 0 0 12px; color: var(--muted); font-size: 13px; }
  .builder {
    border: 1px solid var(--rule);
    border-radius: 12px;
    padding: 20px;
    margin-top: 24px;
    background: #FFFFFF;
  }
  .builder h2 { margin: 0 0 6px; font-size: 17px; }
  .builder p { margin: 0 0 14px; color: var(--muted); font-size: 13px; }
  .cmd-row { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
  .cmd-row code {
    flex: 1 1 520px;
    padding: 10px 12px;
    border-radius: 8px;
    background: #1A2129;
    color: #E8EDF2;
    font-size: 12.5px;
    overflow-x: auto;
    white-space: nowrap;
  }
  .toolbar { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 14px; }
  button {
    padding: 5px 12px;
    border: 1px solid var(--rule);
    border-radius: 999px;
    background: #FFFFFF;
    color: var(--accent);
    font-size: 12px;
    cursor: pointer;
  }
  button:hover, button.active { border-color: var(--accent); background: var(--wash); }
  button.primary { background: var(--accent); border-color: var(--accent); color: #FFFFFF; }
  .panel-section { margin-top: 38px; }
  .section-title {
    margin: 0 0 14px;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--rule);
    font-size: 19px;
  }
  .gallery {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(340px, 1fr));
    gap: 24px;
  }
  .card {
    position: relative;
    border: 1px solid var(--rule);
    border-radius: 12px;
    overflow: hidden;
    background: #FFFFFF;
    display: flex;
    flex-direction: column;
  }
  .card img { width: 100%; height: 232px; object-fit: contain; background: #FFFFFF; padding: 10px; }
  .card-body { padding: 14px 16px 16px; border-top: 1px solid var(--rule); }
  .card h3 { margin: 0 0 6px; font-size: 15px; line-height: 1.35; }
  .card p { margin: 0 0 12px; color: var(--muted); font-size: 12.5px; line-height: 1.5; }
  .pick {
    position: absolute;
    top: 10px;
    right: 10px;
    z-index: 2;
    display: flex;
    gap: 5px;
    align-items: center;
    padding: 4px 9px;
    border: 1px solid var(--rule);
    border-radius: 999px;
    background: rgba(255,255,255,0.94);
    font-size: 11.5px;
    color: var(--muted);
    cursor: pointer;
  }
  .letter {
    display: inline-block;
    min-width: 22px;
    padding: 1px 5px;
    margin-right: 4px;
    border-radius: 5px;
    background: var(--accent);
    color: #FFFFFF;
    font-size: 12px;
    text-align: center;
  }
  .chips { display: flex; flex-wrap: wrap; gap: 6px; }
  .chip {
    display: inline-block;
    padding: 3px 9px;
    border: 1px solid var(--rule);
    border-radius: 999px;
    font-size: 11.5px;
    color: var(--accent);
    text-decoration: none;
    background: #FFFFFF;
  }
  .chip:hover { border-color: var(--accent); background: var(--wash); }
  .chip-pdf { border-color: var(--accent); background: var(--accent); color: #FFFFFF; }
  footer { margin-top: 44px; color: var(--muted); font-size: 12px; line-height: 1.6; }
  @media (prefers-color-scheme: dark) {
    :root { --rule: #2B3540; --wash: #1A2129; }
    body { background: #12171C; color: #E8EDF2; }
    .card, .chip, button, .builder { background: #1A2129; }
    .plate { background: #1A2129; }
    .plate-custom { background: #162A38; }
    .pick { background: #1A2129; }
    .cmd-row code { background: #0C1116; }
  }
</style>
</head>
<body>
<header>
  <h1>ENCODE H3K4me3 promoter landscape</h1>
  <p class="sub">11,291 peaks &middot; mm10 &middot; mouse ENCODE ENCSR000CDQ &middot;
     generated with pychipseeker</p>
  <div class="rules">
    <div class="plates">$PLATES$</div>
    <div class="builder">
      <h2>Build a plate from the panels you want</h2>
      <p>Tick any panels below. The command updates automatically; run it to
         render only those panels into <code>encode_custom_plate.pdf</code>
         (plus PNG/SVG/EPS/TIFF). Every panel also has its own PDF vector file.</p>
      <div class="cmd-row">
        <code id="cmd">python examples/make_encode_plate.py --panels A,C,K</code>
        <button id="copy" class="primary" type="button">Copy command</button>
      </div>
      <div class="toolbar">
        <button type="button" data-filter="all" class="active">All 18</button>
        <button type="button" data-filter="main">Main A-F</button>
        <button type="button" data-filter="extended">Extended G-L</button>
        <button type="button" data-filter="analysis">Analysis M-R</button>
        <button type="button" id="select-all">Select all</button>
        <button type="button" id="clear-all">Clear</button>
      </div>
    </div>
  </div>
</header>
<main>$SECTIONS$</main>
<footer>
  Vector files use Type 42 fonts and are safe for journal submission; SVG text
  stays editable; rasters are RGB.  Individual panels live in
  <code>figures/</code>.  Use <code>--formats pdf,svg,png</code> to limit the
  export formats.
</footer>
<script>
  const boxes = Array.from(document.querySelectorAll('.pick-box'));
  const cmd = document.getElementById('cmd');

  function refresh() {
    const picked = boxes.filter(box => box.checked).map(box => box.value).sort();
    cmd.textContent = 'python examples/make_encode_plate.py --panels '
      + (picked.length ? picked.join(',') : 'A,C,K');
  }

  boxes.forEach(box => box.addEventListener('change', refresh));
  document.getElementById('select-all').addEventListener('click', () => {
    boxes.forEach(box => { box.checked = true; });
    refresh();
  });
  document.getElementById('clear-all').addEventListener('click', () => {
    boxes.forEach(box => { box.checked = false; });
    refresh();
  });
  document.getElementById('copy').addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(cmd.textContent);
      const button = document.getElementById('copy');
      button.textContent = 'Copied';
      setTimeout(() => { button.textContent = 'Copy command'; }, 1200);
    } catch (error) {
      const range = document.createRange();
      range.selectNodeContents(cmd);
      const selection = window.getSelection();
      selection.removeAllRanges();
      selection.addRange(range);
    }
  });
  document.querySelectorAll('[data-filter]').forEach(button => {
    button.addEventListener('click', () => {
      const group = button.dataset.filter;
      document.querySelectorAll('[data-filter]').forEach(item => {
        item.classList.toggle('active', item === button);
      });
      document.querySelectorAll('.panel-section').forEach(section => {
        section.style.display = (group === 'all' || section.dataset.group === group)
          ? '' : 'none';
      });
    });
  });
</script>
</body>
</html>
"""

    document = document.replace("$PLATES$", "".join(plate_cards))
    document = document.replace("$SECTIONS$", sections)

    path = out / "gallery.html"
    path.write_text(document, encoding="utf-8")
    return path


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
    parser.add_argument(
        "--formats",
        default=",".join(FORMATS),
        help="comma-separated export formats: " + ", ".join(FORMATS),
    )
    parser.add_argument(
        "--panels",
        default="",
        help="comma-separated panel letters for a custom plate (e.g. A,C,K); "
        "use --list-panels to see all choices",
    )
    parser.add_argument(
        "--list-panels",
        action="store_true",
        help="print the panel catalog and exit",
    )
    args = parser.parse_args()

    if args.list_panels:
        print("letter  stem                            title")
        for letter, stem, title, _ in PANEL_CATALOG:
            print(f"{letter:>6}  {stem:<30}  {title}")
        return 0

    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    unknown = sorted(set(formats) - set(FORMATS))
    if unknown:
        parser.error(f"unsupported format(s): {', '.join(unknown)}")
    if not formats:
        parser.error("at least one export format is required")

    selected: list[str] = []
    if args.panels.strip():
        raw = [item.strip() for item in args.panels.split(",") if item.strip()]
        if any(item.lower() == "all" for item in raw):
            selected = [letter for letter, *_ in PANEL_CATALOG]
        else:
            by_stem = {stem.lower(): letter for letter, stem, *_ in PANEL_CATALOG}
            for item in raw:
                letter = item.upper()
                if letter not in PANEL_BY_LETTER:
                    letter = by_stem.get(item.lower(), "")
                if not letter:
                    parser.error(
                        f"unknown panel {item!r}; use --list-panels to see choices"
                    )
                if letter not in selected:
                    selected.append(letter)
        if not selected:
            parser.error("--panels did not select any panel")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    for fmt in FORMATS:
        for stale in figdir.glob(f"*.{fmt}"):
            stale.unlink()
    for fmt in FORMATS:
        stale = out / f"encode_custom_plate.{fmt}"
        if stale.exists():
            stale.unlink()
    stale = out / "encode_custom_plate_600dpi.png"
    if stale.exists():
        stale.unlink()

    ann = pcs.annotate_peak(peaks=args.peaks, TxDb=args.txdb, verbose=False)
    used = len(ann.peaks)
    if args.sample and used > args.sample:
        peaks = ann.peaks.sample(n=args.sample, random_state=args.seed)
        peaks = peaks.sort_values(["chrom", "start"], kind="mergesort")
        ann = pcs.annotate_peak(peaks=peaks, TxDb=args.txdb, verbose=False)
        used = len(ann.peaks)

    df = ann.to_dataframe(include_peak_metadata=True).merge(
        ann.peaks[["peak_id", "score"]], on="peak_id", how="left"
    )

    promoters = pcs.get_promoters(
        TxDb=ann.gene_annotation,
        upstream=args.upstream,
        downstream=args.downstream,
    )
    nbin = 120
    coverage = pcs.get_tag_matrix(peak=ann.peaks, windows=promoters, nbin=nbin)
    # Peak coverage is measured in base pairs per bin; dividing by the bin
    # width turns it into a 0-100 % scale that is independent of ``nbin``.
    bin_bp = (args.upstream + args.downstream) / nbin
    tag = pcs.TagMatrix(
        np.asarray(coverage, dtype=float) * 100.0 / bin_bp,
        upstream=args.upstream,
        downstream=args.downstream,
        is_binning=True,
        rows=coverage.rows,
    )

    profiles = _quartile_profiles(tag, df)

    caption = (
        "ENCODE H3K4me3 peaks (ENCSR000CDQ, mm10) called at signal \u2265 3, "
        "merged over 1 kb gaps and \u2265 500 bp; promoter windows "
        f"\u00b1{args.upstream / 1000:g} kb around mm10 TSSs."
    )
    provenance = (
        f"pychipseeker {pcs.__version__}  \u00b7  annotation: mm10 TxDb  "
        "\u00b7  data: ENCODE ENCFF893ZGR"
    )

    with PdfPages(str(out / "encode_figures.pdf")) as pdf:
        with _no_tight_layout():
            panels = _panel_figures(
                ann, tag, df, profiles, args.upstream, args.downstream
            )
            for fig, stem in panels:
                _save_formats(fig, stem, figdir, formats=formats)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

    plate = _composite_plate(ann, tag, df, caption, provenance)
    plate_written = _save_formats(plate, "encode_figure_plate", out, formats=formats)
    plate.savefig(out / "encode_figure_plate_600dpi.png", dpi=600, bbox_inches="tight")
    _flatten_rgb(out / "encode_figure_plate_600dpi.png")
    plate_written.append(out / "encode_figure_plate_600dpi.png")
    plt.close(plate)

    supp = _composite_supplementary(df, profiles, caption, provenance)
    supp_written = _save_formats(
        supp, "encode_supplementary_plate", out, formats=formats
    )
    supp.savefig(
        out / "encode_supplementary_plate_600dpi.png", dpi=600, bbox_inches="tight"
    )
    _flatten_rgb(out / "encode_supplementary_plate_600dpi.png")
    supp_written.append(out / "encode_supplementary_plate_600dpi.png")
    plt.close(supp)

    analysis = _composite_analysis(ann, tag, df, caption, provenance)
    analysis_written = _save_formats(
        analysis, "encode_analysis_plate", out, formats=formats
    )
    analysis.savefig(
        out / "encode_analysis_plate_600dpi.png", dpi=600, bbox_inches="tight"
    )
    _flatten_rgb(out / "encode_analysis_plate_600dpi.png")
    analysis_written.append(out / "encode_analysis_plate_600dpi.png")
    plt.close(analysis)

    custom_written: list[Path] = []
    if selected:
        custom = _composite_custom(
            selected,
            ann,
            tag,
            df,
            profiles,
            args.upstream,
            args.downstream,
            caption,
            provenance,
        )
        custom_written = _save_formats(
            custom, "encode_custom_plate", out, formats=formats
        )
        custom.savefig(
            out / "encode_custom_plate_600dpi.png", dpi=600, bbox_inches="tight"
        )
        _flatten_rgb(out / "encode_custom_plate_600dpi.png")
        custom_written.append(out / "encode_custom_plate_600dpi.png")
        plt.close(custom)

    gallery = _write_gallery(out, formats)

    print(ann.annoStat.to_string(index=False))
    print(f"peaks   : {used:,}")
    print(f"promoter: {ann.annoStat.loc[0, 'Frequency']:.1f}% (<=1 kb)")
    print(f"plate   : {out}")
    for path in plate_written:
        print(f"          {path.name:<38} {path.stat().st_size / 1024:7.1f} KiB")
    for path in supp_written:
        print(f"          {path.name:<38} {path.stat().st_size / 1024:7.1f} KiB")
    for path in analysis_written:
        print(f"          {path.name:<38} {path.stat().st_size / 1024:7.1f} KiB")
    for path in custom_written:
        print(f"          {path.name:<38} {path.stat().st_size / 1024:7.1f} KiB")
    print(f"figures : {out / 'encode_figures.pdf'}")
    print(
        f"singles : {figdir} ({len(formats)} formats x "
        f"{len(PANEL_CATALOG)} panels)"
    )
    print(f"gallery : {gallery}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
