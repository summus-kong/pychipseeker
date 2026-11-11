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
- ``encode_figures.pdf`` -- six pages, one figure each
- ``figures/<panel>.*`` -- individual figures in PNG/PDF/SVG/EPS/TIFF

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


def _panel_figures(ann, tag, df) -> Iterator[tuple[plt.Figure, str]]:
    used = len(ann.peaks)
    yield _figure_pie(ann), "A_annotation_composition"
    yield _figure_width(df), "B_peak_width"
    yield _figure_distance(df), "C_distance_to_tss"
    yield _figure_profile(tag), "D_tag_profile"
    yield _figure_heatmap(tag), "E_tag_heatmap"
    yield _figure_table(ann, used), "F_annotation_stats"


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
    args = parser.parse_args()

    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    unknown = sorted(set(formats) - set(FORMATS))
    if unknown:
        parser.error(f"unsupported format(s): {', '.join(unknown)}")
    if not formats:
        parser.error("at least one export format is required")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    for fmt in FORMATS:
        for stale in figdir.glob(f"*.{fmt}"):
            stale.unlink()

    ann = pcs.annotate_peak(peaks=args.peaks, TxDb=args.txdb, verbose=False)
    used = len(ann.peaks)
    if args.sample and used > args.sample:
        peaks = ann.peaks.sample(n=args.sample, random_state=args.seed)
        peaks = peaks.sort_values(["chrom", "start"], kind="mergesort")
        ann = pcs.annotate_peak(peaks=peaks, TxDb=args.txdb, verbose=False)
        used = len(ann.peaks)

    df = ann.to_dataframe(include_peak_metadata=True)

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
            for fig, stem in _panel_figures(ann, tag, df):
                _save_formats(fig, stem, figdir, formats=formats)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

    plate = _composite_plate(ann, tag, df, caption, provenance)
    plate_written = _save_formats(plate, "encode_figure_plate", out, formats=formats)
    plate.savefig(out / "encode_figure_plate_600dpi.png", dpi=600, bbox_inches="tight")
    _flatten_rgb(out / "encode_figure_plate_600dpi.png")
    plate_written.append(out / "encode_figure_plate_600dpi.png")
    plt.close(plate)

    print(ann.annoStat.to_string(index=False))
    print(f"peaks   : {used:,}")
    print(f"promoter: {ann.annoStat.loc[0, 'Frequency']:.1f}% (<=1 kb)")
    print(f"plate   : {out}")
    for path in plate_written:
        print(f"          {path.name:<38} {path.stat().st_size / 1024:7.1f} KiB")
    print(f"figures : {out / 'encode_figures.pdf'}")
    print(f"singles : {figdir} ({len(formats)} formats x 6 panels)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
