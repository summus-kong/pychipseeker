"""ChIPseeker-style plots implemented with matplotlib.

The functions intentionally mirror the *semantics* of ChIPseeker's plotting
helpers: annotation statistics are the detailed ``annoStat`` categories
(``1st Exon``, ``Other Exon``, promoter bins, ...), ``plot_anno_bar`` uses a
horizontal stacked percentage bar, and ``plot_dist_to_tss`` uses a diverging
5'->3' distance bar rather than a histogram.

Matplotlib cannot reproduce every ggplot2 detail, but labels, categories,
percentages and axis directions follow ChIPseeker.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from .annotation import get_genomic_anno_stat

# ChIPseeker's getCols() palette (Set3-like).
_CHIPSEEKER_COLORS = [
    "#8dd3c7",
    "#ffffb3",
    "#bebada",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
    "#bc80bd",
    "#ccebc5",
    "#ffed6f",
]

_DIST_COLORS = ["#9ecae1", "#3182bd", "#C7A76C", "#86B875", "#39BEB1", "#CD99D8"]


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _as_df(annotation) -> pd.DataFrame:
    if hasattr(annotation, "to_dataframe"):
        df = annotation.to_dataframe(include_peak_metadata=True)
    elif isinstance(annotation, pd.DataFrame):
        df = annotation.copy()
    else:
        raise TypeError("annotation must be a PeakAnnotation or a pandas.DataFrame")

    # Accept ChIPseeker-style column names as input too.
    rename = {
        "distanceToTSS": "distance_to_tss",
        "peakChr": "chrom",
        "peakStart": "start",
        "peakEnd": "end",
        "geneChr": "gene_chrom",
        "geneStart": "gene_start",
        "geneEnd": "gene_end",
        "geneStrand": "gene_strand",
        "geneLength": "gene_length",
        "geneId": "gene_id",
        "transcriptId": "transcript_id",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    return df


def _stat_from_annotation(
    annotation,
    *,
    ignore_promoter_subcategory: bool = False,
    ignore_1st_exon: bool = False,
    ignore_1st_intron: bool = False,
    ignore_downstream: bool = False,
    downstream_distance: Optional[int] = None,
) -> pd.DataFrame:
    if hasattr(annotation, "anno_stat"):
        return annotation.anno_stat(
            ignore_promoter_subcategory=ignore_promoter_subcategory,
            ignore_1st_exon=ignore_1st_exon,
            ignore_1st_intron=ignore_1st_intron,
            ignore_downstream=ignore_downstream,
        )

    if isinstance(annotation, pd.DataFrame):
        df = annotation.copy()
        if {"Feature", "Frequency"}.issubset(df.columns):
            return df
        if "annotation" in df.columns:
            if downstream_distance is None:
                downstream_distance = 300
            return get_genomic_anno_stat(
                df["annotation"],
                downstream_distance=downstream_distance,
                ignore_promoter_subcategory=ignore_promoter_subcategory,
                ignore_1st_exon=ignore_1st_exon,
                ignore_1st_intron=ignore_1st_intron,
                ignore_downstream=ignore_downstream,
            )
    raise TypeError("annotation must provide annotation labels or Feature/Frequency")


def _colors(n: int) -> list[str]:
    if n <= 0:
        return []
    if n <= len(_CHIPSEEKER_COLORS):
        return _CHIPSEEKER_COLORS[:n]
    # Linear interpolation, analogous to colorRampPalette().
    base = np.array(
        [[int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)] for c in _CHIPSEEKER_COLORS],
        dtype=float,
    )
    xp = np.linspace(0, 1, len(base))
    x = np.linspace(0, 1, n)
    rgb = np.column_stack([np.interp(x, xp, base[:, i]) for i in range(3)]).astype(int)
    return ["#%02x%02x%02x" % tuple(v) for v in rgb]


def _style_axes(ax) -> None:
    """Clean up spines/ticks for the ChIPseeker matplotlib look."""

    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=0.8, labelsize=9, labelcolor="#4d4d4d", color="#bfbfbf")
    ax.xaxis.label.set_color("#333333")
    ax.yaxis.label.set_color("#333333")
    for side in ("left", "bottom"):
        if ax.spines.get(side) is not None:
            ax.spines[side].set_color("#cccccc")


def _title(ax, text: str) -> None:
    """Set the standardized figure title."""

    ax.set_title(
        text,
        fontsize=11.5,
        fontweight="semibold",
        color="#1a1a1a",
        pad=12.0,
    )


def _soften(fig) -> None:
    """Apply the shared light-weight typography to a whole figure."""

    if fig._suptitle is not None:
        fig._suptitle.set_fontsize(12.0)
        fig._suptitle.set_fontweight("semibold")
        fig._suptitle.set_color("#1a1a1a")


def plot_anno_pie(
    annotation,
    ax=None,
    title: str = "",
    colors: Optional[Sequence[str]] = None,
    legend: bool = True,
    ndigit: int = 2,
    figsize: tuple[float, float] = (7.0, 4.5),
    **kwargs,
):
    """Pie chart of detailed ChIPseeker annotation categories."""

    plt = _plt()
    stat_names = (
        "ignore_promoter_subcategory",
        "ignore_1st_exon",
        "ignore_1st_intron",
        "ignore_downstream",
    )
    stat_kwargs = {name: kwargs.pop(name) for name in stat_names if name in kwargs}
    stat = _stat_from_annotation(annotation, **stat_kwargs)
    if stat.empty:
        raise ValueError("no annotations to plot")

    values = stat["Frequency"].to_numpy(dtype=float)
    labels = [
        f"{feature} ({round(float(freq), ndigit)}%)"
        for feature, freq in zip(stat["Feature"], values)
    ]
    if colors is None:
        colors = _colors(len(stat))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    autopct = kwargs.pop("autopct", lambda pct: f"{pct:.{ndigit}f}%" if pct >= 3.0 else "")
    pctdistance = kwargs.pop("pctdistance", 0.80)
    wedges, texts, autotexts = ax.pie(
        values,
        labels=None,
        colors=colors,
        startangle=90,
        counterclock=False,
        autopct=autopct,
        pctdistance=pctdistance,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=2.0),
        **kwargs,
    )
    # Percentage labels: white when drawn on top of a dark wedge.
    for autotext, wedge in zip(autotexts, wedges):
        autotext.set_fontsize(9.0)
        autotext.set_fontweight("semibold")
        rgb = wedge.get_facecolor()
        luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
        autotext.set_color("white" if luminance < 0.6 else "#333333")
    # Donut-hole annotation with the number of annotated peaks.
    try:
        n_rows = len(_as_df(annotation))
    except Exception:  # pragma: no cover - defensive; annotation without a table
        n_rows = None
    if n_rows is not None:
        ax.text(
            0.0,
            0.12,
            f"{n_rows}",
            ha="center",
            va="center",
            fontsize=17,
            fontweight="bold",
            color="#333333",
        )
        ax.text(
            0.0,
            -0.16,
            "peaks\nannotated",
            ha="center",
            va="center",
            fontsize=8.5,
            color="#8c8c8c",
        )
    if legend:
        ax.legend(
            wedges,
            labels,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            frameon=False,
            fontsize=9.5,
        )
    _title(ax, title)
    ax.axis("equal")
    fig.tight_layout()
    return fig, ax


def _grouped_stat(
    annotation,
    groupby: str,
    **kwargs,
) -> pd.DataFrame:
    df = _as_df(annotation)
    if groupby not in df.columns:
        raise KeyError(f"groupby column not found: {groupby}")
    if "annotation" not in df.columns:
        raise KeyError("annotation table has no 'annotation' column")

    all_stat = _stat_from_annotation(annotation, **kwargs)
    features = all_stat["Feature"].astype(str).tolist()
    downstream_distance = kwargs.get("downstream_distance")
    if downstream_distance is None:
        downstream_distance = getattr(annotation, "downstream_distance", 300)
    rows: dict[str, pd.Series] = {}
    for group, sub in df.groupby(groupby, dropna=False):
        stat = get_genomic_anno_stat(
            sub["annotation"],
            downstream_distance=downstream_distance,
            ignore_promoter_subcategory=kwargs.get("ignore_promoter_subcategory", False),
            ignore_1st_exon=kwargs.get("ignore_1st_exon", False),
            ignore_1st_intron=kwargs.get("ignore_1st_intron", False),
            ignore_downstream=kwargs.get("ignore_downstream", False),
        )
        rows[str(group)] = (
            stat.set_index("Feature")["Frequency"].reindex(features, fill_value=0.0)
        )
    return pd.DataFrame.from_dict(rows, orient="index")


def plot_anno_bar(
    annotation,
    groupby: Optional[str] = None,
    ax=None,
    title: str = "Feature Distribution",
    xlab: str = "",
    ylab: str = "Percentage(%)",
    figsize: tuple[float, float] = (7.0, 4.5),
    colors: Optional[Sequence[str]] = None,
    **kwargs,
):
    """Horizontal stacked bar chart of detailed annotation categories.

    With ``groupby=None`` a single stacked bar is drawn, like
    ``plotAnnoBar(csAnno)``.  With a column name, one bar per group is drawn,
    like ``plotAnnoBar(list_of_csAnno)``.
    """

    plt = _plt()
    if groupby is None:
        stat = _stat_from_annotation(annotation, **kwargs)
        if stat.empty:
            raise ValueError("no annotations to plot")
        features = stat["Feature"].astype(str).tolist()
        data = pd.DataFrame(
            [stat["Frequency"].to_numpy(dtype=float)],
            index=["All"],
            columns=features,
        )
        single = True
    else:
        data = _grouped_stat(annotation, groupby, **kwargs)
        features = data.columns.astype(str).tolist()
        single = False

    if colors is None:
        colors = _colors(len(features))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    left = np.zeros(len(data), dtype=float)
    y = np.arange(len(data))
    for color, feature in zip(colors, features):
        vals = data[feature].to_numpy(dtype=float)
        ax.barh(
            y,
            vals,
            left=left,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            label=feature,
        )
        left += vals

    # Percentage labels centered in each segment wide enough to hold them.
    for yy in y:
        acc = 0.0
        for color, feature in zip(colors, features):
            val = float(data[feature].iloc[yy])
            if val >= 4.0:
                rgb = tuple(int(color[j : j + 2], 16) / 255 for j in (1, 3, 5))
                luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
                ax.text(
                    acc + val / 2.0,
                    yy,
                    f"{val:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=8.0,
                    color="white" if luminance < 0.62 else "#333333",
                )
            acc += val

    ax.set_xlabel(ylab)
    ax.set_ylabel(xlab)
    ax.set_yticks(y)
    if single:
        ax.set_yticklabels([])
    else:
        ax.set_yticklabels(data.index.tolist())
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        fontsize=9.5,
        reverse=True,
    )
    ax.set_xlim(0, max(100.0, float(np.nanmax(left)) if len(left) else 100.0))
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["left"]].set_visible(False)
    _style_axes(ax)
    _title(ax, title)
    fig.tight_layout()
    return fig, ax


def _distance_breaks(
    distance_breaks: Optional[Sequence[float]],
) -> tuple[np.ndarray, list[str]]:
    if distance_breaks is None:
        breaks = [0, 1000, 3000, 5000, 10000, 100000]
    else:
        breaks = [float(x) for x in distance_breaks]
    breaks = sorted(set(breaks))
    if 0 not in breaks:
        breaks = [0] + breaks
    if not np.isinf(breaks[-1]):
        breaks = breaks + [float("inf")]
    arr = np.asarray(breaks, dtype=float)

    def fmt(value: float) -> str:
        if np.isinf(value):
            return "Inf"
        if value == 0:
            return "0"
        if value >= 1000 and value % 1000 == 0:
            return f"{int(value // 1000)}kb"
        return f"{value / 1000:g}kb"

    labels = [f"{fmt(arr[i])}-{fmt(arr[i + 1])}" for i in range(len(arr) - 2)]
    labels.append(f">{fmt(arr[-2])}")
    return arr, labels


def plot_dist_to_tss(
    annotation,
    groupby: Optional[str] = None,
    ax=None,
    distance_breaks: Optional[Sequence[float]] = None,
    distance_column: str = "distance_to_tss",
    title: str = "Distribution of transcription factor-binding loci relative to TSS",
    xlab: str = "",
    ylab: str = "Binding sites (%) (5'->3')",
    figsize: tuple[float, float] = (7.5, 4.5),
    colors: Optional[Sequence[str]] = None,
    **kwargs,
):
    """Diverging distance-to-TSS bar chart, following ``plotDistToTSS``.

    Positive distances are drawn to the right (downstream) and negative
    distances to the left (upstream).  Zero-distance peaks are split equally
    between the two sides, as in ChIPseeker.
    """

    plt = _plt()
    df = _as_df(annotation)
    if distance_column not in df.columns:
        for candidate in ("distanceToTSS", "distance_to_tss"):
            if candidate in df.columns:
                distance_column = candidate
                break
        else:
            raise KeyError("annotation table has no distance-to-TSS column")
    if groupby is not None and groupby not in df.columns:
        raise KeyError(f"groupby column not found: {groupby}")

    work = pd.DataFrame(
        {
            "distance": pd.to_numeric(df[distance_column], errors="coerce"),
            "group": df[groupby].astype(str) if groupby is not None else "All",
        }
    ).dropna(subset=["distance"])
    if work.empty:
        raise ValueError("no finite distance values to plot")

    breaks, labels = _distance_breaks(distance_breaks)
    work["Feature"] = pd.cut(
        work["distance"].abs(),
        bins=breaks,
        labels=labels,
        include_lowest=True,
        right=True,
    )
    work["sign"] = np.sign(work["distance"]).astype(int)

    # Zero-distance peaks contribute half to each side.
    zero = work[work["sign"] == 0].copy()
    if not zero.empty:
        zero["Frequency"] = 0.5
        zero["sign"] = -1
        work = pd.concat([work, zero], ignore_index=True)
        zero["sign"] = 1
        work = pd.concat([work, zero], ignore_index=True)
    work = work[work["sign"] != 0]

    counts = (
        work.groupby(["group", "Feature", "sign"], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    totals = counts.groupby("group")["count"].transform("sum")
    counts["Frequency"] = counts["count"] / totals * 100

    groups = list(dict.fromkeys(work["group"].astype(str)))
    feature_order = labels
    if colors is None:
        colors = _DIST_COLORS * (len(feature_order) // len(_DIST_COLORS) + 1)
    color_map = dict(zip(feature_order, colors[: len(feature_order)]))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    y = np.arange(len(groups))
    left_neg = {g: 0.0 for g in groups}
    left_pos = {g: 0.0 for g in groups}
    for feature in feature_order:
        for idx, group in enumerate(groups):
            sub = counts[
                (counts["group"].astype(str) == group)
                & (counts["Feature"].astype(str) == str(feature))
            ]
            neg = float(sub.loc[sub["sign"] == -1, "Frequency"].sum())
            pos = float(sub.loc[sub["sign"] == 1, "Frequency"].sum())
            if neg:
                ax.barh(
                    idx,
                    -neg,
                    left=-left_neg[group],
                    color=color_map[feature],
                    edgecolor="white",
                    linewidth=0.3,
                )
                left_neg[group] += neg
            if pos:
                ax.barh(
                    idx,
                    pos,
                    left=left_pos[group],
                    color=color_map[feature],
                    edgecolor="white",
                    linewidth=0.3,
                )
                left_pos[group] += pos

    ax.axvline(0, color="#4d4d4d", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(groups)
    ax.set_xlabel(ylab)
    ax.set_ylabel(xlab)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(length=3, width=0.8, labelsize=9, labelcolor="#4d4d4d", color="#bfbfbf")
    # 5'->3' orientation hints in figure coordinates.
    ax.text(
        0.02,
        1.015,
        "Upstream (5') \u2190",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#8c8c8c",
    )
    ax.text(
        0.98,
        1.015,
        "\u2192 Downstream (3')",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="#8c8c8c",
    )
    _title(ax, title)

    # Symmetric-ish limits and absolute tick labels with TSS at zero.
    max_neg = max(left_neg.values(), default=0.0)
    max_pos = max(left_pos.values(), default=0.0)
    lim = max(np.ceil(max(max_neg, max_pos) / 10.0) * 10.0, 10.0)
    ax.set_xlim(-lim, lim)
    ticks = np.arange(-lim, lim + 1, 10)
    ax.set_xticks(ticks)
    ax.set_xticklabels(["TSS" if t == 0 else f"{abs(int(t))}" for t in ticks])

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=color_map[f]) for f in feature_order
    ]
    ax.legend(
        handles,
        feature_order,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        reverse=True,
    )
    fig.tight_layout()
    return fig, ax


def _coerce_signal_frame(signal: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    """Normalize a bedGraph-like signal table to chrom/start/end/value.

    This is deliberately dependency-light: a four-column tab-delimited file (or
    a ``.bedGraph``/``.bdg`` file, optionally gzipped) can be read by pandas
    without pyBigWig.  ``track``/``browser``/comment lines are ignored.
    """

    columns = ["chrom", "start", "end", "value"]
    if isinstance(signal, (str, Path)):
        df = pd.read_csv(
            signal,
            sep="\t",
            header=None,
            comment="#",
            names=columns,
            usecols=[0, 1, 2, 3],
            dtype={0: str},
            low_memory=False,
        )
    elif isinstance(signal, pd.DataFrame):
        df = signal.copy()
        if not set(columns).issubset(df.columns):
            if len(df.columns) < 4:
                raise ValueError(
                    "signal DataFrame must contain chrom,start,end,value columns"
                )
            first = list(df.columns[:4])
            df = df.rename(
                columns=dict(zip(first, columns))
            )
    else:
        raise TypeError("signal must be a path or a pandas.DataFrame")

    # Drop header/track/browser rows and any malformed records.  bedGraph is
    # 0-based, half-open, so no coordinate conversion is needed here.
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "value"]).copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    df["value"] = df["value"].astype(float)
    df = df.loc[df["end"] > df["start"]]
    if df.empty:
        raise ValueError("signal table contains no valid intervals")
    return df[columns]


def _mean_in_interval(
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    start: float,
    end: float,
) -> float:
    """Length-weighted mean of signal values overlapping ``[start, end)``."""

    # ``starts``/``ends`` are sorted by start and non-overlapping for bedGraph.
    lo = int(np.searchsorted(ends, start, side="right"))
    hi = int(np.searchsorted(starts, end, side="left"))
    if hi <= lo:
        return float("nan")
    overlap_start = np.maximum(starts[lo:hi], start)
    overlap_end = np.minimum(ends[lo:hi], end)
    lengths = overlap_end - overlap_start
    valid = lengths > 0
    if not valid.any():
        return float("nan")
    total = float(np.sum(lengths[valid]))
    if total <= 0:
        return float("nan")
    return float(np.sum(values[lo:hi][valid] * lengths[valid]) / total)


def extract_tss_matrix_from_signal(
    signal: Union[str, Path, pd.DataFrame],
    genes: pd.DataFrame,
    window: int = 2000,
    n_bins: int = 100,
    min_value: Optional[float] = 0.0,
    chrom_sizes: Optional[dict[str, int]] = None,
) -> np.ndarray:
    """Extract a TSS-centered matrix from a bedGraph-like signal table.

    This is the no-pyBigWig alternative to :func:`extract_tss_matrix`.  The
    input may be a bedGraph/bedGraph.gz path or a DataFrame with
    ``chrom,start,end,value`` columns.  Coordinates are 0-based, half-open.
    Bins with no signal are left as ``NaN`` (mirroring pyBigWig's mean stats).

    ``genes`` needs ``chrom``, ``tss`` and ``strand`` columns, such as the
    transcript/gene table exposed by :class:`~pychipseeker.GeneAnnotation`.
    Windows that fall outside a chromosome are skipped.  Chromosome lengths are
    inferred from the signal unless ``chrom_sizes`` is supplied.
    """

    needed = {"chrom", "tss", "strand"}
    missing = needed - set(genes.columns)
    if missing:
        raise ValueError(f"genes table is missing columns: {sorted(missing)}")

    window = int(window)
    n_bins = int(n_bins)
    if window <= 0 or n_bins <= 1:
        raise ValueError("window must be positive and n_bins must be > 1")

    sig = _coerce_signal_frame(signal)
    if chrom_sizes is None:
        chrom_sizes = (
            sig.groupby("chrom", sort=False)["end"].max().astype("int64").to_dict()
        )
    else:
        chrom_sizes = {str(k): int(v) for k, v in chrom_sizes.items()}

    by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, group in sig.groupby("chrom", sort=False):
        group = group.sort_values("start", kind="mergesort")
        by_chrom[str(chrom)] = (
            group["start"].to_numpy(dtype=np.int64),
            group["end"].to_numpy(dtype=np.int64),
            group["value"].to_numpy(dtype=float),
        )

    rows: list[np.ndarray] = []
    for row in genes.itertuples(index=False):
        chrom = str(row.chrom)
        if chrom not in by_chrom:
            continue
        tss = int(row.tss)
        chrom_len = chrom_sizes.get(chrom)
        if chrom_len is None:
            continue
        start = tss - window
        end = tss + window
        if start < 0 or end > chrom_len:
            continue

        starts, ends, values = by_chrom[chrom]
        edges = np.linspace(start, end, n_bins + 1)
        arr = np.asarray(
            [
                _mean_in_interval(starts, ends, values, edges[i], edges[i + 1])
                for i in range(n_bins)
            ],
            dtype=float,
        )
        if np.all(np.isnan(arr)):
            continue
        if min_value is not None:
            arr = np.where(arr < float(min_value), 0.0, arr)
        if str(row.strand) == "-":
            arr = arr[::-1]
        rows.append(arr)

    if not rows:
        raise ValueError("no features could be extracted from the signal track")
    return np.vstack(rows)


def extract_tss_matrix(
    bigwig: Union[str, Path],
    genes: pd.DataFrame,
    window: int = 2000,
    n_bins: int = 100,
    min_value: float = 0.0,
) -> np.ndarray:
    """Extract a features x bins TSS-centered signal matrix from a bigWig file.

    ``pyBigWig`` is a declared dependency of pychipseeker.  A bedGraph table or
    precomputed matrix can be used instead via
    :func:`extract_tss_matrix_from_signal` or ``matrix=...`` in
    ``plot_tss_profile`` / ``plot_tss_heatmap``.
    """

    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - dependency installation issue
        raise ImportError(
            "pyBigWig is required to read .bw/.bigWig files directly. "
            "Reinstall pychipseeker (pyBigWig is a declared dependency), or "
            "use extract_tss_matrix_from_signal() with a bedGraph table or "
            "pass matrix=... to plot_tss_profile/plot_tss_heatmap."
        ) from exc

    needed = {"chrom", "tss", "strand"}
    missing = needed - set(genes.columns)
    if missing:
        raise ValueError(f"genes table is missing columns: {sorted(missing)}")

    window = int(window)
    n_bins = int(n_bins)
    if window <= 0 or n_bins <= 1:
        raise ValueError("window must be positive and n_bins must be > 1")

    bw = pyBigWig.open(str(bigwig))
    if bw is None:
        raise ValueError(f"could not open bigWig: {bigwig}")

    chrom_sizes = bw.chroms()
    rows: list[np.ndarray] = []

    try:
        for row in genes.itertuples(index=False):
            chrom = str(row.chrom)
            tss = int(row.tss)
            strand = str(row.strand)
            chrom_len = chrom_sizes.get(chrom)
            if chrom_len is None:
                continue
            start = tss - window
            end = tss + window
            if start < 0 or end > chrom_len:
                continue

            vals = bw.stats(chrom, start, end, nBins=n_bins, type="mean")
            if vals is None:
                continue
            arr = np.asarray(
                [np.nan if v is None else float(v) for v in vals],
                dtype=float,
            )
            if np.all(np.isnan(arr)):
                continue
            if min_value is not None:
                arr = np.where(arr < min_value, 0.0, arr)
            if strand == "-":
                arr = arr[::-1]
            rows.append(arr)
    finally:
        bw.close()

    if not rows:
        raise ValueError("no features could be extracted from the bigWig track")
    return np.vstack(rows)


def plot_tss_profile(
    matrix: np.ndarray,
    window: int = 2000,
    ax=None,
    title: str = "TSS profile",
    color: str = "#E45756",
    figsize: tuple[float, float] = (6.0, 4.0),
    show_sem: bool = True,
    **kwargs,
):
    """Plot the mean TSS-centered signal profile."""

    plt = _plt()
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("matrix must be 2-dimensional (features x bins)")
    if matrix.shape[0] == 0:
        raise ValueError("matrix has no rows")

    n_bins = matrix.shape[1]
    x = np.linspace(-int(window), int(window), n_bins)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(matrix, axis=0)
        sem = np.nanstd(matrix, axis=0) / np.sqrt(
            np.maximum(np.sum(~np.isnan(matrix), axis=0), 1)
        )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    ax.plot(x, mean, color=color, linewidth=2.0, label="mean", **kwargs)
    # Soft area under the mean curve for a cleaner visual mass.
    ax.fill_between(x, 0.0, mean, color=color, alpha=0.08, linewidth=0, zorder=0.5)
    # Highlight the profile maximum.
    peak_idx = int(np.nanargmax(mean))
    ax.plot(
        x[peak_idx],
        mean[peak_idx],
        "o",
        markersize=4.5,
        color=color,
        markeredgecolor="white",
        markeredgewidth=0.8,
        zorder=6,
    )
    ax.annotate(
        f"{float(mean[peak_idx]):.4g}",
        (float(x[peak_idx]), float(mean[peak_idx])),
        textcoords="offset points",
        xytext=(0, 7),
        ha="center",
        fontsize=8,
        color="#333333",
    )
    if show_sem:
        ax.fill_between(
            x,
            mean - sem,
            mean + sem,
            color=color,
            alpha=0.22,
            linewidth=0,
            label="±SEM",
        )
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.axhline(0, color="#bbbbbb", linewidth=0.8, linestyle="--", zorder=0)
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("Mean signal")
    if show_sem:
        ax.legend(frameon=False, fontsize=9)
    ax.margins(x=0.01)
    _style_axes(ax)
    _title(ax, title)
    fig.tight_layout()
    return fig, ax


def plot_tss_heatmap(
    matrix: np.ndarray,
    window: int = 2000,
    ax=None,
    title: str = "TSS heatmap",
    cmap: str = "viridis",
    figsize: tuple[float, float] = (6.0, 5.0),
    sort_rows: bool = False,
    **kwargs,
):
    """Plot a heatmap of TSS-centered signal.

    With ``sort_rows=True`` rows are ordered by decreasing maximum signal so the
    strongest features sit on top and the map fades downwards.
    """

    plt = _plt()
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("matrix must be 2-dimensional (features x bins)")
    if matrix.shape[0] == 0:
        raise ValueError("matrix has no rows")

    if sort_rows:
        finite = np.nan_to_num(matrix, nan=-np.inf)
        order = np.argsort(-np.nanmax(finite, axis=1), kind="stable")
        matrix = matrix[order]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    masked = np.ma.masked_invalid(matrix)
    from matplotlib import colormaps

    cmap_obj = colormaps[(kwargs.pop("cmap", None) or cmap)].with_extremes(bad="#f5f5f5")
    im = ax.imshow(
        masked,
        aspect="auto",
        cmap=cmap_obj,
        interpolation="nearest",
        extent=(-int(window), int(window), matrix.shape[0], 0),
        **kwargs,
    )
    ax.axvline(0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("Features")
    cbar = fig.colorbar(im, ax=ax, label="Signal", shrink=0.85, pad=0.02)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=8, length=2.5, color="#bfbfbf")
    _style_axes(ax)
    _title(ax, title)
    fig.tight_layout()
    return fig, ax


# ChIPseeker-style aliases.  ``plotProfile``/``plotHeatmap`` in R operate on a
# tag matrix as well, so a direct alias is reasonable here.
plotAnnoPie = plot_anno_pie
plotAnnoBar = plot_anno_bar
plotDistToTSS = plot_dist_to_tss
plotProfile = plot_tss_profile
plotHeatmap = plot_tss_heatmap
