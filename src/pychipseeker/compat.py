"""Compatibility helpers for common ChIPseeker functions.

The implementation is deliberately dependency-light.  It accepts the same
kinds of inputs as the Python core (BED-like paths and pandas DataFrames) and
returns ``(fig, ax)`` pairs from plotting functions, matching the rest of
``pychipseeker``.
"""

from __future__ import annotations

from collections import Counter
from itertools import product
import math
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd

from .annotation import (
    GeneAnnotation,
    PeakAnnotation,
    _normalise_peaks,
    annotate_peak,
    get_genomic_anno_stat,
)
from .io import read_peaks
from .tagmatrix import get_bio_region


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _as_annotation_df(annotation) -> pd.DataFrame:
    if isinstance(annotation, PeakAnnotation):
        return annotation.to_chipseeker()
    if isinstance(annotation, pd.DataFrame):
        return annotation.copy()
    raise TypeError("expected a PeakAnnotation or pandas.DataFrame")


def _rename_to_granges(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename = {
        "peakChr": "seqnames",
        "peakStart": "start",
        "peakEnd": "end",
        "peakStrand": "strand",
        "chrom": "seqnames",
        "chromStart": "start",
        "chromEnd": "end",
        "chr": "seqnames",
    }
    out = out.rename(columns={k: v for k, v in rename.items() if k in out.columns})
    # Only rename chrom/start/end if doing so does not clobber existing
    # annotation coordinates.  peak* columns above have precedence.
    return out


def read_peak_file(
    peak: Union[str, Path, pd.DataFrame],
    as_type: str = "dataframe",
    **kwargs,
) -> pd.DataFrame:
    """Read a peak file; ``as_type='GRanges'`` returns a GRanges-like table."""

    if isinstance(peak, pd.DataFrame):
        df = peak.copy()
    else:
        df = read_peaks(peak, **kwargs)
    if str(as_type).lower() in {"granges", "grange", "gr"}:
        return _rename_to_granges(df)
    return df


def as_granges(annotation) -> pd.DataFrame:
    """Return a GRanges-like table for a ``PeakAnnotation`` or DataFrame."""

    if isinstance(annotation, PeakAnnotation):
        return annotation.as_granges()
    if isinstance(annotation, pd.DataFrame):
        return _rename_to_granges(annotation)
    raise TypeError("expected a PeakAnnotation or pandas.DataFrame")


def get_anno_stat(annotation) -> pd.DataFrame:
    """ChIPseeker-style ``getAnnoStat`` for annotation objects/tables."""

    if isinstance(annotation, PeakAnnotation):
        return annotation.anno_stat()
    if isinstance(annotation, pd.DataFrame):
        if {"Feature", "Frequency"}.issubset(annotation.columns):
            return annotation.copy()
        if "annotation" in annotation.columns:
            return get_genomic_anno_stat(annotation["annotation"])
        # Accept a ChIPseeker-style table with an "annotation" column.
        columns = {str(c).lower(): c for c in annotation.columns}
        if "annotation" in columns:
            return get_genomic_anno_stat(annotation[columns["annotation"]])
    raise TypeError("annotation must provide annotation labels or Feature/Frequency")


def detail_genomic_annotation(annotation) -> pd.DataFrame:
    """Return the detailed per-peak genomic annotation table."""

    if isinstance(annotation, PeakAnnotation):
        return annotation.detail_genomic_annotation.copy()
    if isinstance(annotation, pd.DataFrame):
        return annotation.copy()
    raise TypeError("expected a PeakAnnotation or pandas.DataFrame")


def drop_anno(annotation, distance_to_tss_cutoff: float = 10000):
    """Filter peaks by ``abs(distanceToTSS) < cutoff``.

    This mirrors ``ChIPseeker::dropAnno``: annotation columns are *not*
    removed; peaks far from a TSS are dropped instead.
    """

    cutoff = float(distance_to_tss_cutoff)
    if isinstance(annotation, PeakAnnotation):
        return annotation.drop_anno(cutoff)
    if isinstance(annotation, pd.DataFrame):
        df = annotation.copy()
        distance = None
        for name in ("distance_to_tss", "distanceToTSS"):
            if name in df.columns:
                distance = df[name]
                break
        if distance is None:
            raise ValueError("annotation has no distance_to_tss/distanceToTSS column")
        return df.loc[distance.abs() < cutoff].reset_index(drop=True)
    raise TypeError("expected a PeakAnnotation or pandas.DataFrame")


def annotate_db(
    annotation: PeakAnnotation,
    anno_db,
    columns: Optional[Sequence[str]] = None,
) -> PeakAnnotation:
    """Add OrgDb-like annotation columns to a ``PeakAnnotation``.

    ChIPseeker's ``annoDb`` is an R OrgDb package.  In Python the closest
    dependency-light equivalent is a DataFrame or ``{gene_id: {...}}`` mapping;
    the requested ``columns`` are added by ``gene_id``.
    """

    if not isinstance(annotation, PeakAnnotation):
        raise TypeError("annotation must be a PeakAnnotation")
    df = annotation.annotation.copy()
    if "gene_id" not in df.columns:
        raise ValueError("annotation has no gene_id column")
    key = "gene_id"

    if isinstance(anno_db, pd.DataFrame):
        table = anno_db.copy()
        if key not in table.columns:
            candidate = next(
                (c for c in ("geneId", "gene", "GID", "ENSEMBL") if c in table.columns),
                None,
            )
            if candidate is None:
                candidate = table.columns[0]
            table = table.rename(columns={candidate: key})
        selected = [c for c in (columns or table.columns) if c in table.columns and c != key]
        table = table[[key] + selected].drop_duplicates(key)
        merged = df.merge(table, on=key, how="left", suffixes=("", "_annoDb"))
        annotation.annotation = merged.reset_index(drop=True)
        return annotation

    if isinstance(anno_db, dict):
        records = []
        for gene in df[key].astype(str):
            value = anno_db.get(gene, {})
            if isinstance(value, dict):
                records.append(value)
            else:
                records.append({"SYMBOL": value})
        extra = pd.DataFrame(records)
        selected = [c for c in (columns or extra.columns) if c in extra.columns]
        if selected:
            extra = extra[selected]
        annotation.annotation = pd.concat(
            [df.reset_index(drop=True), extra.reset_index(drop=True)], axis=1
        )
        return annotation

    raise TypeError("annoDb must be a pandas.DataFrame or a mapping")


def annotatePeak(
    peak,
    tssRegion: Sequence[int] = (-3000, 3000),
    TxDb=None,
    level: str = "transcript",
    assignGenomicAnnotation: bool = True,
    genomicAnnotationPriority: Optional[Sequence[str]] = None,
    annoDb=None,
    addFlankGeneInfo: bool = False,
    flankDistance: int = 5000,
    sameStrand: bool = False,
    ignoreOverlap: bool = False,
    ignoreUpstream: bool = False,
    ignoreDownstream: bool = False,
    overlap: str = "TSS",
    verbose: bool = True,
    columns: Sequence[str] = ("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"),
    **kwargs,
) -> PeakAnnotation:
    """ChIPseeker-style ``annotatePeak`` wrapper with camelCase arguments."""

    if "tss_region" in kwargs:
        tssRegion = kwargs.pop("tss_region")
    if "assign_genomic_annotation" in kwargs:
        assignGenomicAnnotation = kwargs.pop("assign_genomic_annotation")
    if "genomic_annotation_priority" in kwargs:
        genomicAnnotationPriority = kwargs.pop("genomic_annotation_priority")
    if "add_flank_gene_info" in kwargs:
        addFlankGeneInfo = kwargs.pop("add_flank_gene_info")
    if "flank_distance" in kwargs:
        flankDistance = kwargs.pop("flank_distance")
    if "same_strand" in kwargs:
        sameStrand = kwargs.pop("same_strand")
    if "ignore_overlap" in kwargs:
        ignoreOverlap = kwargs.pop("ignore_overlap")
    if "ignore_upstream" in kwargs:
        ignoreUpstream = kwargs.pop("ignore_upstream")
    if "ignore_downstream" in kwargs:
        ignoreDownstream = kwargs.pop("ignore_downstream")

    source = TxDb
    if source is None:
        source = kwargs.pop("gtf", None)
    if source is None:
        source = kwargs.pop("annotation", None)
    if source is None:
        raise ValueError("provide TxDb=... or gtf=...")

    priority = (
        tuple(genomicAnnotationPriority)
        if genomicAnnotationPriority is not None
        else (
            "Promoter",
            "5UTR",
            "3UTR",
            "Exon",
            "Intron",
            "Downstream",
            "Intergenic",
        )
    )
    ann = annotate_peak(
        peaks=peak,
        gtf=source,
        tss_region=tssRegion,
        level=level,
        assign_genomic_annotation=assignGenomicAnnotation,
        genomic_annotation_priority=priority,
        add_flank_gene_info=addFlankGeneInfo,
        flank_distance=flankDistance,
        same_strand=sameStrand,
        ignore_overlap=ignoreOverlap,
        ignore_upstream=ignoreUpstream,
        ignore_downstream=ignoreDownstream,
        overlap=overlap,
        verbose=verbose,
        **kwargs,
    )
    if annoDb is not None:
        ann = annotate_db(ann, annoDb, columns=columns)
    return ann


def combine_cs_anno(*annotations) -> PeakAnnotation:
    """Combine several ``PeakAnnotation`` objects (ChIPseeker ``combine_csAnno``)."""
    if len(annotations) == 1 and isinstance(annotations[0], (list, tuple)):
        annotations = tuple(annotations[0])
    if not annotations:
        raise ValueError("at least one annotation is required")
    if not all(isinstance(a, PeakAnnotation) for a in annotations):
        raise TypeError("all inputs must be PeakAnnotation objects")

    first = annotations[0]
    for other in annotations[1:]:
        if tuple(other.tss_region) != tuple(first.tss_region):
            raise ValueError("all annotations must use the same tss_region")
        if other.level != first.level:
            raise ValueError("all annotations must use the same level")
        if bool(other.has_genomic_annotation) != bool(first.has_genomic_annotation):
            raise ValueError("all annotations must agree on has_genomic_annotation")

    peaks = pd.concat([a.peaks for a in annotations], ignore_index=True)
    annotation = pd.concat([a.annotation for a in annotations], ignore_index=True)
    detail = pd.concat(
        [a.detail_genomic_annotation for a in annotations], ignore_index=True
    )
    return PeakAnnotation(
        peaks=peaks,
        annotation=annotation,
        gene_annotation=first.gene_annotation,
        tss_region=first.tss_region,
        level=first.level,
        promoter_bins=first.promoter_bins,
        genomic_annotation_priority=first.genomic_annotation_priority,
        downstream_distance=first.downstream_distance,
        detail_genomic_annotation=detail,
        has_genomic_annotation=first.has_genomic_annotation,
        add_flank_gene_info=first.add_flank_gene_info,
        flank_distance=first.flank_distance,
    )


# ----------------------------------------------------------------------
# set overlap / venn helpers
# ----------------------------------------------------------------------
def _set_elements(value) -> set[str]:
    if isinstance(value, PeakAnnotation):
        if "peak_id" in value.annotation.columns:
            return set(value.annotation["peak_id"].astype(str))
        return set(value.annotation.index.astype(str))
    if isinstance(value, pd.DataFrame):
        df = _normalise_peaks(value)
        return {
            f"{r.chrom}:{int(r.start)}-{int(r.end)}"
            for r in df.itertuples(index=False)
        }
    if isinstance(value, dict):
        return set(map(str, value.keys()))
    return {str(v) for v in value}


def overlap(Sets) -> pd.DataFrame:
    """Calculate exact membership counts for a list/dict of sets.

    The result has one row per membership pattern and a final ``Weight``
    column, matching the shape used by ChIPseeker's ``overlap`` helper.
    """

    if isinstance(Sets, dict):
        names = [str(k) for k in Sets.keys()]
        values = list(Sets.values())
    else:
        values = list(Sets)
        names = [f"Set{i + 1}" for i in range(len(values))]
    if not values:
        raise ValueError("Sets is empty")
    sets = [_set_elements(v) for v in values]
    counts: Counter[tuple[int, ...]] = Counter()
    for element in set().union(*sets):
        membership = tuple(1 if element in s else 0 for s in sets)
        counts[membership] += 1

    rows = []
    for membership in product((0, 1), repeat=len(sets)):
        if not any(membership):
            continue
        row = {name: int(flag) for name, flag in zip(names, membership)}
        row["Weight"] = int(counts.get(membership, 0))
        rows.append(row)
    rows.sort(key=lambda r: (-sum(r[n] for n in names), [r[n] for n in names]))
    return pd.DataFrame(rows, columns=names + ["Weight"])


def vennplot(
    Sets,
    labels: Optional[Sequence[str]] = None,
    title: str = "",
    colors: Optional[Sequence[str]] = None,
    figsize: tuple[float, float] = (5.5, 5.0),
    ax=None,
):
    """Plot a Venn diagram for 1-3 sets, or an overlap bar chart for more."""

    plt = _plt()
    if isinstance(Sets, dict):
        names = [str(k) for k in Sets.keys()]
        values = list(Sets.values())
    else:
        values = list(Sets)
        names = list(labels) if labels is not None else [
            f"Set{i + 1}" for i in range(len(values))
        ]
    if not values:
        raise ValueError("Sets is empty")
    if labels is not None and len(labels) != len(values):
        raise ValueError("labels length must match Sets")
    table = overlap(dict(zip(names, values)))

    def weight(*flags: int) -> int:
        row = table
        for name, flag in zip(names, flags):
            row = row.loc[row[name] == int(flag)]
        return int(row["Weight"].iloc[0]) if len(row) else 0

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    palette = list(colors) if colors is not None else [
        "#8dd3c7",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
        "#fccde5",
    ]
    n = len(names)

    if n <= 3:
        from matplotlib.patches import Circle

        if n == 1:
            centers = [(0.0, 0.0)]
            radius = 0.8
        elif n == 2:
            centers = [(-0.42, 0.0), (0.42, 0.0)]
            radius = 0.68
        else:
            centers = [(-0.5, -0.22), (0.5, -0.22), (0.0, 0.5)]
            radius = 0.62
        for i, (x, y) in enumerate(centers):
            ax.add_patch(
                Circle(
                    (x, y),
                    radius,
                    facecolor=palette[i % len(palette)],
                    edgecolor="white",
                    linewidth=1.2,
                    alpha=0.55,
                )
            )
        # Label positions are deliberately approximate; the numbers are exact.
        if n == 1:
            ax.text(0, 0, str(weight(1)), ha="center", va="center", fontsize=12)
        elif n == 2:
            ax.text(-0.58, 0, str(weight(1, 0)), ha="center", va="center", fontsize=11)
            ax.text(0.58, 0, str(weight(0, 1)), ha="center", va="center", fontsize=11)
            ax.text(0, 0, str(weight(1, 1)), ha="center", va="center", fontsize=11)
        else:
            ax.text(-0.58, -0.28, str(weight(1, 0, 0)), ha="center", va="center", fontsize=10)
            ax.text(0.58, -0.28, str(weight(0, 1, 0)), ha="center", va="center", fontsize=10)
            ax.text(0.0, 0.62, str(weight(0, 0, 1)), ha="center", va="center", fontsize=10)
            ax.text(0.0, -0.4, str(weight(1, 1, 0)), ha="center", va="center", fontsize=10)
            ax.text(-0.28, 0.18, str(weight(1, 0, 1)), ha="center", va="center", fontsize=10)
            ax.text(0.28, 0.18, str(weight(0, 1, 1)), ha="center", va="center", fontsize=10)
            ax.text(0.0, 0.02, str(weight(1, 1, 1)), ha="center", va="center", fontsize=10)
        ax.set_xlim(-1.25, 1.25)
        ax.set_ylim(-1.1, 1.2)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.legend(
            [
                plt.Line2D([0], [0], marker="o", linestyle="", color=palette[i % len(palette)])
                for i in range(n)
            ],
            names,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.02),
            ncol=min(n, 3),
            frameon=False,
        )
    else:
        nonempty = table[table["Weight"] > 0].copy()
        labels_bar = []
        for _, row in nonempty.iterrows():
            labels_bar.append(" & ".join(
                name for name in names if int(row[name]) == 1
            ))
        ax.barh(range(len(nonempty)), nonempty["Weight"], color=palette[0], alpha=0.8)
        ax.set_yticks(range(len(nonempty)))
        ax.set_yticklabels(labels_bar)
        ax.invert_yaxis()
        ax.set_xlabel("Overlap count")

    ax.set_title(title)
    fig.tight_layout()
    return fig, ax


def vennplot_peakfile(
    peakfiles: Sequence[Union[str, Path, pd.DataFrame]],
    labels: Optional[Sequence[str]] = None,
    **kwargs,
):
    """`vennplot` for a list of peak files (peaks compared by coordinates)."""

    sets = []
    for peak in peakfiles:
        df = read_peak_file(peak)
        sets.append(
            {
                f"{r.chrom}:{int(r.start)}-{int(r.end)}"
                for r in df.itertuples(index=False)
            }
        )
    return vennplot(sets, labels=labels, **kwargs)


def vennpie(annotation: PeakAnnotation, r: float = 0.2, figsize=(6.0, 5.0), ax=None):
    """Draw nested pies summarising ChIPseeker's genomic annotation detail."""

    if not isinstance(annotation, PeakAnnotation):
        raise TypeError("vennpie expects a PeakAnnotation")
    detail = annotation.detail_genomic_annotation
    if detail.empty:
        raise ValueError("annotation has no detailGenomicAnnotation")
    distance = annotation.annotation.get("distance_to_tss", pd.Series(dtype=float))
    total = len(annotation.annotation)
    genic = int(detail["genic"].sum()) if "genic" in detail else 0
    intergenic = max(total - genic, 0)
    exon = int(detail["Exon"].sum()) if "Exon" in detail else 0
    intron = int(detail["Intron"].sum()) if "Intron" in detail else 0
    other_genic = max(genic - exon - intron, 0)
    upstream = int(((detail["Promoter"]) & (distance < 0)).sum()) if "Promoter" in detail else 0
    downstream = int(detail["downstream"].sum()) if "downstream" in detail else 0
    distal = int(detail["distal_intergenic"].sum()) if "distal_intergenic" in detail else 0

    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    outer_vals = [genic, intergenic]
    outer_labels = ["Genic", "Intergenic"]
    outer_colors = ["#3182bd", "#fec44f"]
    ax.pie(
        outer_vals,
        radius=1.0,
        colors=outer_colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.28, edgecolor="white"),
    )
    mid_vals = [exon, intron, other_genic]
    mid_labels = ["Exon", "Intron", "Other genic"]
    mid_colors = ["#9ecae1", "#fc9272", "#bebada"]
    ax.pie(
        mid_vals,
        radius=0.70,
        colors=mid_colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.26, edgecolor="white"),
    )
    inner_vals = [upstream, downstream, distal, max(intergenic - upstream - downstream - distal, 0)]
    inner_labels = ["Upstream", "Downstream", "Distal intergenic", "Other intergenic"]
    inner_colors = ["#ffeda0", "#fee0d2", "#d95f0e", "#d9d9d9"]
    ax.pie(
        inner_vals,
        radius=0.42,
        colors=inner_colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.24, edgecolor="white"),
    )
    ax.text(0, 0, f"n={total}", ha="center", va="center", fontsize=9)
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=c)
        for c in outer_colors + mid_colors + inner_colors
    ]
    labels = [
        f"{label} ({value})"
        for label, value in zip(
            outer_labels + mid_labels + inner_labels,
            outer_vals + mid_vals + inner_vals,
        )
    ]
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax.set_aspect("equal")
    ax.set_title("Genomic annotation")
    fig.tight_layout()
    return fig, ax


# ----------------------------------------------------------------------
# coverage / randomization / overlap enrichment
# ----------------------------------------------------------------------
def _chrom_sort_key(chrom: str):
    text = str(chrom)
    head = text.rstrip("0123456789")
    tail = text[len(head):]
    return (head, int(tail) if tail.isdigit() else 0, text)


def covplot(
    peak: Union[str, Path, pd.DataFrame, Sequence],
    weightCol: Optional[str] = None,
    xlab: str = "Genomic Coordinates",
    ylab: str = "Coverage",
    title: str = "",
    chrs: Optional[Sequence[str]] = None,
    xlim: Optional[Sequence[int]] = None,
    lower: Optional[float] = None,
    fill_color: Optional[str] = None,
    figsize: Optional[tuple[float, float]] = None,
):
    """Plot peak coverage along chromosomes (ChIPseeker ``covplot``)."""

    if isinstance(peak, (list, tuple)):
        frames = [_normalise_peaks(item) for item in peak]
    else:
        frames = [_normalise_peaks(peak)]
    all_chroms = sorted(
        set().union(*(set(df["chrom"].astype(str)) for df in frames)),
        key=_chrom_sort_key,
    )
    if chrs is not None:
        selected = [str(c) for c in chrs]
        all_chroms = [c for c in all_chroms if c in selected]
    if not all_chroms:
        raise ValueError("no chromosomes to plot")

    plt = _plt()
    n = len(all_chroms)
    ncols = 1 if n == 1 else 2
    nrows = int(np.ceil(n / ncols))
    figsize = figsize or (10.0, max(2.2 * nrows, 2.5))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=False)
    axes_flat = axes.ravel()
    palette = [
        fill_color or "#3182bd",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
    ]
    for ax in axes_flat[n:]:
        ax.axis("off")

    for idx, chrom in enumerate(all_chroms):
        ax = axes_flat[idx]
        for fi, df in enumerate(frames):
            sub = df[df["chrom"].astype(str) == chrom].copy()
            if sub.empty:
                continue
            weights = (
                pd.to_numeric(sub[weightCol], errors="coerce").fillna(0.0).to_numpy()
                if weightCol is not None and weightCol in sub.columns
                else np.ones(len(sub), dtype=float)
            )
            events: list[tuple[int, float]] = []
            for row, weight in zip(sub.itertuples(index=False), weights):
                events.append((int(row.start), float(weight)))
                events.append((int(row.end), -float(weight)))
            events.sort()
            xs: list[int] = []
            ys: list[float] = []
            current = 0.0
            last_x = None
            for x, delta in events:
                if last_x is not None and x != last_x:
                    xs.append(last_x)
                    ys.append(current)
                current += delta
                last_x = x
                xs.append(x)
                ys.append(current)
            if xlim is not None and len(xlim) == 2:
                lo, hi = int(xlim[0]), int(xlim[1])
                xs_arr = np.asarray(xs)
                ys_arr = np.asarray(ys)
                keep = (xs_arr >= lo) & (xs_arr <= hi)
                xs = xs_arr[keep].tolist()
                ys = ys_arr[keep].tolist()
            if lower is not None:
                ys = [y if y >= float(lower) else 0.0 for y in ys]
            color = palette[fi % len(palette)]
            ax.fill_between(xs, 0, ys, step="post", color=color, alpha=0.65, linewidth=0)
            ax.plot(xs, ys, drawstyle="steps-post", color=color, linewidth=0.8)
        ax.set_xlim(*(xlim if xlim is not None else (
            int(min(df["start"].min() for df in frames if len(df))),
            int(max(df["end"].max() for df in frames if len(df))),
        )))
        ax.set_title(chrom)
        ax.set_xlabel(xlab if idx >= n - ncols else "")
        ax.set_ylabel(ylab)
        ax.set_facecolor("#EBEBEB")
        ax.grid(color="white", linewidth=1.1)
        ax.set_axisbelow(True)
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(length=0, labelsize=9, labelcolor="#4D4D4D")
        ax.xaxis.label.set_color("#333333")
        ax.yaxis.label.set_color("#333333")
    if len(frames) > 1:
        handles = [
            plt.Line2D([0], [0], color=palette[i % len(palette)], linewidth=3)
            for i in range(len(frames))
        ]
        axes_flat[0].legend(handles, [f"Set{i + 1}" for i in range(len(frames))], frameon=False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes_flat if n > 1 else axes_flat[0]


def _chrom_sizes_from_peaks(peaks: pd.DataFrame) -> dict[str, int]:
    if peaks.empty:
        return {}
    out = peaks.groupby("chrom")["end"].max().astype(int).to_dict()
    return {str(k): max(int(v), 1) for k, v in out.items()}


def shuffle(
    peak: Union[str, Path, pd.DataFrame],
    chrom_sizes: Optional[Union[dict, pd.DataFrame]] = None,
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    n: Optional[int] = None,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """Randomly relocate peaks within chromosomes, preserving width."""

    peaks = _normalise_peaks(peak)
    if n is not None:
        peaks = peaks.sample(n=int(n), replace=len(peaks) < int(n), random_state=seed).reset_index(drop=True)
    if chrom_sizes is None and TxDb is not None:
        model = (
            TxDb
            if isinstance(TxDb, GeneAnnotation)
            else GeneAnnotation.from_gtf(TxDb)
        )
        chrom_sizes = model.genes.groupby("chrom")["end"].max().astype(int).to_dict()
    if chrom_sizes is None:
        chrom_sizes = _chrom_sizes_from_peaks(peaks)
    if isinstance(chrom_sizes, pd.DataFrame):
        size_col = "size" if "size" in chrom_sizes.columns else "end"
        chrom_sizes = dict(zip(chrom_sizes["chrom"].astype(str), chrom_sizes[size_col]))
    sizes = {str(k): int(v) for k, v in dict(chrom_sizes).items()}

    rng = np.random.default_rng(seed)
    out = peaks.copy()
    for i, row in out.iterrows():
        chrom = str(row["chrom"])
        width = int(row["end"]) - int(row["start"])
        chrom_len = max(int(sizes.get(chrom, 0)), width + 1)
        new_start = int(rng.integers(0, chrom_len - width + 1))
        out.at[i, "start"] = new_start
        out.at[i, "end"] = new_start + width
    return out


def _count_query_overlaps(query: pd.DataFrame, target: pd.DataFrame) -> int:
    if query.empty or target.empty:
        return 0
    target_by_chrom = {
        str(chrom): (
            group["start"].to_numpy(dtype=np.int64),
            group["end"].to_numpy(dtype=np.int64),
        )
        for chrom, group in target.groupby("chrom", sort=False)
    }
    count = 0
    for row in query.itertuples(index=False):
        arrays = target_by_chrom.get(str(row.chrom))
        if arrays is None:
            continue
        starts, ends = arrays
        order = np.argsort(starts, kind="mergesort")
        starts = starts[order]
        ends = ends[order]
        lo = int(np.searchsorted(starts, int(row.end), side="left"))
        if lo <= 0:
            continue
        # Any target starting before query.end and ending after query.start.
        candidate_ends = ends[:lo]
        if np.any(candidate_ends > int(row.start)):
            count += 1
    return count


def _p_adjust(pvalues: Sequence[float], method: str = "BH") -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    method = str(method).upper()
    if method in {"NONE", "NOP"}:
        return p
    n = len(p)
    if n == 0:
        return p
    if method == "BONFERRONI":
        return np.minimum(p * n, 1.0)
    if method != "BH":
        raise ValueError("pAdjustMethod must be one of BH, bonferroni or none")
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def enrich_peak_overlap(
    queryPeak: Union[str, Path, pd.DataFrame],
    targetPeak,
    chrom_sizes: Optional[Union[dict, pd.DataFrame]] = None,
    nShuffle: int = 1000,
    pAdjustMethod: str = "BH",
    seed: Optional[int] = None,
    pool: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """Permutation test for overlap between query and target peak sets."""

    query = _normalise_peaks(queryPeak)
    targets = (
        list(targetPeak)
        if isinstance(targetPeak, (list, tuple))
        else [targetPeak]
    )
    if chrom_sizes is None:
        chrom_sizes = _chrom_sizes_from_peaks(query)
    if isinstance(chrom_sizes, pd.DataFrame):
        size_col = "size" if "size" in chrom_sizes.columns else "end"
        chrom_sizes = dict(zip(chrom_sizes["chrom"].astype(str), chrom_sizes[size_col]))
    sizes = {str(k): int(v) for k, v in dict(chrom_sizes).items()}
    rng = np.random.default_rng(seed)
    rows = []
    for target in targets:
        tgt = _normalise_peaks(target)
        observed = _count_query_overlaps(query, tgt)
        pvalue = float("nan")
        if nShuffle and nShuffle > 0 and sizes:
            exceed = 0
            for _ in range(int(nShuffle)):
                shuffled = shuffle(query, chrom_sizes=sizes, seed=int(rng.integers(0, 2**31 - 1)))
                if _count_query_overlaps(shuffled, tgt) >= observed:
                    exceed += 1
            pvalue = (exceed + 1.0) / (int(nShuffle) + 1.0)
        rows.append(
            {
                "qSample": "query",
                "tSample": "target",
                "qLen": len(query),
                "tLen": len(tgt),
                "N_OL": observed,
                "pvalue": pvalue,
            }
        )
    out = pd.DataFrame(rows)
    out["p.adjust"] = _p_adjust(out["pvalue"].fillna(1.0).to_numpy(), pAdjustMethod)
    return out


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """P(X >= k) for a hypergeometric distribution."""

    if k <= 0:
        return 1.0
    upper = min(n, N)
    if k > upper:
        return 0.0
    terms = []
    for i in range(int(k), int(upper) + 1):
        terms.append(
            math.exp(
                _log_comb(n, i)
                + _log_comb(M - n, N - i)
                - _log_comb(M, N)
            )
        )
    return float(min(max(math.fsum(terms), 0.0), 1.0))


def _annotation_gene_ids(
    peak: Union[str, Path, pd.DataFrame],
    *,
    TxDb=None,
    gtf=None,
    annotation=None,
    level: str = "gene",
    tss_region: Sequence[int] = (-3000, 3000),
) -> tuple[set[str], int]:
    source = TxDb if TxDb is not None else (gtf if gtf is not None else annotation)
    if source is None:
        raise ValueError("provide TxDb/gtf/annotation for enrich_anno_overlap")
    ann = annotate_peak(
        peak,
        source,
        tss_region=tss_region,
        level=level if level in {"gene", "transcript"} else "gene",
        assign_genomic_annotation=False,
        drop_unmapped=False,
    )
    df = ann.annotation
    gene_col = "gene_id" if "gene_id" in df.columns else "geneId"
    if gene_col not in df.columns:
        raise ValueError("annotation result has no gene_id column")
    total = len(ann.gene_annotation.genes)
    return set(df[gene_col].dropna().astype(str)), total


def enrich_anno_overlap(
    queryPeak: Union[str, Path, pd.DataFrame],
    targetPeak,
    TxDb=None,
    gtf=None,
    annotation=None,
    level: str = "gene",
    tss_region: Sequence[int] = (-3000, 3000),
    pAdjustMethod: str = "BH",
    **kwargs,
) -> pd.DataFrame:
    """Hypergeometric enrichment of genes annotated to two peak sets."""

    query_genes, total_genes = _annotation_gene_ids(
        queryPeak,
        TxDb=TxDb,
        gtf=gtf,
        annotation=annotation,
        level=level,
        tss_region=tss_region,
    )
    targets = (
        list(targetPeak) if isinstance(targetPeak, (list, tuple)) else [targetPeak]
    )
    rows = []
    for target in targets:
        target_genes, _ = _annotation_gene_ids(
            target,
            TxDb=TxDb,
            gtf=gtf,
            annotation=annotation,
            level=level,
            tss_region=tss_region,
        )
        overlap_genes = query_genes & target_genes
        pvalue = _hypergeom_sf(
            len(overlap_genes),
            total_genes,
            len(query_genes),
            len(target_genes),
        )
        rows.append(
            {
                "qSample": "query",
                "tSample": "target",
                "qLen": len(query_genes),
                "tLen": len(target_genes),
                "N_OL": len(overlap_genes),
                "pvalue": pvalue,
            }
        )
    out = pd.DataFrame(rows)
    out["p.adjust"] = _p_adjust(out["pvalue"].to_numpy(), pAdjustMethod)
    return out


def upsetplot(
    annotation: PeakAnnotation,
    order_by: str = "freq",
    vennpie: bool = False,
    n_intersections: int = 20,
    figsize: tuple[float, float] = (7.0, 4.0),
    ax=None,
):
    """Approximate ChIPseeker's UpSet plot for detailed annotation flags."""

    if not isinstance(annotation, PeakAnnotation):
        raise TypeError("upsetplot expects a PeakAnnotation")
    detail = annotation.detail_genomic_annotation
    if detail.empty:
        raise ValueError("annotation has no detailGenomicAnnotation")
    bool_cols = [
        c
        for c in detail.columns
        if pd.api.types.is_bool_dtype(detail[c])
        or set(pd.unique(detail[c].dropna())).issubset({0, 1, True, False})
    ]
    if not bool_cols:
        raise ValueError("no boolean annotation columns to plot")
    flags = detail[bool_cols].fillna(False).astype(bool)
    counts = Counter(
        tuple(bool_cols[i] for i in np.flatnonzero(row.to_numpy()))
        for _, row in flags.iterrows()
    )
    items = list(counts.items())
    if str(order_by).lower() in {"degree", "deg"}:
        items.sort(key=lambda kv: (-len(kv[0]), -kv[1]))
    else:
        items.sort(key=lambda kv: (-kv[1], -len(kv[0])))
    items = items[: int(n_intersections)]

    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    labels = [
        " & ".join(parts) if parts else "None"
        for parts, _ in items
    ]
    values = [value for _, value in items]
    y = np.arange(len(items))
    ax.barh(y, values, color="#3182bd", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Number of peaks")
    ax.set_title("Annotation intersections")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig, ax


def readPeakFile(peak, **kwargs) -> pd.DataFrame:
    """ChIPseeker-style ``readPeakFile`` with ``as='GRanges'`` support."""

    as_type = kwargs.pop(
        "as_type",
        kwargs.pop("as_", kwargs.pop("as", "dataframe")),
    )
    return read_peak_file(peak, as_type=as_type, **kwargs)


# ----------------------------------------------------------------------
# ChIPseeker GEO helpers (bundled gsminfo snapshot)
# ----------------------------------------------------------------------
# ChIPseeker ships a ``gsminfo`` table (~10^5 GEO ChIP-seq samples collected
# from GEOmetadb) and derives getGEOspecies / getGEOgenomeVersion /
# getGEOInfo / downloadGEObedFiles / downloadGSMbedFiles from it.  The same
# snapshot is bundled here as ``pychipseeker/data/gsminfo.csv.gz`` with the
# original column order (data last updated upstream: Mar 03, 2015, extended
# by subsequent upstream releases).  A custom snapshot can be supplied with
# ``gsminfo=<DataFrame or path>`` to all of these functions.

_GEO_SIMPLIFY_COLUMNS = (
    "series_id",
    "gsm",
    "organism",
    "title",
    "supplementary_file",
    "genomeVersion",
    "pubmed_id",
)

# Columns R keeps verbatim; letter-only replacements of "_ch1" would be
# applied by the upstream pipeline, so the bundled table already carries
# the clean names.


def _load_gsminfo(gsminfo=None) -> pd.DataFrame:
    """Load the bundled (or user supplied) GEO ``gsminfo`` snapshot."""

    from importlib.resources import as_file, files

    if gsminfo is None:
        with as_file(files("pychipseeker").joinpath("data", "gsminfo.csv.gz")) as path:
            return pd.read_csv(
                path,
                compression="gzip",
                dtype={"series_id": "string", "gsm": "string", "genomeVersion": "string"},
            )
    if isinstance(gsminfo, pd.DataFrame):
        return gsminfo.copy()
    if isinstance(gsminfo, (str, Path)):
        path = Path(gsminfo)
        if not path.exists():
            raise FileNotFoundError(path)
        compression = "gzip" if str(path).endswith(".gz") else "infer"
        return pd.read_csv(
            path,
            compression=compression,
            dtype={"series_id": "string", "gsm": "string", "genomeVersion": "string"},
        )
    raise TypeError("gsminfo must be None, a pandas.DataFrame, or a path")


def _geo_https_url(url: str) -> str:
    """Rewrite ``ftp://`` GEO URLs to ``https://`` (upstream issue #254)."""

    out = str(url)
    if out.startswith("ftp://"):
        out = "https://" + out[len("ftp://"):]
    return out


def _download_geo_internal(info: pd.DataFrame, destDir=None) -> list[Path]:
    """Download GEO supplementary files; skip missing/failed URLs silently."""

    if destDir is None:
        destDir = Path.cwd()
    dest_dir = Path(destDir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    import urllib.request

    downloaded: list[Path] = []
    for url in info["supplementary_file"]:
        if pd.isna(url) or not str(url).strip():
            continue
        remote = _geo_https_url(url)
        destfile = dest_dir / remote.replace("\\", "/").split("/")[-1]
        if destfile.exists():
            continue
        try:
            with urllib.request.urlopen(remote, timeout=120) as response:
                destfile.write_bytes(response.read())
            downloaded.append(destfile)
        except Exception as exc:  # mirror R: message(url, ..., ", skip")
            print(f"{remote}: {exc}, skip", flush=True)
    return downloaded


def getGEOspecies(gsminfo=None) -> pd.DataFrame:
    """Species statistics collected from the bundled GEO snapshot.

    Mirrors ChIPseeker's ``getGEOspecies``: an alphabetically sorted count of
    the ``organism`` column.
    """

    species = _load_gsminfo(gsminfo)["organism"]
    counts = species.astype("string").fillna("NA").value_counts().sort_index()
    return pd.DataFrame({"species": counts.index.astype(str), "Freq": counts.to_numpy()})


def getGEOgenomeVersion(gsminfo=None) -> pd.DataFrame:
    """Genome-version statistics collected from GEO ChIP-seq samples."""

    info = _load_gsminfo(gsminfo)
    gv = info[["organism", "genomeVersion"]].copy()
    gv["genomeVersion"] = gv["genomeVersion"].astype("string").fillna("NA").astype(str)
    counts = gv["genomeVersion"].value_counts().sort_index()
    res = pd.DataFrame({"genomeVersion": counts.index.astype(str), "Freq": counts.to_numpy()})
    unique_gv = gv.assign(organism=gv["organism"].astype("string").fillna("NA").astype(str))
    unique_gv = unique_gv.drop_duplicates()
    merged = unique_gv.merge(res, on="genomeVersion", how="right")
    return merged[["organism", "genomeVersion", "Freq"]].sort_values(
        ["genomeVersion", "organism"], ignore_index=True
    )


def getGEOInfo(genome, simplify: bool = True, gsminfo=None) -> pd.DataFrame:
    """Subset the GEO snapshot by genome version (ChIPseeker ``getGEOInfo``)."""

    info = _load_gsminfo(gsminfo)
    res = info[info["genomeVersion"].astype(str) == str(genome)].copy()
    res = res.reset_index(drop=True)
    if simplify:
        columns = [c for c in _GEO_SIMPLIFY_COLUMNS if c in res.columns]
        res = res[columns]
    return res


def downloadGEObedFiles(genome, destDir=None, gsminfo=None) -> list[Path]:
    """Download all BED supplementary files for ``genome`` from GEO.

    Returns the list of destination paths created (R returns ``NULL``).
    Existing files are not re-downloaded, and unreachable URLs are skipped
    with a message, exactly like the R implementation.
    """

    info = getGEOInfo(genome, simplify=True, gsminfo=gsminfo)
    return _download_geo_internal(info, destDir)


def downloadGSMbedFiles(GSM, destDir=None, gsminfo=None) -> list[Path]:
    """Download BED supplementary files for the given GSM accessions."""

    info = _load_gsminfo(gsminfo)
    if isinstance(GSM, (str, Path)) or not hasattr(GSM, "__iter__"):
        GSM = [GSM]
    wanted = {str(g) for g in GSM}
    subset = info[info["gsm"].astype(str).isin(wanted)]
    return _download_geo_internal(subset, destDir)


def get_sample_files(directory=None) -> dict[str, Path]:
    """Paths of the bundled sample peak files (ChIPseeker ``getSampleFiles``).

    Returns ``{sample_name: path}``, with names derived from the file names
    using the same regular expression as the R implementation
    (``GSM\\d+_(\\w+)_\\w+_.*`` -> ``\\1``).

    ``directory=...`` can point at a folder of ``GSM<number>_<name>_*.bed.gz``
    files instead of the bundled one.
    """

    import re
    from importlib.resources import as_file, files

    if directory is None:
        with as_file(files("pychipseeker").joinpath("data", "GEO_sample_data")) as dir_path:
            result: dict[str, Path] = {}
            for fname in sorted(p.name for p in dir_path.iterdir() if p.is_file()):
                label = re.sub(r"GSM\d+_(\w+_\w+)_.*", r"\1", fname)
                label = re.sub(r"_Chip.+", "", label)
                result[label] = dir_path / fname
            return result
    dir_path = Path(directory)
    if not dir_path.is_dir():
        raise FileNotFoundError(f"sample data directory not found: {dir_path}")
    result = {}
    for fname in sorted(p.name for p in dir_path.iterdir() if p.is_file()):
        label = re.sub(r"GSM\d+_(\w+_\w+)_.*", r"\1", fname)
        label = re.sub(r"_Chip.+", "", label)
        result[label] = dir_path / fname
    return result


getSampleFiles = get_sample_files
get_geo_species = getGEOspecies
get_geo_genome_version = getGEOgenomeVersion
get_geo_info = getGEOInfo
download_geo_bed_files = downloadGEObedFiles
download_gsm_bed_files = downloadGSMbedFiles


# ChIPseeker spelling aliases.  Python identifiers cannot contain a dot, so
# ``as.GRanges`` and ``vennplot.peakfile`` are exposed as ``as_granges`` /
# ``vennplot_peakfile`` (and, for interactive users, as function attributes).
asGRanges = as_granges
getAnnoStat = get_anno_stat
dropAnno = drop_anno
combine_csAnno = combine_cs_anno
vennplot.peakfile = vennplot_peakfile
enrichPeakOverlap = enrich_peak_overlap
enrichAnnoOverlap = enrich_anno_overlap
detailGenomicAnnotation = detail_genomic_annotation
