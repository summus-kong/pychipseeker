"""High-level summaries for annotated ChIP-seq peak sets.

These helpers sit on top of :class:`~pychipseeker.annotation.PeakAnnotation`
and return plain DataFrames, so the same numbers can feed a table, a plot, or
a report.  All functions are dependency-light (NumPy and pandas only).

The summaries are deliberately opinionated about the questions people ask
first:

* ``annotation_summary`` -- how many peaks fall in each category, and how wide
  they are;
* ``peaks_per_gene`` / ``top_genes`` -- which genes carry the most peaks;
* ``gene_type_summary`` -- which gene biotypes the peaks are assigned to;
* ``distance_summary`` -- how close the peaks are to transcription start
  sites;
* ``signal_groups`` / ``tss_enrichment_by_quartile`` -- whether stronger peaks
  produce sharper TSS profiles.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from .annotation import PeakAnnotation

__all__ = [
    "annotation_summary",
    "distance_summary",
    "gene_type_summary",
    "peaks_per_gene",
    "signal_groups",
    "top_genes",
    "tss_enrichment_by_quartile",
]


def _annotation_frame(annotation: Union[PeakAnnotation, pd.DataFrame]) -> pd.DataFrame:
    """Return a peak table with annotation metadata and a numeric width."""

    if isinstance(annotation, PeakAnnotation):
        frame = annotation.to_dataframe(include_peak_metadata=True)
        if "score" not in frame.columns and "score" in annotation.peaks.columns:
            frame = frame.merge(
                annotation.peaks[["peak_id", "score"]],
                on="peak_id",
                how="left",
            )
    elif isinstance(annotation, pd.DataFrame):
        frame = annotation.copy()
    else:
        raise TypeError("annotation must be a PeakAnnotation or a pandas.DataFrame")

    required = {"chrom", "start", "end", "annotation"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"annotation table is missing columns: {sorted(missing)}")

    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame.dropna(subset=["start", "end"]).copy()
    frame["width"] = (frame["end"] - frame["start"]).astype(float)
    frame = frame[frame["width"] > 0].copy()
    if frame.empty:
        raise ValueError("annotation table has no valid peaks")
    return frame


def annotation_summary(
    annotation: Union[PeakAnnotation, pd.DataFrame],
) -> pd.DataFrame:
    """Summarise peak count, width and signal for each annotation category.

    Returns one row per detailed ChIPseeker label (for example
    ``Promoter (<=1kb)`` or ``5' UTR``), sorted by decreasing peak count.
    ``percent`` is the share of all annotated peaks, not of the category.
    """

    frame = _annotation_frame(annotation)
    frame["annotation"] = frame["annotation"].astype(str)
    grouped = frame.groupby("annotation", sort=False)
    out = grouped.agg(
        count=("width", "size"),
        median_width_bp=("width", "median"),
        mean_width_bp=("width", "mean"),
        total_bp=("width", "sum"),
    ).reset_index()
    out["percent"] = out["count"] / out["count"].sum() * 100.0
    if "score" in frame.columns:
        medians = grouped["score"].median()
        out["median_signal"] = out["annotation"].map(medians).astype(float)
    out = out.sort_values(
        ["count", "annotation"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)
    columns = [
        "annotation",
        "count",
        "percent",
        "median_width_bp",
        "mean_width_bp",
        "total_bp",
    ]
    if "median_signal" in out.columns:
        columns.append("median_signal")
    return out[columns]


def peaks_per_gene(
    annotation: Union[PeakAnnotation, pd.DataFrame],
) -> pd.DataFrame:
    """Return one row per gene with its peak count and peak-level summaries."""

    frame = _annotation_frame(annotation)
    if "gene_id" in frame.columns:
        keys = ["gene_id"]
    elif "gene_name" in frame.columns:
        keys = ["gene_name"]
    else:
        raise ValueError("annotation table has neither gene_id nor gene_name")

    if "gene_name" not in frame.columns:
        frame["gene_name"] = frame[keys[0]].astype(str)
    frame["gene_id"] = frame.get("gene_id", frame[keys[0]]).astype(str)
    frame["gene_name"] = frame["gene_name"].fillna(frame["gene_id"]).astype(str)

    grouped = frame.groupby("gene_id", sort=False)
    out = grouped.agg(
        gene_name=("gene_name", "first"),
        n_peaks=("width", "size"),
        total_bp=("width", "sum"),
        median_width_bp=("width", "median"),
    ).reset_index()
    if "score" in frame.columns:
        out["median_signal"] = out["gene_id"].map(grouped["score"].median())
    if "distance_to_tss" in frame.columns:
        out["median_distance_to_tss"] = out["gene_id"].map(
            grouped["distance_to_tss"].median()
        )

    sort_columns = ["n_peaks"]
    ascending = [False]
    if "median_signal" in out.columns:
        sort_columns.append("median_signal")
        ascending.append(False)
    sort_columns.append("gene_name")
    ascending.append(True)
    out = out.sort_values(
        sort_columns, ascending=ascending, kind="stable"
    ).reset_index(drop=True)
    return out


def top_genes(
    annotation: Union[PeakAnnotation, pd.DataFrame],
    n: int = 20,
) -> pd.DataFrame:
    """Return the ``n`` genes carrying the most annotated peaks."""

    if n <= 0:
        raise ValueError("n must be positive")
    out = peaks_per_gene(annotation).head(int(n)).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1, dtype=int))
    return out


def gene_type_summary(
    annotation: Union[PeakAnnotation, pd.DataFrame],
    top: int = 8,
) -> pd.DataFrame:
    """Summarise peaks by gene biotype, folding the tail into ``Other``.

    ``percent`` always refers to the complete peak set, so the returned rows
    sum to approximately 100 even when categories are folded.
    """

    if top <= 0:
        raise ValueError("top must be positive")
    frame = _annotation_frame(annotation)
    if "gene_type" not in frame.columns:
        raise ValueError("annotation table has no gene_type column")

    labels = frame["gene_type"].fillna("unknown").astype(str)
    counts = labels.value_counts()
    total = int(counts.sum())
    keep = counts.head(int(top))
    rows = [
        {
            "gene_type": str(name),
            "count": int(value),
            "percent": float(value) / total * 100.0,
        }
        for name, value in keep.items()
    ]
    rest = total - int(keep.sum())
    if rest > 0:
        rows.append(
            {
                "gene_type": "Other",
                "count": int(rest),
                "percent": float(rest) / total * 100.0,
            }
        )
    return pd.DataFrame(rows, columns=["gene_type", "count", "percent"])


def distance_summary(
    annotation: Union[PeakAnnotation, pd.DataFrame],
    thresholds: Sequence[int] = (1000, 10000, 100000),
) -> pd.DataFrame:
    """Return distance-to-TSS metrics as a two-column ``Metric``/``Value`` table."""

    frame = _annotation_frame(annotation)
    if "distance_to_tss" in frame.columns:
        distance = pd.to_numeric(frame["distance_to_tss"], errors="coerce")
    elif "nearest_distance" in frame.columns:
        distance = pd.to_numeric(frame["nearest_distance"], errors="coerce")
    else:
        raise ValueError("annotation table has no distance_to_tss column")

    distance = np.abs(distance.to_numpy(dtype=float))
    distance = distance[np.isfinite(distance)]
    if distance.size == 0:
        raise ValueError("annotation table has no finite distance values")

    rows: list[tuple[str, object]] = [
        ("Peaks", int(distance.size)),
        ("Median distance (bp)", float(np.median(distance))),
        ("Mean distance (bp)", float(distance.mean())),
        ("90th percentile (bp)", float(np.percentile(distance, 90))),
        ("Max distance (bp)", float(distance.max())),
        ("At TSS (distance = 0)", float(np.mean(distance == 0) * 100.0)),
    ]
    for cutoff in thresholds:
        cutoff = int(cutoff)
        label = f"Within {_human_bp(cutoff)}"
        rows.append((label, float(np.mean(distance <= cutoff) * 100.0)))
    return pd.DataFrame(rows, columns=["Metric", "Value"])


def _human_bp(value: int) -> str:
    if value >= 1_000_000 and value % 1_000_000 == 0:
        return f"{value // 1_000_000} Mb"
    if value >= 1_000 and value % 1_000 == 0:
        return f"{value // 1_000} kb"
    return f"{value} bp"


def signal_groups(
    tag,
    annotation: Union[PeakAnnotation, pd.DataFrame],
    n_groups: int = 4,
) -> pd.DataFrame:
    """Assign tag-matrix rows to balanced groups by peak signal.

    Each promoter inherits the strongest peak score annotated to its gene.
    When row metadata or scores are unavailable, the row mean coverage is used
    instead, so the helper still works for synthetic matrices.

    Returns one row per tag-matrix row with ``row``, ``group``, ``label`` and
    ``score`` columns.  Groups are balanced by rank, which avoids empty bins
    when many scores are tied.
    """

    if n_groups < 2:
        raise ValueError("n_groups must be at least 2")
    matrix = np.asarray(tag, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("tag must be a 2-dimensional matrix")
    n_rows = int(matrix.shape[0])
    if n_rows == 0:
        raise ValueError("tag has no rows")
    n_groups = min(int(n_groups), n_rows)

    rows = getattr(tag, "rows", None)
    score: Optional[np.ndarray] = None

    if rows is not None and len(rows) == n_rows and "gene_id" in rows.columns:
        frame = _annotation_frame(annotation)
        if "score" in frame.columns:
            per_gene = frame.groupby("gene_id")["score"].max()
            score = rows["gene_id"].astype(str).map(per_gene).to_numpy(dtype=float)
            if not np.isfinite(score).any():
                score = None

    if score is None:
        score = np.nanmean(matrix, axis=1)
    score = np.asarray(score, dtype=float)
    finite = np.isfinite(score)
    if not finite.any():
        raise ValueError("no finite signal values to group")
    # Genes without an assigned peak have no signal evidence; treating them as
    # zero keeps them in the lowest group instead of injecting a large tied
    # block at the median.
    score = np.where(finite, score, 0.0)

    order = np.argsort(score, kind="stable")
    ranks = np.empty(n_rows, dtype=int)
    ranks[order] = np.arange(n_rows)
    groups = np.minimum((ranks * int(n_groups)) // n_rows, int(n_groups) - 1)

    labels = []
    for index in range(int(n_groups)):
        if index == 0:
            labels.append(f"Q{index + 1} lowest")
        elif index == int(n_groups) - 1:
            labels.append(f"Q{index + 1} highest")
        else:
            labels.append(f"Q{index + 1}")
    return pd.DataFrame(
        {
            "row": np.arange(n_rows, dtype=int),
            "group": groups.astype(int),
            "label": [labels[g] for g in groups],
            "score": score,
        }
    )


def tss_enrichment_by_quartile(
    tag,
    annotation: Union[PeakAnnotation, pd.DataFrame],
    n_groups: int = 4,
    *,
    center_fraction: float = 0.05,
    flank_fraction: float = 0.10,
) -> pd.DataFrame:
    """Summarise TSS enrichment for signal-ranked groups of tag-matrix rows.

    Enrichment is the mean coverage in a central window around the anchor
    divided by the mean coverage in the outer flanks.  A value of 1 means the
    centre looks like the flanks; larger values mean a sharper TSS peak.
    """

    if not 0.0 < center_fraction < 1.0:
        raise ValueError("center_fraction must be in (0, 1)")
    if not 0.0 < flank_fraction < 1.0:
        raise ValueError("flank_fraction must be in (0, 1)")

    matrix = np.asarray(tag, dtype=float)
    groups = signal_groups(tag, annotation, n_groups=n_groups)
    n_bins = int(matrix.shape[1])
    upstream = getattr(tag, "upstream", None)
    downstream = getattr(tag, "downstream", None)
    if upstream is not None and downstream is not None and upstream + downstream > 0:
        center = n_bins * float(upstream) / float(upstream + downstream)
    else:
        center = (n_bins - 1) / 2.0

    half_center = max(1, int(round(n_bins * center_fraction / 2.0)))
    flank_width = max(1, int(round(n_bins * flank_fraction / 2.0)))
    center_lo = max(0, int(round(center)) - half_center)
    center_hi = min(n_bins, int(round(center)) + half_center)

    rows = []
    for index in range(int(groups["group"].max()) + 1):
        mask = groups["group"].to_numpy() == index
        block = matrix[mask]
        if block.size == 0:
            continue
        profile = np.nanmean(block, axis=0)
        center_mean = float(np.nanmean(profile[center_lo:center_hi]))
        flank_parts = []
        if flank_width > 0:
            flank_parts.append(profile[:flank_width])
            flank_parts.append(profile[-flank_width:])
        flank_mean = float(np.nanmean(np.concatenate(flank_parts)))
        enrichment = center_mean / flank_mean if flank_mean > 0 else np.nan
        rows.append(
            {
                "quartile": str(groups.loc[mask, "label"].iloc[0]),
                "n": int(mask.sum()),
                "median_signal": float(np.median(groups.loc[mask, "score"])),
                "tss_enrichment": float(enrichment),
                "max_coverage": float(np.nanmax(profile)),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "quartile",
            "n",
            "median_signal",
            "tss_enrichment",
            "max_coverage",
        ],
    )
