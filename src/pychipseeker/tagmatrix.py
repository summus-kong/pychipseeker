"""ChIPseeker-style tag matrix, profile and heatmap helpers.

ChIPseeker's ``getTagMatrix()`` builds a matrix from **peak coverage** over a
set of genomic windows (usually promoters).  This module implements that core
behaviour and also exposes a bigWig-backed matrix using ``pyBigWig`` for users
who want to use an external coverage track instead of the peak intervals
themselves.

Coordinates are 0-based, half-open throughout this module.  ``TagMatrix`` is a
``numpy.ndarray`` subclass carrying the ChIPseeker-style metadata
(``upstream``, ``downstream``, ``type``, ``by``, ``label`` and
``is_binning``) needed by the plotting functions.
"""

from __future__ import annotations

from pathlib import Path
from statistics import NormalDist
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd

from .annotation import GeneAnnotation, _normalise_peaks
from .io import read_peaks


class TagMatrix(np.ndarray):
    """A float matrix with ChIPseeker-style region metadata."""

    def __new__(
        cls,
        matrix: np.ndarray,
        *,
        upstream: Optional[int] = None,
        downstream: Optional[int] = None,
        type: Optional[str] = None,
        by: Optional[str] = None,
        label: Optional[str] = None,
        is_binning: bool = False,
        rows: Optional[pd.DataFrame] = None,
    ) -> "TagMatrix":
        obj = np.asarray(matrix, dtype=float).view(cls)
        obj.upstream = upstream
        obj.downstream = downstream
        obj.type = type
        obj.by = by
        obj.label = label
        obj.is_binning = bool(is_binning)
        obj.rows = rows
        return obj

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return
        self.upstream = getattr(obj, "upstream", None)
        self.downstream = getattr(obj, "downstream", None)
        self.type = getattr(obj, "type", None)
        self.by = getattr(obj, "by", None)
        self.label = getattr(obj, "label", None)
        self.is_binning = getattr(obj, "is_binning", False)
        self.rows = getattr(obj, "rows", None)

    @property
    def matrix(self) -> np.ndarray:
        return np.asarray(self)

    def __array_wrap__(self, out_arr, context=None, return_scalar=False):
        # Keep ordinary numpy ufuncs returning plain ndarrays; metadata is only
        # meaningful for the original matrix, not for a derived mean/std.
        if out_arr.ndim == 0:
            return out_arr[()]
        return np.asarray(out_arr)


def _resolve_annotation(
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
) -> GeneAnnotation:
    source = annotation if annotation is not None else gtf
    if source is None:
        source = TxDb
    if source is None:
        raise ValueError("provide an annotation/gtf or a GeneAnnotation object")
    if isinstance(source, GeneAnnotation):
        return source
    if isinstance(source, (str, Path, pd.DataFrame)):
        return GeneAnnotation.from_gtf(source)
    raise TypeError("annotation must be a GTF path, DataFrame or GeneAnnotation")


def _feature_table(annotation: GeneAnnotation, by: str) -> pd.DataFrame:
    by_lower = str(by).lower()
    if by_lower in {"gene", "transcript"}:
        return annotation.features(by_lower).copy()
    if by_lower in {"exon", "intron"}:
        return annotation.category_table(by_lower.capitalize()).copy()
    if by_lower in {"5utr", "five_prime_utr", "five-utr"}:
        return annotation.category_table("5UTR").copy()
    if by_lower in {"3utr", "three_prime_utr", "three-utr"}:
        return annotation.category_table("3UTR").copy()
    if by_lower == "utr":
        five = annotation.category_table("5UTR").copy()
        three = annotation.category_table("3UTR").copy()
        five["utr_type"] = "5UTR"
        three["utr_type"] = "3UTR"
        return pd.concat([five, three], ignore_index=True)
    raise ValueError(
        "by must be one of gene, transcript, exon, intron, 3UTR, 5UTR or UTR"
    )


def _normalise_region_type(value: str) -> str:
    text = str(value).lower().replace("-", "_")
    if text in {"start_site", "start", "tss"}:
        return "start_site"
    if text in {"end_site", "end", "tts"}:
        return "end_site"
    if text == "body":
        return "body"
    raise ValueError("type must be 'start_site', 'end_site' or 'body'")


def _region_label(region_type: str, by: str) -> str:
    if region_type == "start_site":
        return "TSS" if str(by).lower() in {"gene", "transcript"} else "Start"
    if region_type == "end_site":
        return "TES" if str(by).lower() in {"gene", "transcript"} else "End"
    return f"{by} body"


def _make_region_table(
    table: pd.DataFrame,
    *,
    region_type: str,
    upstream: int,
    downstream: int,
    by: str,
) -> pd.DataFrame:
    if table.empty:
        out = pd.DataFrame(columns=["chrom", "start", "end", "strand"])
        out.attrs.update(
            type=region_type,
            by=by,
            label=_region_label(region_type, by),
            upstream=int(upstream),
            downstream=int(downstream),
        )
        return out

    df = table.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end"]).copy()
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    if "strand" not in df.columns:
        df["strand"] = "*"
    df["strand"] = df["strand"].astype(str)
    df.loc[~df["strand"].isin(["+", "-"]), "strand"] = "*"

    if region_type == "body":
        # For bodies, upstream/downstream extend the annotated interval.
        start = df["start"].to_numpy(dtype=np.int64) - int(upstream)
        end = df["end"].to_numpy(dtype=np.int64) + int(downstream)
    else:
        if region_type == "start_site":
            # 0-based TSS base: gene start on '+', gene end-1 on '-'.
            tss = df["tss"].to_numpy(dtype=np.int64) if "tss" in df.columns else np.where(
                df["strand"].to_numpy() == "-",
                df["end"].to_numpy(dtype=np.int64) - 1,
                df["start"].to_numpy(dtype=np.int64),
            )
            plus = df["strand"].to_numpy() == "+"
            minus = df["strand"].to_numpy() == "-"
            start = np.where(plus, tss - int(upstream), df["start"].to_numpy())
            end = np.where(plus, tss + int(downstream), df["end"].to_numpy())
            # Minus-strand: upstream is to the right, downstream to the left.
            start = np.where(minus, tss + 1 - int(downstream), start)
            end = np.where(minus, tss + 1 + int(upstream), end)
            # Unknown strand: treat like '+'.
            unknown = ~(plus | minus)
            start = np.where(unknown, tss - int(upstream), start)
            end = np.where(unknown, tss + int(downstream), end)
        else:  # end_site
            plus = df["strand"].to_numpy() == "+"
            minus = df["strand"].to_numpy() == "-"
            start = np.where(plus, df["end"].to_numpy() - int(upstream), df["start"].to_numpy())
            end = np.where(plus, df["end"].to_numpy() + int(downstream), df["end"].to_numpy())
            start = np.where(minus, df["start"].to_numpy() + 1 - int(downstream), start)
            end = np.where(minus, df["start"].to_numpy() + 1 + int(upstream), end)
            unknown = ~(plus | minus)
            start = np.where(unknown, df["end"].to_numpy() - int(upstream), start)
            end = np.where(unknown, df["end"].to_numpy() + int(downstream), end)

    df["start"] = start.astype(np.int64)
    df["end"] = end.astype(np.int64)
    df = df[df["end"] > df["start"]].copy()
    df = df.drop_duplicates(["chrom", "start", "end", "strand"]).reset_index(drop=True)
    df.attrs.update(
        type=region_type,
        by=by,
        label=_region_label(region_type, by),
        upstream=int(upstream),
        downstream=int(downstream),
    )
    return df


def get_bio_region(
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    upstream: int = 1000,
    downstream: int = 1000,
    by: str = "gene",
    type: str = "start_site",
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """Return ChIPseeker-style promoter/body regions from a GTF annotation.

    ``TxDb`` is accepted for API familiarity; in Python it can be a
    :class:`~pychipseeker.GeneAnnotation`, a GTF/GFF path, or a feature table.
    """

    model = _resolve_annotation(annotation=annotation, gtf=gtf, TxDb=TxDb)
    region_type = _normalise_region_type(type)
    table = _feature_table(model, by)
    return _make_region_table(
        table,
        region_type=region_type,
        upstream=int(upstream),
        downstream=int(downstream),
        by=by,
    )


def get_promoters(
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    upstream: int = 1000,
    downstream: int = 1000,
    by: str = "gene",
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """ChIPseeker-compatible promoter regions (``start_site`` bioregions)."""

    return get_bio_region(
        TxDb=TxDb,
        upstream=upstream,
        downstream=downstream,
        by=by,
        type="start_site",
        annotation=annotation,
        gtf=gtf,
    )


def make_bio_region_from_granges(
    gr: pd.DataFrame,
    by: str,
    type: str,
    upstream: int = 1000,
    downstream: int = 1000,
) -> pd.DataFrame:
    """Make a bioregion from an arbitrary interval table.

    The Python analogue of ChIPseeker's ``GRanges`` input is a DataFrame with
    ``chrom``, ``start``, ``end`` and optionally ``strand``/``tss`` columns.
    """

    if not isinstance(gr, pd.DataFrame):
        raise TypeError("gr must be a pandas.DataFrame with genomic intervals")
    needed = {"chrom", "start", "end"}
    missing = needed - set(gr.columns)
    if missing:
        raise ValueError(f"gr is missing columns: {sorted(missing)}")
    return _make_region_table(
        gr,
        region_type=_normalise_region_type(type),
        upstream=int(upstream),
        downstream=int(downstream),
        by=by,
    )


def _window_frame(windows: Union[pd.DataFrame, str, Path]) -> pd.DataFrame:
    if isinstance(windows, (str, Path)):
        windows = pd.read_csv(windows, sep="\t")
    if not isinstance(windows, pd.DataFrame):
        raise TypeError("windows must be a DataFrame or a tab-delimited path")
    needed = {"chrom", "start", "end"}
    missing = needed - set(windows.columns)
    if missing:
        raise ValueError(f"windows is missing columns: {sorted(missing)}")
    out = windows.copy()
    out["chrom"] = out["chrom"].astype(str)
    out["start"] = pd.to_numeric(out["start"], errors="raise").astype(np.int64)
    out["end"] = pd.to_numeric(out["end"], errors="raise").astype(np.int64)
    if "strand" not in out.columns:
        out["strand"] = "*"
    out["strand"] = out["strand"].astype(str)
    out.loc[~out["strand"].isin(["+", "-"]), "strand"] = "*"
    return out


def _peak_arrays(peaks: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, group in peaks.groupby("chrom", sort=False):
        group = group.sort_values("start", kind="mergesort")
        starts = group["start"].to_numpy(dtype=np.int64)
        ends = group["end"].to_numpy(dtype=np.int64)
        weights = (
            group["_weight"].to_numpy(dtype=float)
            if "_weight" in group.columns
            else np.ones(len(group), dtype=float)
        )
        out[str(chrom)] = (starts, ends, weights)
    return out


def _coverage_matrix(
    peaks: pd.DataFrame,
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
) -> tuple[np.ndarray, pd.DataFrame, bool]:
    if windows.empty:
        raise ValueError("windows is empty")
    widths = (windows["end"] - windows["start"]).to_numpy(dtype=np.int64)
    if n_bins is None and len(np.unique(widths)) != 1:
        raise ValueError(
            "windows must have equal widths when nbin is not supplied; "
            "use nbin=... for variable-width regions"
        )

    by_chrom = _peak_arrays(peaks)
    rows: list[np.ndarray] = []
    kept: list[pd.Series] = []
    is_binning = n_bins is not None

    for _, window in windows.iterrows():
        chrom = str(window["chrom"])
        start = int(window["start"])
        end = int(window["end"])
        width = end - start
        if width <= 0:
            continue
        bins = int(n_bins) if n_bins is not None else int(width)
        if bins <= 1:
            raise ValueError("nbin must be > 1")
        if chrom not in by_chrom:
            continue
        pstarts, pends, pweights = by_chrom[chrom]
        lo = int(np.searchsorted(pends, start, side="right"))
        hi = int(np.searchsorted(pstarts, end, side="left"))
        if hi <= lo:
            continue

        edges = np.linspace(start, end, bins + 1)
        values = np.zeros(bins, dtype=float)
        for j in range(lo, hi):
            a = max(start, int(pstarts[j]))
            b = min(end, int(pends[j]))
            if b <= a:
                continue
            weight = float(pweights[j])
            first = int(np.searchsorted(edges, a, side="right") - 1)
            first = max(0, min(first, bins - 1))
            last = int(np.searchsorted(edges, b, side="left") - 1)
            last = max(0, min(last, bins - 1))
            for bi in range(first, last + 1):
                left = max(a, edges[bi])
                right = min(b, edges[bi + 1])
                if right > left:
                    values[bi] += weight * float(right - left)

        if np.all(values == 0):
            continue
        if str(window["strand"]) == "-" and not ignore_strand:
            values = values[::-1]
        rows.append(values)
        kept.append(window)

    if not rows:
        raise ValueError("no windows overlap any peak")
    matrix = np.vstack(rows)
    row_frame = pd.DataFrame(kept).reset_index(drop=True)
    return matrix, row_frame, is_binning


def _bigwig_matrix(
    bigwig: Union[str, Path],
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
    value_type: str = "sum",
) -> tuple[np.ndarray, pd.DataFrame, bool]:
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - dependency installation issue
        raise ImportError(
            "pyBigWig is required for bigWig-backed tag matrices; "
            "reinstall pychipseeker (pyBigWig is a declared dependency)"
        ) from exc

    if str(value_type) not in {"sum", "mean", "max", "min", "coverage", "std"}:
        raise ValueError("unsupported bigWig value_type")

    bw = pyBigWig.open(str(bigwig))
    if bw is None:
        raise ValueError(f"could not open bigWig: {bigwig}")
    chrom_sizes = bw.chroms()
    rows: list[np.ndarray] = []
    kept: list[pd.Series] = []
    is_binning = n_bins is not None
    try:
        for _, window in windows.iterrows():
            chrom = str(window["chrom"])
            start = int(window["start"])
            end = int(window["end"])
            chrom_len = chrom_sizes.get(chrom)
            if chrom_len is None or start < 0 or end > chrom_len or end <= start:
                continue
            bins = int(n_bins) if n_bins is not None else int(end - start)
            if bins <= 1:
                raise ValueError("nbin must be > 1")
            vals = bw.stats(chrom, start, end, nBins=bins, type=value_type)
            if vals is None:
                continue
            values = np.asarray(
                [np.nan if v is None else float(v) for v in vals], dtype=float
            )
            values = np.nan_to_num(values, nan=0.0)
            if str(window["strand"]) == "-" and not ignore_strand:
                values = values[::-1]
            rows.append(values)
            kept.append(window)
    finally:
        bw.close()

    if not rows:
        raise ValueError("no windows overlap the bigWig track")
    return np.vstack(rows), pd.DataFrame(kept).reset_index(drop=True), is_binning


def get_tag_matrix(
    peak: Optional[Union[str, Path, pd.DataFrame]] = None,
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
    windows: Optional[Union[pd.DataFrame, str, Path]] = None,
    type: str = "start_site",
    by: str = "gene",
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    weightCol: Optional[str] = None,
    nbin: Optional[int] = None,
    verbose: bool = False,
    ignore_strand: bool = False,
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    bigwig: Optional[Union[str, Path]] = None,
    value_type: str = "sum",
) -> Union[TagMatrix, list[TagMatrix]]:
    """Build a ChIPseeker-style tag matrix.

    With ``bigwig=...`` (or a ``peak`` path ending in ``.bw``/``.bigWig``),
    the matrix is extracted from a bigWig track using ``pyBigWig``.  Otherwise
    the matrix is built from peak coverage, exactly like ChIPseeker's
    ``getTagMatrix()``.

    If ``windows`` is a list/tuple of DataFrames, a list of matrices is
    returned, mirroring ChIPseeker's multi-window behaviour.

    If ``windows`` is omitted, a bioregion is built from ``TxDb``/``gtf``/
    ``annotation`` using ``upstream``, ``downstream``, ``by`` and ``type``.
    """

    if bigwig is None and isinstance(peak, (str, Path)):
        suffix = str(peak).lower()
        if suffix.endswith((".bw", ".bigwig")):
            bigwig = peak
            peak = None

    if isinstance(windows, (list, tuple)):
        if not windows:
            raise ValueError("windows is empty")
        matrices = []
        for i, one_window in enumerate(windows):
            matrix = get_tag_matrix(
                peak=peak,
                upstream=upstream,
                downstream=downstream,
                windows=one_window,
                type=type,
                by=by,
                TxDb=TxDb,
                weightCol=weightCol,
                nbin=nbin,
                verbose=verbose,
                ignore_strand=ignore_strand,
                annotation=annotation,
                gtf=gtf,
                bigwig=bigwig,
                value_type=value_type,
            )
            matrix.label = (
                getattr(one_window, "attrs", {}).get("label")
                or getattr(one_window, "attrs", {}).get("type")
                or f"Window{i + 1}"
            )
            matrices.append(matrix)
        return matrices

    if windows is None:
        if upstream is None:
            upstream = 3000
        if downstream is None:
            downstream = 3000
        windows_df = get_bio_region(
            TxDb=TxDb,
            upstream=int(upstream),
            downstream=int(downstream),
            by=by,
            type=type,
            annotation=annotation,
            gtf=gtf,
        )
    else:
        windows_df = _window_frame(windows)
        upstream = int(upstream) if upstream is not None else windows_df.attrs.get("upstream")
        downstream = (
            int(downstream) if downstream is not None else windows_df.attrs.get("downstream")
        )
        if by == "gene" and windows_df.attrs.get("by"):
            by = str(windows_df.attrs["by"])
        if type == "start_site" and windows_df.attrs.get("type"):
            type = str(windows_df.attrs["type"])

    if bigwig is not None:
        matrix, kept, is_binning = _bigwig_matrix(
            bigwig,
            windows_df,
            n_bins=nbin,
            ignore_strand=ignore_strand,
            value_type=value_type,
        )
    else:
        if peak is None:
            raise ValueError("provide peak intervals or bigwig=...")
        peaks = _normalise_peaks(peak)
        if weightCol is not None:
            if weightCol not in peaks.columns:
                # Preserve extra BED columns through _normalise_peaks; common
                # weight columns are named in the BED extra columns.
                raise KeyError(f"weight column not found: {weightCol}")
            peaks["_weight"] = pd.to_numeric(peaks[weightCol], errors="coerce").fillna(0.0)
        matrix, kept, is_binning = _coverage_matrix(
            peaks,
            windows_df,
            n_bins=nbin,
            ignore_strand=ignore_strand,
        )

    if verbose:
        print(f">> tag matrix: {matrix.shape[0]} regions x {matrix.shape[1]} bins")
    return TagMatrix(
        matrix,
        upstream=upstream,
        downstream=downstream,
        type=type,
        by=by,
        label=windows_df.attrs.get("label"),
        is_binning=is_binning,
        rows=kept,
    )


def get_tag_matrix_from_bigwig(
    bigwig: Union[str, Path],
    windows: Union[pd.DataFrame, str, Path],
    nbin: Optional[int] = None,
    ignore_strand: bool = False,
    value_type: str = "sum",
) -> TagMatrix:
    """Convenience wrapper for a bigWig-backed :class:`TagMatrix`."""

    windows_df = _window_frame(windows)
    matrix, kept, is_binning = _bigwig_matrix(
        bigwig,
        windows_df,
        n_bins=nbin,
        ignore_strand=ignore_strand,
        value_type=value_type,
    )
    return TagMatrix(
        matrix,
        upstream=windows_df.attrs.get("upstream"),
        downstream=windows_df.attrs.get("downstream"),
        type=windows_df.attrs.get("type"),
        by=windows_df.attrs.get("by"),
        label=windows_df.attrs.get("label"),
        is_binning=is_binning,
        rows=kept,
    )


def example_signal_path() -> Path:
    """Path to the bundled ENCODE H3K4me3 demo bigWig.

    The bundled file is a compact subset of ENCODE ``ENCFF893ZGR``
    (H3K4me3 fold change over control, *Mus musculus* mm10, experiment
    ``ENCSR000CDQ``): one 20 Mb window of chr5 re-binned to 25 bp means.
    See ``example_signal_provenance()`` for full details.
    """
    from importlib.resources import as_file, files

    with as_file(
        files("pychipseeker").joinpath("data", "ENCODE_demo_signal.bw")
    ) as path:
        return Path(path)


def example_signal_provenance() -> dict:
    """Provenance record for :func:`example_signal_path`."""
    import json
    from importlib.resources import as_file, files

    with as_file(
        files("pychipseeker").joinpath("data", "ENCODE_demo_signal.json")
    ) as path:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)


def _as_matrix(tag_matrix) -> tuple[np.ndarray, Optional[int], Optional[int]]:
    if isinstance(tag_matrix, TagMatrix):
        return np.asarray(tag_matrix, dtype=float), tag_matrix.upstream, tag_matrix.downstream
    arr = np.asarray(tag_matrix, dtype=float)
    if arr.ndim != 2:
        raise ValueError("tagMatrix must be 2-dimensional")
    return arr, getattr(tag_matrix, "upstream", None), getattr(tag_matrix, "downstream", None)


def _x_coordinates(
    n: int, upstream: Optional[int], downstream: Optional[int]
) -> np.ndarray:
    if upstream is not None and downstream is not None:
        return np.linspace(-int(upstream), int(downstream), n)
    return np.arange(n, dtype=float)


def _plt():
    import matplotlib.pyplot as plt

    return plt


def plot_avg_prof(
    tagMatrix,
    xlim: Optional[Sequence[int]] = None,
    xlab: str = "Genomic Region (5'->3')",
    ylab: str = "Peak Count Frequency",
    conf: Optional[float] = 0.95,
    facet: str = "none",
    free_y: bool = True,
    origin_label: str = "TSS",
    ax=None,
    figsize: tuple[float, float] = (6.0, 4.0),
    color: Optional[str] = None,
    label: Optional[str] = None,
    **kwargs,
):
    """Plot the mean tag-matrix profile with an optional confidence band.

    This follows ``ChIPseeker::plotAvgProf`` semantics.  The confidence band is
    a percentile bootstrap of the row-wise mean (ChIPseeker uses a bootstrap
    confidence interval as well); exact interval values are not intended to be
    byte-identical.
    """

    plt = _plt()
    arr, upstream, downstream = _as_matrix(tagMatrix)
    if arr.shape[0] == 0:
        raise ValueError("tagMatrix has no rows")
    if xlim is not None:
        if len(xlim) != 2:
            raise ValueError("xlim must be two values")
        upstream, downstream = abs(int(xlim[0])), int(xlim[1])

    x = _x_coordinates(arr.shape[1], upstream, downstream)
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(arr, axis=0)
        n = np.sum(~np.isnan(arr), axis=0)
        sem = np.nanstd(arr, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))

    title = kwargs.pop("title", "Peak profile")
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    from .plotting import _title as _title_style  # deferred: avoid import cycles

    line = ax.plot(
        x,
        mean,
        color=color or "#E45756",
        linewidth=2.0,
        label=label,
        **kwargs,
    )[0]
    # Soft area under the mean curve for a cleaner visual mass.
    ax.fill_between(x, 0.0, mean, color=line.get_color(), alpha=0.08, linewidth=0, zorder=0.5)
    # Highlight the profile maximum.
    peak_idx = int(np.nanargmax(mean))
    ax.plot(
        x[peak_idx],
        mean[peak_idx],
        "o",
        markersize=4.5,
        color=line.get_color(),
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

    if conf is not None:
        conf = float(conf)
        if not 0 < conf < 1:
            raise ValueError("conf must be between 0 and 1")
        if conf > 0:
            z = NormalDist().inv_cdf(0.5 + conf / 2.0)
            ax.fill_between(
                x,
                mean - z * sem,
                mean + z * sem,
                color=line.get_color(),
                alpha=0.22,
                linewidth=0,
            )

    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    # ChIPseeker draws a dashed baseline at zero signal.
    ax.axhline(0, color="#bbbbbb", linewidth=0.8, linestyle="--", zorder=0)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if label:
        ax.legend(frameon=False, fontsize=9)
    _title_style(ax, title)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=0.8, labelsize=9, labelcolor="#4d4d4d", color="#bfbfbf")
    fig.tight_layout()
    return fig, ax


def plot_peak_prof2(
    peak: Union[str, Path, pd.DataFrame],
    upstream: int = 3000,
    downstream: int = 3000,
    conf: Optional[float] = 0.95,
    by: str = "gene",
    type: str = "start_site",
    weightCol: Optional[str] = None,
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    xlab: str = "Genomic Region (5'->3')",
    ylab: str = "Peak Count Frequency",
    facet: str = "none",
    free_y: bool = True,
    verbose: bool = False,
    nbin: Optional[int] = None,
    ignore_strand: bool = False,
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    **kwargs,
):
    """Build a tag matrix from peaks and plot its average profile."""

    tag_matrix = get_tag_matrix(
        peak=peak,
        upstream=upstream,
        downstream=downstream,
        type=type,
        by=by,
        TxDb=TxDb,
        weightCol=weightCol,
        nbin=nbin,
        verbose=verbose,
        ignore_strand=ignore_strand,
        annotation=annotation,
        gtf=gtf,
    )
    return plot_avg_prof(
        tag_matrix,
        xlab=xlab,
        ylab=ylab,
        conf=conf,
        facet=facet,
        free_y=free_y,
        **kwargs,
    )


def plot_peak_prof(
    tagMatrix=None,
    peak: Optional[Union[str, Path, pd.DataFrame]] = None,
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
    conf: Optional[float] = 0.95,
    by: str = "gene",
    type: str = "start_site",
    windows_name: Optional[Sequence[str]] = None,
    weightCol: Optional[str] = None,
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    xlab: str = "Genomic Region (5'->3')",
    ylab: str = "Peak Count Frequency",
    facet: str = "row",
    free_y: bool = True,
    verbose: bool = False,
    nbin: Optional[int] = None,
    ignore_strand: bool = False,
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    **kwargs,
):
    """Plot a profile from a tag matrix, or build one from peaks."""

    if tagMatrix is not None:
        if isinstance(tagMatrix, (list, tuple)):
            return plot_multi_prof(
                tagMatrix,
                xlab=xlab,
                ylab=ylab,
                conf=conf,
                facet=facet,
                free_y=free_y,
                **kwargs,
            )
        return plot_avg_prof(
            tagMatrix,
            xlab=xlab,
            ylab=ylab,
            conf=conf,
            facet=facet,
            free_y=free_y,
            **kwargs,
        )

    if peak is None:
        raise ValueError("provide tagMatrix or peak")
    return plot_peak_prof2(
        peak=peak,
        upstream=3000 if upstream is None else upstream,
        downstream=3000 if downstream is None else downstream,
        conf=conf,
        by=by,
        type=type,
        weightCol=weightCol,
        TxDb=TxDb,
        xlab=xlab,
        ylab=ylab,
        facet=facet,
        free_y=free_y,
        verbose=verbose,
        nbin=nbin,
        ignore_strand=ignore_strand,
        annotation=annotation,
        gtf=gtf,
        **kwargs,
    )


def plot_multi_prof(
    tagMatrices: Sequence,
    xlab: str = "Genomic Region (5'->3')",
    ylab: str = "Peak Count Frequency",
    conf: Optional[float] = 0.95,
    facet: str = "row",
    free_y: bool = True,
    **kwargs,
):
    """Plot several tag matrices, one line (or facet) per matrix."""

    plt = _plt()
    items = list(tagMatrices)
    if not items:
        raise ValueError("tagMatrices is empty")
    if facet in {"row", "column"}:
        fig, axes = plt.subplots(
            len(items),
            1 if facet == "row" else len(items),
            figsize=(6.0, max(2.0, 2.2 * len(items))),
            squeeze=False,
            sharex=True,
        )
        axes_flat = axes.ravel()
        for i, matrix in enumerate(items):
            arr, upstream, downstream = _as_matrix(matrix)
            label = getattr(matrix, "label", None) or f"Set{i + 1}"
            plot_avg_prof(
                matrix,
                xlab=xlab if i == len(items) - 1 else "",
                ylab=ylab,
                conf=conf,
                free_y=free_y,
                ax=axes_flat[i],
                label=label,
                **kwargs,
            )
        fig.tight_layout()
        return fig, axes_flat[0]

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for i, matrix in enumerate(items):
        label = getattr(matrix, "label", None) or f"Set{i + 1}"
        plot_avg_prof(
            matrix,
            xlab=xlab,
            ylab=ylab,
            conf=conf,
            free_y=free_y,
            ax=ax,
            label=label,
            **kwargs,
        )
    ax.legend(frameon=False)
    return fig, ax


def tag_heatmap(
    tagMatrix,
    xlab: str = "Genomic Region (5'->3')",
    ylab: str = "",
    title: Optional[str] = None,
    palette: Optional[str] = None,
    ax=None,
    figsize: tuple[float, float] = (6.0, 4.0),
    **kwargs,
):
    """Heatmap of a tag matrix, centered on the TSS/region anchor."""

    plt = _plt()
    arr, upstream, downstream = _as_matrix(tagMatrix)
    if arr.shape[0] == 0:
        raise ValueError("tagMatrix has no rows")
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # ChIPseeker's heatmap is row-sorted by signal strength.  Strongest rows sit
    # on top and the map fades downwards; sorting also makes the visual pattern
    # stable across input orders.
    order = np.argsort(-np.nanmean(arr, axis=1), kind="stable")
    data = arr[order]
    extent = (
        -int(upstream) if upstream is not None else 0,
        int(downstream) if downstream is not None else data.shape[1],
        data.shape[0],
        0,
    )
    from matplotlib import colormaps
    from .plotting import _title as _title_style  # deferred: avoid import cycles

    cmap_obj = colormaps[(kwargs.pop("cmap", None) or palette or "Reds")].with_extremes(
        bad="#f5f5f5"
    )
    im = ax.imshow(
        np.ma.masked_invalid(data),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap_obj,
        extent=extent,
        **kwargs,
    )
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab or "Peaks")
    cbar = fig.colorbar(im, ax=ax, label="Signal", shrink=0.85, pad=0.02)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=8, length=2.5, color="#bfbfbf")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=0.8, labelsize=9, labelcolor="#4d4d4d", color="#bfbfbf")
    _title_style(ax, title or "Tag heatmap")
    fig.tight_layout()
    return fig, ax


def peak_heatmap(
    peak: Union[str, Path, pd.DataFrame],
    weightCol: Optional[str] = None,
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    upstream: int = 1000,
    downstream: int = 1000,
    xlab: str = "",
    ylab: str = "",
    title: Optional[str] = None,
    palette: Optional[str] = None,
    verbose: bool = False,
    by: str = "gene",
    type: str = "start_site",
    nbin: Optional[int] = None,
    ignore_strand: bool = False,
    windows: Optional[Union[pd.DataFrame, str, Path]] = None,
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    **kwargs,
):
    """Build a tag matrix from peaks and draw a heatmap."""

    if windows is None:
        windows = get_bio_region(
            TxDb=TxDb,
            upstream=upstream,
            downstream=downstream,
            by=by,
            type=type,
            annotation=annotation,
            gtf=gtf,
        )
    tag_matrix = get_tag_matrix(
        peak=peak,
        windows=windows,
        weightCol=weightCol,
        nbin=nbin,
        verbose=verbose,
        ignore_strand=ignore_strand,
        upstream=upstream,
        downstream=downstream,
        type=type,
        by=by,
    )
    return tag_heatmap(
        tag_matrix,
        xlab=xlab,
        ylab=ylab,
        title=title,
        palette=palette,
        **kwargs,
    )


def peak_profile_heatmap(
    peak: Union[str, Path, pd.DataFrame],
    upstream: int = 3000,
    downstream: int = 3000,
    by: str = "gene",
    type: str = "start_site",
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    weightCol: Optional[str] = None,
    nbin: Optional[int] = None,
    conf: Optional[float] = 0.95,
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    **kwargs,
):
    """Combine average profile and heatmap in one figure."""

    tag_matrix = get_tag_matrix(
        peak=peak,
        upstream=upstream,
        downstream=downstream,
        by=by,
        type=type,
        TxDb=TxDb,
        weightCol=weightCol,
        nbin=nbin,
        annotation=annotation,
        gtf=gtf,
    )
    plt = _plt()
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(6.5, 6.0),
        gridspec_kw={"height_ratios": [1, 2]},
        sharex=False,
    )
    plot_avg_prof(tag_matrix, conf=conf, ax=axes[0], **kwargs)
    tag_heatmap(tag_matrix, ax=axes[1], **kwargs)
    fig.tight_layout()
    return fig, axes


def plot_peak_prof_multi_windows(
    peak: Union[str, Path, pd.DataFrame],
    windows: Sequence[pd.DataFrame],
    windows_name: Optional[Sequence[str]] = None,
    weightCol: Optional[str] = None,
    nbin: Optional[int] = None,
    conf: Optional[float] = 0.95,
    ignore_strand: bool = False,
    **kwargs,
):
    """Plot profiles for several bioregions/windows on the same peak set."""

    matrices = []
    for i, window in enumerate(windows):
        matrix = get_tag_matrix(
            peak=peak,
            windows=window,
            weightCol=weightCol,
            nbin=nbin,
            ignore_strand=ignore_strand,
        )
        matrix.label = (
            windows_name[i]
            if windows_name is not None and i < len(windows_name)
            else getattr(window, "attrs", {}).get("label", f"Window{i + 1}")
        )
        matrices.append(matrix)
    return plot_multi_prof(matrices, conf=conf, **kwargs)


def peak_heatmap_multiple_sets(
    peaks: Sequence[Union[str, Path, pd.DataFrame]],
    labels: Optional[Sequence[str]] = None,
    upstream: int = 1000,
    downstream: int = 1000,
    by: str = "gene",
    type: str = "start_site",
    TxDb: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    weightCol: Optional[str] = None,
    nbin: Optional[int] = None,
    ignore_strand: bool = False,
    *,
    annotation: Optional[Union[GeneAnnotation, str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    **kwargs,
):
    """Heatmaps for several peak sets over the same bioregion."""

    windows = get_bio_region(
        TxDb=TxDb,
        upstream=upstream,
        downstream=downstream,
        by=by,
        type=type,
        annotation=annotation,
        gtf=gtf,
    )
    plt = _plt()
    n = len(peaks)
    fig, axes = plt.subplots(
        1,
        n,
        figsize=(max(3.0 * n, 4.0), 4.0),
        squeeze=False,
    )
    for i, peak in enumerate(peaks):
        matrix = get_tag_matrix(
            peak=peak,
            windows=windows,
            weightCol=weightCol,
            nbin=nbin,
            ignore_strand=ignore_strand,
        )
        title = labels[i] if labels is not None and i < len(labels) else f"Set{i + 1}"
        tag_heatmap(matrix, title=title, ax=axes[0, i], **kwargs)
    fig.tight_layout()
    return fig, axes[0]


# ChIPseeker spelling aliases.
getTagMatrix = get_tag_matrix
getPromoters = get_promoters
getBioRegion = get_bio_region
makeBioRegionFromGranges = make_bio_region_from_granges
plotAvgProf = plot_avg_prof
plotAvgProf2 = plot_peak_prof2
plotPeakProf = plot_peak_prof
plotPeakProf2 = plot_peak_prof2
plotPeakProf_MultiWindows = plot_peak_prof_multi_windows
plotMultiProf = plot_multi_prof
tagHeatmap = tag_heatmap
peakHeatmap = peak_heatmap
peakHeatmap_multiple_Sets = peak_heatmap_multiple_sets
peak_Profile_Heatmap = peak_profile_heatmap
