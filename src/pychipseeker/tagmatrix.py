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
from typing import List, Optional, Sequence, Union, cast

import numpy as np
import pandas as pd

from ._bigwig import (
    VECTORISED_VALUE_TYPES,
    WINDOW_MODE,
    BlockJob,
    IntervalIntegrator,
    plan_blocks,
    run_blocks,
    stats_matrix,
    window_matrix,
)
from .annotation import GeneAnnotation, _normalise_peaks
from .genomes import is_assembly_name
from .parallel import (
    chunk_bounds,
    effective_jobs,
    parallel_map,
    process_start_method,
)
from .seqlevels import harmonize_frame


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
    if isinstance(source, str) and is_assembly_name(source):
        from .genomes import load_genome  # deferred: keep import graph flat

        return load_genome(source)
    from .genomes import is_genome_spec, load_annotation

    if isinstance(source, (str, Path)) and is_genome_spec(str(source)):
        return load_annotation(str(source))
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


def _peak_arrays(
    peaks: pd.DataFrame,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Return ``(starts, ends, weights)`` arrays per chromosome.

    The arrays are left in input order: the coverage builder turns them into a
    piecewise-constant track and does not need pre-sorted intervals.
    """

    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, group in peaks.groupby("chrom", sort=False):
        starts = group["start"].to_numpy(dtype=np.int64)
        ends = group["end"].to_numpy(dtype=np.int64)
        weights = (
            group["_weight"].to_numpy(dtype=float)
            if "_weight" in group.columns
            else np.ones(len(group), dtype=float)
        )
        out[str(chrom)] = (starts, ends, weights)
    return out


_CoverageTrack = tuple[np.ndarray, np.ndarray, np.ndarray]


def _coverage_track(
    starts: np.ndarray,
    ends: np.ndarray,
    weights: np.ndarray,
) -> Optional[_CoverageTrack]:
    """Turn peak intervals into a coverage track and its running integral.

    ``levels[k]`` is the summed peak weight on ``[coords[k], coords[k + 1])``
    and ``integral[k]`` is the integral of the coverage from ``coords[0]`` to
    ``coords[k]``.  A bin's coverage is then a difference of two integrals, so
    every window can be evaluated with a few vectorised lookups instead of a
    Python loop over overlapping peaks and bins.
    """

    valid = ends > starts
    if not np.any(valid):
        return None
    starts = starts[valid]
    ends = ends[valid]
    weights = weights[valid]
    if starts.size == 0:
        return None
    coords, inverse = np.unique(
        np.concatenate([starts, ends]), return_inverse=True
    )
    if coords.size < 2:
        return None
    changes = np.concatenate([weights, -weights])
    deltas = np.bincount(
        inverse, weights=changes, minlength=coords.size
    ).astype(float, copy=False)
    # Track how many peaks are active in each segment so rounding residues in
    # the cumulative weight cannot leak phantom coverage into peak-free gaps.
    active_changes = np.concatenate(
        [
            np.ones(starts.size, dtype=np.int64),
            -np.ones(starts.size, dtype=np.int64),
        ]
    )
    active = np.cumsum(
        np.bincount(inverse, weights=active_changes, minlength=coords.size).astype(
            np.int64
        )
    )[:-1]
    levels = np.cumsum(deltas)[:-1]
    levels[active == 0] = 0.0
    integral = np.zeros(coords.size, dtype=float)
    integral[1:] = np.cumsum(levels * np.diff(coords))
    return coords, levels, integral


def _integral_at(track: _CoverageTrack, queries: np.ndarray) -> np.ndarray:
    """Evaluate the running coverage integral at arbitrary genomic positions."""

    coords, levels, integral = track
    position = np.asarray(queries, dtype=float)
    idx = np.searchsorted(coords, position, side="right") - 1
    np.clip(idx, 0, coords.size - 2, out=idx)
    values = integral[idx] + levels[idx] * (position - coords[idx])
    values[position <= coords[0]] = 0.0
    values[position >= coords[-1]] = integral[-1]
    return values


# Cap the number of coverage values computed in one vectorised batch so that
# per-base windows (nbin=None) do not allocate an enormous edge matrix.
_COVERAGE_BLOCK_ELEMENTS = 2_000_000


def _coverage_rows(
    by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Build peak-coverage rows for one block of windows.

    Each chromosome is reduced once to a piecewise-constant coverage track;
    windows are then evaluated as differences of the track's cumulative
    integral.  This keeps the per-window cost proportional to the number of
    bins and independent of the number of peaks.
    """

    empty_width = int(n_bins) if n_bins is not None else 0
    empty = np.empty((0, empty_width), dtype=float)
    if windows.empty:
        return empty, pd.DataFrame(columns=windows.columns)

    chroms = windows["chrom"].astype(str).to_numpy()
    starts = windows["start"].to_numpy(dtype=np.int64)
    ends = windows["end"].to_numpy(dtype=np.int64)
    strands = windows["strand"].astype(str).to_numpy()
    widths = ends - starts
    positive = widths > 0

    if n_bins is None:
        if np.any(widths[positive] <= 1):
            raise ValueError("nbin must be > 1")
        bins_per_row = widths
    else:
        if int(n_bins) <= 1 and bool(np.any(positive)):
            raise ValueError("nbin must be > 1")
        bins_per_row = np.full(len(windows), int(n_bins), dtype=np.int64)

    valid = positive & (bins_per_row > 1)

    # Group by chromosome *and* bin count: one integral lookup serves a whole
    # batch of windows with the same shape.
    groups: dict[tuple[str, int], list[int]] = {}
    for index in np.flatnonzero(valid):
        key = (chroms[index], int(bins_per_row[index]))
        groups.setdefault(key, []).append(int(index))

    tracks: dict[str, Optional[_CoverageTrack]] = {}
    blocks: list[np.ndarray] = []
    kept_indices: list[np.ndarray] = []
    for (chrom, bins), indices in groups.items():
        if chrom not in by_chrom:
            continue
        if chrom not in tracks:
            tracks[chrom] = _coverage_track(*by_chrom[chrom])
        track = tracks[chrom]
        if track is None:
            continue

        index_array = np.asarray(indices, dtype=np.int64)
        # Keep the edge batch at a predictable size; per-base windows can span
        # tens of kilobases, so the chunk size shrinks with the bin count.
        chunk = max(1, min(4096, _COVERAGE_BLOCK_ELEMENTS // (bins + 1)))
        fractions = np.linspace(0.0, 1.0, bins + 1)
        for offset in range(0, index_array.size, chunk):
            block = index_array[offset : offset + chunk]
            left = starts[block].astype(float)[:, None]
            right = ends[block].astype(float)[:, None]
            edges = left + (right - left) * fractions
            values = np.diff(_integral_at(track, edges), axis=1)
            keep = np.any(values != 0.0, axis=1)
            if not np.any(keep):
                continue
            blocks.append(values[keep])
            kept_indices.append(block[keep])

    if not blocks:
        return empty, pd.DataFrame(columns=windows.columns)

    matrix = np.vstack(blocks)
    order = np.concatenate(kept_indices)
    # Groups are visited per chromosome; restore the caller's window order.
    permutation = np.argsort(order, kind="stable")
    matrix = matrix[permutation]
    order = order[permutation]

    if not ignore_strand:
        minus = strands[order] == "-"
        if np.any(minus):
            matrix[minus] = matrix[minus][:, ::-1]
    return matrix, windows.iloc[order].reset_index(drop=True)


_COVERAGE_MIN_WINDOWS = 50_000
# Spawning workers on Windows/macOS costs roughly a second per process, so a
# parallel coverage matrix only pays off for very large window sets there.
_COVERAGE_SPAWN_MIN_WINDOWS = 400_000


def _coverage_parallel_floor() -> int:
    if process_start_method() == "spawn":
        return _COVERAGE_SPAWN_MIN_WINDOWS
    return _COVERAGE_MIN_WINDOWS


_COVERAGE_WORKER: Optional[tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], Optional[int], bool]] = None


def _coverage_worker_init(
    by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    n_bins: Optional[int],
    ignore_strand: bool,
) -> None:
    global _COVERAGE_WORKER
    _COVERAGE_WORKER = (by_chrom, n_bins, bool(ignore_strand))


def _coverage_worker(windows: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    if _COVERAGE_WORKER is None:  # pragma: no cover - pool guard
        raise RuntimeError("coverage matrix worker was not initialised")
    by_chrom, n_bins, ignore_strand = _COVERAGE_WORKER
    return _coverage_rows(
        by_chrom,
        windows,
        n_bins=n_bins,
        ignore_strand=ignore_strand,
    )


def _coverage_matrix(
    peaks: pd.DataFrame,
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
    n_jobs: Optional[int] = None,
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
    workers = effective_jobs(n_jobs, tasks=len(windows))
    if (
        workers > 1
        and len(windows) >= _coverage_parallel_floor() * workers
    ):
        chunks = [
            windows.iloc[start:stop].reset_index(drop=True)
            for start, stop in chunk_bounds(len(windows), workers * 4)
        ]
        results = parallel_map(
            _coverage_worker,
            chunks,
            n_jobs=workers,
            initializer=_coverage_worker_init,
            initargs=(by_chrom, n_bins, bool(ignore_strand)),
        )
        matrices = [block for block, _ in results if block.size]
        kept_frames = [frame for _, frame in results if len(frame)]
        if not matrices:
            raise ValueError("no windows overlap any peak")
        matrix = np.vstack(matrices)
        row_frame = pd.concat(kept_frames, ignore_index=True)
        return matrix, row_frame, n_bins is not None

    matrix, row_frame = _coverage_rows(
        by_chrom,
        windows,
        n_bins=n_bins,
        ignore_strand=ignore_strand,
    )
    if matrix.size == 0:
        raise ValueError("no windows overlap any peak")
    return matrix, row_frame, n_bins is not None


_BIGWIG_MATRIX_MIN_WINDOWS = 4000


def _bigwig_matrix(
    bigwig: Union[str, Path],
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
    value_type: str = "sum",
    n_jobs: Optional[int] = None,
) -> tuple[np.ndarray, pd.DataFrame, bool]:
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - dependency installation issue
        raise ImportError(
            "pyBigWig is required for bigWig-backed tag matrices; install with "
            "`pip install pychipseeker[bigwig]` (on Windows use "
            "`conda install -c conda-forge pybigwig`)"
        ) from exc

    if str(value_type) not in {"sum", "mean", "max", "min", "coverage", "std"}:
        raise ValueError(
            "unsupported bigWig value_type; choose sum, mean, max, min, "
            "coverage or std"
        )

    bw = pyBigWig.open(str(bigwig))
    if bw is None:
        raise ValueError(f"could not open bigWig: {bigwig}")
    try:
        chrom_sizes = bw.chroms()
        windows = windows.reset_index(drop=True)
        if windows.empty:
            raise ValueError("windows is empty")

        # Match Ensembl-style target names (1/MT) to UCSC query names
        # (chr1/chrM), and vice versa.
        windows, _chrom_map = harmonize_frame(
            windows, chrom_sizes.keys(), column="chrom", copy=True
        )
        starts = pd.to_numeric(windows["start"], errors="coerce").to_numpy(
            dtype=np.int64
        )
        ends = pd.to_numeric(windows["end"], errors="coerce").to_numpy(
            dtype=np.int64
        )
        chroms = windows["chrom"].astype(str).to_numpy()
        strands = windows["strand"].astype(str).to_numpy()
        chrom_lengths = np.array(
            [chrom_sizes.get(chrom, -1) for chrom in chroms], dtype=np.int64
        )
        widths = ends - starts
        is_binning = n_bins is not None
        bins_per_row = np.full(len(windows), -1, dtype=np.int64)
        if n_bins is None:
            if len(np.unique(widths)) != 1:
                raise ValueError(
                    "windows must have equal widths when nbin is not supplied; "
                    "use nbin=... for variable-width regions"
                )
            bins_per_row[:] = widths
        else:
            bins_per_row[:] = int(n_bins)
        if np.any(bins_per_row[widths > 0] <= 1):
            raise ValueError("nbin must be > 1")

        valid = (
            (chrom_lengths > 0)
            & (starts >= 0)
            & (ends <= chrom_lengths)
            & (ends > starts)
            & (bins_per_row > 1)
        )
        rows_idx = np.flatnonzero(valid)
        if rows_idx.size == 0:
            raise ValueError("no windows overlap the bigWig track")
        unique_bins = np.unique(bins_per_row[rows_idx])
        if unique_bins.size != 1:
            raise ValueError(
                "windows must have equal widths when nbin is not supplied"
            )
        n_bins_eff = int(unique_bins[0])
        matrix = np.zeros((rows_idx.size, n_bins_eff), dtype=float)

        intervals_fn = getattr(bw, "intervals", None)
        workers = effective_jobs(n_jobs, tasks=len(rows_idx))
        if (
            workers > 1
            and len(rows_idx) >= _BIGWIG_MATRIX_MIN_WINDOWS * workers
            and intervals_fn is not None
        ):
            tasks = plan_blocks(
                chroms,
                starts,
                starts,
                ends,
                rows_idx,
                n_blocks=workers * 8,
            )
            job = BlockJob(
                kind=WINDOW_MODE,
                bigwig=str(bigwig),
                n_bins=n_bins_eff,
                value_type=str(value_type),
                ignore_strand=ignore_strand,
                starts=starts,
                ends=ends,
                strands=strands,
            )
            for block_rows, block in run_blocks(job, tasks, workers):
                if block is not None:
                    matrix[np.searchsorted(rows_idx, block_rows)] = block
        else:
            vectorised = intervals_fn is not None and str(
                value_type
            ) in VECTORISED_VALUE_TYPES
            if vectorised:
                # Read each chromosome once and integrate the intervals at the
                # bin edges.  This replaces one bw.stats() call per window,
                # which dominates the runtime for transcript-scale matrices.
                for chrom in pd.unique(chroms[rows_idx]):
                    index = rows_idx[chroms[rows_idx] == chrom]
                    if intervals_fn is None:
                        continue
                    intervals = intervals_fn(str(chrom))
                    if not intervals:
                        continue
                    integrator = IntervalIntegrator(intervals)
                    arr = window_matrix(
                        integrator,
                        starts[index],
                        ends[index],
                        n_bins_eff,
                        str(value_type),
                    )
                    matrix[np.searchsorted(rows_idx, index)] = np.nan_to_num(
                        arr, nan=0.0
                    )
            else:
                # max/min/std need per-bin extrema and are not expressible as
                # a simple integral.  Keep the exact bw.stats fallback for
                # them, and for minimal pyBigWig stand-ins without intervals().
                for out_row, idx in enumerate(rows_idx):
                    matrix[out_row] = stats_matrix(
                        bw,
                        str(chroms[idx]),
                        starts[idx : idx + 1],
                        ends[idx : idx + 1],
                        n_bins_eff,
                        str(value_type),
                    )[0]

        if not ignore_strand:
            minus = strands[rows_idx] == "-"
            if minus.any():
                matrix[minus] = matrix[minus][:, ::-1]
        kept = windows.iloc[rows_idx].reset_index(drop=True)
        return matrix, kept, is_binning
    finally:
        bw.close()


def _bin_signal(values: np.ndarray, n_bins: Optional[int], value_type: str) -> np.ndarray:
    """Aggregate per-base coverage into ``n_bins`` bins."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("coverage array is empty")
    if n_bins is None or int(n_bins) == values.size:
        return values
    n_bins = int(n_bins)
    if n_bins <= 1:
        raise ValueError("nbin must be > 1")
    edges = (np.arange(n_bins + 1, dtype=np.int64) * values.size) // n_bins
    out = np.zeros(n_bins, dtype=float)
    for i in range(n_bins):
        chunk = values[edges[i] : edges[i + 1]]
        if chunk.size == 0:
            continue
        if value_type == "sum":
            out[i] = float(chunk.sum())
        elif value_type == "max":
            out[i] = float(chunk.max())
        elif value_type == "coverage":
            out[i] = float((chunk > 0).mean())
        else:  # mean
            out[i] = float(chunk.mean())
    return out


_BAM_MATRIX_MIN_WINDOWS = 500
_BAM_MATRIX_WORKER: Optional[tuple[str, Optional[int], bool, str, int, bool]] = None


def _bam_matrix_rows(
    bam: Union[str, Path],
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
    value_type: str,
    min_mapping_quality: int,
    exclude_duplicates: bool,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Extract coverage rows for one block of windows."""

    try:
        import pysam  # type: ignore
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "pysam is required for BAM-backed tag matrices; install with "
            "`pip install pychipseeker[bam]` (on Windows use "
            "`conda install -c bioconda pysam`)"
        ) from exc

    handle = pysam.AlignmentFile(str(bam), "rb")
    try:
        references = list(getattr(handle, "references", []) or [])
        windows = windows.reset_index(drop=True)
        if references:
            windows, _ = harmonize_frame(
                windows, references, column="chrom", copy=True
            )
        strands = windows["strand"].astype(str).to_numpy()
        rows: list[np.ndarray] = []
        kept: list[pd.Series] = []
        for row_idx, window in windows.iterrows():
            chrom = str(window["chrom"])
            start = int(window["start"])
            end = int(window["end"])
            if end <= start or start < 0 or chrom not in references:
                continue
            try:
                coverage = handle.count_coverage(
                    chrom,
                    start,
                    end,
                    quality_threshold=0,
                    read_callback=lambda read: not (
                        read.is_unmapped
                        or read.is_secondary
                        or read.is_supplementary
                        or (exclude_duplicates and read.is_duplicate)
                        or int(read.mapping_quality) < int(min_mapping_quality)
                    ),
                )
            except (ValueError, KeyError):
                continue
            if not coverage:
                continue
            values = np.sum(np.asarray(coverage, dtype=float), axis=0)
            if values.size == 0:
                continue
            row = _bin_signal(values, n_bins, value_type)
            if str(strands[row_idx]) == "-" and not ignore_strand:
                row = row[::-1]
            rows.append(row)
            kept.append(window)

        if not rows:
            width = int(n_bins) if n_bins is not None else 0
            return np.empty((0, width), dtype=float), pd.DataFrame(columns=windows.columns)
        return np.vstack(rows), pd.DataFrame(kept).reset_index(drop=True)
    finally:
        handle.close()


def _bam_matrix_worker_init(
    bam: Union[str, Path],
    n_bins: Optional[int],
    ignore_strand: bool,
    value_type: str,
    min_mapping_quality: int,
    exclude_duplicates: bool,
) -> None:
    global _BAM_MATRIX_WORKER
    _BAM_MATRIX_WORKER = (
        str(bam),
        n_bins,
        bool(ignore_strand),
        value_type,
        int(min_mapping_quality),
        bool(exclude_duplicates),
    )


def _bam_matrix_worker(
    windows: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame]:
    if _BAM_MATRIX_WORKER is None:  # pragma: no cover - pool guard
        raise RuntimeError("BAM matrix worker was not initialised")
    bam, n_bins, ignore_strand, value_type, min_mapq, exclude_duplicates = (
        _BAM_MATRIX_WORKER
    )
    return _bam_matrix_rows(
        bam,
        windows,
        n_bins=n_bins,
        ignore_strand=ignore_strand,
        value_type=value_type,
        min_mapping_quality=min_mapq,
        exclude_duplicates=exclude_duplicates,
    )


def _bam_matrix(
    bam: Union[str, Path],
    windows: pd.DataFrame,
    *,
    n_bins: Optional[int],
    ignore_strand: bool,
    value_type: str = "mean",
    min_mapping_quality: int = 0,
    exclude_duplicates: bool = True,
    n_jobs: Optional[int] = None,
) -> tuple[np.ndarray, pd.DataFrame, bool]:
    """Extract a read-coverage matrix from a BAM/CRAM file.

    Coverage is counted per base from alignments passing the mapping-quality
    and duplicate filters, then binned.  ``value_type`` accepts ``mean``
    (default), ``sum``, ``max`` or ``coverage``.  ``n_jobs`` splits the
    windows across processes; each worker reads its own chromosome regions
    from the indexed BAM.
    """

    value_type = str(value_type)
    if value_type not in {"mean", "sum", "max", "coverage"}:
        raise ValueError("BAM value_type must be mean, sum, max or coverage")

    # Fail early with the friendly optional-dependency message.
    try:
        import pysam  # type: ignore  # noqa: F401
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "pysam is required for BAM-backed tag matrices; install with "
            "`pip install pychipseeker[bam]` (on Windows use "
            "`conda install -c bioconda pysam`)"
        ) from exc

    windows = windows.reset_index(drop=True)
    workers = effective_jobs(n_jobs, tasks=len(windows))
    if workers <= 1 or len(windows) < _BAM_MATRIX_MIN_WINDOWS:
        matrix, kept = _bam_matrix_rows(
            bam,
            windows,
            n_bins=n_bins,
            ignore_strand=ignore_strand,
            value_type=value_type,
            min_mapping_quality=min_mapping_quality,
            exclude_duplicates=exclude_duplicates,
        )
    else:
        chunks = [
            windows.iloc[start:stop].reset_index(drop=True)
            for start, stop in chunk_bounds(len(windows), workers * 4)
        ]
        results = parallel_map(
            _bam_matrix_worker,
            chunks,
            n_jobs=workers,
            initializer=_bam_matrix_worker_init,
            initargs=(
                bam,
                n_bins,
                ignore_strand,
                value_type,
                min_mapping_quality,
                exclude_duplicates,
            ),
        )
        matrices = [chunk for chunk, _ in results if chunk.size]
        kept_frames = [frame for _, frame in results if len(frame)]
        if not matrices:
            raise ValueError("no windows overlap the BAM file")
        matrix = np.vstack(matrices)
        kept = pd.concat(kept_frames, ignore_index=True)

    if matrix.size == 0:
        raise ValueError("no windows overlap the BAM file")
    lengths = {row.size for row in matrix}
    if len(lengths) != 1:
        raise ValueError(
            "windows must have equal widths when nbin is not supplied"
        )
    return matrix, kept.reset_index(drop=True), n_bins is not None


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
    bam: Optional[Union[str, Path]] = None,
    value_type: str = "sum",
    min_mapping_quality: int = 0,
    exclude_duplicates: bool = True,
    n_jobs: Optional[int] = None,
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

    ``n_jobs`` splits bigWig chromosome blocks or BAM windows across processes;
    peak-coverage matrices are vectorised per chromosome and only fall back to
    processes for very large window sets.  ``-1`` uses every core and results
    match the sequential run.
    """

    if bigwig is None and isinstance(peak, (str, Path)):
        suffix = str(peak).lower()
        if suffix.endswith((".bw", ".bigwig")):
            bigwig = peak
            peak = None
        elif suffix.endswith((".bam", ".cram")):
            bam = peak
            peak = None

    if isinstance(windows, (list, tuple)):
        if not windows:
            raise ValueError("windows is empty")
        matrices: List[TagMatrix] = []
        for i, one_window in enumerate(windows):
            tag_matrix = cast(
                TagMatrix,
                get_tag_matrix(
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
                    bam=bam,
                    value_type=value_type,
                    min_mapping_quality=min_mapping_quality,
                    exclude_duplicates=exclude_duplicates,
                    n_jobs=n_jobs,
                ),
            )
            tag_matrix.label = (
                getattr(one_window, "attrs", {}).get("label")
                or getattr(one_window, "attrs", {}).get("type")
                or f"Window{i + 1}"
            )
            matrices.append(tag_matrix)
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

    if bigwig is not None and bam is not None:
        raise ValueError("supply only one of bigwig=... or bam=...")
    if bigwig is not None:
        matrix, kept, is_binning = _bigwig_matrix(
            bigwig,
            windows_df,
            n_bins=nbin,
            ignore_strand=ignore_strand,
            value_type=value_type,
            n_jobs=n_jobs,
        )
    elif bam is not None:
        matrix, kept, is_binning = _bam_matrix(
            bam,
            windows_df,
            n_bins=nbin,
            ignore_strand=ignore_strand,
            value_type=value_type,
            min_mapping_quality=min_mapping_quality,
            exclude_duplicates=exclude_duplicates,
            n_jobs=n_jobs,
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
            n_jobs=n_jobs,
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
    n_jobs: Optional[int] = None,
) -> TagMatrix:
    """Convenience wrapper for a bigWig-backed :class:`TagMatrix`."""

    windows_df = _window_frame(windows)
    matrix, kept, is_binning = _bigwig_matrix(
        bigwig,
        windows_df,
        n_bins=nbin,
        ignore_strand=ignore_strand,
        value_type=value_type,
        n_jobs=n_jobs,
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


def get_tag_matrix_from_bam(
    bam: Union[str, Path],
    windows: Union[pd.DataFrame, str, Path],
    nbin: Optional[int] = None,
    ignore_strand: bool = False,
    value_type: str = "mean",
    min_mapping_quality: int = 0,
    exclude_duplicates: bool = True,
    n_jobs: Optional[int] = None,
) -> TagMatrix:
    """Convenience wrapper for a BAM-backed :class:`TagMatrix`."""

    windows_df = _window_frame(windows)
    matrix, kept, is_binning = _bam_matrix(
        bam,
        windows_df,
        n_bins=nbin,
        ignore_strand=ignore_strand,
        value_type=value_type,
        min_mapping_quality=min_mapping_quality,
        exclude_duplicates=exclude_duplicates,
        n_jobs=n_jobs,
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
        files("pychipseeker").joinpath("data").joinpath("ENCODE_demo_signal.bw")
    ) as path:
        return Path(path)


def example_signal_provenance() -> dict:
    """Provenance record for :func:`example_signal_path`."""
    import json
    from importlib.resources import as_file, files

    with as_file(
        files("pychipseeker").joinpath("data").joinpath("ENCODE_demo_signal.json")
    ) as path:
        with open(path, encoding="utf-8") as fh:
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

    from .plotting import _style_axes as _gg_axes  # deferred: avoid import cycles

    line = ax.plot(
        x,
        mean,
        color=color or "#E45756",
        linewidth=1.5,
        label=label,
        **kwargs,
    )[0]

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
                alpha=0.2,
                linewidth=0,
            )

    ax.axvline(0, color="#555555", linewidth=0.8, linestyle="--")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if label:
        ax.legend(frameon=False)
    ax.set_title(title)
    _gg_axes(ax, grid_axis="both")
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
    """Heatmap of a tag matrix, centered on the TSS/region anchor.

    Rows are sorted by mean signal so the most enriched regions are drawn at
    the top of the map.
    """

    plt = _plt()
    arr, upstream, downstream = _as_matrix(tagMatrix)
    if arr.shape[0] == 0:
        raise ValueError("tagMatrix has no rows")
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Row-sort by signal strength, strongest first.  This matches
    # ``plot_tss_heatmap(sort_rows=True)`` and the convention used by genome
    # browsers and deepTools, where the most enriched rows sit at the top.
    # Sorting also makes the visual pattern stable across input orders.
    order = np.argsort(-np.nanmean(arr, axis=1), kind="stable")
    data = arr[order]
    extent = (
        -int(upstream) if upstream is not None else 0,
        int(downstream) if downstream is not None else data.shape[1],
        data.shape[0],
        0,
    )
    from matplotlib import colormaps

    cmap_obj = colormaps[(kwargs.pop("cmap", None) or palette or "Reds")].with_extremes(
        bad="#EBEBEB"
    )
    im = ax.imshow(
        np.ma.masked_invalid(data),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap_obj,
        extent=extent,
        **kwargs,
    )
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab or "Peaks")
    ax.set_title(title or "Tag heatmap")
    fig.colorbar(im, ax=ax, shrink=0.8)
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

    matrices: List[TagMatrix] = []
    for i, window in enumerate(windows):
        tag_matrix = cast(
            TagMatrix,
            get_tag_matrix(
                peak=peak,
                windows=window,
                weightCol=weightCol,
                nbin=nbin,
                ignore_strand=ignore_strand,
            ),
        )
        tag_matrix.label = (
            windows_name[i]
            if windows_name is not None and i < len(windows_name)
            else getattr(window, "attrs", {}).get("label", f"Window{i + 1}")
        )
        matrices.append(tag_matrix)
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
