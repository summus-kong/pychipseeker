"""Shared bigWig interval integration and block-parallel extraction.

Both the TSS-profile helper in :mod:`pychipseeker.plotting` and the
bigWig-backed tag matrix in :mod:`pychipseeker.tagmatrix` read a bigWig track
once per chromosome and integrate the intervals at the bin edges.  This module
holds that shared math plus a small block-parallel runner.

``pyBigWig`` handles are not picklable, so every worker opens its own handle in
the pool initializer.  Work units are genomic blocks (a chromosome segment),
which keeps reads sequential on disk while still balancing dense and sparse
regions across cores.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

TSS_MODE = "tss"
WINDOW_MODE = "window"

VECTORISED_VALUE_TYPES = frozenset({"sum", "mean", "coverage"})


class IntervalIntegrator:
    """Cumulative length/value integrals over sorted bigWig intervals."""

    __slots__ = (
        "starts",
        "ends",
        "values",
        "lengths",
        "cum_len",
        "cum_val",
        "count",
    )

    def __init__(self, intervals: Sequence[tuple[int, int, float]]) -> None:
        table = np.asarray(intervals, dtype=np.float64)
        if table.ndim != 2 or table.shape[0] == 0 or table.shape[1] < 3:
            raise ValueError("intervals must be a non-empty sequence of records")
        order = np.argsort(table[:, 0], kind="mergesort")
        table = table[order]
        self.starts = table[:, 0].astype(np.int64)
        self.ends = table[:, 1].astype(np.int64)
        self.values = table[:, 2]
        self.lengths = (self.ends - self.starts).astype(np.float64)
        self.cum_len = np.concatenate([[0.0], np.cumsum(self.lengths)])
        self.cum_val = np.concatenate(
            [[0.0], np.cumsum(self.values * self.lengths)]
        )
        self.count = len(self.starts)

    def integrate(
        self, points: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return covered bases and integrated signal before each point."""

        index = np.searchsorted(self.ends, points, side="right")
        cursor = np.minimum(index, self.count - 1)
        part = np.where(
            index < self.count,
            np.clip(
                points - self.starts[cursor], 0.0, self.lengths[cursor]
            ),
            0.0,
        )
        return (
            self.cum_len[index] + part,
            self.cum_val[index] + part * self.values[cursor],
        )


def tss_matrix(
    integrator: IntervalIntegrator,
    tss: np.ndarray,
    up: int,
    down: int,
    n_bins: int,
    min_value: Optional[float] = None,
) -> np.ndarray:
    """Bin a TSS-centred matrix, keeping NaN for uncovered bins."""

    offsets = (np.arange(n_bins + 1, dtype=np.int64) * (int(up) + int(down))) // int(
        n_bins
    )
    edges = (np.asarray(tss, dtype=np.int64) - int(up))[:, None] + offsets[None, :]
    covered, integral = integrator.integrate(edges.ravel())
    covered = covered.reshape(-1, n_bins + 1)
    integral = integral.reshape(-1, n_bins + 1)
    bin_len = np.diff(covered, axis=1)
    bin_val = np.diff(integral, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        arr = np.where(bin_len > 0, bin_val / bin_len, np.nan)
    if min_value is not None:
        arr = np.where(arr < float(min_value), 0.0, arr)
    return arr


def window_matrix(
    integrator: IntervalIntegrator,
    starts: np.ndarray,
    ends: np.ndarray,
    n_bins: int,
    value_type: str,
) -> np.ndarray:
    """Bin variable-width windows using ChIPseeker's tag-matrix semantics."""

    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    widths = ends - starts
    offsets = (
        np.arange(n_bins + 1, dtype=np.int64)[None, :] * widths[:, None]
    ) // int(n_bins)
    edges = starts[:, None] + offsets
    covered, integral = integrator.integrate(edges.ravel())
    covered = covered.reshape(-1, n_bins + 1)
    integral = integral.reshape(-1, n_bins + 1)
    bin_len = np.diff(covered, axis=1)
    bin_val = np.diff(integral, axis=1)
    if value_type == "sum":
        return bin_val
    if value_type == "mean":
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(bin_len > 0, bin_val / bin_len, np.nan)
    if value_type == "coverage":
        bin_width = np.diff(edges, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(bin_width > 0, bin_len / bin_width, np.nan)
    raise ValueError(f"unsupported value_type for interval integration: {value_type}")


def stats_matrix(
    handle,
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    n_bins: int,
    value_type: str,
) -> np.ndarray:
    """Per-window ``bw.stats`` fallback for extrema value types."""

    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    out = np.zeros((len(starts), int(n_bins)), dtype=float)
    for i, (start, end) in enumerate(zip(starts, ends)):
        values = handle.stats(
            str(chrom),
            int(start),
            int(end),
            nBins=int(n_bins),
            type=str(value_type),
        )
        if values is None:
            continue
        out[i] = np.nan_to_num(
            np.asarray(
                [np.nan if v is None else float(v) for v in values],
                dtype=float,
            ),
            nan=0.0,
        )
    return out


@dataclass
class BlockJob:
    """Immutable configuration pickled once into every worker."""

    kind: str
    bigwig: str
    n_bins: int
    value_type: str = "mean"
    up: int = 0
    down: int = 0
    min_value: Optional[float] = None
    ignore_strand: bool = False
    dtype: str = "float64"
    tss: Optional[np.ndarray] = None
    starts: Optional[np.ndarray] = None
    ends: Optional[np.ndarray] = None
    strands: Optional[np.ndarray] = None


def plan_blocks(
    chroms: np.ndarray,
    order_key: np.ndarray,
    region_start: np.ndarray,
    region_end: np.ndarray,
    rows: np.ndarray,
    *,
    n_blocks: int,
    min_block: int = 256,
) -> list[tuple[str, np.ndarray, int, int]]:
    """Split ``rows`` into chromosome-local genomic blocks.

    Blocks are ordered from the densest chromosome to the sparsest so the
    process pool starts the largest units first and finishes more evenly.
    """

    chroms = np.asarray(chroms)
    order_key = np.asarray(order_key)
    region_start = np.asarray(region_start, dtype=np.int64)
    region_end = np.asarray(region_end, dtype=np.int64)
    rows = np.asarray(rows, dtype=np.int64)
    if rows.size == 0:
        return []

    names, counts = np.unique(chroms[rows], return_counts=True)
    dense_first = np.argsort(-counts, kind="mergesort")
    per_block = max(int(min_block), -(-int(rows.size) // max(1, int(n_blocks))))

    tasks: list[tuple[str, np.ndarray, int, int]] = []
    for k in dense_first:
        chrom = str(names[k])
        block_rows = rows[chroms[rows] == chrom]
        block_rows = block_rows[
            np.argsort(order_key[block_rows], kind="mergesort")
        ]
        for start in range(0, len(block_rows), per_block):
            block = block_rows[start : start + per_block]
            tasks.append(
                (
                    chrom,
                    block,
                    int(region_start[block].min()),
                    int(region_end[block].max()),
                )
            )
    return tasks


_WORKER_JOB: Optional[BlockJob] = None
_WORKER_HANDLE = None


def _worker_init(job: BlockJob) -> None:
    global _WORKER_JOB, _WORKER_HANDLE
    import pyBigWig

    _WORKER_JOB = job
    _WORKER_HANDLE = pyBigWig.open(job.bigwig)


def _worker_close() -> None:
    global _WORKER_HANDLE
    if _WORKER_HANDLE is not None:
        try:
            _WORKER_HANDLE.close()
        finally:
            _WORKER_HANDLE = None


def _worker_run(
    task: tuple[str, np.ndarray, int, int],
) -> tuple[np.ndarray, Optional[np.ndarray]]:
    """Extract one genomic block; returns ``(rows, matrix_or_None)``."""

    if _WORKER_JOB is None or _WORKER_HANDLE is None:  # pragma: no cover - guard
        raise RuntimeError("bigWig worker was not initialised")

    job = _WORKER_JOB
    chrom, rows, region_start, region_end = task
    intervals = _WORKER_HANDLE.intervals(
        str(chrom), int(region_start), int(region_end)
    )
    if not intervals:
        return rows, None

    integrator = IntervalIntegrator(intervals)
    if job.kind == TSS_MODE:
        if job.tss is None or job.strands is None:  # pragma: no cover - guard
            raise RuntimeError("TSS blocks require tss/strands arrays")
        arr = tss_matrix(
            integrator,
            job.tss[rows],
            job.up,
            job.down,
            job.n_bins,
            job.min_value,
        )
        minus = job.strands[rows] == "-"
        if minus.any():
            arr[minus] = arr[minus][:, ::-1]
        return rows, arr.astype(job.dtype, copy=False)

    if job.starts is None or job.ends is None or job.strands is None:
        raise RuntimeError("window blocks require starts/ends/strands arrays")
    if str(job.value_type) in VECTORISED_VALUE_TYPES:
        arr = window_matrix(
            integrator,
            job.starts[rows],
            job.ends[rows],
            job.n_bins,
            str(job.value_type),
        )
        arr = np.nan_to_num(arr, nan=0.0)
    else:
        arr = stats_matrix(
            _WORKER_HANDLE,
            chrom,
            job.starts[rows],
            job.ends[rows],
            job.n_bins,
            str(job.value_type),
        )
    # Strand orientation is applied by the caller once the full matrix has
    # been assembled, so blocks must stay in genomic order here.
    return rows, arr.astype(job.dtype, copy=False)


def run_blocks(
    job: BlockJob,
    tasks: Sequence[tuple[str, np.ndarray, int, int]],
    n_jobs: Optional[int],
) -> list[tuple[np.ndarray, Optional[np.ndarray]]]:
    """Run block extraction sequentially or in a process pool."""

    from .parallel import effective_jobs, parallel_map

    workers = effective_jobs(n_jobs, tasks=len(tasks))
    if workers <= 1 or len(tasks) <= 1:
        _worker_init(job)
        try:
            return [_worker_run(task) for task in tasks]
        finally:
            _worker_close()
    return parallel_map(
        _worker_run,
        list(tasks),
        n_jobs=workers,
        initializer=_worker_init,
        initargs=(job,),
    )
