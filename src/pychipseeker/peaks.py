"""Peak calling helpers.

The threshold caller turns a signal bigWig (for example an ENCODE
fold-change track) into BED-like peaks without requiring MACS2.  It is meant
for cases where only a signal track is available; when a BAM file and MACS2/3
are installed, :func:`call_peaks_macs` provides the conventional caller.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd

from .io import read_peaks, write_table
from .parallel import effective_jobs, parallel_map


def _open_bigwig(path: Union[str, Path]):
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "pyBigWig is required to call peaks from a bigWig; install with "
            "`pip install pychipseeker[bigwig]` (on Windows use "
            "`conda install -c conda-forge pybigwig`)"
        ) from exc
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise ValueError(f"could not open bigWig: {path}")
    return handle


def _call_peaks_chrom(
    bw,
    chrom: str,
    *,
    threshold: float,
    max_gap: int,
    min_width: int,
) -> list[dict[str, object]]:
    """Call peaks for one chromosome from an open bigWig handle."""

    intervals = bw.intervals(chrom)
    if intervals is None or len(intervals) == 0:
        return []

    table = np.asarray(intervals, dtype=np.float64)
    if table.ndim != 2 or table.shape[0] == 0:
        return []
    starts = table[:, 0].astype(np.int64)
    ends = table[:, 1].astype(np.int64)
    values = table[:, 2]
    keep = np.isfinite(values) & (values >= threshold) & (ends > starts)
    if not keep.any():
        return []
    starts = starts[keep]
    ends = ends[keep]
    values = values[keep]
    order = np.argsort(starts, kind="mergesort")
    starts = starts[order]
    ends = ends[order]
    values = values[order]

    # A new run starts when the next interval is farther than max_gap from the
    # maximum end seen so far.  Using the cumulative maximum also handles
    # nested/overlapping bigWig intervals.
    cumulative_end = np.maximum.accumulate(ends)
    breaks = np.flatnonzero(starts[1:] > cumulative_end[:-1] + int(max_gap)) + 1
    group_starts = np.concatenate(([0], breaks))
    group_ends = np.concatenate((breaks, [len(starts)]))
    run_start = starts[group_starts]
    run_end = np.maximum.reduceat(ends, group_starts)
    run_sum = np.add.reduceat(values * (ends - starts), group_starts)
    run_max = np.maximum.reduceat(values, group_starts)
    widths = run_end - run_start

    rows: list[dict[str, object]] = []
    for group, (a, b) in enumerate(zip(group_starts, group_ends)):
        width = int(widths[group])
        if width < min_width:
            continue
        local = int(np.argmax(values[a:b]))
        run_peak = int(starts[a + local])
        rows.append(
            {
                "chrom": str(chrom),
                "start": int(run_start[group]),
                "end": int(run_end[group]),
                "score": round(float(run_max[group]), 6),
                "strand": ".",
                "signalValue": float(run_sum[group]) / max(width, 1),
                "peak": run_peak - int(run_start[group]),
            }
        )
    return rows


_CALLPEAK_WORKER: Optional[tuple[str, float, int, int]] = None


def _callpeak_worker_init(
    bigwig: Union[str, Path], threshold: float, max_gap: int, min_width: int
) -> None:
    global _CALLPEAK_WORKER
    _CALLPEAK_WORKER = (str(bigwig), float(threshold), int(max_gap), int(min_width))


def _callpeak_worker(chrom: str) -> list[dict[str, object]]:
    if _CALLPEAK_WORKER is None:  # pragma: no cover - pool initialiser guard
        raise RuntimeError("peak-calling worker was not initialised")
    bigwig, threshold, max_gap, min_width = _CALLPEAK_WORKER
    handle = _open_bigwig(bigwig)
    try:
        return _call_peaks_chrom(
            handle,
            chrom,
            threshold=threshold,
            max_gap=max_gap,
            min_width=min_width,
        )
    finally:
        handle.close()


def call_peaks(
    bigwig: Union[str, Path],
    *,
    threshold: float = 3.0,
    max_gap: int = 1000,
    min_width: int = 500,
    chroms: Optional[Sequence[str]] = None,
    sample: Optional[int] = None,
    seed: Optional[int] = None,
    n_jobs: Optional[int] = None,
) -> pd.DataFrame:
    """Call peaks as runs of signal above ``threshold``.

    Consecutive intervals are merged when the gap between them is at most
    ``max_gap`` bases; merged runs shorter than ``min_width`` are discarded.
    Returns a BED6+ table with ``signalValue`` (mean signal over the peak) and
    ``peak`` (offset of the highest base from the peak start).

    ``n_jobs`` calls chromosomes in parallel processes; ``-1`` uses every
    available core.  Each worker opens the bigWig read-only, so memory stays
    flat even for whole-genome tracks.
    """

    threshold = float(threshold)
    max_gap = int(max_gap)
    min_width = int(min_width)
    if max_gap < 0 or min_width < 0:
        raise ValueError("max_gap and min_width must be non-negative")

    probe = _open_bigwig(bigwig)
    try:
        available = list(probe.chroms().keys())
    finally:
        probe.close()
    selected = available if chroms is None else [str(c) for c in chroms]

    workers = effective_jobs(n_jobs, tasks=len(selected))
    if workers > 1 and len(selected) > 1:
        chunks = parallel_map(
            _callpeak_worker,
            selected,
            n_jobs=workers,
            initializer=_callpeak_worker_init,
            initargs=(bigwig, threshold, max_gap, min_width),
        )
        rows = [row for chunk in chunks for row in chunk]
    else:
        bw = _open_bigwig(bigwig)
        try:
            rows = [
                row
                for chrom in selected
                for row in _call_peaks_chrom(
                    bw,
                    chrom,
                    threshold=threshold,
                    max_gap=max_gap,
                    min_width=min_width,
                )
            ]
        finally:
            bw.close()

    if not rows:
        raise ValueError(
            "no peaks passed the threshold; lower --threshold or --min-width"
        )
    peaks = pd.DataFrame.from_records(rows)
    if sample is not None and len(peaks) > int(sample):
        peaks = peaks.sample(n=int(sample), random_state=seed).sort_index()
        peaks = peaks.reset_index(drop=True)
    peaks["name"] = [f"Peak_{i + 1}" for i in range(len(peaks))]
    return peaks.reset_index(drop=True)


def _overlaps_any(
    peaks: pd.DataFrame,
    blacklist: pd.DataFrame,
) -> np.ndarray:
    """Boolean mask of peaks overlapping any blacklist interval."""

    result = np.zeros(len(peaks), dtype=bool)
    required = {"chrom", "start", "end"}
    if not required.issubset(blacklist.columns):
        raise ValueError("blacklist must contain chrom/start/end columns")
    bl = blacklist.copy()
    bl["start"] = pd.to_numeric(bl["start"], errors="coerce")
    bl["end"] = pd.to_numeric(bl["end"], errors="coerce")
    bl = bl.dropna(subset=["start", "end"])
    bl["chrom"] = bl["chrom"].astype(str)
    for chrom, group in peaks.groupby("chrom", sort=False):
        sub = bl[bl["chrom"] == str(chrom)].sort_values("start")
        if sub.empty:
            continue
        b_start = sub["start"].to_numpy(dtype=np.int64)
        b_end = sub["end"].to_numpy(dtype=np.int64)
        cumulative_end = np.maximum.accumulate(b_end)
        p_start = group["start"].to_numpy(dtype=np.int64)
        p_end = group["end"].to_numpy(dtype=np.int64)
        lo = np.searchsorted(cumulative_end, p_start, side="right")
        hi = np.searchsorted(b_start, p_end, side="left")
        overlap = hi > lo
        result[group.index.to_numpy()] = overlap
    return result


def filter_peaks(
    peaks: Union[str, Path, pd.DataFrame],
    *,
    blacklist: Optional[Union[str, Path, pd.DataFrame]] = None,
    drop_mito: bool = False,
    min_width: Optional[int] = None,
    max_width: Optional[int] = None,
    min_score: Optional[float] = None,
    chroms: Optional[Sequence[str]] = None,
    exclude_chroms: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Filter a peak table by blacklist, chromosome, width and score.

    Blacklist intervals are removed when they overlap a peak by at least one
    base.  Mitochondrial chromosomes are recognised in both UCSC (``chrM``)
    and Ensembl (``MT``) spellings.
    """

    if isinstance(peaks, (str, Path)):
        frame = read_peaks(peaks)
    elif isinstance(peaks, pd.DataFrame):
        frame = peaks.copy()
    else:
        raise TypeError("peaks must be a path or a pandas.DataFrame")
    required = {"chrom", "start", "end"}
    if not required.issubset(frame.columns):
        raise ValueError(f"peak table is missing columns: {sorted(required - set(frame.columns))}")

    frame = frame.copy()
    frame["chrom"] = frame["chrom"].astype(str)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame.dropna(subset=["start", "end"])
    frame = frame[frame["end"] > frame["start"]].copy()
    frame["start"] = frame["start"].astype(np.int64)
    frame["end"] = frame["end"].astype(np.int64)

    if blacklist is not None:
        if isinstance(blacklist, (str, Path)):
            bl = read_peaks(blacklist)
        elif isinstance(blacklist, pd.DataFrame):
            bl = blacklist
        else:
            raise TypeError("blacklist must be a path or a pandas.DataFrame")
        from .seqlevels import harmonize_frame

        bl, _ = harmonize_frame(
            bl, frame["chrom"].unique(), column="chrom", copy=True
        )
        frame = frame.loc[~_overlaps_any(frame, bl)].copy()

    if drop_mito:
        from .seqlevels import canonical_chrom

        frame = frame[
            frame["chrom"].map(canonical_chrom) != "m"
        ].copy()

    if chroms is not None:
        from .seqlevels import canonical_chrom

        keep_chroms = {canonical_chrom(value) for value in chroms}
        frame = frame[frame["chrom"].map(canonical_chrom).isin(keep_chroms)].copy()
    if exclude_chroms is not None:
        from .seqlevels import canonical_chrom

        drop_chroms = {canonical_chrom(value) for value in exclude_chroms}
        frame = frame[~frame["chrom"].map(canonical_chrom).isin(drop_chroms)].copy()

    widths = frame["end"] - frame["start"]
    if min_width is not None:
        frame = frame[widths >= int(min_width)].copy()
        widths = frame["end"] - frame["start"]
    if max_width is not None:
        frame = frame[widths <= int(max_width)].copy()
    if min_score is not None:
        if "score" not in frame.columns:
            raise ValueError("min_score requires a score column")
        scores = pd.to_numeric(frame["score"], errors="coerce")
        frame = frame[scores >= float(min_score)].copy()

    return frame.reset_index(drop=True)


def call_peaks_macs(
    treatment: Union[str, Path],
    *,
    control: Optional[Union[str, Path]] = None,
    genome: Optional[str] = None,
    outdir: Union[str, Path] = "macs_output",
    backend: str = "macs2",
    extra_args: Iterable[str] = (),
    name: str = "pychipseeker",
) -> pd.DataFrame:
    """Run MACS2/MACS3 (if installed) and return the called narrowPeak table.

    ``backend`` is ``"macs2"`` or ``"macs3"``.  A genome size (``hs``, ``mm``,
    ``1.0e9``, ...) is passed through to ``-g`` when supplied.
    """

    executable = shutil.which(backend)
    if executable is None:
        raise ImportError(
            f"{backend} was not found on PATH; install MACS2/MACS3 or use "
            "call_peaks() for a threshold caller"
        )
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    command = [
        executable,
        "callpeak",
        "-t",
        str(treatment),
        "-n",
        str(name),
        "--outdir",
        str(outdir),
    ]
    if control is not None:
        command += ["-c", str(control)]
    if genome:
        command += ["-g", str(genome)]
    command += list(extra_args)
    subprocess.run(command, check=True)

    candidates = [
        outdir / f"{name}_peaks.narrowPeak",
        outdir / f"{name}_peaks.broadPeak",
    ]
    for candidate in candidates:
        if candidate.exists():
            return read_peaks(candidate)
    raise FileNotFoundError(
        f"{backend} finished but no narrowPeak/broadPeak file was found in {outdir}"
    )


def write_peaks(peaks: pd.DataFrame, path: Union[str, Path]) -> None:
    """Write a peak table as BED6+ (``chrom start end name score strand``)."""

    columns = ["chrom", "start", "end", "name", "score", "strand"]
    extra = [column for column in peaks.columns if column not in columns]
    write_table(peaks[columns + extra], path)
