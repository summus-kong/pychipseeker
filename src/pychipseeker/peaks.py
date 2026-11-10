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


def _open_bigwig(path: Union[str, Path]):
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "pyBigWig is required to call peaks from a bigWig; install with "
            "`pip install pychipseeker[bigwig]`"
        ) from exc
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise ValueError(f"could not open bigWig: {path}")
    return handle


def call_peaks(
    bigwig: Union[str, Path],
    *,
    threshold: float = 3.0,
    max_gap: int = 1000,
    min_width: int = 500,
    chroms: Optional[Sequence[str]] = None,
    sample: Optional[int] = None,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """Call peaks as runs of signal above ``threshold``.

    Consecutive intervals are merged when the gap between them is at most
    ``max_gap`` bases; merged runs shorter than ``min_width`` are discarded.
    Returns a BED6+ table with ``signalValue`` (mean signal over the peak) and
    ``peak`` (offset of the highest base from the peak start).
    """

    threshold = float(threshold)
    max_gap = int(max_gap)
    min_width = int(min_width)
    if max_gap < 0 or min_width < 0:
        raise ValueError("max_gap and min_width must be non-negative")

    bw = _open_bigwig(bigwig)
    rows: list[dict[str, object]] = []
    try:
        available = list(bw.chroms().keys())
        selected = available if chroms is None else [str(c) for c in chroms]
        for chrom in selected:
            intervals = bw.intervals(chrom)
            if not intervals:
                continue
            starts, ends, values = [], [], []
            for start, end, value in intervals:
                if value is None or not np.isfinite(value) or value < threshold:
                    continue
                starts.append(int(start))
                ends.append(int(end))
                values.append(float(value))
            if not starts:
                continue

            run_start = starts[0]
            run_end = ends[0]
            run_sum = values[0] * (ends[0] - starts[0])
            run_max = values[0]
            run_peak = starts[0]
            for start, end, value in zip(starts[1:], ends[1:], values[1:]):
                if start - run_end <= max_gap:
                    if value > run_max:
                        run_max = value
                        run_peak = start
                    run_sum += value * (end - start)
                    run_end = max(run_end, end)
                    continue

                width = run_end - run_start
                if width >= min_width:
                    rows.append(
                        {
                            "chrom": str(chrom),
                            "start": run_start,
                            "end": run_end,
                            "name": f"Peak_{len(rows) + 1}",
                            "score": round(run_max, 6),
                            "strand": ".",
                            "signalValue": run_sum / max(width, 1),
                            "peak": run_peak - run_start,
                        }
                    )
                run_start, run_end = start, end
                run_sum = value * (end - start)
                run_max = value
                run_peak = start

            width = run_end - run_start
            if width >= min_width:
                rows.append(
                    {
                        "chrom": str(chrom),
                        "start": run_start,
                        "end": run_end,
                        "name": f"Peak_{len(rows) + 1}",
                        "score": round(run_max, 6),
                        "strand": ".",
                        "signalValue": run_sum / max(width, 1),
                        "peak": run_peak - run_start,
                    }
                )
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
