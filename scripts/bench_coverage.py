"""Synthetic scaling benchmark for the peak-coverage tag-matrix builder.

Run from the repository root::

    python scripts/bench_coverage.py
    python scripts/bench_coverage.py --peaks 100000 --windows 100000
    python scripts/bench_coverage.py --windows 200000 --jobs 1 4

The script generates random peak and window intervals on a small synthetic
genome, then times ``get_tag_matrix(peak=..., windows=..., nbin=...)`` without
reading any annotation file.  It is intended for tracking the vectorised
coverage path in ``docs/performance.md``.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import pychipseeker as pcs  # noqa: E402


def synthetic_intervals(
    n_peaks: int,
    n_windows: int,
    *,
    n_chroms: int = 20,
    window: int = 6_000,
    seed: int = 20261111,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    sizes = np.linspace(30_000_000, 60_000_000, n_chroms, dtype=np.int64)
    names = np.array([f"chr{i + 1}" for i in range(n_chroms)])

    peak_chrom = rng.integers(0, n_chroms, n_peaks)
    peak_start = (rng.random(n_peaks) * (sizes[peak_chrom] - 2_000)).astype(np.int64)
    peaks = pd.DataFrame(
        {
            "chrom": names[peak_chrom],
            "start": peak_start,
            "end": peak_start + rng.integers(200, 1_000, n_peaks),
            "strand": np.where(rng.random(n_peaks) < 0.5, "-", "+"),
        }
    )

    win_chrom = rng.integers(0, n_chroms, n_windows)
    win_start = (rng.random(n_windows) * (sizes[win_chrom] - window - 1)).astype(
        np.int64
    )
    windows = pd.DataFrame(
        {
            "chrom": names[win_chrom],
            "start": win_start,
            "end": win_start + window,
            "strand": np.where(rng.random(n_windows) < 0.5, "-", "+"),
        }
    )
    return peaks, windows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", type=int, default=100_000)
    parser.add_argument("--windows", type=int, default=100_000)
    parser.add_argument("--bins", type=int, default=100)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument(
        "--jobs",
        type=int,
        nargs="+",
        default=[1],
        help="worker counts to benchmark (-1 = all cores)",
    )
    args = parser.parse_args(argv)

    peaks, windows = synthetic_intervals(
        args.peaks, args.windows, seed=args.seed
    )
    print(f"pychipseeker {pcs.__version__}")
    print(
        f"{args.peaks:,} peaks, {args.windows:,} windows "
        f"x {args.bins} bins"
    )

    for jobs in args.jobs:
        best = float("inf")
        shape = None
        for _ in range(max(1, args.repeat)):
            start = time.perf_counter()
            matrix = pcs.get_tag_matrix(
                peak=peaks,
                windows=windows,
                nbin=args.bins,
                n_jobs=jobs,
            )
            best = min(best, time.perf_counter() - start)
            shape = matrix.shape
        rate = args.windows / best if best else float("inf")
        print(
            f"n_jobs={jobs:<3} {best:8.2f} s   {rate:10,.0f} windows/s   "
            f"shape={shape}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
