"""Extract H3K4me3 peaks from a real ENCODE signal bigWig.

Run::

    python examples/call_encode_peaks.py \
        --bigwig C:/path/ENCFF893ZGR.bigWig \
        --out examples/encode_peaks.bed

The ENCODE track is a per-base signal track, not a peak call.  Peaks are
defined as maximal runs of bases with ``score >= --threshold``, bridged
across gaps of at most ``--max-gap`` bp, and kept when at least
``--min-width`` bp wide.  ``score`` is the run's maximum signal, ``name``
is ``Peak_<n>`` and the strand is ``.``.

This is a *threshold caller*, not MACS2: it is documented as such in the
figure caption and in the accompanying report.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def call_peaks(
    bigwig: str | Path,
    threshold: float = 3.0,
    max_gap: int = 1000,
    min_width: int = 500,
    chroms: list[str] | None = None,
):
    """Return a BED6-like ``DataFrame`` of threshold-called peaks."""

    import pyBigWig
    import pandas as pd

    bw = pyBigWig.open(str(bigwig))
    try:
        sizes = bw.chroms()
        if chroms is None:
            canonical = {f"chr{i}" for i in range(1, 20)} | {"chrX", "chrY"}
            chroms = [c for c in sizes if c in canonical]
        records = []
        for chrom in chroms:
            intervals = bw.intervals(chrom) or []
            if not intervals:
                continue
            starts = np.array([x[0] for x in intervals], dtype=np.int64)
            ends = np.array([x[1] for x in intervals], dtype=np.int64)
            values = np.array([x[2] for x in intervals], dtype=float)
            mask = values >= threshold
            hot_s, hot_e, hot_v = starts[mask], ends[mask], values[mask]
            if hot_s.size == 0:
                continue

            breaks = np.flatnonzero(hot_s[1:] - hot_e[:-1] > max_gap) + 1
            run_s = np.r_[hot_s[0], hot_s[breaks]]
            run_e = np.r_[hot_e[breaks - 1], hot_e[-1]]
            widths = run_e - run_s
            keep = widths >= min_width
            if not keep.any():
                continue
            peak_max = np.maximum.reduceat(hot_v, np.r_[0, breaks])
            for s, e, v in zip(run_s[keep], run_e[keep], peak_max[keep]):
                records.append((chrom, int(s), int(e), float(v)))
    finally:
        bw.close()

    df = pd.DataFrame(records, columns=["chrom", "start", "end", "score"])
    if df.empty:
        return df
    df = df.sort_values(["chrom", "start"], kind="mergesort").reset_index(drop=True)
    df.insert(3, "name", [f"Peak_{i + 1}" for i in range(len(df))])
    df["strand"] = "."
    return df


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bigwig", required=True)
    parser.add_argument("--out", default="examples/encode_peaks.bed")
    parser.add_argument("--threshold", type=float, default=3.0)
    parser.add_argument("--max-gap", type=int, default=1000)
    parser.add_argument("--min-width", type=int, default=500)
    parser.add_argument(
        "--sample",
        type=int,
        default=0,
        help="keep only this many peaks (random, fixed seed); 0 keeps all",
    )
    parser.add_argument("--seed", type=int, default=20261026)
    args = parser.parse_args()

    df = call_peaks(
        args.bigwig,
        threshold=args.threshold,
        max_gap=args.max_gap,
        min_width=args.min_width,
    )
    print(f"called {len(df)} peaks at signal >= {args.threshold}")
    if args.sample and len(df) > args.sample:
        df = (
            df.sample(n=args.sample, random_state=args.seed)
            .sort_values(["chrom", "start"], kind="mergesort")
            .reset_index(drop=True)
        )
        df["name"] = [f"Peak_{i + 1}" for i in range(len(df))]
        print(f"sampled {len(df)} peaks (seed {args.seed})")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df[["chrom", "start", "end", "name", "score", "strand"]].to_csv(
        out, sep="\t", header=False, index=False
    )
    widths = (df["end"] - df["start"]).astype(float)
    print(f"width: median {widths.median():.0f} bp, mean {widths.mean():.0f} bp")
    print(f"score : median {df['score'].median():.1f}, max {df['score'].max():.1f}")
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KiB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
