"""Small benchmark for annotation, tag matrices and bigWig extraction.

Run from the repository root::

    python scripts/bench.py
    python scripts/bench.py --txdb mm10 --peaks examples/encode_peaks.bed
    python scripts/bench.py --txdb mm10 --bigwig signal.bw --sample 2000
    python scripts/bench.py --jobs -1 --repeat 3

The script intentionally measures end-to-end calls, not micro-operations, so
the numbers can be compared with the figures in ``docs/performance.md``.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

import pychipseeker as pcs  # noqa: E402


def _timed(label: str, func) -> float:
    start = time.perf_counter()
    value = func()
    elapsed = time.perf_counter() - start
    print(f"{label:<44} {elapsed:8.2f} s")
    return elapsed, value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--txdb", default="mm10", help="assembly name")
    parser.add_argument(
        "--peaks",
        default=str(ROOT / "examples" / "encode_peaks.bed"),
        help="peak file",
    )
    parser.add_argument("--repeat", type=int, default=2)
    parser.add_argument("--bigwig", default=None, help="optional signal bigWig")
    parser.add_argument(
        "--sample",
        type=int,
        default=3000,
        help="transcripts sampled for the bigWig matrix (0 = all)",
    )
    parser.add_argument("--window", type=int, default=2000)
    parser.add_argument("--bins", type=int, default=100)
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="worker processes for annotation/tag matrices (-1 = all cores)",
    )
    parser.add_argument(
        "--tagmatrix",
        action="store_true",
        help="also time peak-coverage tag matrices over promoters",
    )
    args = parser.parse_args(argv)

    print(f"pychipseeker {pcs.__version__}")
    _, model = _timed(f"load_genome({args.txdb!r})", lambda: pcs.load_genome(args.txdb))

    annotation = None
    for index in range(1, max(1, args.repeat) + 1):
        _, annotation = _timed(
            f"annotate_peak #{index} (jobs={args.jobs})",
            lambda: pcs.annotate_peak(
                peaks=args.peaks,
                TxDb=args.txdb,
                verbose=False,
                n_jobs=args.jobs,
            ),
        )
    if annotation is not None:
        print(f"{'peaks annotated':<44} {len(annotation.peaks):8,}")

    if args.tagmatrix:
        promoters = pcs.get_promoters(TxDb=model, upstream=3000, downstream=3000)
        _, matrix = _timed(
            f"get_tag_matrix(peak) ({len(promoters):,} windows)",
            lambda: pcs.get_tag_matrix(
                peak=args.peaks,
                windows=promoters,
                nbin=100,
                n_jobs=args.jobs,
            ),
        )
        print(f"{'tag matrix shape':<44} {matrix.shape}")

    if args.bigwig:
        transcripts = model.features("transcript")
        if args.sample and args.sample < len(transcripts):
            transcripts = transcripts.sample(n=args.sample, random_state=20261026)
        _, matrix = _timed(
            f"extract_tss_matrix({len(transcripts):,} TSS, jobs={args.jobs})",
            lambda: pcs.extract_tss_matrix(
                args.bigwig,
                transcripts,
                window=args.window,
                n_bins=args.bins,
                n_jobs=args.jobs,
            ),
        )
        print(f"{'matrix shape':<44} {matrix.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
