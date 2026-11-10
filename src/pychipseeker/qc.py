"""Quality-control metrics for ChIP-seq peak sets and signal tracks.

The helpers mirror the checks ENCODE recommends before interpreting peaks:

* peak-count and width summaries, mitochondrial/blacklist fractions;
* FRiP (fraction of reads in peaks) for bigWig coverage or BAM alignments;
* TSS enrichment, the ratio of signal at transcription start sites to the
  flanking background.

They are dependency-light: bigWig support uses the optional ``pyBigWig``
package and BAM support uses the optional ``pysam`` package, matching the rest
of pychipseeker.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd

from .io import read_peaks
from .seqlevels import canonical_chrom, harmonize_frame


def _coerce_peaks(peaks: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    if isinstance(peaks, (str, Path)):
        frame = read_peaks(peaks)
    elif isinstance(peaks, pd.DataFrame):
        frame = peaks.copy()
    else:
        raise TypeError("peaks must be a path or a pandas.DataFrame")
    required = {"chrom", "start", "end"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"peak table is missing columns: {sorted(missing)}")
    frame["chrom"] = frame["chrom"].astype(str)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame.dropna(subset=["start", "end"])
    frame = frame[frame["end"] > frame["start"]].copy()
    frame["start"] = frame["start"].astype(np.int64)
    frame["end"] = frame["end"].astype(np.int64)
    return frame.reset_index(drop=True)


def _chrom_sizes_from_genome(genome: Union[str, Path]) -> Dict[str, int]:
    from .genomes import load_genome

    annotation = load_genome(genome)
    sizes: Dict[str, int] = {}
    for table in (
        annotation.genes,
        annotation.transcripts,
        annotation.exons,
    ):
        if table is None or table.empty:
            continue
        for chrom, value in table.groupby("chrom", sort=False)["end"].max().items():
            sizes[str(chrom)] = max(sizes.get(str(chrom), 0), int(value))
    return sizes


def _merge_intervals(
    peaks: pd.DataFrame,
    chrom_sizes: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """Merge overlapping/book-ended peak intervals per chromosome."""

    rows: list[tuple[str, int, int]] = []
    for chrom, group in peaks.groupby("chrom", sort=False):
        starts = group["start"].to_numpy(dtype=np.int64)
        ends = group["end"].to_numpy(dtype=np.int64)
        if chrom_sizes is not None:
            size = int(chrom_sizes.get(str(chrom), 0))
            if size > 0:
                keep = (starts < size) & (ends > 0)
                starts = np.clip(starts[keep], 0, size)
                ends = np.clip(ends[keep], 0, size)
                good = ends > starts
                starts = starts[good]
                ends = ends[good]
        if starts.size == 0:
            continue
        order = np.argsort(starts, kind="mergesort")
        starts = starts[order]
        ends = ends[order]
        run_start = int(starts[0])
        run_end = int(ends[0])
        for start, end in zip(starts[1:], ends[1:]):
            start, end = int(start), int(end)
            if start <= run_end:
                run_end = max(run_end, end)
                continue
            rows.append((str(chrom), run_start, run_end))
            run_start, run_end = start, end
        rows.append((str(chrom), run_start, run_end))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def peak_qc(
    peaks: Union[str, Path, pd.DataFrame],
    *,
    blacklist: Optional[Union[str, Path, pd.DataFrame]] = None,
    chrom_sizes: Optional[Dict[str, int]] = None,
    genome: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """Return a long-form table of peak quality metrics.

    Columns are ``Metric`` and ``Value``.  Fractions are reported in the
    0-1 range; widths and coordinates are in base pairs.
    """

    frame = _coerce_peaks(peaks)
    if frame.empty:
        raise ValueError("peak table is empty")
    if chrom_sizes is None and genome is not None:
        chrom_sizes = _chrom_sizes_from_genome(genome)
    if chrom_sizes is not None:
        frame, _ = harmonize_frame(frame, chrom_sizes.keys(), column="chrom", copy=False)

    widths = (frame["end"] - frame["start"]).to_numpy(dtype=np.int64)
    n_peaks = len(frame)
    metrics: list[tuple[str, object]] = [
        ("Peaks", n_peaks),
        ("Chromosomes", int(frame["chrom"].nunique())),
        ("Total width (bp)", int(widths.sum())),
        ("Mean width (bp)", float(widths.mean())),
        ("Median width (bp)", float(np.median(widths))),
        ("Min width (bp)", int(widths.min())),
        ("Max width (bp)", int(widths.max())),
        (
            "Mitochondrial peaks",
            int((frame["chrom"].map(canonical_chrom) == "m").sum()),
        ),
        (
            "Mitochondrial fraction",
            float((frame["chrom"].map(canonical_chrom) == "m").mean()),
        ),
        (
            "Duplicate fraction",
            float(frame.duplicated(subset=["chrom", "start", "end"]).mean()),
        ),
    ]

    if chrom_sizes:
        total = sum(int(size) for size in chrom_sizes.values())
        if total > 0:
            metrics.append(("Peaks per Mb", n_peaks / (total / 1_000_000.0)))

    if blacklist is not None:
        from .peaks import _overlaps_any

        if isinstance(blacklist, (str, Path)):
            blacklist_frame = read_peaks(blacklist)
        elif isinstance(blacklist, pd.DataFrame):
            blacklist_frame = blacklist
        else:
            raise TypeError("blacklist must be a path or a pandas.DataFrame")
        blacklist_frame, _ = harmonize_frame(
            blacklist_frame, frame["chrom"].unique(), column="chrom", copy=True
        )
        hits = _overlaps_any(frame, blacklist_frame)
        metrics.append(("Blacklisted peaks", int(hits.sum())))
        metrics.append(("Blacklisted fraction", float(hits.mean())))

    return pd.DataFrame(metrics, columns=["Metric", "Value"])


def _bigwig_handle(path: Union[str, Path]):
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "pyBigWig is required for bigWig QC; install with "
            "`pip install pychipseeker[bigwig]`"
        ) from exc
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise ValueError(f"could not open bigWig: {path}")
    return handle


def _bam_handle(path: Union[str, Path]):
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "pysam is required for BAM/CRAM QC; install with "
            "`pip install pychipseeker[bam]`"
        ) from exc
    return pysam.AlignmentFile(str(path), "rb")


def _frip_bigwig(
    bigwig: Union[str, Path],
    peaks: pd.DataFrame,
    *,
    chrom_sizes: Optional[Dict[str, int]] = None,
) -> float:
    bw = _bigwig_handle(bigwig)
    try:
        sizes = dict(bw.chroms())
        if chrom_sizes:
            sizes = {str(k): int(v) for k, v in chrom_sizes.items()}
        peaks, _ = harmonize_frame(peaks, sizes.keys(), column="chrom", copy=True)
        merged = _merge_intervals(peaks, sizes)
        total = 0.0
        for chrom, size in sizes.items():
            value = bw.stats(str(chrom), 0, int(size), nBins=1, type="sum")
            if value and value[0] is not None:
                total += float(value[0])
        in_peaks = 0.0
        for row in merged.itertuples(index=False):
            size = sizes.get(str(row.chrom))
            if size is None:
                continue
            start = max(0, int(row.start))
            end = min(int(size), int(row.end))
            if end <= start:
                continue
            value = bw.stats(
                str(row.chrom), start, end, nBins=1, type="sum"
            )
            if value and value[0] is not None:
                in_peaks += float(value[0])
    finally:
        bw.close()
    if total <= 0:
        raise ValueError("bigWig contains no signal")
    return float(in_peaks / total)


def _frip_bam(
    bam: Union[str, Path],
    peaks: pd.DataFrame,
    *,
    exclude_duplicates: bool = True,
    min_mapping_quality: int = 0,
) -> float:
    handle = _bam_handle(bam)
    merged = _merge_intervals(peaks)
    references = list(getattr(handle, "references", []) or [])
    if references:
        merged, _ = harmonize_frame(
            merged, references, column="chrom", copy=False
        )
    by_chrom: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, group in merged.groupby("chrom", sort=False):
        starts = group["start"].to_numpy(dtype=np.int64)
        ends = group["end"].to_numpy(dtype=np.int64)
        order = np.argsort(starts, kind="mergesort")
        by_chrom[str(chrom)] = (starts[order], ends[order])

    total = 0
    in_peaks = 0
    try:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if exclude_duplicates and read.is_duplicate:
                continue
            if int(read.mapping_quality) < int(min_mapping_quality):
                continue
            chrom = read.reference_name
            if chrom is None:
                continue
            start = int(read.reference_start)
            end = read.reference_end
            if end is None:
                end = start + max(1, int(read.query_length or 1))
            end = int(end)
            total += 1
            if str(chrom) not in by_chrom:
                continue
            peak_starts, peak_ends = by_chrom[str(chrom)]
            index = int(np.searchsorted(peak_ends, start, side="right"))
            if index < len(peak_starts) and peak_starts[index] < end:
                in_peaks += 1
    finally:
        handle.close()
    if total <= 0:
        raise ValueError("BAM contains no mapped reads")
    return float(in_peaks / total)


def frip(
    peaks: Union[str, Path, pd.DataFrame],
    *,
    bigwig: Optional[Union[str, Path]] = None,
    bam: Optional[Union[str, Path]] = None,
    chrom_sizes: Optional[Dict[str, int]] = None,
    exclude_duplicates: bool = True,
    min_mapping_quality: int = 0,
) -> float:
    """Fraction of signal/reads falling inside peaks.

    Supply exactly one of ``bigwig`` (signal sum) or ``bam`` (read counts).
    Peaks are merged first so overlapping intervals are not double-counted.
    """

    if (bigwig is None) == (bam is None):
        raise ValueError("supply exactly one of bigwig=... or bam=...")
    frame = _coerce_peaks(peaks)
    if frame.empty:
        raise ValueError("peak table is empty")
    if bigwig is not None:
        return _frip_bigwig(bigwig, frame, chrom_sizes=chrom_sizes)
    return _frip_bam(
        bam,
        frame,
        exclude_duplicates=exclude_duplicates,
        min_mapping_quality=min_mapping_quality,
    )


def _bam_tss_matrix(
    bam: Union[str, Path],
    genes: pd.DataFrame,
    *,
    upstream: int,
    downstream: int,
    n_bins: int,
) -> np.ndarray:
    handle = _bam_handle(bam)
    matrix = np.full((len(genes), n_bins), np.nan, dtype=np.float64)
    try:
        for row_idx, row in enumerate(genes.itertuples(index=False)):
            chrom = str(row.chrom)
            tss = int(row.tss)
            start = max(0, tss - int(upstream))
            end = tss + int(downstream)
            if end <= start:
                continue
            try:
                coverage = handle.count_coverage(
                    chrom, start, end, quality_threshold=0
                )
            except (ValueError, KeyError):
                continue
            if not coverage:
                continue
            values = np.sum(
                np.asarray(coverage, dtype=np.float64), axis=0
            )
            if values.size == 0 or not np.isfinite(values).any():
                continue
            if values.size == n_bins:
                binned = values
            else:
                edges = (
                    np.arange(n_bins + 1, dtype=np.int64) * values.size
                ) // n_bins
                binned = np.array(
                    [
                        values[edges[i] : edges[i + 1]].mean()
                        if edges[i + 1] > edges[i]
                        else 0.0
                        for i in range(n_bins)
                    ],
                    dtype=np.float64,
                )
            if not np.isfinite(binned).all():
                binned = np.nan_to_num(binned, nan=0.0)
            if str(row.strand) == "-":
                binned = binned[::-1]
            matrix[row_idx] = binned
    finally:
        handle.close()
    return matrix


def _profile_score(
    matrix: np.ndarray,
    *,
    upstream: int,
    downstream: int,
    center_window: int,
    flank_window: Optional[int],
) -> dict[str, float]:
    if matrix.ndim != 2 or matrix.shape[1] < 2:
        raise ValueError("matrix must be features x bins")
    valid = np.isfinite(matrix)
    counts = valid.sum(axis=0)
    totals = np.where(valid, matrix, 0.0).sum(axis=0)
    n_bins = matrix.shape[1]
    profile = np.divide(
        totals,
        counts,
        out=np.full(n_bins, np.nan),
        where=counts > 0,
    )
    x = np.linspace(-int(upstream), int(downstream), n_bins)
    center = np.abs(x) <= int(center_window)
    if flank_window is None:
        flank_threshold = max(0, min(int(upstream), int(downstream)) - int(center_window))
    else:
        flank_threshold = int(flank_window)
    flank = np.abs(x) >= flank_threshold
    if not center.any() or not flank.any():
        raise ValueError("window is too small for TSS enrichment")

    n_features = int(np.sum(~np.all(np.isnan(matrix), axis=1)))
    if n_features == 0:
        return {
            "TSS enrichment": float("nan"),
            "TSS center signal": float("nan"),
            "TSS flank signal": float("nan"),
            "TSS features": 0,
        }

    def _mean(values: np.ndarray) -> float:
        finite = values[np.isfinite(values)]
        # Bins without coverage contribute zero signal rather than NaN.
        return float(finite.mean()) if finite.size else 0.0

    center_signal = _mean(profile[center])
    flank_signal = _mean(profile[flank])
    # A pseudocount keeps the ratio finite when the flanks have no coverage,
    # which is common for compact peak-only tracks.
    pseudocount = 1.0
    score = (center_signal + pseudocount) / (flank_signal + pseudocount)
    return {
        "TSS enrichment": float(score),
        "TSS center signal": center_signal,
        "TSS flank signal": flank_signal,
        "TSS features": n_features,
    }


def tss_enrichment(
    genes: pd.DataFrame,
    *,
    bigwig: Optional[Union[str, Path]] = None,
    bam: Optional[Union[str, Path]] = None,
    upstream: int = 2000,
    downstream: int = 2000,
    n_bins: int = 100,
    sample: Optional[int] = 2000,
    seed: int = 20261026,
    center_window: int = 100,
    flank_window: Optional[int] = None,
) -> dict[str, float]:
    """Compute the ENCODE-style TSS enrichment score.

    The score is the mean signal within ``+/- center_window`` bp of the TSS
    divided by the mean signal in the outer flanks of the window.  Supply
    exactly one of ``bigwig`` or ``bam``; BAM coverage is counted per base with
    ``pysam``.  A reproducible subsample of ``sample`` features is used unless
    ``sample=None``.
    """

    if (bigwig is None) == (bam is None):
        raise ValueError("supply exactly one of bigwig=... or bam=...")
    required = {"chrom", "tss", "strand"}
    missing = required - set(genes.columns)
    if missing:
        raise ValueError(f"genes table is missing columns: {sorted(missing)}")

    frame = genes.reset_index(drop=True).copy()
    if sample is not None and int(sample) > 0 and int(sample) < len(frame):
        frame = frame.sample(n=int(sample), random_state=int(seed)).reset_index(
            drop=True
        )
    frame["chrom"] = frame["chrom"].astype(str)
    frame["strand"] = frame["strand"].astype(str)
    frame["tss"] = pd.to_numeric(frame["tss"], errors="coerce")
    frame = frame.dropna(subset=["tss"]).reset_index(drop=True)
    if frame.empty:
        raise ValueError("genes table has no valid TSS coordinates")

    if bigwig is not None:
        from .plotting import extract_tss_matrix

        matrix = extract_tss_matrix(
            bigwig,
            frame,
            n_bins=int(n_bins),
            upstream=int(upstream),
            downstream=int(downstream),
            drop_empty=False,
            return_mask=False,
        )
    else:
        handle = _bam_handle(bam)
        try:
            references = list(getattr(handle, "references", []) or [])
        finally:
            handle.close()
        if references:
            frame, _ = harmonize_frame(
                frame, references, column="chrom", copy=False
            )
        matrix = _bam_tss_matrix(
            bam,
            frame,
            upstream=int(upstream),
            downstream=int(downstream),
            n_bins=int(n_bins),
        )
    return _profile_score(
        np.asarray(matrix, dtype=float),
        upstream=int(upstream),
        downstream=int(downstream),
        center_window=int(center_window),
        flank_window=flank_window,
    )


def qc_summary(
    peaks: Union[str, Path, pd.DataFrame],
    *,
    bigwig: Optional[Union[str, Path]] = None,
    bam: Optional[Union[str, Path]] = None,
    genes: Optional[pd.DataFrame] = None,
    blacklist: Optional[Union[str, Path, pd.DataFrame]] = None,
    chrom_sizes: Optional[Dict[str, int]] = None,
    genome: Optional[Union[str, Path]] = None,
    tss_upstream: int = 2000,
    tss_downstream: int = 2000,
    tss_sample: Optional[int] = 2000,
    seed: int = 20261026,
    exclude_duplicates: bool = True,
) -> pd.DataFrame:
    """Combine peak metrics, FRiP and TSS enrichment in one long table."""

    table = peak_qc(
        peaks,
        blacklist=blacklist,
        chrom_sizes=chrom_sizes,
        genome=genome,
    )
    rows: list[tuple[str, object]] = table[["Metric", "Value"]].values.tolist()
    if bigwig is not None or bam is not None:
        rows.append(
            (
                "FRiP",
                frip(
                    peaks,
                    bigwig=bigwig,
                    bam=bam,
                    chrom_sizes=chrom_sizes,
                    exclude_duplicates=exclude_duplicates,
                ),
            )
        )
    if genes is not None and (bigwig is not None or bam is not None):
        enrichment = tss_enrichment(
            genes,
            bigwig=bigwig,
            bam=bam,
            upstream=int(tss_upstream),
            downstream=int(tss_downstream),
            sample=tss_sample,
            seed=int(seed),
        )
        rows.extend(enrichment.items())
    return pd.DataFrame(rows, columns=["Metric", "Value"])
