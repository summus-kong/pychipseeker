"""Small interval-index helpers used by the annotation engine.

The implementation is deliberately dependency-light.  It uses a fixed-bin
index for overlap queries and sorted TSS arrays for nearest-TSS queries.  This
is fast enough for typical peak sets and avoids requiring pyranges/pybedtools.

The ``follow``/``precede`` helpers use the same 5'->3' orientation convention
as Bioconductor's ``GenomicRanges`` methods: on a minus-strand query,
"downstream" means decreasing genomic coordinates.  Peaks or features with an
unknown strand (``*`` or ``.``) are treated as unstranded.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Optional, Tuple, cast

import numpy as np
import pandas as pd


def _normalise_strand(value: object) -> str:
    text = str(value)
    if text in {"+", "-"}:
        return text
    return "*"


def _normalise_strands(values: Iterable[object]) -> np.ndarray:
    """Vectorised :func:`_normalise_strand` for a column of values."""

    if hasattr(values, "to_numpy"):
        arr = np.asarray(cast(Any, values).astype(str), dtype="U1")
    else:
        arr = np.asarray([str(v) for v in values], dtype="U1")
    out = np.full(arr.shape, "*", dtype="U1")
    out[arr == "+"] = "+"
    out[arr == "-"] = "-"
    return out


def _strand_compatible(
    query: str,
    subject: str,
    same_strand: bool,
) -> bool:
    if not same_strand:
        return True
    query = _normalise_strand(query)
    subject = _normalise_strand(subject)
    if query == "*" or subject == "*":
        return True
    return query == subject


def _min_row_for_value(
    vals: np.ndarray, rows: np.ndarray, value: int
) -> int:
    """Smallest row index among entries equal to ``value``.

    ``vals`` is sorted but may contain ties; GenomicRanges tie-breaking uses
    the first feature in subject order, i.e. the smallest row index.
    """

    lo = int(np.searchsorted(vals, value, side="left"))
    hi = int(np.searchsorted(vals, value, side="right"))
    if hi - lo == 1:
        return int(rows[lo])
    return int(rows[lo:hi].min())


def _first_greater(
    vals: np.ndarray, rows: np.ndarray, threshold: int
) -> Optional[int]:
    """Row index of the smallest value strictly above ``threshold``."""

    i = int(np.searchsorted(vals, threshold, side="right"))
    if i >= len(vals):
        return None
    return _min_row_for_value(vals, rows, int(vals[i]))


def _last_smaller(
    vals: np.ndarray, rows: np.ndarray, threshold: int
) -> Optional[int]:
    """Row index of the largest value strictly below ``threshold``."""

    j = int(np.searchsorted(vals, threshold, side="left")) - 1
    if j < 0:
        return None
    return _min_row_for_value(vals, rows, int(vals[j]))


class BinnedIntervalIndex:
    """A simple bin-based interval index for overlap queries.

    Parameters
    ----------
    bin_size:
        Bin size in base pairs.  One megabase works well for genes; 100 kb is
        a good compromise for exons/UTRs.
    """

    def __init__(self, bin_size: int = 1_000_000) -> None:
        if int(bin_size) <= 0:
            raise ValueError("bin_size must be positive")
        self.bin_size = int(bin_size)
        self._bins: Dict[str, Dict[int, np.ndarray]] = {}
        self._starts: Optional[np.ndarray] = None
        self._ends: Optional[np.ndarray] = None
        self._strands: Optional[np.ndarray] = None

    def build(
        self,
        df: pd.DataFrame,
        chrom_col: str = "chrom",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: str = "strand",
    ) -> "BinnedIntervalIndex":
        """Index a DataFrame of intervals.

        The DataFrame is expected to have 0-based, half-open coordinates.
        """

        if len(df) == 0:
            self._starts = np.empty(0, dtype=np.int64)
            self._ends = np.empty(0, dtype=np.int64)
            self._strands = np.empty(0, dtype="U1")
            self._bins = {}
            return self

        missing = {chrom_col, start_col, end_col} - set(df.columns)
        if missing:
            raise ValueError(f"missing columns for interval index: {sorted(missing)}")

        starts = df[start_col].to_numpy(dtype=np.int64, copy=False)
        ends = df[end_col].to_numpy(dtype=np.int64, copy=False)
        if np.any(ends <= starts):
            bad = int(np.flatnonzero(ends <= starts)[0])
            raise ValueError(
                f"invalid interval at row {bad}: end ({ends[bad]}) <= start ({starts[bad]})"
            )

        self._starts = starts
        self._ends = ends
        if strand_col in df.columns:
            self._strands = _normalise_strands(df[strand_col])
        else:
            self._strands = np.full(len(df), "*", dtype="U1")
        self._bins = {}

        chroms = df[chrom_col].astype(str).to_numpy()
        chrom_codes, chrom_names = pd.factorize(chroms, sort=False)

        # Vectorised bin construction.  Most intervals touch one or two bins;
        # long genes can span several, so expand the (row, bin) pairs first and
        # group them with a single lexsort instead of a Python loop.
        bin_starts = starts // self.bin_size
        bin_ends = (ends - 1) // self.bin_size
        counts = (bin_ends - bin_starts + 1).astype(np.int64)
        total = int(counts.sum())
        row_ids = np.repeat(np.arange(len(df), dtype=np.int64), counts)
        first_bins = np.repeat(bin_starts, counts)
        offsets = np.arange(total, dtype=np.int64) - np.repeat(
            np.cumsum(counts) - counts, counts
        )
        bin_ids = first_bins + offsets

        order = np.lexsort((row_ids, bin_ids, chrom_codes[row_ids]))
        rows_sorted = row_ids[order]
        bins_sorted = bin_ids[order]
        codes_sorted = chrom_codes[rows_sorted]

        change = np.empty(total, dtype=bool)
        change[0] = True
        np.not_equal(bins_sorted[1:], bins_sorted[:-1], out=change[1:])
        change[1:] |= codes_sorted[1:] != codes_sorted[:-1]
        starts_at = np.flatnonzero(change)
        stops_at = np.append(starts_at[1:], total)

        for lo, hi in zip(starts_at, stops_at):
            chrom = chrom_names[int(codes_sorted[lo])]
            bin_id = int(bins_sorted[lo])
            self._bins.setdefault(str(chrom), {})[bin_id] = rows_sorted[lo:hi]

        return self

    def query(
        self,
        chrom: str,
        start: int,
        end: int,
        strand: Optional[str] = None,
        same_strand: bool = False,
    ) -> np.ndarray:
        """Return row indices overlapping ``[start, end)``.

        If ``strand`` and ``same_strand`` are supplied, only features compatible
        with the query strand are returned.  A feature/query with unknown
        strand (``*`` or ``.``) is compatible with everything.
        """

        if self._starts is None or self._ends is None:
            raise RuntimeError("BinnedIntervalIndex.build() must be called first")
        if int(end) <= int(start):
            return np.empty(0, dtype=np.int64)

        chrom_bins = self._bins.get(str(chrom))
        if not chrom_bins:
            return np.empty(0, dtype=np.int64)

        b0 = int(start) // self.bin_size
        b1 = (int(end) - 1) // self.bin_size
        parts = [
            chrom_bins[b] for b in range(b0, b1 + 1) if b in chrom_bins
        ]
        if not parts:
            return np.empty(0, dtype=np.int64)

        if len(parts) == 1:
            idx = parts[0]
        else:
            idx = np.unique(np.concatenate(parts))
        mask = (self._starts[idx] < int(end)) & (self._ends[idx] > int(start))
        if same_strand and strand is not None and self._strands is not None:
            q = _normalise_strand(strand)
            s = self._strands[idx]
            mask &= (q == "*") | (s == "*") | (q == s)
        return np.sort(idx[mask])


class _SortedFeatureIndex:
    """Base class for sorted start/TSS feature indexes."""

    def __init__(self, features: pd.DataFrame, key: str) -> None:
        self._features = features.reset_index(drop=True)
        self._key = key
        self._by_key: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

        chroms = self._features["chrom"].astype(str).to_numpy()
        strands = _normalise_strands(self._features["strand"])
        rows = np.arange(len(self._features), dtype=np.int64)
        values = self._features[key].to_numpy(dtype=np.int64, copy=False)

        # Build the full per-chromosome index and the three strand views once.
        for chrom in pd.unique(chroms):
            chrom_mask = chroms == chrom
            for strand_key in ("all", "+", "-", "*"):
                if strand_key == "all":
                    mask = chrom_mask
                else:
                    mask = chrom_mask & (strands == strand_key)
                if not mask.any():
                    continue
                vals = values[mask]
                rws = rows[mask]
                order = np.argsort(vals, kind="mergesort")
                self._by_key[(str(chrom), strand_key)] = (
                    vals[order],
                    rws[order],
                    strands[mask][order],
                )

    def _arrays(
        self,
        chrom: str,
        strand: str,
        same_strand: bool,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        chrom = str(chrom)
        if not same_strand:
            return self._by_key.get((chrom, "all"))
        q = _normalise_strand(strand)
        if q == "*":
            return self._by_key.get((chrom, "all"))
        if q == "+":
            entries = [
                self._by_key.get((chrom, "+")),
                self._by_key.get((chrom, "*")),
            ]
        else:
            entries = [
                self._by_key.get((chrom, "-")),
                self._by_key.get((chrom, "*")),
            ]
        entries = [e for e in entries if e is not None]
        if not entries:
            return None
        if len(entries) == 1:
            return entries[0]
        present = [
            cast(Tuple[np.ndarray, np.ndarray, np.ndarray], entry) for entry in entries
        ]
        vals = np.concatenate([e[0] for e in present])
        rows = np.concatenate([e[1] for e in present])
        strands = np.concatenate([e[2] for e in present])
        order = np.argsort(vals, kind="mergesort")
        return vals[order], rows[order], strands[order]

class TSSIndex(_SortedFeatureIndex):
    """Sorted transcription-start-site index for nearest-TSS queries."""

    def __init__(self, features: pd.DataFrame) -> None:
        if "chrom" not in features or "tss" not in features:
            raise ValueError("features must contain 'chrom' and 'tss' columns")
        super().__init__(features, key="tss")

    def nearest(self, chrom: str, pos: int) -> Optional[Tuple[int, int]]:
        """Return ``(feature_row_index, pos - tss)`` for the nearest TSS.

        ``pos - tss`` is the raw signed delta; the caller must apply strand
        orientation to obtain the biological signed distance.
        """

        entry = self._by_key.get((str(chrom), "all"))
        if entry is None:
            return None

        arr, rows, _ = entry
        if len(arr) == 0:
            return None

        i = int(np.searchsorted(arr, int(pos), side="left"))
        candidates: list[int] = []
        if i > 0:
            candidates.append(i - 1)
        if i < len(arr):
            candidates.append(i)

        best = min(
            candidates,
            key=lambda j: (abs(int(arr[j]) - int(pos)), int(arr[j]), int(rows[j])),
        )
        return int(rows[best]), int(pos) - int(arr[best])

    def within(
        self,
        chrom: str,
        start: int,
        end: int,
        strand: Optional[str] = None,
        same_strand: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return feature row indices and TSS values in ``[start, end)``."""

        entry = self._arrays(chrom, strand or "*", same_strand)
        if entry is None:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)

        arr, rows, strands = entry
        lo = int(np.searchsorted(arr, int(start), side="left"))
        hi = int(np.searchsorted(arr, int(end), side="left"))
        arr = arr[lo:hi]
        rows = rows[lo:hi]
        if same_strand and strand is not None:
            q = _normalise_strand(strand)
            mask = (q == "*") | (strands[lo:hi] == "*") | (strands[lo:hi] == q)
            rows = rows[mask]
            arr = arr[mask]
        return rows, arr

    def follow_precede(
        self,
        chrom: str,
        peak_start: int,
        peak_end: int,
        strand: str = "*",
        same_strand: bool = False,
    ) -> Tuple[Optional[int], Optional[int]]:
        """Return nearest downstream and upstream feature row indices.

        ``follow`` is downstream in the query's 5'->3' orientation and
        ``precede`` is upstream.  The peak end is interpreted as an inclusive
        coordinate (BED end - 1).
        """

        entry = self._arrays(chrom, strand, same_strand)
        if entry is None:
            return None, None

        vals, rows, _ = entry
        if len(vals) == 0:
            return None, None

        q = _normalise_strand(strand) if same_strand else "*"
        end_inclusive = int(peak_end) - 1
        start = int(peak_start)

        # IRanges semantics: ``follow`` is the subject interval directly
        # *followed by* the query, i.e. the upstream neighbour; ``precede`` is
        # the downstream neighbour.  On a minus-strand query the 5'->3'
        # direction is reversed.  The TSS array is sorted, so binary search
        # replaces the previous full-chromosome masks.
        if q == "-":
            follow_idx = _first_greater(vals, rows, end_inclusive)
            precede_idx = _last_smaller(vals, rows, start)
        else:
            follow_idx = _last_smaller(vals, rows, start)
            precede_idx = _first_greater(vals, rows, end_inclusive)
        return follow_idx, precede_idx

    def overlapping(
        self,
        chrom: str,
        start: int,
        end: int,
        strand: str = "*",
        same_strand: bool = False,
    ) -> np.ndarray:
        """Return feature rows whose TSS point lies in ``[start, end)``."""

        rows, _ = self.within(chrom, start, end, strand=strand, same_strand=same_strand)
        return np.sort(rows)


class SpanIndex(_SortedFeatureIndex):
    """Sorted full-feature index used for downstream/flanking queries."""

    def __init__(self, features: pd.DataFrame) -> None:
        missing = {"chrom", "start", "end", "strand"} - set(features.columns)
        if missing:
            raise ValueError(f"features are missing columns: {sorted(missing)}")
        self._spans = features.reset_index(drop=True)
        self._starts = self._spans["start"].to_numpy(dtype=np.int64, copy=False)
        self._ends = self._spans["end"].to_numpy(dtype=np.int64, copy=False)
        self._by_key: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        self._by_end: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

        chroms = self._spans["chrom"].astype(str).to_numpy()
        strands = _normalise_strands(self._spans["strand"])
        rows = np.arange(len(self._spans), dtype=np.int64)

        for chrom in pd.unique(chroms):
            chrom_mask = chroms == chrom
            for strand_key in ("all", "+", "-", "*"):
                if strand_key == "all":
                    mask = chrom_mask
                else:
                    mask = chrom_mask & (strands == strand_key)
                if not mask.any():
                    continue
                starts = self._starts[mask]
                ends = self._ends[mask]
                rws = rows[mask]
                order = np.argsort(starts, kind="mergesort")
                self._by_key[(str(chrom), strand_key)] = (
                    starts[order],
                    rws[order],
                    strands[mask][order],
                )
                order = np.argsort(ends, kind="mergesort")
                self._by_end[(str(chrom), strand_key)] = (
                    ends[order],
                    rws[order],
                    strands[mask][order],
                )

    def _view(
        self,
        chrom: str,
        strand: str,
        same_strand: bool,
        by_end: bool,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        chrom = str(chrom)
        if not same_strand:
            return (self._by_end if by_end else self._by_key).get((chrom, "all"))
        q = _normalise_strand(strand)
        if q == "*":
            return (self._by_end if by_end else self._by_key).get((chrom, "all"))
        keys = ("+", "*") if q == "+" else ("-", "*")
        source = self._by_end if by_end else self._by_key
        entries = [source.get((chrom, key)) for key in keys]
        entries = [e for e in entries if e is not None]
        if not entries:
            return None
        if len(entries) == 1:
            return entries[0]
        present = [
            cast(Tuple[np.ndarray, np.ndarray, np.ndarray], entry) for entry in entries
        ]
        vals = np.concatenate([e[0] for e in present])
        rows = np.concatenate([e[1] for e in present])
        strands = np.concatenate([e[2] for e in present])
        order = np.argsort(vals, kind="mergesort")
        return vals[order], rows[order], strands[order]

    def follow_precede(
        self,
        chrom: str,
        peak_start: int,
        peak_end: int,
        strand: str = "*",
        same_strand: bool = False,
    ) -> Tuple[Optional[int], Optional[int]]:
        """Return nearest upstream (``follow``) and downstream (``precede``)
        feature spans.

        A feature must lie completely upstream/downstream of the peak, matching
        Bioconductor's ``follow``/``precede`` behaviour on full features.  This
        is used by ChIPseeker's downstream annotation, which looks at the gene
        immediately upstream of a peak (the peak lies downstream of that gene).
        """

        q = _normalise_strand(strand) if same_strand else "*"
        start = int(peak_start)
        end_inclusive = int(peak_end) - 1

        # ``_by_key`` is sorted by start and ``_by_end`` by end, so both
        # nearest-neighbour lookups are binary searches rather than full
        # chromosome masks.
        start_entry = self._view(chrom, strand, same_strand, by_end=False)
        end_entry = self._view(chrom, strand, same_strand, by_end=True)
        follow_idx = None
        precede_idx = None
        if q == "-":
            # Upstream is to the right; downstream is to the left.
            if start_entry is not None:
                vals, rows, _ = start_entry
                follow_idx = _first_greater(vals, rows, end_inclusive)
            if end_entry is not None:
                vals, rows, _ = end_entry
                precede_idx = _last_smaller(vals, rows, start)
            return follow_idx, precede_idx

        if end_entry is not None:
            vals, rows, _ = end_entry
            follow_idx = _last_smaller(vals, rows, start)
        if start_entry is not None:
            vals, rows, _ = start_entry
            precede_idx = _first_greater(vals, rows, end_inclusive)
        return follow_idx, precede_idx

    def overlapping(
        self,
        chrom: str,
        start: int,
        end: int,
        strand: str = "*",
        same_strand: bool = False,
    ) -> np.ndarray:
        """Return feature rows overlapping ``[start, end)``."""

        start_view = self._view(chrom, strand, same_strand, by_end=False)
        end_view = self._view(chrom, strand, same_strand, by_end=True)
        if start_view is None or end_view is None:
            return np.empty(0, dtype=np.int64)
        starts, start_rows, _ = start_view
        ends, end_rows, _ = end_view
        n_rows = len(start_rows)

        # Features starting before ``end`` and features ending after ``start``
        # must both be filtered.  Scan whichever side is shorter.
        before_end = int(np.searchsorted(starts, int(end), side="left"))
        after_start = int(np.searchsorted(ends, int(start), side="right"))
        if before_end <= n_rows - after_start:
            candidates = start_rows[:before_end]
            mask = self._ends[candidates] > int(start)
        else:
            candidates = end_rows[after_start:]
            mask = self._starts[candidates] < int(end)
        return np.sort(candidates[mask])
