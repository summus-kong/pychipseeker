"""Interval indexes: binning, strand handling and nearest-feature queries."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pychipseeker.intervals import (
    BinnedIntervalIndex,
    SpanIndex,
    TSSIndex,
    _first_greater,
    _last_smaller,
    _min_row_for_value,
    _normalise_strand,
    _normalise_strands,
    _strand_compatible,
)


def _features(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand"])


def test_strand_helpers_normalise_and_compare() -> None:
    assert _normalise_strand("+") == "+"
    assert _normalise_strand("-") == "-"
    assert _normalise_strand(".") == "*"
    assert _normalise_strand(None) == "*"
    assert _normalise_strand(1) == "*"

    values = pd.Series(["+", "-", ".", None, "x"])
    np.testing.assert_array_equal(
        _normalise_strands(values), np.asarray(["+", "-", "*", "*", "*"])
    )
    np.testing.assert_array_equal(
        _normalise_strands(["+", "-", "."]), np.asarray(["+", "-", "*"])
    )

    assert _strand_compatible("+", "-", same_strand=False)
    assert _strand_compatible("+", "+", same_strand=True)
    assert _strand_compatible("+", "*", same_strand=True)
    assert _strand_compatible("*", "-", same_strand=True)
    assert not _strand_compatible("+", "-", same_strand=True)


def test_min_row_and_neighbour_helpers() -> None:
    vals = np.asarray([10, 20, 20, 30])
    rows = np.asarray([9, 7, 5, 1])
    assert _min_row_for_value(vals, rows, 20) == 5
    assert _first_greater(vals, rows, 5) == 9
    assert _first_greater(vals, rows, 10) == 5
    assert _first_greater(vals, rows, 30) is None
    assert _last_smaller(vals, rows, 30) == 5
    assert _last_smaller(vals, rows, 20) == 9
    assert _last_smaller(vals, rows, 10) is None


def test_binned_index_validation_and_empty_build() -> None:
    with pytest.raises(ValueError, match="bin_size"):
        BinnedIntervalIndex(0)

    index = BinnedIntervalIndex(10)
    with pytest.raises(RuntimeError, match="build"):
        index.query("chr1", 0, 10)

    empty = index.build(pd.DataFrame({"chrom": [], "start": [], "end": []}))
    assert empty.query("chr1", 0, 10).tolist() == []

    frame = pd.DataFrame({"chrom": ["chr1"], "start": [0]})
    with pytest.raises(ValueError, match="missing columns"):
        BinnedIntervalIndex(10).build(frame)
    with pytest.raises(ValueError, match="invalid interval"):
        BinnedIntervalIndex(10).build(
            pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [10]})
        )


def test_binned_index_overlap_across_bins_and_chromosomes() -> None:
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "start": [0, 5, 15, 0],
            "end": [3, 25, 18, 10],
        }
    )
    index = BinnedIntervalIndex(10).build(frame)
    # The long interval 5-25 spans bins 0, 1 and 2 and must appear once.
    assert index.query("chr1", 20, 30).tolist() == [1]
    assert index.query("chr1", 0, 100).tolist() == [0, 1, 2]
    assert index.query("chr1", 3, 6).tolist() == [1]
    assert index.query("chr1", 100, 200).tolist() == []
    assert index.query("chr2", 0, 10).tolist() == [3]
    assert index.query("chrZ", 0, 10).tolist() == []
    assert index.query("chr1", 10, 10).tolist() == []

    # A query starting in an empty bin but spanning an indexed one.
    assert index.query("chr1", 2, 16).tolist() == [0, 1, 2]


def test_binned_index_strand_filtering() -> None:
    frame = _features(
        [
            ("chr1", 0, 10, "+"),
            ("chr1", 20, 30, "-"),
            ("chr1", 40, 50, "."),
        ]
    )
    index = BinnedIntervalIndex(100).build(frame)
    assert index.query("chr1", 0, 100, strand="+", same_strand=True).tolist() == [0, 2]
    assert index.query("chr1", 0, 100, strand="-", same_strand=True).tolist() == [1, 2]
    assert index.query("chr1", 0, 100, strand="*", same_strand=True).tolist() == [
        0,
        1,
        2,
    ]
    assert index.query("chr1", 0, 100, strand="+", same_strand=False).tolist() == [
        0,
        1,
        2,
    ]
    # strand=None disables filtering even when same_strand=True.
    assert index.query("chr1", 0, 100, same_strand=True).tolist() == [0, 1, 2]


def test_tss_index_validation_empty_and_nearest() -> None:
    with pytest.raises(ValueError, match="chrom.*tss"):
        TSSIndex(pd.DataFrame({"chrom": ["chr1"], "start": [0]}))

    empty = TSSIndex(pd.DataFrame({"chrom": [], "tss": [], "strand": []}))
    assert empty.nearest("chr1", 100) is None
    assert empty.within("chr1", 0, 100)[0].tolist() == []
    assert empty.follow_precede("chr1", 0, 100) == (None, None)
    assert empty.overlapping("chr1", 0, 100).tolist() == []

    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "tss": [100, 200, 200, 50],
            "strand": ["+", "-", "+", "+"],
        }
    )
    index = TSSIndex(frame)
    assert index.nearest("chr1", 150) == (0, 50)
    assert index.nearest("chr1", 100) == (0, 0)
    assert index.nearest("chr1", 1000) == (2, 800)
    assert index.nearest("chr1", 0) == (0, -100)
    assert index.nearest("chrZ", 100) is None


def test_tss_index_within_and_strand_views() -> None:
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "tss": [10, 20, 30, 40],
            "strand": ["+", "-", "*", "+"],
        }
    )
    index = TSSIndex(frame)
    rows, values = index.within("chr1", 10, 40)
    assert rows.tolist() == [0, 1, 2] and values.tolist() == [10, 20, 30]

    rows, _ = index.within("chr1", 0, 50, strand="+", same_strand=True)
    assert rows.tolist() == [0, 2, 3]
    rows, _ = index.within("chr1", 0, 50, strand="-", same_strand=True)
    assert rows.tolist() == [1, 2]
    rows, _ = index.within("chr1", 0, 50, strand="*", same_strand=True)
    assert rows.tolist() == [0, 1, 2, 3]
    assert index.within("chrZ", 0, 50)[0].tolist() == []

    minus_only = TSSIndex(
        pd.DataFrame({"chrom": ["chr1"], "tss": [10], "strand": ["-"]})
    )
    assert (
        minus_only.within("chr1", 0, 50, strand="+", same_strand=True)[0].tolist() == []
    )

    assert index.overlapping("chr1", 20, 41).tolist() == [1, 2, 3]
    assert index.overlapping("chr1", 10, 10).tolist() == []


def test_tss_index_follow_precede_orientation_and_boundaries() -> None:
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "tss": [50, 100, 200, 300],
            "strand": ["+", "+", "+", "-"],
        }
    )
    index = TSSIndex(frame)

    # Unstranded: follow = last TSS before the peak start, precede = first
    # TSS after the inclusive peak end.
    assert index.follow_precede("chr1", 100, 200) == (0, 2)
    assert index.follow_precede("chr1", 50, 50) == (None, 0)
    assert index.follow_precede("chr1", 0, 1000) == (None, None)
    assert index.follow_precede("chrZ", 0, 100) == (None, None)

    # Minus strand flips 5'->3' direction; only the minus-strand feature is
    # eligible when same_strand=True.
    assert index.follow_precede("chr1", 100, 200, strand="-", same_strand=True) == (
        3,
        None,
    )
    # Unknown strand is compatible with everything.
    assert index.follow_precede("chr1", 100, 200, strand="*", same_strand=True) == (
        0,
        2,
    )


def test_span_index_validation_and_empty() -> None:
    with pytest.raises(ValueError, match="missing columns"):
        SpanIndex(pd.DataFrame({"chrom": ["chr1"], "start": [0]}))

    empty = SpanIndex(pd.DataFrame({"chrom": [], "start": [], "end": [], "strand": []}))
    assert empty.overlapping("chr1", 0, 10).tolist() == []
    assert empty.follow_precede("chr1", 0, 10) == (None, None)


def test_span_index_overlap_both_candidate_scans() -> None:
    frame = _features(
        [
            ("chr1", 0, 10, "+"),
            ("chr1", 20, 30, "-"),
            ("chr1", 40, 50, "*"),
            ("chr2", 0, 100, "+"),
        ]
    )
    index = SpanIndex(frame)
    assert index.overlapping("chr1", 5, 25).tolist() == [0, 1]
    assert index.overlapping("chr1", 10, 20).tolist() == []
    assert index.overlapping("chr1", 0, 100).tolist() == [0, 1, 2]
    assert index.overlapping("chr2", 50, 60).tolist() == [3]
    assert index.overlapping("chrZ", 0, 10).tolist() == []

    # Many features start before the query end, few end before its start:
    # this exercises the end-sorted scan.
    dense = _features([("chr1", i, i + 1, "+") for i in range(100)])
    dense_index = SpanIndex(dense)
    assert dense_index.overlapping("chr1", 99, 1000).tolist() == [99]


def test_span_index_strand_views_and_overlap() -> None:
    frame = _features(
        [
            ("chr1", 0, 10, "+"),
            ("chr1", 20, 30, "-"),
            ("chr1", 40, 50, "*"),
            ("chr2", 0, 10, "+"),
        ]
    )
    index = SpanIndex(frame)
    assert index.overlapping("chr1", 0, 100, strand="+", same_strand=True).tolist() == [
        0,
        2,
    ]
    assert index.overlapping("chr1", 0, 100, strand="-", same_strand=True).tolist() == [
        1,
        2,
    ]
    # Only the unknown-strand feature remains for a minus query on chr2.
    assert (
        index.overlapping("chr2", 0, 100, strand="-", same_strand=True).tolist() == []
    )
    assert index.overlapping("chr1", 0, 100, strand="*", same_strand=True).tolist() == [
        0,
        1,
        2,
    ]
    assert index.overlapping("chr1", 0, 100, same_strand=True).tolist() == [0, 1, 2]


def test_span_index_follow_precede_orientation_and_boundaries() -> None:
    frame = _features(
        [
            ("chr1", 0, 10, "+"),
            ("chr1", 20, 30, "+"),
            ("chr1", 40, 50, "+"),
            ("chr1", 60, 70, "-"),
        ]
    )
    index = SpanIndex(frame)

    # Unstranded: follow = feature ending before the peak start,
    # precede = feature starting after the inclusive peak end.
    assert index.follow_precede("chr1", 20, 30) == (0, 2)
    assert index.follow_precede("chr1", 0, 1000) == (None, None)
    assert index.follow_precede("chr1", 0, 10) == (None, 1)
    assert index.follow_precede("chrZ", 0, 10) == (None, None)

    # Minus strand flips the direction; only the minus-strand feature is
    # eligible when same_strand=True.
    assert index.follow_precede("chr1", 20, 30, strand="-", same_strand=True) == (
        3,
        None,
    )
    assert index.follow_precede("chr1", 80, 90, strand="-", same_strand=True) == (
        None,
        3,
    )
