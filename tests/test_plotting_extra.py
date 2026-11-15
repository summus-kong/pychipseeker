"""Branch-focused tests for plotting helpers and TSS extraction backends."""

from __future__ import annotations

import sys
import types

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pychipseeker.plotting as plotting
from pychipseeker import annotate_peak
from pychipseeker.plotting import (
    _as_df,
    _coerce_signal_frame,
    _colors,
    _distance_breaks,
    _grouped_stat,
    _mean_in_interval,
    _resolve_tss_window,
    _stat_from_annotation,
    extract_tss_matrix,
    extract_tss_matrix_from_signal,
    plot_anno_bar,
    plot_anno_pie,
    plot_dist_to_tss,
    plot_tss_heatmap,
    plot_tss_profile,
)


@pytest.fixture
def annotation(examples_dir):
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        upstream=10,
        downstream=10,
    )


@pytest.fixture
def frame(annotation) -> pd.DataFrame:
    data = annotation.to_dataframe(include_peak_metadata=True)
    data["group"] = np.where(np.arange(len(data)) % 2 == 0, "A", "B")
    return data


def test_as_df_rejects_unknown_types_and_renames_columns():
    with pytest.raises(TypeError, match="PeakAnnotation"):
        _as_df([1, 2, 3])

    renamed = _as_df(
        pd.DataFrame(
            {
                "peakChr": ["chr1"],
                "peakStart": [0],
                "peakEnd": [10],
                "distanceToTSS": [5],
            }
        )
    )
    assert {"chrom", "start", "end", "distance_to_tss"} <= set(renamed.columns)


def test_stat_from_annotation_variants():
    direct = pd.DataFrame({"Feature": ["Promoter"], "Frequency": [100.0]})
    assert _stat_from_annotation(direct).equals(direct)

    labels = pd.DataFrame(
        {"annotation": ["Promoter", "Exon", "Distal Intergenic"]}
    )
    stat = _stat_from_annotation(labels, downstream_distance=100)
    assert not stat.empty

    with pytest.raises(TypeError, match="Feature/Frequency"):
        _stat_from_annotation(pd.DataFrame({"x": [1]}))


def test_color_palette_handles_empty_and_large_requests():
    assert _colors(0) == []
    assert _colors(3) == plotting._CHIPSEEKER_COLORS[:3]

    colors = _colors(20)
    assert len(colors) == 20
    assert all(color.startswith("#") and len(color) == 7 for color in colors)
    assert colors[0] == plotting._CHIPSEEKER_COLORS[0]


def test_pie_empty_and_external_axes(annotation):
    with pytest.raises(ValueError, match="no annotations"):
        plot_anno_pie(pd.DataFrame({"Feature": [], "Frequency": []}))

    fig, ax = plt.subplots()
    returned_fig, returned_ax = plot_anno_pie(annotation, ax=ax, title="demo")
    assert returned_fig is fig
    assert returned_ax is ax
    plt.close(fig)


def test_bar_grouped_and_errors(frame):
    fig, ax = plt.subplots()
    plot_anno_bar(frame, groupby="group", ax=ax)
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["A", "B"]
    plt.close(fig)

    with pytest.raises(ValueError, match="no annotations"):
        plot_anno_bar(pd.DataFrame({"Feature": [], "Frequency": []}))
    with pytest.raises(KeyError, match="groupby column"):
        plot_anno_bar(frame, groupby="missing")
    with pytest.raises(KeyError, match="groupby column"):
        _grouped_stat(frame, "missing")
    with pytest.raises(KeyError, match="annotation"):
        _grouped_stat(pd.DataFrame({"group": ["A"]}), "group")


def test_distance_breaks_custom_and_formatting():
    arr, labels = _distance_breaks(None)
    assert arr[0] == 0
    assert np.isinf(arr[-1])
    assert len(labels) == len(arr) - 1

    arr, labels = _distance_breaks([1000, 0, 500])
    assert arr.tolist() == [0.0, 500.0, 1000.0, np.inf]
    assert labels[-1] == ">1kb"

    arr, _ = _distance_breaks([-5, 2500, float("inf")])
    assert arr[0] == 0
    assert np.isinf(arr[-1])

    _, labels = _distance_breaks([1500])
    assert labels[-1] == ">1.5kb"


def test_dist_to_tss_grouped_custom_breaks_and_errors(frame):
    fig, ax = plt.subplots()
    plot_dist_to_tss(
        frame,
        groupby="group",
        ax=ax,
        distance_breaks=[0, 100, 1000],
    )
    assert set(tick.get_text() for tick in ax.get_yticklabels()) == {"A", "B"}
    plt.close(fig)

    with pytest.raises(KeyError, match="distance-to-TSS"):
        plot_dist_to_tss(pd.DataFrame({"chrom": ["chr1"]}))
    with pytest.raises(KeyError, match="groupby column"):
        plot_dist_to_tss(frame, groupby="missing")
    with pytest.raises(ValueError, match="no finite distance"):
        plot_dist_to_tss(pd.DataFrame({"distance_to_tss": [np.nan]}))


def test_grouped_stat_defaults_missing_downstream_distance(annotation):
    frame = annotation.to_dataframe(include_peak_metadata=True)
    frame["group"] = "A"

    class Stub:
        downstream_distance = None

        def to_dataframe(self, include_peak_metadata=True):
            return frame

        def anno_stat(self, **kwargs):
            return pd.DataFrame({"Feature": ["Promoter"], "Frequency": [100.0]})

    stat = _grouped_stat(Stub(), "group")
    assert not stat.empty
    assert stat.index.tolist() == ["A"]


def test_dist_to_tss_falls_back_to_known_distance_columns(frame):
    fig, ax = plot_dist_to_tss(frame, distance_column="missing")
    assert fig is not None
    plt.close(fig)


def test_coerce_signal_frame_variants():
    renamed = _coerce_signal_frame(
        pd.DataFrame({"c": ["chr1"], "s": [0], "e": [10], "v": [2.0]})
    )
    assert renamed.columns.tolist() == ["chrom", "start", "end", "value"]
    assert renamed.loc[0, "value"] == 2.0

    with pytest.raises(ValueError, match="must contain"):
        _coerce_signal_frame(pd.DataFrame({"chrom": ["chr1"], "start": [0]}))
    with pytest.raises(TypeError, match="path or a pandas.DataFrame"):
        _coerce_signal_frame(["chr1", 0, 10, 1])
    with pytest.raises(ValueError, match="no valid intervals"):
        _coerce_signal_frame(
            pd.DataFrame(
                {"chrom": ["chr1"], "start": [10], "end": [10], "value": [1.0]}
            )
        )


def test_mean_in_interval_without_overlap():
    starts = np.array([10], dtype=np.int64)
    ends = np.array([10], dtype=np.int64)
    values = np.array([1.0])
    assert np.isnan(_mean_in_interval(starts, ends, values, 0, 5))
    # Overlapping search indices can still describe a zero-length interval.
    assert np.isnan(_mean_in_interval(starts, ends, values, 0, 20))


def test_extract_signal_validation_and_empty_intervals():
    genes = pd.DataFrame({"chrom": ["chr1"], "tss": [50], "strand": ["+"]})
    signal = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [100], "value": [2.0]}
    )
    with pytest.raises(ValueError, match="missing columns"):
        extract_tss_matrix_from_signal(signal, pd.DataFrame({"chrom": ["chr1"]}))
    with pytest.raises(ValueError, match="n_bins must be > 1"):
        extract_tss_matrix_from_signal(signal, genes, n_bins=1)
    with pytest.raises(ValueError, match="upstream/downstream"):
        _resolve_tss_window(0, 0, 0)

    matrix = extract_tss_matrix_from_signal(
        signal,
        genes,
        window=20,
        n_bins=4,
        chrom_sizes={"chr1": 1000},
        min_value=None,
    )
    assert matrix.shape == (1, 4)
    np.testing.assert_allclose(matrix[0], 2.0)

    # Unknown chromosome and missing chromosome length both leave NaNs.
    with pytest.raises(ValueError, match="no features could be extracted"):
        extract_tss_matrix_from_signal(
            signal,
            pd.DataFrame({"chrom": ["chr2"], "tss": [50], "strand": ["+"]}),
            window=20,
            n_bins=4,
            drop_empty=True,
        )
    with pytest.raises(ValueError, match="no features could be extracted"):
        extract_tss_matrix_from_signal(
            signal,
            genes,
            window=20,
            n_bins=4,
            chrom_sizes={},
            drop_empty=True,
        )

    # A valid chromosome with no signal in the window is also empty.
    gap = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [10], "value": [1.0]}
    )
    with pytest.raises(ValueError, match="no features could be extracted"):
        extract_tss_matrix_from_signal(
            gap,
            genes,
            window=20,
            n_bins=4,
            chrom_sizes={"chr1": 1000},
            drop_empty=True,
        )


class _IntervalBigWig:
    def __init__(self, path):
        self.path = path
        self.closed = False

    def chroms(self):
        return {"chr1": 1000, "chr2": 1000}

    def intervals(self, chrom):
        return [(0, 1000, 1.0)] if chrom == "chr1" else None

    def stats(self, chrom, start, end, nBins=1, type="mean"):
        return [1.0] * int(nBins)

    def close(self):
        self.closed = True


class _StatsOnlyBigWig:
    def __init__(self, path):
        self.path = path
        self.closed = False

    def chroms(self):
        return {"chr1": 1000}

    def stats(self, chrom, start, end, nBins=1, type="mean"):
        if start == 100:
            return None
        return [2.0] * int(nBins)

    def close(self):
        self.closed = True


def _install_pybigwig(monkeypatch, factory) -> None:
    fake = types.ModuleType("pyBigWig")
    fake.open = factory
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)


def test_extract_tss_matrix_validation(monkeypatch):
    _install_pybigwig(monkeypatch, _IntervalBigWig)
    genes = pd.DataFrame({"chrom": ["chr1"], "tss": [100], "strand": ["+"]})
    with pytest.raises(ValueError, match="missing columns"):
        extract_tss_matrix("x.bw", pd.DataFrame({"chrom": ["chr1"]}))
    with pytest.raises(ValueError, match="n_bins must be > 1"):
        extract_tss_matrix("x.bw", genes, n_bins=1)

    far = pd.DataFrame({"chrom": ["chr1"], "tss": [5000], "strand": ["+"]})
    with pytest.raises(ValueError, match="no features"):
        extract_tss_matrix("x.bw", far, window=100, n_bins=4)


def test_extract_tss_matrix_unopenable(monkeypatch):
    _install_pybigwig(monkeypatch, lambda path: None)
    genes = pd.DataFrame({"chrom": ["chr1"], "tss": [100], "strand": ["+"]})
    with pytest.raises(ValueError, match="could not open bigWig"):
        extract_tss_matrix("x.bw", genes, window=10, n_bins=2)


def test_extract_tss_matrix_interval_path_and_masks(monkeypatch):
    _install_pybigwig(monkeypatch, _IntervalBigWig)
    genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chr1"],
            "tss": [200, 200, 900],
            "strand": ["+", "+", "-"],
        }
    )
    matrix, mask = extract_tss_matrix(
        "x.bw",
        genes,
        window=100,
        n_bins=4,
        min_value=None,
        return_mask=True,
    )
    assert matrix.shape == (3, 4)
    np.testing.assert_allclose(matrix[0], 1.0)
    assert np.isnan(matrix[1]).all()
    np.testing.assert_allclose(matrix[2], 1.0)
    assert mask.tolist() == [True, False, True]

    dropped = extract_tss_matrix(
        "x.bw", genes, window=100, n_bins=4, min_value=None, drop_empty=True
    )
    assert dropped.shape == (2, 4)

    empty = pd.DataFrame(
        {"chrom": ["chr2", "chr2"], "tss": [200, 300], "strand": ["+", "+"]}
    )
    with pytest.raises(ValueError, match="no features"):
        extract_tss_matrix(
            "x.bw", empty, window=100, n_bins=4, drop_empty=True
        )


def test_extract_tss_matrix_parallel_assignment(monkeypatch):
    _install_pybigwig(monkeypatch, _IntervalBigWig)
    monkeypatch.setattr(plotting, "_TSS_MATRIX_MIN_FEATURES", 1)

    def fake_run_blocks(job, tasks, workers):
        rows = np.concatenate([task[1] for task in tasks])
        yield rows, np.full((rows.size, job.n_bins), 5.0)

    monkeypatch.setattr(plotting, "run_blocks", fake_run_blocks)
    genes = pd.DataFrame(
        {"chrom": ["chr1", "chr1"], "tss": [200, 400], "strand": ["+", "-"]}
    )
    matrix = extract_tss_matrix(
        "x.bw", genes, window=100, n_bins=4, n_jobs=2
    )
    np.testing.assert_allclose(matrix, 5.0)


def test_extract_tss_matrix_stats_fallback(monkeypatch):
    _install_pybigwig(monkeypatch, _StatsOnlyBigWig)
    genes = pd.DataFrame(
        {"chrom": ["chr1", "chr1"], "tss": [200, 300], "strand": ["+", "-"]}
    )
    matrix, mask = extract_tss_matrix(
        "x.bw",
        genes,
        window=100,
        n_bins=4,
        return_mask=True,
    )
    assert np.isnan(matrix[0]).all()
    np.testing.assert_allclose(matrix[1], 2.0)
    assert mask.tolist() == [False, True]

    dropped = extract_tss_matrix(
        "x.bw", genes, window=100, n_bins=4, min_value=None, drop_empty=True
    )
    assert dropped.shape == (1, 4)


def test_tss_plot_validation_and_external_axes():
    matrix = np.ones((3, 5), dtype=float)
    with pytest.raises(ValueError, match="2-dimensional"):
        plot_tss_profile(np.ones(5))
    with pytest.raises(ValueError, match="no rows"):
        plot_tss_profile(np.empty((0, 5)))
    with pytest.raises(ValueError, match="2-dimensional"):
        plot_tss_heatmap(np.ones(5))
    with pytest.raises(ValueError, match="no rows"):
        plot_tss_heatmap(np.empty((0, 5)))

    fig, ax = plt.subplots()
    returned_fig, returned_ax = plot_tss_profile(matrix, ax=ax, show_sem=False)
    assert returned_fig is fig
    assert returned_ax is ax
    plt.close(fig)

    fig, ax = plt.subplots()
    returned_fig, returned_ax = plot_tss_heatmap(matrix, ax=ax)
    assert returned_fig is fig
    assert returned_ax is ax
    plt.close(fig)


def test_tss_heatmap_sorts_rows():
    matrix = np.zeros((3, 5), dtype=float)
    matrix[:, 2] = [1.0, 3.0, 2.0]
    fig, ax = plot_tss_heatmap(matrix, sort_rows=True)
    drawn = np.asarray(ax.images[0].get_array(), dtype=float)
    assert drawn[0, 2] == 3.0
    assert drawn[-1, 2] == 1.0
    plt.close(fig)
