"""Validation, dispatch and optional-backend branches for tag-matrix helpers.

These tests concentrate on the small guard clauses that keep the public
functions predictable: region construction, window normalisation, explicit
error messages, bigWig/BAM dispatch and the plotting preconditions.  The
numeric behaviour of the coverage builder is already covered by
``tests/test_tagmatrix.py``.
"""

from __future__ import annotations

import sys
import types

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pychipseeker.tagmatrix as tagmatrix
from pychipseeker.annotation import GeneAnnotation, read_gtf
from pychipseeker.tagmatrix import (
    TagMatrix,
    _as_matrix,
    _bam_matrix,
    _bam_matrix_rows,
    _bigwig_matrix,
    _bin_signal,
    _coverage_matrix,
    _coverage_parallel_floor,
    _coverage_rows,
    _coverage_track,
    _coverage_worker,
    _coverage_worker_init,
    _feature_table,
    _make_region_table,
    _normalise_region_type,
    _peak_arrays,
    _region_label,
    _resolve_annotation,
    _window_frame,
    _x_coordinates,
    get_bio_region,
    get_tag_matrix,
    make_bio_region_from_granges,
    plot_avg_prof,
    plot_multi_prof,
    plot_peak_prof,
    tag_heatmap,
)


@pytest.fixture
def model(examples_dir) -> GeneAnnotation:
    return GeneAnnotation.from_gtf(examples_dir / "example.gtf")


@pytest.fixture
def windows() -> pd.DataFrame:
    return pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10], "strand": ["+"]})


def test_tag_matrix_metadata_helpers():
    matrix = TagMatrix(np.zeros((2, 3)), upstream=10, downstream=5)
    assert isinstance(matrix.matrix, np.ndarray)
    assert not isinstance(matrix.matrix, TagMatrix)
    assert matrix.__array_finalize__(None) is None

    with pytest.raises(ValueError, match="2-dimensional"):
        _as_matrix(np.zeros(3))
    assert _x_coordinates(3, None, None).tolist() == [0.0, 1.0, 2.0]
    assert _x_coordinates(3, 2, 4).tolist() == [-2.0, 1.0, 4.0]


def test_resolve_annotation_from_paths_and_frames(examples_dir, monkeypatch):
    with pytest.raises(ValueError, match="provide an annotation"):
        _resolve_annotation()

    from_path = _resolve_annotation(gtf=examples_dir / "example.gtf")
    assert isinstance(from_path, GeneAnnotation)

    from_frame = _resolve_annotation(annotation=read_gtf(examples_dir / "example.gtf"))
    assert len(from_frame.genes) == len(from_path.genes)

    with pytest.raises(TypeError, match="annotation must be"):
        _resolve_annotation(annotation=object())

    import pychipseeker.genomes as genomes

    sentinel = object()
    monkeypatch.setattr(genomes, "load_genome", lambda name: sentinel)
    assert _resolve_annotation(annotation="hg38") is sentinel

    spec = "ensembl:test"
    monkeypatch.setattr(genomes, "is_genome_spec", lambda value: str(value) == spec)
    monkeypatch.setattr(genomes, "load_annotation", lambda value: sentinel)
    assert _resolve_annotation(annotation=spec) is sentinel


def test_feature_table_aliases(model):
    assert not _feature_table(model, "exon").empty
    assert not _feature_table(model, "intron").empty
    for by in ("5UTR", "five_prime_utr", "five-utr"):
        assert "utr_type" not in _feature_table(model, by).columns
    for by in ("3UTR", "three_prime_utr", "three-utr"):
        assert "utr_type" not in _feature_table(model, by).columns

    utr = _feature_table(model, "utr")
    assert "utr_type" in utr.columns

    with pytest.raises(ValueError, match="by must be one of"):
        _feature_table(model, "promoter")


def test_region_type_and_label_helpers():
    assert _normalise_region_type("start") == "start_site"
    for value in ("end_site", "end", "tts", "END"):
        assert _normalise_region_type(value) == "end_site"
    assert _normalise_region_type("body") == "body"
    with pytest.raises(ValueError, match="type must be"):
        _normalise_region_type("middle")

    assert _region_label("end_site", "gene") == "TES"
    assert _region_label("end_site", "exon") == "End"
    assert _region_label("body", "exon") == "exon body"


def test_make_region_table_empty_and_missing_strand():
    empty = _make_region_table(
        pd.DataFrame(columns=["chrom", "start", "end"]),
        region_type="start_site",
        upstream=10,
        downstream=20,
        by="gene",
    )
    assert empty.empty
    assert empty.attrs["label"] == "TSS"
    assert empty.attrs["upstream"] == 10

    body = _make_region_table(
        pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [20]}),
        region_type="body",
        upstream=0,
        downstream=0,
        by="gene",
    )
    assert body["strand"].tolist() == ["*"]


def test_make_region_table_end_site_coordinates():
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [100, 300, 500],
            "end": [200, 400, 600],
            "strand": ["+", "-", "."],
        }
    )
    out = _make_region_table(
        frame,
        region_type="end_site",
        upstream=10,
        downstream=5,
        by="gene",
    ).set_index("strand")

    assert (out.loc["+", "start"], out.loc["+", "end"]) == (190, 205)
    assert (out.loc["-", "start"], out.loc["-", "end"]) == (296, 311)
    assert (out.loc["*", "start"], out.loc["*", "end"]) == (590, 605)
    assert out.attrs["label"] == "TES"


def test_make_bio_region_from_granges_validation():
    with pytest.raises(TypeError, match="DataFrame"):
        make_bio_region_from_granges([1, 2], by="gene", type="start_site")
    with pytest.raises(ValueError, match="missing columns"):
        make_bio_region_from_granges(
            pd.DataFrame({"chrom": ["chr1"]}), by="gene", type="start_site"
        )

    frame = pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [200], "strand": ["+"]}
    )
    out = make_bio_region_from_granges(
        frame, by="gene", type="end_site", upstream=5, downstream=5
    )
    assert (out.loc[0, "start"], out.loc[0, "end"]) == (195, 205)
    assert out.attrs["label"] == "TES"


def test_get_bio_region_uses_resolved_annotation(model):
    body = get_bio_region(annotation=model, type="body", upstream=0, downstream=0)
    assert body.attrs["label"] == "gene body"
    with pytest.raises(ValueError, match="type must be"):
        get_bio_region(annotation=model, type="middle")


def test_window_frame_path_and_validation(tmp_path):
    frame = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    path = tmp_path / "windows.tsv"
    frame.to_csv(path, sep="\t", index=False)

    loaded = _window_frame(path)
    assert loaded["strand"].tolist() == ["*"]
    assert _window_frame(frame)["strand"].tolist() == ["*"]

    from_file = _window_frame(
        pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10], "strand": ["?"]})
    )
    assert from_file["strand"].tolist() == ["*"]

    with pytest.raises(TypeError, match="DataFrame or a tab-delimited path"):
        _window_frame(["chr1", 0, 10])
    with pytest.raises(ValueError, match="windows is missing columns"):
        _window_frame(pd.DataFrame({"chrom": ["chr1"]}))


def test_coverage_track_and_row_edge_cases(windows):
    assert (
        _coverage_track(
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=float),
        )
        is None
    )
    assert (
        _coverage_track(
            np.array([10], dtype=np.int64),
            np.array([10], dtype=np.int64),
            np.array([1.0]),
        )
        is None
    )

    empty = pd.DataFrame({"chrom": [], "start": [], "end": [], "strand": []})
    matrix, kept = _coverage_rows({}, empty, n_bins=4, ignore_strand=False)
    assert matrix.shape == (0, 4)
    assert kept.empty

    one_base = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [1], "strand": ["+"]}
    )
    with pytest.raises(ValueError, match="nbin must be > 1"):
        _coverage_rows({}, one_base, n_bins=None, ignore_strand=False)
    with pytest.raises(ValueError, match="nbin must be > 1"):
        _coverage_rows({}, one_base, n_bins=1, ignore_strand=False)

    # Peaks with zero-length intervals produce no coverage track at all.
    degenerate = pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [10]})
    matrix, kept = _coverage_rows(
        _peak_arrays(degenerate), windows, n_bins=4, ignore_strand=False
    )
    assert matrix.size == 0
    assert kept.empty

    # A valid track outside the window leaves no non-zero row.
    elsewhere = pd.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200]})
    far_window = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [10], "strand": ["+"]}
    )
    matrix, kept = _coverage_rows(
        _peak_arrays(elsewhere), far_window, n_bins=4, ignore_strand=False
    )
    assert matrix.size == 0
    assert kept.empty


def test_coverage_matrix_validation(windows):
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    empty = windows.iloc[0:0]
    with pytest.raises(ValueError, match="windows is empty"):
        _coverage_matrix(peaks, empty, n_bins=4, ignore_strand=False)

    unequal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 10],
            "end": [100, 200],
            "strand": ["+", "+"],
        }
    )
    with pytest.raises(ValueError, match="equal widths"):
        _coverage_matrix(peaks, unequal, n_bins=None, ignore_strand=False)

    far = pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [110], "strand": ["+"]}
    )
    with pytest.raises(ValueError, match="no windows overlap"):
        _coverage_matrix(peaks, far, n_bins=4, ignore_strand=False)


def test_coverage_matrix_parallel_without_matches(monkeypatch):
    monkeypatch.setattr(tagmatrix, "_coverage_parallel_floor", lambda: 1)

    def fake_parallel_map(func, items, *, n_jobs, initializer=None, initargs=()):
        return [
            (np.empty((0, 4)), pd.DataFrame(columns=["chrom", "start", "end"]))
            for _ in items
        ]

    monkeypatch.setattr(tagmatrix, "parallel_map", fake_parallel_map)
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    windows = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [1000, 2000],
            "end": [1010, 2010],
            "strand": ["+", "+"],
        }
    )
    with pytest.raises(ValueError, match="no windows overlap any peak"):
        _coverage_matrix(peaks, windows, n_bins=4, ignore_strand=False, n_jobs=2)


def test_coverage_parallel_floor_and_worker(monkeypatch, windows):
    monkeypatch.setattr(tagmatrix, "process_start_method", lambda: "spawn")
    assert _coverage_parallel_floor() == tagmatrix._COVERAGE_SPAWN_MIN_WINDOWS
    monkeypatch.setattr(tagmatrix, "process_start_method", lambda: "fork")
    assert _coverage_parallel_floor() == tagmatrix._COVERAGE_MIN_WINDOWS

    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    _coverage_worker_init(_peak_arrays(peaks), 2, False)
    try:
        matrix, kept = _coverage_worker(windows)
    finally:
        tagmatrix._COVERAGE_WORKER = None
    assert matrix.shape == (1, 2)
    assert kept["start"].tolist() == [0]


def _install_fake_pybigwig(monkeypatch, opener):
    fake = types.ModuleType("pyBigWig")
    fake.open = opener
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)


class _FakeBigWig:
    def __init__(self, path):
        self.path = path
        self.closed = False

    def chroms(self):
        return {"chr1": 1000, "chr2": 500}

    def stats(self, chrom, start, end, nBins=1, type="mean"):
        return [1.0] * int(nBins)

    def intervals(self, chrom):
        if chrom == "chr2":
            return None
        return [(0, 1000, 2.0)]

    def close(self):
        self.closed = True


def test_bigwig_matrix_validation(monkeypatch, windows):
    _install_fake_pybigwig(monkeypatch, _FakeBigWig)
    with pytest.raises(ValueError, match="unsupported bigWig value_type"):
        _bigwig_matrix(
            "x.bw", windows, n_bins=2, ignore_strand=False, value_type="median"
        )
    with pytest.raises(ValueError, match="windows is empty"):
        _bigwig_matrix("x.bw", windows.iloc[0:0], n_bins=2, ignore_strand=False)

    unequal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 0],
            "end": [10, 20],
            "strand": ["+", "+"],
        }
    )
    with pytest.raises(ValueError, match="equal widths"):
        _bigwig_matrix("x.bw", unequal, n_bins=None, ignore_strand=False)
    with pytest.raises(ValueError, match="nbin must be > 1"):
        _bigwig_matrix("x.bw", windows, n_bins=1, ignore_strand=False)

    outside = pd.DataFrame(
        {"chrom": ["chr1"], "start": [990], "end": [1200], "strand": ["+"]}
    )
    with pytest.raises(ValueError, match="no windows overlap"):
        _bigwig_matrix("x.bw", outside, n_bins=2, ignore_strand=False)


def test_bigwig_matrix_unopenable(monkeypatch, windows):
    _install_fake_pybigwig(monkeypatch, lambda path: None)
    with pytest.raises(ValueError, match="could not open bigWig"):
        _bigwig_matrix("x.bw", windows, n_bins=2, ignore_strand=False)


def test_bigwig_matrix_vectorised_and_missing_intervals(monkeypatch):
    _install_fake_pybigwig(monkeypatch, _FakeBigWig)
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [10, 10],
            "strand": ["+", "-"],
        }
    )
    matrix, kept, is_binning = _bigwig_matrix(
        "x.bw", frame, n_bins=2, ignore_strand=False, value_type="mean"
    )
    assert is_binning is True
    assert kept["chrom"].tolist() == ["chr1", "chr2"]
    np.testing.assert_allclose(matrix[0], 2.0)
    np.testing.assert_allclose(matrix[1], 0.0)


def test_bigwig_matrix_without_explicit_bins(monkeypatch):
    _install_fake_pybigwig(monkeypatch, _FakeBigWig)
    frame = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [10], "strand": ["+"]}
    )
    matrix, kept, is_binning = _bigwig_matrix(
        "x.bw", frame, n_bins=None, ignore_strand=False, value_type="mean"
    )
    assert is_binning is False
    assert matrix.shape == (1, 10)
    assert len(kept) == 1
    np.testing.assert_allclose(matrix, 2.0)


def test_bigwig_matrix_parallel_assignment(monkeypatch):
    _install_fake_pybigwig(monkeypatch, _FakeBigWig)
    monkeypatch.setattr(tagmatrix, "_BIGWIG_MATRIX_MIN_WINDOWS", 1)

    def fake_run_blocks(job, tasks, workers):
        rows = np.concatenate([task[1] for task in tasks])
        yield rows, np.full((rows.size, job.n_bins), 3.0)

    monkeypatch.setattr(tagmatrix, "run_blocks", fake_run_blocks)
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 10],
            "end": [10, 20],
            "strand": ["+", "-"],
        }
    )
    matrix, kept, _ = _bigwig_matrix(
        "x.bw", frame, n_bins=2, ignore_strand=False, n_jobs=2
    )
    assert matrix.shape == (2, 2)
    np.testing.assert_allclose(matrix, 3.0)
    assert len(kept) == 2


def test_bin_signal_value_types_and_validation():
    values = np.array([0.0, 2.0, 4.0, 0.0])
    with pytest.raises(ValueError, match="coverage array is empty"):
        _bin_signal(np.array([]), 2, "mean")

    np.testing.assert_allclose(_bin_signal(values, None, "mean"), values)
    np.testing.assert_allclose(_bin_signal(values, len(values), "sum"), values)
    with pytest.raises(ValueError, match="nbin must be > 1"):
        _bin_signal(values, 1, "mean")

    np.testing.assert_allclose(_bin_signal(values, 2, "sum"), [2.0, 4.0])
    np.testing.assert_allclose(_bin_signal(values, 2, "max"), [2.0, 4.0])
    np.testing.assert_allclose(_bin_signal(values, 2, "coverage"), [0.5, 0.5])
    np.testing.assert_allclose(
        _bin_signal(values, 5, "mean"), [0.0, 0.0, 2.0, 4.0, 0.0]
    )


class _FakeAlignmentFile:
    references = ["chr1", "chr2"]

    def __init__(self, path, mode="rb"):
        self.path = path
        self.mode = mode
        self.closed = False

    def count_coverage(
        self, chrom, start, end, quality_threshold=0, read_callback=None
    ):
        width = end - start
        if start == 5:
            raise ValueError("bad region")
        if chrom == "chr2":
            return ()
        if start == 7:
            return tuple(np.zeros(0) for _ in range(4))
        values = np.ones(width)
        return (values, np.zeros(width), np.zeros(width), np.zeros(width))

    def close(self):
        self.closed = True


def _install_fake_pysam(monkeypatch):
    fake = types.ModuleType("pysam")
    fake.AlignmentFile = _FakeAlignmentFile
    monkeypatch.setitem(sys.modules, "pysam", fake)


def test_bam_matrix_rows_filters_bad_windows(monkeypatch):
    _install_fake_pysam(monkeypatch)
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chrZ", "chr1"],
            "start": [0, 5, 7, 0, 0, 0],
            "end": [10, 10, 10, 10, 10, 10],
            "strand": ["+", "+", "+", "+", "+", "-"],
        }
    )
    matrix, kept = _bam_matrix_rows(
        "x.bam",
        frame,
        n_bins=5,
        ignore_strand=False,
        value_type="mean",
        min_mapping_quality=0,
        exclude_duplicates=True,
    )
    assert kept["start"].tolist() == [0, 0]
    assert kept["strand"].tolist() == ["+", "-"]
    assert matrix.shape == (2, 5)
    np.testing.assert_allclose(matrix, 1.0)


def test_bam_matrix_rows_without_usable_windows(monkeypatch):
    _install_fake_pysam(monkeypatch)
    frame = pd.DataFrame(
        {
            "chrom": ["chr3", "chr1"],
            "start": [-1, 20],
            "end": [10, 10],
            "strand": ["+", "+"],
        }
    )
    matrix, kept = _bam_matrix_rows(
        "x.bam",
        frame,
        n_bins=4,
        ignore_strand=False,
        value_type="mean",
        min_mapping_quality=0,
        exclude_duplicates=True,
    )
    assert matrix.shape == (0, 4)
    assert kept.empty


def test_bam_matrix_validation(monkeypatch, windows):
    _install_fake_pysam(monkeypatch)
    with pytest.raises(ValueError, match="BAM value_type"):
        _bam_matrix(
            "x.bam", windows, n_bins=2, ignore_strand=False, value_type="median"
        )
    with pytest.raises(ValueError, match="no windows overlap"):
        _bam_matrix("x.bam", windows.iloc[0:0], n_bins=2, ignore_strand=False)

    unequal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 0],
            "end": [10, 6],
            "strand": ["+", "+"],
        }
    )
    with pytest.raises(ValueError, match="equal widths"):
        _bam_matrix("x.bam", unequal, n_bins=None, ignore_strand=False)


def test_bam_matrix_parallel_dispatch(monkeypatch):
    _install_fake_pysam(monkeypatch)
    monkeypatch.setattr(tagmatrix, "_BAM_MATRIX_MIN_WINDOWS", 1)

    def fake_parallel_map(func, items, *, n_jobs, initializer=None, initargs=()):
        if initializer is not None:
            initializer(*initargs)
        return [func(item) for item in items]

    monkeypatch.setattr(tagmatrix, "parallel_map", fake_parallel_map)
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 0],
            "end": [10, 10],
            "strand": ["+", "-"],
        }
    )
    matrix, kept, is_binning = _bam_matrix(
        "x.bam", frame, n_bins=2, ignore_strand=False, n_jobs=2
    )
    assert matrix.shape == (2, 2)
    assert kept["strand"].tolist() == ["+", "-"]
    assert is_binning is True

    unusable = pd.DataFrame(
        {
            "chrom": ["chr3", "chr3"],
            "start": [0, 10],
            "end": [10, 20],
            "strand": ["+", "+"],
        }
    )
    with pytest.raises(ValueError, match="no windows overlap"):
        _bam_matrix("x.bam", unusable, n_bins=2, ignore_strand=False, n_jobs=2)


def test_get_tag_matrix_detects_signal_paths(monkeypatch, windows):
    _install_fake_pybigwig(monkeypatch, _FakeBigWig)
    from_bw = get_tag_matrix(
        peak="signal.bw", windows=windows, nbin=2, value_type="mean"
    )
    assert from_bw.shape == (1, 2)
    np.testing.assert_allclose(from_bw, 2.0)

    _install_fake_pysam(monkeypatch)
    for suffix in ("bam", "cram"):
        from_alignment = get_tag_matrix(
            peak=f"reads.{suffix}", windows=windows, nbin=2, value_type="mean"
        )
        assert from_alignment.shape == (1, 2)
        np.testing.assert_allclose(from_alignment, 1.0)


def test_get_tag_matrix_dispatch_errors(windows):
    with pytest.raises(ValueError, match="windows is empty"):
        get_tag_matrix(peak=None, windows=[], nbin=2)

    with pytest.raises(ValueError, match="supply only one"):
        get_tag_matrix(windows=windows, bigwig="x.bw", bam="x.bam", nbin=2)

    with pytest.raises(ValueError, match="provide peak intervals"):
        get_tag_matrix(windows=windows, nbin=2)

    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    with pytest.raises(KeyError, match="weight column not found"):
        get_tag_matrix(peak=peaks, windows=windows, weightCol="missing", nbin=2)


def test_get_tag_matrix_defaults_and_verbose(model, capsys):
    peaks = pd.DataFrame(
        {"chrom": ["chr1", "chr2"], "start": [1000, 5000], "end": [1100, 5100]}
    )
    matrix = get_tag_matrix(peak=peaks, annotation=model, nbin=5, verbose=True)
    assert matrix.upstream == 3000
    assert matrix.downstream == 3000
    assert "tag matrix:" in capsys.readouterr().out


def test_plot_avg_prof_uses_xlim():
    matrix = np.ones((2, 5), dtype=float)
    fig, ax = plot_avg_prof(matrix, xlim=[-2, 3])
    x = ax.lines[0].get_xdata()
    assert x[0] == -2
    assert x[-1] == 3
    plt.close(fig)


def test_plot_peak_prof_builds_from_peaks(model):
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [1000, 6900],
            "end": [1100, 7000],
        }
    )
    fig, ax = plot_peak_prof(
        peak=peaks,
        annotation=model,
        upstream=100,
        downstream=100,
        nbin=10,
    )
    assert "Genomic Region" in ax.get_xlabel()
    plt.close(fig)


def test_plot_preconditions_and_list_dispatch():
    matrix = np.ones((2, 5), dtype=float)
    with pytest.raises(ValueError, match="tagMatrix has no rows"):
        plot_avg_prof(np.zeros((0, 5)))
    with pytest.raises(ValueError, match="tagMatrix has no rows"):
        tag_heatmap(np.zeros((0, 5)))
    with pytest.raises(ValueError, match="xlim must be two values"):
        plot_avg_prof(matrix, xlim=[1])
    with pytest.raises(ValueError, match="conf must be between"):
        plot_avg_prof(matrix, conf=1.5)

    fig, axes = plot_peak_prof([matrix, matrix], facet="row")
    assert axes is not None
    plt.close(fig)

    with pytest.raises(ValueError, match="provide tagMatrix or peak"):
        plot_peak_prof()
    with pytest.raises(ValueError, match="tagMatrices is empty"):
        plot_multi_prof([])
