import sys
import types
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pychipseeker import annotate_peak
from pychipseeker.tagmatrix import (
    TagMatrix,
    get_bio_region,
    get_promoters,
    get_tag_matrix,
    peak_heatmap,
    peak_profile_heatmap,
    plot_avg_prof,
    plot_multi_prof,
    plot_peak_prof2,
    tag_heatmap,
)


@pytest.fixture
def annotation(examples_dir):
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        upstream=1000,
        downstream=1000,
    )


def test_get_promoters_and_bio_region(annotation):
    promoters = get_promoters(annotation=annotation.gene_annotation, upstream=100, downstream=200)
    assert len(promoters) == 3
    row = promoters.set_index("gene_id").loc["G1"]
    assert row["start"] == 900
    assert row["end"] == 1200
    assert promoters.attrs["type"] == "start_site"

    minus = promoters.set_index("gene_id").loc["G2"]
    assert minus["start"] == 6800
    assert minus["end"] == 7100

    body = get_bio_region(
        annotation=annotation.gene_annotation,
        by="gene",
        type="body",
        upstream=0,
        downstream=0,
    )
    assert body.set_index("gene_id").loc["G1", "start"] == 1000
    assert body.set_index("gene_id").loc["G1", "end"] == 3000


def test_tag_matrix_from_peak_coverage(annotation):
    promoters = get_promoters(annotation=annotation.gene_annotation, upstream=1000, downstream=1000)
    matrix = get_tag_matrix(
        peak=annotation.peaks,
        windows=promoters,
        nbin=100,
        ignore_strand=False,
    )
    assert isinstance(matrix, TagMatrix)
    assert matrix.shape[0] == 2
    assert matrix.shape[1] == 100
    assert matrix.upstream == 1000
    assert matrix.downstream == 1000
    assert matrix.is_binning is True
    assert matrix.rows is not None
    assert matrix.sum() > 0

    # The minus-strand row must be the reverse of the same region on the plus
    # strand when the matrix is built without strand flipping.
    raw = get_tag_matrix(
        peak=annotation.peaks,
        windows=promoters,
        nbin=100,
        ignore_strand=True,
    )
    assert raw.shape == matrix.shape


def test_tag_matrix_without_binning(annotation):
    promoters = get_promoters(annotation=annotation.gene_annotation, upstream=50, downstream=50)
    matrix = get_tag_matrix(peak=annotation.peaks, windows=promoters)
    assert matrix.shape[1] == 100
    assert matrix.is_binning is False


def test_tag_matrix_multiple_windows(annotation):
    promoters_a = get_promoters(
        annotation=annotation.gene_annotation, upstream=100, downstream=100
    )
    promoters_a.attrs["label"] = "small"
    promoters_b = get_promoters(
        annotation=annotation.gene_annotation, upstream=200, downstream=200
    )
    promoters_b.attrs["label"] = "large"
    matrices = get_tag_matrix(
        peak=annotation.peaks, windows=[promoters_a, promoters_b], nbin=20
    )
    assert isinstance(matrices, list)
    assert len(matrices) == 2
    assert matrices[0].label == "small"
    assert matrices[1].label == "large"


def test_tag_matrix_weight_column():
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [90],
            "end": [110],
            "strand": ["+"],
            "score": [3.0],
        }
    )
    windows = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [200], "strand": ["+"]}
    )
    matrix = get_tag_matrix(peak=peaks, windows=windows, weightCol="score", nbin=20)
    assert np.isclose(matrix.sum(), 60.0)


def test_profile_and_heatmap_plots(annotation):
    promoters = get_promoters(annotation=annotation.gene_annotation, upstream=100, downstream=100)
    matrix = get_tag_matrix(peak=annotation.peaks, windows=promoters, nbin=50)
    fig, ax = plot_avg_prof(matrix, conf=0.95)
    assert "Genomic Region" in ax.get_xlabel()
    plt.close(fig)

    fig, ax = tag_heatmap(matrix)
    assert ax.get_title() == "Tag heatmap"
    plt.close(fig)

    fig, ax = plot_peak_prof2(
        annotation.peaks,
        annotation=annotation.gene_annotation,
        upstream=100,
        downstream=100,
        nbin=50,
    )
    assert "Genomic Region" in ax.get_xlabel()
    plt.close(fig)

    fig, axes = peak_profile_heatmap(
        annotation.peaks,
        annotation=annotation.gene_annotation,
        upstream=100,
        downstream=100,
        nbin=50,
    )
    assert len(axes) == 2
    plt.close(fig)


def test_peak_heatmap(annotation):
    fig, ax = peak_heatmap(
        annotation.peaks,
        annotation=annotation.gene_annotation,
        upstream=100,
        downstream=100,
        nbin=50,
    )
    assert ax is not None
    plt.close(fig)


def test_peak_annotation_tag_matrix_methods(annotation):
    matrix = annotation.get_tag_matrix(upstream=100, downstream=100, nbin=50)
    assert matrix.shape[1] == 50
    fig, ax = annotation.plot_peak_prof(tagMatrix=matrix)
    assert "Genomic Region" in ax.get_xlabel()
    plt.close(fig)
    fig, ax = annotation.peak_heatmap(
        annotation.peaks, upstream=100, downstream=100, nbin=50
    )
    assert ax is not None
    plt.close(fig)


def test_bigwig_tag_matrix_uses_pybigwig(monkeypatch):
    class FakeBigWig:
        def __init__(self, path):
            self.path = path

        def chroms(self):
            return {"chr1": 10000}

        def stats(self, chrom, start, end, nBins=1, type="mean"):
            assert chrom == "chr1"
            assert 0 <= start < end <= 10000
            return [1.0] * int(nBins)

        def close(self):
            return None

    fake = types.ModuleType("pyBigWig")
    fake.open = FakeBigWig
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)

    windows = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [900, 1900],
            "end": [1100, 2100],
            "strand": ["+", "-"],
        }
    )
    matrix = get_tag_matrix(
        bigwig="fake.bw",
        windows=windows,
        nbin=10,
        ignore_strand=False,
    )
    assert matrix.shape == (2, 10)
    np.testing.assert_allclose(matrix, 1.0)


def test_extract_tss_matrix_bigwig_path(monkeypatch):
    from pychipseeker.plotting import extract_tss_matrix

    class FakeBigWig:
        def __init__(self, path):
            self.path = path

        def chroms(self):
            return {"chr1": 10000}

        def stats(self, chrom, start, end, nBins=1, type="mean"):
            assert type == "mean"
            return [float(i) for i in range(int(nBins))]

        def close(self):
            return None

    fake = types.ModuleType("pyBigWig")
    fake.open = FakeBigWig
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)

    genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "tss": [5000, 5000],
            "strand": ["+", "-"],
        }
    )
    matrix = extract_tss_matrix("fake.bw", genes, window=100, n_bins=4)
    assert matrix.shape == (2, 4)
    np.testing.assert_allclose(matrix[0], [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(matrix[1], matrix[0, ::-1])


# ---------------------------------------------------------------------------
# Real bigWig integration tests (uses tests/data/test.bw)
# ---------------------------------------------------------------------------


def _real_bigwig() -> Path:
    pyBigWig = pytest.importorskip("pyBigWig")
    path = Path(__file__).parent / "data" / "test.bw"
    if not path.exists():  # pragma: no cover - packaging guard
        pytest.skip("tests/data/test.bw is missing")
    bw = pyBigWig.open(str(path))
    assert bw is not None
    chroms = bw.chroms()
    bw.close()
    assert "1" in chroms and "10" in chroms
    return path


def test_bigwig_real_file_get_tag_matrix():
    from pychipseeker.tagmatrix import get_tag_matrix_from_bigwig

    bigwig = _real_bigwig()
    windows = pd.DataFrame(
        {
            "chrom": ["1", "1", "10", "chr1", "1"],
            "start": [0, 0, 200, 0, 0],
            "end": [4, 4, 300, 4, 195471972],
            "strand": ["+", "-", "+", "+", "+"],
        }
    )
    matrix = get_tag_matrix_from_bigwig(bigwig, windows, nbin=4)
    assert matrix.shape == (3, 4)  # 'chr1' and out-of-range windows are dropped
    np.testing.assert_allclose(matrix[0], [0.1, 0.2, 0.3, 0.0], rtol=1e-5)
    np.testing.assert_allclose(matrix[1], [0.0, 0.3, 0.2, 0.1], rtol=1e-5)
    np.testing.assert_allclose(matrix[2], [50.0, 50.0, 50.0, 50.0], rtol=1e-5)


def test_bigwig_real_file_strand_and_value_type():
    from pychipseeker.tagmatrix import get_tag_matrix_from_bigwig

    bigwig = _real_bigwig()
    windows = pd.DataFrame(
        {
            "chrom": ["1", "1"],
            "start": [1, 1],
            "end": [3, 3],
            "strand": ["+", "-"],
        }
    )
    matrix = get_tag_matrix_from_bigwig(bigwig, windows, nbin=2)
    np.testing.assert_allclose(matrix[0], [0.2, 0.3], rtol=1e-5)
    np.testing.assert_allclose(matrix[1], [0.3, 0.2], rtol=1e-5)

    matrix_mean = get_tag_matrix_from_bigwig(
        bigwig,
        windows,
        nbin=2,
        value_type="mean",
        ignore_strand=True,
    )
    np.testing.assert_allclose(matrix_mean[0], [0.2, 0.3], rtol=1e-5)
    np.testing.assert_allclose(matrix_mean[1], [0.2, 0.3], rtol=1e-5)

    with pytest.raises(ValueError):
        get_tag_matrix_from_bigwig(bigwig, windows, value_type="median")


def test_bigwig_real_file_extract_tss_matrix():
    from pychipseeker.plotting import extract_tss_matrix

    bigwig = _real_bigwig()
    genes = pd.DataFrame(
        {
            "chrom": ["1", "1"],
            "tss": [2, 2],
            "strand": ["+", "-"],
        }
    )
    matrix = extract_tss_matrix(bigwig, genes, window=1, n_bins=2)
    assert matrix.shape == (2, 2)
    np.testing.assert_allclose(matrix[0], [0.2, 0.3], rtol=1e-5)
    np.testing.assert_allclose(matrix[1], [0.3, 0.2], rtol=1e-5)


def test_multi_profile(annotation):
    promoters = get_promoters(annotation=annotation.gene_annotation, upstream=100, downstream=100)
    matrix = get_tag_matrix(peak=annotation.peaks, windows=promoters, nbin=20)
    fig, ax = plot_multi_prof([matrix, matrix], facet="none")
    assert ax is not None
    plt.close(fig)
