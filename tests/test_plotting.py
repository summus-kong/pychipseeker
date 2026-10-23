import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pychipseeker import annotate_peak
from pychipseeker.plotting import (
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


def test_plot_anno_pie(annotation):
    fig, ax = plot_anno_pie(annotation)
    assert fig is not None
    assert ax is not None
    plt.close(fig)


def test_plot_anno_bar(annotation):
    fig, ax = plot_anno_bar(annotation)
    assert fig is not None
    assert ax.get_xlabel() == "Percentage(%)"
    plt.close(fig)


def test_plot_dist_to_tss(annotation):
    fig, ax = plot_dist_to_tss(annotation)
    assert fig is not None
    assert "Binding sites" in ax.get_xlabel()
    plt.close(fig)


def test_plot_tss_profile_and_heatmap():
    rng = np.random.default_rng(0)
    matrix = rng.normal(size=(20, 50))
    fig, ax = plot_tss_profile(matrix, window=1000)
    assert fig is not None
    assert "Distance to TSS" in ax.get_xlabel()
    plt.close(fig)

    fig, ax = plot_tss_heatmap(matrix, window=1000)
    assert fig is not None
    plt.close(fig)


def test_extract_tss_matrix_from_signal():
    signal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "start": [0, 100, 200, 300],
            "end": [100, 200, 300, 400],
            "value": [1.0, 2.0, 3.0, 4.0],
        }
    )
    genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "tss": [200, 200],
            "strand": ["+", "-"],
        }
    )
    matrix = extract_tss_matrix_from_signal(
        signal, genes, window=100, n_bins=4, min_value=None
    )
    assert matrix.shape == (2, 4)
    # Bins are [100,150), [150,200), [200,250), [250,300).
    np.testing.assert_allclose(matrix[0], [2.0, 2.0, 3.0, 3.0])
    np.testing.assert_allclose(matrix[1], matrix[0, ::-1])


def test_extract_tss_matrix_from_signal_file(tmp_path):
    path = tmp_path / "signal.bedgraph"
    path.write_text(
        "# a comment\n"
        "track type=bedGraph\n"
        "chr1\t0\t100\t1\n"
        "chr1\t100\t200\t2\n"
        "chr1\t200\t300\t3\n"
        "chr1\t300\t400\t4\n",
        encoding="utf-8",
    )
    genes = pd.DataFrame(
        {"chrom": ["chr1"], "tss": [200], "strand": ["+"]}
    )
    matrix = extract_tss_matrix_from_signal(
        path, genes, window=100, n_bins=4, min_value=None
    )
    np.testing.assert_allclose(matrix[0], [2.0, 2.0, 3.0, 3.0])


def test_annotation_tss_plots_accept_matrix(annotation):
    matrix = np.ones((5, 20), dtype=float)
    fig, ax = annotation.plot_tss_profile(matrix=matrix, window=500)
    assert fig is not None
    assert "Distance to TSS" in ax.get_xlabel()
    plt.close(fig)

    fig, ax = annotation.plot_tss_heatmap(matrix=matrix, window=500)
    assert fig is not None
    plt.close(fig)


def test_annotation_tss_plots_accept_signal(annotation):
    signal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [10000, 20000],
            "value": [1.0, 2.0],
        }
    )
    fig, ax = annotation.plot_tss_profile(signal=signal, window=10, n_bins=5)
    assert fig is not None
    assert "Distance to TSS" in ax.get_xlabel()
    plt.close(fig)

    fig, ax = annotation.plot_tss_heatmap(signal=signal, window=10, n_bins=5)
    assert fig is not None
    plt.close(fig)


def test_extract_tss_matrix_optional_dependency(tmp_path):
    genes = pd.DataFrame({"chrom": ["chr1"], "tss": [1000], "strand": ["+"]})
    try:
        import pyBigWig  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="pyBigWig"):
            extract_tss_matrix(tmp_path / "missing.bw", genes)
    else:
        pytest.skip("pyBigWig is installed; optional import path not exercised")
