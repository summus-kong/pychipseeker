"""Process-pool helpers and n_jobs parity checks."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import pychipseeker as pcs
from pychipseeker.annotation import GeneAnnotation
from pychipseeker.parallel import chunk_bounds, effective_jobs
from pychipseeker.tagmatrix import get_tag_matrix_from_bam


def test_effective_jobs() -> None:
    assert effective_jobs(None, tasks=10) == 1
    assert effective_jobs(1, tasks=10) == 1
    assert effective_jobs(4, tasks=10) == 4
    assert effective_jobs(-1, tasks=2) == 2
    assert effective_jobs(8, tasks=3) == 3
    assert effective_jobs(-1, tasks=0) == 1
    with pytest.raises(ValueError, match="n_jobs"):
        effective_jobs(0)
    with pytest.raises(ValueError, match="integer"):
        effective_jobs("two")  # type: ignore[arg-type]


def test_chunk_bounds_cover_every_row() -> None:
    bounds = chunk_bounds(10, 3)
    assert bounds[0][0] == 0
    assert bounds[-1][1] == 10
    assert [stop - start for start, stop in bounds] == [4, 4, 2]
    assert chunk_bounds(0, 4) == []
    assert chunk_bounds(3, 10) == [(0, 1), (1, 2), (2, 3)]


def test_annotate_peak_parallel_matches_sequential(examples_dir: Path) -> None:
    annotation = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    rng = np.random.default_rng(7)
    n_peaks = 2500
    peaks = pd.DataFrame(
        {
            "chrom": np.repeat("chr1", n_peaks),
            "start": rng.integers(900, 9000, n_peaks),
        }
    )
    peaks["end"] = peaks["start"] + rng.integers(100, 600, n_peaks)

    sequential = pcs.annotate_peak(peaks, TxDb=annotation, n_jobs=1)
    parallel = pcs.annotate_peak(peaks, TxDb=annotation, n_jobs=2)

    pd.testing.assert_frame_equal(
        sequential.to_dataframe(), parallel.to_dataframe()
    )
    pd.testing.assert_frame_equal(sequential.anno_stat(), parallel.anno_stat())


def test_call_peaks_parallel_matches_sequential() -> None:
    pytest.importorskip("pyBigWig")
    path = Path(__file__).parent / "data" / "test.bw"
    if not path.exists():  # pragma: no cover - packaging guard
        pytest.skip("tests/data/test.bw is missing")

    options = dict(threshold=0.05, min_width=1, max_gap=0)
    sequential = pcs.call_peaks(path, n_jobs=1, **options)
    parallel = pcs.call_peaks(path, n_jobs=2, **options)
    pd.testing.assert_frame_equal(sequential, parallel)


def test_bam_tag_matrix_parallel_matches_sequential(small_bam: Path) -> None:
    n_windows = 600
    windows = pd.DataFrame(
        {
            "chrom": np.repeat("chr1", n_windows),
            "start": np.repeat(0, n_windows),
            "end": np.repeat(200, n_windows),
            "strand": np.repeat("+", n_windows),
        }
    )
    sequential = get_tag_matrix_from_bam(small_bam, windows, nbin=4, n_jobs=1)
    parallel = get_tag_matrix_from_bam(small_bam, windows, nbin=4, n_jobs=2)
    np.testing.assert_allclose(sequential, parallel)
    assert sequential.shape == (n_windows, 4)


def test_bam_tss_enrichment_parallel_matches_sequential(small_bam: Path) -> None:
    from pychipseeker.qc import tss_enrichment

    genes = pd.DataFrame(
        {
            "chrom": np.repeat("chr1", 600),
            "tss": np.repeat(110, 600),
            "strand": np.repeat("+", 600),
        }
    )
    options = dict(
        bam=small_bam,
        upstream=100,
        downstream=100,
        n_bins=4,
        sample=None,
    )
    sequential = tss_enrichment(genes, n_jobs=1, **options)
    parallel = tss_enrichment(genes, n_jobs=2, **options)
    for key, value in sequential.items():
        assert parallel[key] == pytest.approx(value)
