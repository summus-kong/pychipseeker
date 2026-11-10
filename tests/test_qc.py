"""Peak/QC metrics, FRiP and BAM-backed matrices."""

from __future__ import annotations

import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pychipseeker.qc import frip, peak_qc, tss_enrichment
from pychipseeker.tagmatrix import example_signal_path, get_tag_matrix_from_bam


def _real_bigwig() -> Path:
    pytest.importorskip("pyBigWig")
    path = Path(__file__).parent / "data" / "test.bw"
    if not path.exists():  # pragma: no cover - packaging guard
        pytest.skip("tests/data/test.bw is missing")
    return path


def test_peak_qc_blacklist_and_mito() -> None:
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chrM"],
            "start": [100, 100, 300, 0],
            "end": [200, 200, 500, 100],
        }
    )
    blacklist = pd.DataFrame({"chrom": ["1"], "start": [150], "end": [250]})
    table = peak_qc(
        peaks,
        blacklist=blacklist,
        chrom_sizes={"chr1": 10_000, "chrM": 16_000},
    )
    metrics = table.set_index("Metric")["Value"]
    assert metrics["Peaks"] == 4
    assert metrics["Chromosomes"] == 2
    assert metrics["Mean width (bp)"] == pytest.approx(125.0)
    assert metrics["Mitochondrial peaks"] == 1
    assert metrics["Mitochondrial fraction"] == pytest.approx(0.25)
    assert metrics["Duplicate fraction"] == pytest.approx(0.25)
    assert metrics["Blacklisted peaks"] == 2
    assert metrics["Blacklisted fraction"] == pytest.approx(0.5)
    assert metrics["Peaks per Mb"] == pytest.approx(4 / (26_000 / 1e6))


def test_frip_bigwig_harmonises_chromosomes() -> None:
    pyBigWig = pytest.importorskip("pyBigWig")
    bigwig = _real_bigwig()
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [4]})
    value = frip(peaks, bigwig=bigwig)
    bw = pyBigWig.open(str(bigwig))
    try:
        total = sum(
            (bw.stats(chrom, 0, size, nBins=1, type="sum") or [0])[0] or 0.0
            for chrom, size in bw.chroms().items()
        )
        inside = bw.stats("1", 0, 4, nBins=1, type="sum")[0]
    finally:
        bw.close()
    assert value == pytest.approx(inside / total)


def test_tss_enrichment_uses_center_and_flanks() -> None:
    pytest.importorskip("pyBigWig")
    genes = pd.DataFrame(
        {"chrom": ["chr5"], "tss": [121_712_100], "strand": ["+"]}
    )
    result = tss_enrichment(
        genes,
        bigwig=example_signal_path(),
        upstream=2000,
        downstream=2000,
        n_bins=100,
    )
    assert result["TSS features"] == 1
    assert np.isfinite(result["TSS enrichment"])
    assert result["TSS enrichment"] > 1


def test_frip_bam_excludes_duplicates(monkeypatch: pytest.MonkeyPatch) -> None:
    reads = [
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr1",
            reference_start=100,
            reference_end=150,
            query_length=50,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr1",
            reference_start=120,
            reference_end=170,
            query_length=50,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr1",
            reference_start=400,
            reference_end=450,
            query_length=50,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=True,
            mapping_quality=60,
            reference_name="chr1",
            reference_start=130,
            reference_end=180,
            query_length=50,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr2",
            reference_start=100,
            reference_end=150,
            query_length=50,
        ),
    ]

    class Handle:
        references = ["chr1"]

        def fetch(self, until_eof=True):
            return iter(reads)

        def close(self):
            return None

    fake = types.ModuleType("pysam")
    fake.AlignmentFile = lambda *args, **kwargs: Handle()
    monkeypatch.setitem(sys.modules, "pysam", fake)

    peaks = pd.DataFrame({"chrom": ["1"], "start": [100], "end": [200]})
    value = frip(peaks, bam="fake.bam")
    # chr2 still contributes to the denominator even though it has no peaks.
    assert value == pytest.approx(2 / 4)
    with pytest.raises(ValueError, match="exactly one"):
        frip(peaks)


def test_bam_tag_matrix_binning_and_strand(monkeypatch: pytest.MonkeyPatch) -> None:
    class Handle:
        references = ["chr1"]

        def count_coverage(
            self, chrom, start, end, quality_threshold=0, read_callback="all"
        ):
            values = np.zeros(end - start, dtype=float)
            half = (end - start) // 2
            values[:half] = 1.0
            values[half:] = 3.0
            return (
                values,
                np.zeros_like(values),
                np.zeros_like(values),
                np.zeros_like(values),
            )

        def close(self):
            return None

    fake = types.ModuleType("pysam")
    fake.AlignmentFile = lambda *args, **kwargs: Handle()
    monkeypatch.setitem(sys.modules, "pysam", fake)

    windows = pd.DataFrame(
        {
            "chrom": ["1", "1"],
            "start": [0, 0],
            "end": [8, 8],
            "strand": ["+", "-"],
        }
    )
    matrix = get_tag_matrix_from_bam("fake.bam", windows, nbin=4, value_type="mean")
    np.testing.assert_allclose(matrix[0], [1.0, 1.0, 3.0, 3.0])
    np.testing.assert_allclose(matrix[1], [3.0, 3.0, 1.0, 1.0])

    coverage = get_tag_matrix_from_bam(
        "fake.bam", windows.iloc[[0]], nbin=2, value_type="coverage"
    )
    np.testing.assert_allclose(coverage[0], [1.0, 1.0])


def test_frip_bam_real_file(small_bam: Path) -> None:
    """The indexed counting path must agree with the read-level definition."""

    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [90], "end": [200]})
    assert frip(peaks, bam=small_bam) == pytest.approx(2 / 6)
    # MAPQ 10 drops the low-quality read from the denominator.
    assert frip(peaks, bam=small_bam, min_mapping_quality=10) == pytest.approx(2 / 5)
    # Keeping duplicates adds one read both inside and outside the peak.
    assert frip(
        peaks, bam=small_bam, exclude_duplicates=False
    ) == pytest.approx(3 / 7)


def test_frip_bam_indexed_counting_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the C-level count() branch without a real BAM file."""

    reads = [
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr1",
            reference_start=100,
            reference_end=150,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr1",
            reference_start=400,
            reference_end=450,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name="chr2",
            reference_start=10,
            reference_end=60,
        ),
    ]

    class Handle:
        references = ["chr1", "chr2"]

        def has_index(self) -> bool:
            return True

        def count(
            self,
            contig=None,
            start=None,
            stop=None,
            until_eof: bool = False,
            read_callback="all",
        ) -> int:
            if until_eof:
                selected = reads
            else:
                selected = [
                    read
                    for read in reads
                    if read.reference_name == contig
                    and read.reference_start < int(stop)
                    and read.reference_end > int(start)
                ]
            return sum(1 for read in selected if read_callback(read))

        def close(self) -> None:
            return None

    fake = types.ModuleType("pysam")
    fake.AlignmentFile = lambda *args, **kwargs: Handle()
    monkeypatch.setitem(sys.modules, "pysam", fake)

    peaks = pd.DataFrame({"chrom": ["1"], "start": [90], "end": [200]})
    assert frip(peaks, bam="fake.bam") == pytest.approx(1 / 3)


def test_tss_enrichment_real_bam(small_bam: Path) -> None:
    genes = pd.DataFrame({"chrom": ["chr1"], "tss": [110], "strand": ["+"]})
    result = tss_enrichment(
        genes,
        bam=small_bam,
        upstream=100,
        downstream=100,
        n_bins=4,
        sample=None,
        center_window=50,
    )
    assert result["TSS features"] == 1
    assert np.isfinite(result["TSS enrichment"])
    assert result["TSS enrichment"] > 0


def test_bam_tag_matrix_real_file(small_bam: Path) -> None:
    windows = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 0],
            "end": [200, 200],
            "strand": ["+", "-"],
        }
    )
    matrix = get_tag_matrix_from_bam(
        small_bam, windows, nbin=4, value_type="mean"
    )
    assert matrix.shape == (2, 4)
    # The reverse-strand row is the forward profile flipped.
    np.testing.assert_allclose(matrix[1], matrix[0][::-1])
