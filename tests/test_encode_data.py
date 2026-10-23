"""Tests for the bundled ENCODE H3K4me3 demo signal bigWig."""

import json
import sys

import pytest

from pychipseeker.tagmatrix import example_signal_path, example_signal_provenance


def _bundle_path():
    try:
        return example_signal_path()
    except FileNotFoundError:  # pragma: no cover - package-data fallback
        from pathlib import Path

        path = Path(__file__).resolve().parent / "src" / "pychipseeker" / "data" / "ENCODE_demo_signal.bw"
        return path if path.exists() else None


def _pybigwig():
    pytest.importorskip("pyBigWig")
    import pyBigWig

    return pyBigWig


def test_example_signal_path_exists():
    path = example_signal_path()
    assert path.exists()
    assert path.stat().st_size < 1_500_000, "bundled demo signal should stay tiny"


def test_bundle_chroms_and_header():
    pyBigWig = _pybigwig()
    path = _bundle_path()
    if path is None:
        pytest.skip("bundled ENCODE demo signal not present")
    bw = pyBigWig.open(str(path))
    try:
        assert dict(bw.chroms()) == {"chr5": 151834684}
        header = bw.header()
        assert header["nLevels"] == 0
        assert header["nBasesCovered"] > 0
        assert header["minVal"] is not None
    finally:
        bw.close()


def test_bundle_signal_values():
    pyBigWig = _pybigwig()
    path = _bundle_path()
    if path is None:
        pytest.skip("bundled ENCODE demo signal not present")
    bw = pyBigWig.open(str(path))
    try:
        intervals = bw.intervals("chr5", 0, 151834684)
        if intervals is None:
            pytest.fail("no intervals found in bundled demo signal")
        assert len(intervals) > 10_000
        starts = [interval[0] for interval in intervals]
        assert starts == sorted(starts)
        # interval tuples are (start, end, value)
        assert all(interval[1] - interval[0] == 25 for interval in intervals)
        assert any(interval[2] >= 1.0 for interval in intervals)
        mean = bw.stats(
            "chr5",
            intervals[0][0],
            intervals[-1][1],
            nBins=1,
            type="mean",
        )[0]
        assert mean and mean > 0
    finally:
        bw.close()


def test_provenance_record():
    prov = example_signal_provenance()
    assert prov["source_file"] == "ENCFF893ZGR"
    assert prov["experiment"] == "ENCSR000CDQ"
    assert prov["assembly"] == "mm10"
    assert prov["target"].startswith("H3K4me3")
    subset = prov["subset"]
    assert subset["chromosome"] == "chr5"
    assert subset["bin_bp"] == 25
    assert subset["retained_bins"] > 0
    assert prov["bundle_sha256"]


def test_provenance_sha256_matches_bundle():
    import hashlib
    from pathlib import Path

    prov = example_signal_provenance()
    path = _bundle_path()
    if path is None:
        pytest.skip("bundled ENCODE demo signal not present")
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    assert digest == prov["bundle_sha256"]
