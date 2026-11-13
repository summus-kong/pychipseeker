"""Threshold peak calling, filtering and the MACS wrapper."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pychipseeker.peaks import (
    _call_peaks_chrom,
    _callpeak_worker,
    _callpeak_worker_init,
    _overlaps_any,
    call_peaks,
    call_peaks_macs,
    filter_peaks,
    write_peaks,
)


class FakeBigWig:
    def __init__(self, intervals, chroms=None):
        self._intervals = intervals
        self._chroms = chroms or {"chr1": 1000}
        self.closed = False

    def chroms(self):
        return dict(self._chroms)

    def intervals(self, chrom):
        value = self._intervals.get(str(chrom)) if isinstance(self._intervals, dict) else self._intervals
        return value

    def close(self):
        self.closed = True


def test_call_peaks_chrom_merges_runs_and_filters_width() -> None:
    bw = FakeBigWig(
        [
            (0, 10, 5.0),
            (20, 30, 2.0),  # gap 10 -> same run
            (100, 120, 9.0),  # gap 70 -> new run
            (110, 130, 1.0),  # nested interval, cumulative end handling
        ]
    )
    rows = _call_peaks_chrom(bw, "chr1", threshold=1.0, max_gap=10, min_width=10)
    assert [(row["start"], row["end"]) for row in rows] == [(0, 30), (100, 130)]
    assert rows[0]["score"] == pytest.approx(5.0)
    assert rows[0]["signalValue"] == pytest.approx((10 * 5 + 10 * 2) / 30)
    assert rows[0]["peak"] == 0
    assert rows[1]["peak"] == 0

    narrow = _call_peaks_chrom(bw, "chr1", threshold=1.0, max_gap=10, min_width=40)
    assert narrow == []


def test_call_peaks_chrom_returns_empty_for_missing_or_weak_signal() -> None:
    assert _call_peaks_chrom(FakeBigWig({}), "chr1", threshold=1, max_gap=0, min_width=0) == []
    assert (
        _call_peaks_chrom(
            FakeBigWig(np.empty((0, 3))),
            "chr1",
            threshold=1,
            max_gap=0,
            min_width=0,
        )
        == []
    )
    assert (
        _call_peaks_chrom(
            FakeBigWig([(0, 10, 0.1), (20, 30, np.nan)]),
            "chr1",
            threshold=1.0,
            max_gap=0,
            min_width=0,
        )
        == []
    )


def test_call_peaks_real_bigwig_sample_and_errors() -> None:
    bigwig = Path(__file__).parent / "data" / "test.bw"
    pytest.importorskip("pyBigWig")
    peaks = call_peaks(
        bigwig,
        threshold=0.5,
        max_gap=10,
        min_width=5,
        chroms=["1"],
        sample=1,
        seed=3,
        n_jobs=1,
    )
    assert list(peaks.columns)[:3] == ["chrom", "start", "end"]
    assert len(peaks) == 1
    assert peaks.iloc[0]["name"].startswith("Peak_")
    assert peaks.iloc[0]["end"] - peaks.iloc[0]["start"] >= 5

    with pytest.raises(ValueError, match="no peaks"):
        call_peaks(bigwig, threshold=10.0, chroms=["1"], n_jobs=1)
    with pytest.raises(ValueError, match="non-negative"):
        call_peaks(bigwig, max_gap=-1)
    with pytest.raises(ValueError, match="non-negative"):
        call_peaks(bigwig, min_width=-1)


def test_callpeak_worker_uses_initialiser(monkeypatch: pytest.MonkeyPatch) -> None:
    bigwig = Path(__file__).parent / "data" / "test.bw"
    pytest.importorskip("pyBigWig")
    _callpeak_worker_init(bigwig, 0.5, 10, 5)
    rows = _callpeak_worker("1")
    assert len(rows) == 1

    monkeypatch.setattr("pychipseeker.peaks._CALLPEAK_WORKER", None)
    with pytest.raises(RuntimeError, match="not initialised"):
        _callpeak_worker("1")


def test_overlaps_any_handles_nan_and_missing_columns() -> None:
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "start": [0, 50, 100, 0],
            "end": [10, 60, 110, 10],
        }
    )
    blacklist = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [5, None, 90],
            "end": [15, 100, 105],
        }
    )
    assert _overlaps_any(peaks, blacklist).tolist() == [True, False, True, False]

    with pytest.raises(ValueError, match="chrom/start/end"):
        _overlaps_any(peaks, pd.DataFrame({"chrom": ["chr1"]}))


def test_filter_peaks_dataframe_and_path_filters() -> None:
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "1", "chrM", "MT", "chr2", "chr3"],
            "start": [0, 20, 0, 10, 0, 0],
            "end": [10, 40, 10, 20, 50, 10],
            "score": [1, 2, 3, 4, 5, 6],
        }
    )
    blacklist = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [5]})

    filtered = filter_peaks(peaks, blacklist=blacklist, drop_mito=True)
    assert "chr1" not in set(filtered["chrom"])
    assert "chrM" not in set(filtered["chrom"])
    assert "MT" not in set(filtered["chrom"])

    filtered = filter_peaks(
        peaks,
        chroms=["chr1"],
        min_width=10,
        max_width=30,
        min_score=2,
    )
    assert filtered["chrom"].tolist() == ["1"]
    assert (filtered["end"] - filtered["start"]).between(10, 30).all()

    filtered = filter_peaks(peaks, exclude_chroms=["chrM", "MT"])
    assert set(filtered["chrom"]) == {"chr1", "1", "chr2", "chr3"}

    path = Path(__file__).resolve().parents[1] / "examples" / "example_peaks.bed"
    loaded = filter_peaks(path, chroms=["chr1"], min_width=10)
    assert not loaded.empty


def test_filter_peaks_type_and_column_errors() -> None:
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    with pytest.raises(TypeError, match="peaks must be"):
        filter_peaks([1, 2, 3])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="missing columns"):
        filter_peaks(pd.DataFrame({"chrom": ["chr1"]}))
    with pytest.raises(TypeError, match="blacklist must be"):
        filter_peaks(peaks, blacklist=[("chr1", 0, 10)])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="score"):
        filter_peaks(peaks, min_score=1)


def test_write_peaks_orders_columns(tmp_path: Path) -> None:
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [10],
            "signalValue": [1.5],
            "name": ["p1"],
            "score": [2],
            "strand": ["+"],
        }
    )
    out = tmp_path / "peaks.bed"
    write_peaks(peaks, out)
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert header[:6] == ["chrom", "start", "end", "name", "score", "strand"]
    assert header[-1] == "signalValue"


def test_call_peaks_macs_missing_backend(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("pychipseeker.peaks.shutil.which", lambda _name: None)
    with pytest.raises(ImportError, match="not found on PATH"):
        call_peaks_macs("treatment.bam")


def test_call_peaks_macs_invokes_backend(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls = []

    def fake_run(command, check):
        calls.append((command, check))
        outdir = Path(command[command.index("--outdir") + 1])
        name = command[command.index("-n") + 1]
        bed = outdir / f"{name}_peaks.narrowPeak"
        bed.write_text(
            "chr1\t10\t20\tpeak1\t100\t.\t2.0\t3.0\t4.0\t5\n", encoding="utf-8"
        )

    monkeypatch.setattr("pychipseeker.peaks.shutil.which", lambda name: name)
    monkeypatch.setattr("pychipseeker.peaks.subprocess.run", fake_run)

    peaks = call_peaks_macs(
        "treatment.bam",
        control="control.bam",
        genome="mm",
        outdir=tmp_path,
        backend="macs3",
        extra_args=("--nomodel",),
        name="run",
    )
    assert peaks["chrom"].tolist() == ["chr1"]
    command, check = calls[0]
    assert command[0] == "macs3"
    assert "-c" in command and "control.bam" in command
    assert "-g" in command and "mm" in command
    assert "--nomodel" in command
    assert check is True


def test_call_peaks_macs_missing_output(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr("pychipseeker.peaks.shutil.which", lambda name: name)
    monkeypatch.setattr("pychipseeker.peaks.subprocess.run", lambda *a, **k: None)
    with pytest.raises(FileNotFoundError, match="no narrowPeak"):
        call_peaks_macs("treatment.bam", outdir=tmp_path)
