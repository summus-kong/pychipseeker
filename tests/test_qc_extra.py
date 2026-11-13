"""Edge paths for peak QC, FRiP and TSS enrichment."""

from __future__ import annotations

import sys
import types
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pychipseeker import qc
from pychipseeker.qc import (
    _bam_tss_matrix,
    _bam_tss_rows,
    _bam_tss_worker,
    _bam_tss_worker_init,
    _bigwig_handle,
    _chrom_sizes_from_genome,
    _coerce_peaks,
    _frip_bam,
    _frip_bigwig,
    _merge_intervals,
    _profile_score,
    _read_keep_callback,
    frip,
    peak_qc,
    qc_summary,
    tss_enrichment,
)


def test_coerce_peaks_type_columns_and_empty_path(examples_dir) -> None:
    with pytest.raises(TypeError, match="peaks must be"):
        _coerce_peaks([1, 2, 3])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="missing columns"):
        _coerce_peaks(pd.DataFrame({"chrom": ["chr1"]}))

    loaded = _coerce_peaks(examples_dir / "example_peaks.bed")
    assert loaded["start"].dtype.kind == "i"

    zero = pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [10]})
    with pytest.raises(ValueError, match="peak table is empty"):
        peak_qc(zero)


def test_chrom_sizes_from_genome_uses_max_end(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = SimpleNamespace(
        genes=pd.DataFrame({"chrom": ["chr1", "chr1"], "end": [100, 250]}),
        transcripts=pd.DataFrame({"chrom": ["chr2"], "end": [50]}),
        exons=None,
    )
    monkeypatch.setattr("pychipseeker.genomes.load_genome", lambda name: fake)
    assert _chrom_sizes_from_genome("demo") == {"chr1": 250, "chr2": 50}


def test_peak_qc_genome_blacklist_path_and_bad_type(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    monkeypatch.setattr(
        qc,
        "_chrom_sizes_from_genome",
        lambda genome: {"chr1": 1000, "chr2": 1000},
    )
    peaks = pd.DataFrame(
        {"chrom": ["1", "chr2"], "start": [0, 0], "end": [100, 100]}
    )
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t50\t60\n", encoding="utf-8")
    table = peak_qc(peaks, genome="demo", blacklist=blacklist)
    metrics = table.set_index("Metric")["Value"]
    assert metrics["Peaks"] == 2
    assert metrics["Peaks per Mb"] == pytest.approx(2 / (2000 / 1e6))
    assert metrics["Blacklisted peaks"] == 1

    with pytest.raises(TypeError, match="blacklist must be"):
        peak_qc(peaks, blacklist=["chr1"])  # type: ignore[arg-type]


def test_merge_intervals_clips_and_drops_out_of_range() -> None:
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "start": [0, 50, 900, 10],
            "end": [60, 70, 1100, 20],
        }
    )
    merged = _merge_intervals(peaks, chrom_sizes={"chr1": 1000})
    rows = [(row.chrom, row.start, row.end) for row in merged.itertuples(index=False)]
    assert rows == [
        ("chr1", 0, 70),
        ("chr1", 900, 1000),
        ("chr2", 10, 20),
    ]

    outside = pd.DataFrame({"chrom": ["chr1"], "start": [5000], "end": [6000]})
    assert _merge_intervals(outside, chrom_sizes={"chr1": 1000}).empty


def test_frip_bigwig_overrides_sizes_and_rejects_empty_signal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Handle:
        def chroms(self):
            return {"1": 100, "2": 100}

        def stats(self, chrom, start, end, nBins, type):
            if start == 0 and end == 100:
                return [10.0]
            return [4.0]

        def close(self):
            return None

    monkeypatch.setattr(qc, "_bigwig_handle", lambda path: Handle())
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [20]})
    assert _frip_bigwig("fake.bw", peaks) == pytest.approx(4 / 20)

    sizes = {"1": 100, "2": 100}
    assert _frip_bigwig("fake.bw", peaks, chrom_sizes=sizes) == pytest.approx(4 / 20)

    class Empty(Handle):
        def stats(self, chrom, start, end, nBins, type):
            return [None]

    monkeypatch.setattr(qc, "_bigwig_handle", lambda path: Empty())
    with pytest.raises(ValueError, match="no signal"):
        _frip_bigwig("fake.bw", peaks)


def test_bigwig_handle_none(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = types.ModuleType("pyBigWig")
    fake.open = lambda path: None
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)
    with pytest.raises(ValueError, match="could not open bigWig"):
        _bigwig_handle("missing.bw")


def test_read_keep_callback_filters() -> None:
    keep = _read_keep_callback(True, 10)

    def read(**kwargs):
        defaults = dict(
            is_unmapped=False,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
        )
        defaults.update(kwargs)
        return SimpleNamespace(**defaults)

    assert keep(read())
    assert not keep(read(is_unmapped=True))
    assert not keep(read(is_secondary=True))
    assert not keep(read(is_supplementary=True))
    assert not keep(read(is_duplicate=True))
    assert not keep(read(mapping_quality=5))
    assert _read_keep_callback(False, 0)(read(is_duplicate=True))


def test_frip_bam_streaming_edge_cases(monkeypatch: pytest.MonkeyPatch) -> None:
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
            reference_start=400,
            reference_end=450,
            query_length=50,
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
            query_length=50,
        ),
        SimpleNamespace(
            is_unmapped=True,
            is_secondary=False,
            is_supplementary=False,
            is_duplicate=False,
            mapping_quality=60,
            reference_name=None,
            reference_start=0,
            reference_end=None,
            query_length=None,
        ),
        SimpleNamespace(
            is_unmapped=False,
            is_secondary=True,
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
            is_supplementary=True,
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
            is_duplicate=True,
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
            mapping_quality=5,
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
            reference_name=None,
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
            reference_start=100,
            reference_end=None,
            query_length=None,
        ),
    ]

    class Handle:
        references = ["chr1", "chr2"]

        def has_index(self):
            return False

        def fetch(self, until_eof=True):
            return iter(reads)

        def close(self):
            return None

    monkeypatch.setattr(qc, "_bam_handle", lambda path: Handle())
    peaks = pd.DataFrame({"chrom": ["1", "1"], "start": [90, 120], "end": [200, 210]})
    # Two merged peaks collapse to one interval; reads 1 and 10 overlap it.
    assert _frip_bam("fake.bam", peaks, min_mapping_quality=10) == pytest.approx(
        2 / 4
    )

    class NoReads(Handle):
        def fetch(self, until_eof=True):
            return iter(reads[3:4])

        def has_index(self):
            return False

    monkeypatch.setattr(qc, "_bam_handle", lambda path: NoReads())
    with pytest.raises(ValueError, match="no mapped reads"):
        _frip_bam("fake.bam", peaks)


def test_bam_tss_rows_binning_strand_and_skips(monkeypatch: pytest.MonkeyPatch) -> None:
    closed = []

    class Handle:
        def count_coverage(self, chrom, start, end, quality_threshold=0):
            if chrom == "bad":
                raise ValueError("bad")
            if chrom == "bad2":
                raise KeyError("bad2")
            if chrom == "empty":
                return ()
            if chrom == "nan":
                return (
                    np.full(4, np.nan),
                    np.zeros(4),
                    np.zeros(4),
                    np.zeros(4),
                )
            if chrom == "zero":
                return (np.array([]),) * 4
            values = np.asarray([0.0, 1.0, 2.0, 3.0])
            return (values, np.zeros(4), np.zeros(4), np.zeros(4))

        def close(self):
            closed.append(True)

    genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "bad", "bad2", "empty", "zero"],
            "tss": [100, 100, 100, 100, 100, 100],
            "strand": ["+", "-", "+", "+", "+", "+"],
        }
    )
    matrix = _bam_tss_rows(
        "fake.bam", genes, upstream=2, downstream=2, n_bins=2, handle=Handle()
    )
    np.testing.assert_allclose(matrix[0], [0.5, 2.5])
    np.testing.assert_allclose(matrix[1], [2.5, 0.5])
    assert np.isnan(matrix[2:]).all()
    assert closed == []

    direct = _bam_tss_rows(
        "fake.bam",
        genes.iloc[[0]],
        upstream=2,
        downstream=2,
        n_bins=4,
        handle=Handle(),
    )
    np.testing.assert_allclose(direct[0], [0.0, 1.0, 2.0, 3.0])

    skipped = _bam_tss_rows(
        "fake.bam",
        genes.iloc[[0]],
        upstream=0,
        downstream=0,
        n_bins=2,
        handle=Handle(),
    )
    assert np.isnan(skipped).all()

    def fake_handle(path):
        return Handle()

    monkeypatch.setattr(qc, "_bam_handle", fake_handle)
    owned = _bam_tss_rows(
        "fake.bam", genes.iloc[[0]], upstream=2, downstream=2, n_bins=2
    )
    assert owned.shape == (1, 2)
    assert closed == [True]


def test_bam_tss_worker_guard_and_arguments(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qc, "_BAM_TSS_WORKER", None)
    with pytest.raises(RuntimeError, match="not initialised"):
        _bam_tss_worker(pd.DataFrame())

    calls = []

    def fake_rows(bam, genes, *, upstream, downstream, n_bins):
        calls.append((bam, len(genes), upstream, downstream, n_bins))
        return np.ones((len(genes), n_bins))

    monkeypatch.setattr(qc, "_bam_tss_rows", fake_rows)
    _bam_tss_worker_init("fake.bam", 10, 20, 3)
    out = _bam_tss_worker(pd.DataFrame({"chrom": ["chr1", "chr1"]}))
    assert out.shape == (2, 3)
    assert calls == [("fake.bam", 2, 10, 20, 3)]


def test_bam_tss_matrix_sequential_and_parallel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = []

    def fake_rows(bam, genes, *, upstream, downstream, n_bins):
        calls.append(len(genes))
        return np.ones((len(genes), n_bins))

    monkeypatch.setattr(qc, "_bam_tss_rows", fake_rows)
    monkeypatch.setattr(qc, "effective_jobs", lambda n_jobs, tasks: 1)
    small = pd.DataFrame({"chrom": ["chr1"] * 4})
    out = _bam_tss_matrix(
        "fake.bam", small, upstream=1, downstream=1, n_bins=2, n_jobs=1
    )
    assert out.shape == (4, 2)
    assert calls == [4]

    calls.clear()
    monkeypatch.setattr(qc, "effective_jobs", lambda n_jobs, tasks: 2)

    def fake_map(func, items, n_jobs, initializer, initargs):
        initializer(*initargs)
        return [func(item) for item in items]

    monkeypatch.setattr(qc, "parallel_map", fake_map)
    large = pd.DataFrame({"chrom": ["chr1"] * 600})
    out = _bam_tss_matrix(
        "fake.bam", large, upstream=1, downstream=1, n_bins=2, n_jobs=2
    )
    assert out.shape == (600, 2)
    assert len(calls) > 1
    assert sum(calls) == 600


def test_profile_score_edges_and_windows() -> None:
    with pytest.raises(ValueError, match="features x bins"):
        _profile_score(
            np.ones(5),
            upstream=10,
            downstream=10,
            center_window=2,
            flank_window=None,
        )
    with pytest.raises(ValueError, match="features x bins"):
        _profile_score(
            np.ones((2, 1)),
            upstream=10,
            downstream=10,
            center_window=2,
            flank_window=None,
        )

    all_nan = np.full((2, 10), np.nan)
    result = _profile_score(
        all_nan, upstream=100, downstream=100, center_window=25, flank_window=None
    )
    assert result["TSS features"] == 0
    assert np.isnan(result["TSS enrichment"])

    with pytest.raises(ValueError, match="too small"):
        _profile_score(
            np.ones((2, 10)),
            upstream=10,
            downstream=10,
            center_window=1,
            flank_window=100,
        )

    flat = _profile_score(
        np.ones((2, 10)),
        upstream=100,
        downstream=100,
        center_window=25,
        flank_window=50,
    )
    assert flat["TSS enrichment"] == pytest.approx(1.0)
    assert flat["TSS features"] == 2

    centered = np.zeros((2, 10))
    centered[:, 4:6] = 3.0
    scored = _profile_score(
        centered, upstream=100, downstream=100, center_window=25, flank_window=None
    )
    assert scored["TSS center signal"] > scored["TSS flank signal"]
    assert scored["TSS enrichment"] > 1


def test_tss_enrichment_errors_and_sampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    genes = pd.DataFrame(
        {"chrom": ["chr1"] * 5, "tss": [100] * 5, "strand": ["+"] * 5}
    )
    with pytest.raises(ValueError, match="exactly one"):
        tss_enrichment(genes)
    with pytest.raises(ValueError, match="exactly one"):
        tss_enrichment(genes, bigwig="a.bw", bam="a.bam")
    with pytest.raises(ValueError, match="missing columns"):
        tss_enrichment(pd.DataFrame({"chrom": ["chr1"]}), bigwig="a.bw")
    with pytest.raises(ValueError, match="no valid TSS"):
        tss_enrichment(
            pd.DataFrame({"chrom": ["chr1"], "tss": [None], "strand": ["+"]}),
            bigwig="a.bw",
        )

    calls = []

    def fake_matrix(bigwig, frame, **kwargs):
        calls.append((bigwig, len(frame), kwargs))
        return np.ones((len(frame), 10))

    monkeypatch.setattr("pychipseeker.plotting.extract_tss_matrix", fake_matrix)
    result = tss_enrichment(
        genes,
        bigwig="fake.bw",
        upstream=100,
        downstream=100,
        n_bins=10,
        sample=2,
        seed=1,
        center_window=25,
    )
    assert result["TSS features"] == 2
    assert calls[0][1] == 2
    assert calls[0][2]["n_bins"] == 10
    assert calls[0][2]["upstream"] == 100


def test_tss_enrichment_bam_harmonises_and_scores(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Handle:
        references = ["chr1"]

        def close(self):
            return None

    monkeypatch.setattr(qc, "_bam_handle", lambda path: Handle())
    captured = {}

    def fake_matrix(bam, frame, *, upstream, downstream, n_bins, n_jobs):
        captured["frame"] = frame.copy()
        captured["args"] = (bam, upstream, downstream, n_bins)
        return np.zeros((len(frame), 10))

    monkeypatch.setattr(qc, "_bam_tss_matrix", fake_matrix)
    genes = pd.DataFrame(
        {"chrom": ["1", "1"], "tss": [100, 200], "strand": ["+", "-"]}
    )
    result = tss_enrichment(
        genes, bam="fake.bam", upstream=10, downstream=10, n_bins=10, sample=None
    )
    assert captured["frame"]["chrom"].tolist() == ["chr1", "chr1"]
    assert captured["args"] == ("fake.bam", 10, 10, 10)
    assert result["TSS enrichment"] == pytest.approx(1.0)
    assert result["TSS features"] == 2


def test_frip_requires_exactly_one_and_qc_summary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    with pytest.raises(ValueError, match="exactly one"):
        frip(peaks)
    with pytest.raises(ValueError, match="exactly one"):
        frip(peaks, bigwig="a.bw", bam="a.bam")
    with pytest.raises(ValueError, match="peak table is empty"):
        frip(
            pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [0]}),
            bigwig="a.bw",
        )

    monkeypatch.setattr(
        qc,
        "peak_qc",
        lambda *a, **k: pd.DataFrame({"Metric": ["Peaks"], "Value": [1]}),
    )
    monkeypatch.setattr(qc, "frip", lambda *a, **k: 0.5)
    monkeypatch.setattr(
        qc,
        "tss_enrichment",
        lambda *a, **k: {"TSS enrichment": 2.0, "TSS features": 3},
    )
    genes = pd.DataFrame({"chrom": ["chr1"], "tss": [5], "strand": ["+"]})
    summary = qc_summary(peaks, bigwig="fake.bw", genes=genes)
    metrics = summary.set_index("Metric")["Value"]
    assert metrics["FRiP"] == 0.5
    assert metrics["TSS enrichment"] == 2.0

    plain = qc_summary(peaks)
    assert list(plain["Metric"]) == ["Peaks"]
