"""Tests for the report builder, peak caller, cache helpers and UpSet plot."""

from __future__ import annotations

import json
import sys
import types

import numpy as np
import pandas as pd
import pytest

import pychipseeker as pcs
from pychipseeker.io import read_peaks, write_table
from pychipseeker.peaks import call_peaks, write_peaks


def test_call_peaks_threshold_merge_and_write(tmp_path, monkeypatch):
    class Handle:
        def __init__(self, path):
            self.path = path

        def chroms(self):
            return {"chr1": 5000, "chr2": 1000}

        def intervals(self, chrom):
            if chrom == "chr1":
                return [
                    (100, 200, 5.0),
                    (250, 300, 4.0),  # 50 bp gap -> merged
                    (1000, 1200, 6.0),
                    (3000, 3200, 2.0),  # below threshold
                ]
            return [(0, 200, 7.0)]

        def close(self):
            return None

    fake = types.ModuleType("pyBigWig")
    fake.open = Handle
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)

    peaks = call_peaks(
        "fake.bw", threshold=3.0, max_gap=100, min_width=150
    )
    assert len(peaks) == 3
    first = peaks.iloc[0]
    assert (first["chrom"], first["start"], first["end"]) == ("chr1", 100, 300)
    assert first["peak"] in {0, 150}  # highest base can be at either interval
    assert first["signalValue"] == pytest.approx((5 * 100 + 4 * 50) / 200)

    out = tmp_path / "called.bed"
    write_peaks(peaks, out)
    back = read_peaks(out)
    assert len(back) == 3
    assert list(back.columns[:6]) == [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]


def test_report_builds_multipage_pdf(tmp_path, examples_dir):
    output = tmp_path / "report.pdf"
    path = pcs.build_report(
        peaks=examples_dir / "example_peaks.bed",
        gtf=examples_dir / "example.gtf",
        output=output,
        sample=10,
        include_tag_matrix=False,
    )
    assert path == output
    assert output.exists()
    assert output.stat().st_size > 10_000
    assert output.read_bytes().startswith(b"%PDF")


def test_cache_helpers_list_and_clear(tmp_path, monkeypatch):
    from pychipseeker import genomes as genomes_mod

    monkeypatch.setattr(genomes_mod, "cache_dir", lambda: tmp_path)
    (tmp_path / "genome_dm6.npz").write_bytes(b"npz")
    (tmp_path / "genome_dm6.json").write_text(
        json.dumps({"provider": "ensembl", "release": "110"}), encoding="utf-8"
    )
    (tmp_path / "genome_dm6.gtf.gz").write_bytes(b"gtf")
    gtf_dir = tmp_path / "gtf"
    gtf_dir.mkdir()
    (gtf_dir / "spec_abc.gtf").write_text("chr1\n", encoding="utf-8")

    rows = {row["assembly"]: row for row in pcs.cached_genomes()}
    assert rows["dm6"]["source"] == "cache"
    assert rows["dm6"]["provider"] == "ensembl"

    removed = pcs.clear_genome_cache("dm6")
    assert (tmp_path / "genome_dm6.npz") in removed
    assert not (tmp_path / "genome_dm6.npz").exists()
    assert not (tmp_path / "genome_dm6.gtf.gz").exists()

    pcs.clear_genome_cache()
    assert not (gtf_dir / "spec_abc.gtf").exists()


def test_upsetplot_uses_dot_matrix(examples_dir):
    import matplotlib.pyplot as plt

    ann = pcs.annotate_peak(
        peaks=examples_dir / "example_peaks.bed",
        gtf=examples_dir / "example.gtf",
        verbose=False,
    )
    fig, ax = pcs.upsetplot(ann, n_intersections=5)
    assert len(fig.axes) >= 3
    assert ax is not None
    plt.close(fig)


def test_extract_tss_matrix_intervals_and_strand(monkeypatch):
    from pychipseeker.plotting import extract_tss_matrix

    class FakeBigWig:
        def __init__(self, path):
            self.path = path

        def chroms(self):
            return {"chr1": 1000}

        def intervals(self, chrom):
            assert chrom == "chr1"
            return [
                (0, 100, 1.0),
                (100, 200, 2.0),
                (200, 300, 3.0),
                (300, 400, 4.0),
            ]

        def close(self):
            return None

    fake = types.ModuleType("pyBigWig")
    fake.open = FakeBigWig
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)

    genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "tss": [250, 250],
            "strand": ["+", "-"],
        }
    )
    matrix = extract_tss_matrix("fake.bw", genes, window=100, n_bins=4)
    np.testing.assert_allclose(matrix[0], [2.0, 3.0, 3.0, 4.0])
    np.testing.assert_allclose(matrix[1], matrix[0, ::-1])


def test_write_table_csv(tmp_path):
    df = pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]})
    path = tmp_path / "peaks.csv"
    write_table(df, path)
    assert read_peaks(path).loc[0, "end"] == 2
