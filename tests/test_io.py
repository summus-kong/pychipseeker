import gzip
from pathlib import Path

import pandas as pd
import pytest

from pychipseeker.io import read_gtf, read_peaks


def test_read_gtf_coordinates_and_features(examples_dir):
    df = read_gtf(examples_dir / "example.gtf")
    gene = df[df["gene_id"] == "G1"].iloc[0]
    assert gene["feature"] == "gene"
    assert gene["start"] == 1000  # GTF 1001 -> 0-based
    assert gene["end"] == 3000
    assert gene["strand"] == "+"
    assert gene["gene_name"] == "GeneA"

    assert "five_prime_utr" in set(df["feature"])
    assert "three_prime_utr" in set(df["feature"])


def test_read_gtf_gzip(tmp_path, examples_dir):
    src = examples_dir / "example.gtf"
    dst = tmp_path / "example.gtf.gz"
    with src.open("rb") as f_in, gzip.open(dst, "wb") as f_out:
        f_out.write(f_in.read())

    df = read_gtf(dst)
    assert len(df) > 0
    assert "G1" in set(df["gene_id"].dropna())


def test_read_peaks_bed(examples_dir):
    df = read_peaks(examples_dir / "example_peaks.bed")
    assert list(df.columns[:4]) == ["chrom", "start", "end", "name"]
    assert df.loc[0, "chrom"] == "chr1"
    assert df.loc[0, "start"] == 950
    assert df.loc[0, "end"] == 1050


def test_read_peaks_bed_with_extra_columns(tmp_path):
    path = tmp_path / "peaks.bed"
    path.write_text("chr1\t10\t20\tp1\t0\t+\tfoo\tbar\n", encoding="utf-8")
    df = read_peaks(path)
    assert list(df.columns[:6]) == ["chrom", "start", "end", "name", "score", "strand"]
    assert "extra_1" in df.columns
    assert df.loc[0, "extra_1"] == "foo"


def test_read_peaks_tsv(tmp_path):
    path = tmp_path / "peaks.tsv"
    pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [200], "name": ["p1"]}
    ).to_csv(path, sep="\t", index=False)
    df = read_peaks(path)
    assert df.loc[0, "start"] == 100
    assert df.loc[0, "name"] == "p1"


def test_read_peaks_headerless_tsv(tmp_path):
    path = tmp_path / "peaks.tsv"
    path.write_text("chr1\t10\t20\nchr1\t30\t40\n", encoding="utf-8")
    df = read_peaks(path)
    assert len(df) == 2
    assert list(df.columns[:3]) == ["chrom", "start", "end"]


def test_read_peaks_invalid(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("chrom\tfoo\tbar\nchr1\t1\t2\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_peaks(path)
