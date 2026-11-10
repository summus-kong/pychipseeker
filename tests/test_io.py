import gzip

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


def test_read_gff3_fasta_and_multi_parent(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t101\t500\t.\t+\t.\tID=gene:g1;Name=Gene1\n"
        "chr1\t.\tmRNA\t101\t500\t.\t+\t.\tID=transcript:t1;Parent=gene:g1\n"
        "chr1\t.\tmRNA\t101\t500\t.\t+\t.\tID=transcript:t2;Parent=gene:g1\n"
        "chr1\t.\texon\t101\t200\t.\t+\t.\tParent=transcript:t1,transcript:t2\n"
        "##FASTA\n"
        ">chr1\n"
        "ACGTACGT\n",
        encoding="utf-8",
    )
    df = read_gtf(path)
    exons = df[df["feature"] == "exon"]
    assert len(exons) == 2
    assert set(exons["transcript_id"]) == {"t1", "t2"}
    assert set(exons["gene_id"]) == {"g1"}
    # FASTA sequence rows must not become features.
    assert "ACGTACGT" not in set(df["chrom"])


def test_read_gtf_feature_filter_and_chunks(tmp_path, examples_dir):
    full = read_gtf(examples_dir / "example.gtf")
    exons = read_gtf(examples_dir / "example.gtf", features="exon")
    assert set(exons["feature"]) == {"exon"}
    assert len(exons) == (full["feature"] == "exon").sum()

    chunked = read_gtf(examples_dir / "example.gtf", chunksize=2)
    assert len(chunked) == len(full)
    pd.testing.assert_frame_equal(
        chunked.reset_index(drop=True), full.reset_index(drop=True)
    )


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


def test_read_narrowpeak_columns(tmp_path):
    path = tmp_path / "peaks.narrowPeak"
    path.write_text(
        "chr1\t100\t200\tp1\t500\t+\t10.5\t8.1\t7.2\t50\n"
        "chr1\t300\t450\tp2\t400\t-\t9.5\t6.1\t5.2\t75\n",
        encoding="utf-8",
    )
    df = read_peaks(path)
    assert list(df.columns) == [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
    ]
    assert df.loc[0, "signalValue"] == 10.5
    assert df.loc[1, "peak"] == 75


def test_read_bed_header_and_track_lines(tmp_path):
    path = tmp_path / "peaks.bed"
    path.write_text(
        "track name=test\n"
        "browser position chr1\n"
        "#chrom\tstart\tend\tname\tscore\tstrand\n"
        "chr1\t5\t9\tn\t0\t+\n",
        encoding="utf-8",
    )
    df = read_peaks(path)
    assert list(df.columns[:6]) == [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]
    assert df.loc[0, "start"] == 5


def test_read_bigbed_with_pybigtools(monkeypatch, tmp_path):
    import sys
    import types

    from pychipseeker import io as io_mod

    path = tmp_path / "peaks.bb"
    path.write_bytes(b"fake bigBed")

    class Handle:
        def chroms(self):
            return {"chr1": 1000}

        def records(self, chrom, start, end):
            assert chrom == "chr1"
            return [(100, 200, "p1\t500\t+"), (300, 400, "p2\t100\t-")]

        def close(self):
            return None

    fake = types.ModuleType("pybigtools")
    fake.open = lambda _path: Handle()
    monkeypatch.setitem(sys.modules, "pybigtools", fake)
    monkeypatch.setattr(io_mod.shutil, "which", lambda _name: None)

    df = read_peaks(path)
    assert list(df.columns[:6]) == [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]
    assert df.loc[1, "start"] == 300
    assert df.loc[1, "name"] == "p2"


def test_read_bam_with_pysam(monkeypatch, tmp_path):
    import sys
    import types

    path = tmp_path / "reads.bam"
    path.write_bytes(b"fake bam")

    class Read:
        is_unmapped = False
        reference_start = 100
        reference_end = 150
        query_name = "read1"
        mapping_quality = 60
        is_reverse = False
        reference_name = "chr1"

    class Bam:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def fetch(self, until_eof=True):
            assert until_eof
            return [Read()]

    fake = types.ModuleType("pysam")
    fake.AlignmentFile = lambda _path, _mode: Bam()
    monkeypatch.setitem(sys.modules, "pysam", fake)

    df = read_peaks(path)
    assert df.loc[0, "chrom"] == "chr1"
    assert df.loc[0, "start"] == 100
    assert df.loc[0, "end"] == 150
    assert df.loc[0, "name"] == "read1"
    assert df.loc[0, "strand"] == "+"


def test_bigbed_requires_backend(monkeypatch, tmp_path):
    import sys

    from pychipseeker import io as io_mod

    path = tmp_path / "peaks.bb"
    path.write_bytes(b"fake bigBed")
    monkeypatch.setattr(io_mod.shutil, "which", lambda _name: None)
    monkeypatch.setitem(sys.modules, "pybigtools", None)
    with pytest.raises(ImportError, match="bigBedToBed"):
        read_peaks(path)
