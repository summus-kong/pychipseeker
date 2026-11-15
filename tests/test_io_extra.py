"""Branch coverage for the GTF/peak readers in :mod:`pychipseeker.io`."""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

from pychipseeker import io as io_mod
from pychipseeker.io import (
    _bed_frame,
    _bed_header_info,
    _coerce_peaks,
    _GtfTextReader,
    _is_bed_header,
    _open_text,
    _peak_variant,
    read_gtf,
    read_peaks,
    write_table,
)


def test_gtf_text_reader_read_and_readline() -> None:
    reader = _GtfTextReader(io.StringIO("alpha\nbeta\n##FASTA\nseq\n"))
    assert reader.readable() is True
    assert reader.readline() == "alpha\n"
    assert reader.readline() == "beta\n"
    assert reader.readline() == ""

    reader = _GtfTextReader(io.StringIO("alpha\nbeta\n##FASTA\nseq\n"))
    assert reader.read(size=3) == "alp"
    assert reader.read() == "ha\nbeta\n"
    assert reader.read() == ""

    reader = _GtfTextReader(io.StringIO("no trailing newline"))
    assert reader.readline() == "no trailing newline"
    assert reader.readline() == ""


def test_read_gtf_errors_and_empty_chunks(tmp_path: Path, examples_dir: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_gtf(tmp_path / "missing.gtf")
    with pytest.raises(ValueError, match="chunksize"):
        read_gtf(examples_dir / "example.gtf", chunksize=0)

    invalid = tmp_path / "invalid.gtf"
    invalid.write_text(
        'chr1\t.\tgene\tbad\t100\t.\t+\t.\tgene_id "G1";\n',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no records found"):
        read_gtf(invalid)

    mixed = tmp_path / "mixed.gtf"
    mixed.write_text(
        'chr1\t.\tgene\tbad\t100\t.\t+\t.\tgene_id "G0";\n'
        'chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id "G1";\n',
        encoding="utf-8",
    )
    parsed = read_gtf(mixed, chunksize=1)
    assert set(parsed["gene_id"]) == {"G1"}

    with pytest.raises(ValueError, match="no records matching features"):
        read_gtf(examples_dir / "example.gtf", features="intron")


def test_open_text_compression_suffixes(tmp_path: Path) -> None:
    payload = b"chr1\t1\t2\n"
    for suffix, opener in (
        (".gz", gzip.open),
        (".bz2", bz2.open),
        (".xz", lzma.open),
    ):
        path = tmp_path / f"peaks{suffix}"
        with opener(path, "wb") as handle:
            handle.write(payload)
        with _open_text(path) as handle:
            assert handle.read() == "chr1\t1\t2\n"

    plain = tmp_path / "peaks.bed"
    plain.write_text("chr1\t1\t2\n", encoding="utf-8")
    with _open_text(plain) as handle:
        assert handle.read() == "chr1\t1\t2\n"


def test_bed_header_detection_branches(tmp_path: Path) -> None:
    assert _is_bed_header(["chr1", "start"]) is False

    blank = tmp_path / "blank.bed"
    blank.write_text("\nchr1\t1\t2\n", encoding="utf-8")
    assert _bed_header_info(blank) == (None, 1)

    comment = tmp_path / "comment.bed"
    comment.write_text("# a comment\nchr1\t1\t2\n", encoding="utf-8")
    assert _bed_header_info(comment) == (None, 1)

    only_comments = tmp_path / "comments.bed"
    only_comments.write_text("# one\n# two\n", encoding="utf-8")
    assert _bed_header_info(only_comments) == (None, 0)


def test_coerce_peaks_and_bed_frame_errors() -> None:
    positional = _coerce_peaks(pd.DataFrame({"a": ["chr1"], "b": [1], "c": [2]}))
    assert list(positional.columns[:3]) == ["chrom", "start", "end"]

    with pytest.raises(ValueError, match="missing columns"):
        _coerce_peaks(pd.DataFrame({"a": [1], "b": [2]}))
    with pytest.raises(ValueError, match="no peak records"):
        _bed_frame([])

    frame = _bed_frame([["chr1", 1, 2, "name"]])
    assert list(frame.columns[:4]) == ["chrom", "start", "end", "name"]

    identified = _coerce_peaks(
        pd.DataFrame(
            {"chrom": ["chr1"], "start": [1], "end": [2], "peakId": ["p1"]}
        )
    )
    assert identified.loc[0, "peak_id"] == "p1"


def test_peak_variant_suffixes() -> None:
    assert _peak_variant("peaks.narrowPeak.gz") == "narrowpeak"
    assert _peak_variant("peaks.broadPeak.bz2") == "broadpeak"
    assert _peak_variant("peaks.gappedPeak.xz") == "gappedpeak"
    assert _peak_variant("peaks.bed") == "bed"
    assert _peak_variant("peaks.txt") is None


def test_read_bigbed_with_converter(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    path = tmp_path / "peaks.bb"
    path.write_bytes(b"fake bigBed")

    def fake_run(command, check, stdout, stderr):
        assert check is True
        Path(command[2]).write_text(
            "chr1\t10\t20\tp1\t500\t+\n", encoding="utf-8"
        )

    monkeypatch.setattr(io_mod.shutil, "which", lambda name: "bigBedToBed")
    monkeypatch.setattr(io_mod.subprocess, "run", fake_run)
    df = read_peaks(path)
    assert df.loc[0, "name"] == "p1"
    assert read_peaks(path, fmt="bb").loc[0, "name"] == "p1"


def test_read_bigbed_record_shapes(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    path = tmp_path / "peaks.bb"
    path.write_bytes(b"fake bigBed")

    class Handle:
        def chroms(self):
            return {"chr1": 1000, "chr2": 1000}

        def records(self, chrom, start, end):
            if chrom == "chr2":
                return None
            return [
                (100, 200, "p1\t500\t+"),
                (300, 400, 700, "+"),
            ]

        def close(self):
            return None

    fake = types.ModuleType("pybigtools")
    fake.open = lambda _path: Handle()
    monkeypatch.setitem(sys.modules, "pybigtools", fake)
    monkeypatch.setattr(io_mod.shutil, "which", lambda _name: None)

    df = read_peaks(path)
    assert len(df) == 2
    assert df.loc[0, "name"] == "p1"
    assert df.loc[1, "name"] == 700


def test_read_bam_filters_and_errors(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    path = tmp_path / "reads.bam"
    path.write_bytes(b"fake bam")
    monkeypatch.setitem(sys.modules, "pysam", None)
    with pytest.raises(ImportError, match="pysam"):
        read_peaks(path)

    class Read:
        def __init__(
            self,
            *,
            unmapped: bool = False,
            end=None,
            query_length=50,
        ) -> None:
            self.is_unmapped = unmapped
            self.reference_start = 100
            self.reference_end = end
            self.query_length = query_length
            self.query_name = "read"
            self.mapping_quality = 60
            self.is_reverse = False
            self.reference_name = "chr1"

    class Bam:
        def __init__(self, reads) -> None:
            self._reads = reads

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def fetch(self, until_eof=True):
            return self._reads

    fake = types.ModuleType("pysam")
    fake.AlignmentFile = lambda _path, _mode: Bam(
        [Read(unmapped=True), Read(end=None, query_length=None)]
    )
    monkeypatch.setitem(sys.modules, "pysam", fake)

    df = read_peaks(path, fmt="bam")
    assert len(df) == 1
    assert df.loc[0, "end"] == 101

    fake.AlignmentFile = lambda _path, _mode: Bam([Read(unmapped=True)])
    with pytest.raises(ValueError, match="no aligned reads"):
        read_peaks(path, fmt="bam")


def test_read_peaks_dispatch_and_errors(tmp_path: Path) -> None:
    missing = tmp_path / "missing.bed"
    with pytest.raises(FileNotFoundError):
        read_peaks(missing)

    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t10\t20\tp1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="fmt must be one of"):
        read_peaks(bed, fmt="xlsx")

    narrow = tmp_path / "peaks.narrowPeak"
    narrow.write_text(
        "chr1\t10\t20\tp1\t500\t+\t1\t2\t3\t4\n", encoding="utf-8"
    )
    explicit = read_peaks(narrow, fmt="narrowpeak")
    assert explicit.loc[0, "qValue"] == 3
    assert read_peaks(narrow, fmt="bed").loc[0, "name"] == "p1"

    tsv = tmp_path / "peaks.tsv"
    pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]}).to_csv(
        tsv, sep="\t", index=False
    )
    assert read_peaks(tsv, fmt="tsv").loc[0, "start"] == 1

    csv = tmp_path / "peaks.csv"
    pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]}).to_csv(
        csv, index=False
    )
    assert read_peaks(csv, fmt="csv").loc[0, "start"] == 1
    assert read_peaks(csv).loc[0, "start"] == 1

    short = tmp_path / "short.bed"
    short.write_text("chr1\t10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="at least 3 columns"):
        read_peaks(short)

    too_short = tmp_path / "tiny.tsv"
    too_short.write_text("chr1\t10\nchr1\t20\n", encoding="utf-8")
    with pytest.raises(ValueError, match="chrom/start/end"):
        read_peaks(too_short)


def test_read_peaks_header_mapping_branches(tmp_path: Path) -> None:
    header = tmp_path / "mapped.bed"
    header.write_text(
        "#chrom\tstart\tend\tsignalValue\tpValue\tqValue\tpeak\n"
        "chr1\t10\t20\t9.5\t1.0\t2.0\t15\n",
        encoding="utf-8",
    )
    mapped = read_peaks(header)
    assert mapped.loc[0, "signalValue"] == 9.5
    assert mapped.loc[0, "pValue"] == 1.0
    assert mapped.loc[0, "qValue"] == 2.0
    assert mapped.loc[0, "name"] == 15

    empty_header = tmp_path / "empty-header.bed"
    empty_header.write_text(
        "chrom\tstart\tend\t\tcustom\n"
        "chr1\t10\t20\tignored\tvalue\n",
        encoding="utf-8",
    )
    empty_mapped = read_peaks(empty_header)
    assert empty_mapped.loc[0, "custom"] == "value"
    assert any(column.startswith("extra_") for column in empty_mapped.columns)

    extra = tmp_path / "extra.bed"
    extra.write_text(
        "chrom\tstart\tend\tname\tname\tcustom\t\n"
        "chr1\t10\t20\tp1\tp2\tvalue\tx\n",
        encoding="utf-8",
    )
    extras = read_peaks(extra)
    assert extras.loc[0, "name"] == "p1"
    assert extras.loc[0, "custom"] == "value"
    assert extras.columns[4].startswith("extra_")
    assert extras.columns[6].startswith("extra_")

    short_header = tmp_path / "short-header.bed"
    short_header.write_text(
        "chrom\tstart\tend\tname\n"
        "chr1\t10\t20\tp1\t500\t+\n",
        encoding="utf-8",
    )
    padded = read_peaks(short_header)
    assert padded.loc[0, "score"] == 500
    assert padded.loc[0, "strand"] == "+"


def test_write_table_suffixes(tmp_path: Path) -> None:
    df = pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]})
    tsv = tmp_path / "nested" / "out.tsv"
    csv = tmp_path / "out.csv"
    write_table(df, tsv)
    write_table(df, csv)
    assert tsv.read_text(encoding="utf-8").splitlines()[0] == "chrom\tstart\tend"
    assert csv.read_text(encoding="utf-8").splitlines()[0] == "chrom,start,end"
