"""Tests for scripts/build_gsminfo.py (offline gsminfo rebuild)."""

import csv
import gzip
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_gsminfo  # noqa: E402


def _write_ucsc_release(path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["species", "ucsc_version", "release_date", "release_name"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "species": "Homo sapiens",
                    "ucsc_version": "hg38",
                    "release_date": "Dec. 2013",
                    "release_name": "GRCh38",
                },
                {
                    "species": "Homo sapiens",
                    "ucsc_version": "hg19",
                    "release_date": "Feb. 2009",
                    "release_name": "GRCh37 (hg19)",
                },
                {
                    "species": "Mus musculus",
                    "ucsc_version": "mm9",
                    "release_date": "Jul. 2009",
                    "release_name": "Build 37",
                },
            ]
        )


@pytest.fixture
def geometadb(tmp_path):
    """A miniature GEOmetadb.sqlite mirroring the upstream query surface."""

    path = tmp_path / "GEOmetadb.sqlite"
    con = sqlite3.connect(path)
    cur = con.cursor()
    cur.execute("CREATE TABLE gpl (gpl TEXT, technology TEXT)")
    cur.execute(
        "CREATE TABLE gsm ("
        "gsm TEXT, series_id TEXT, gpl TEXT, organism_ch1 TEXT, title TEXT,"
        " characteristics_ch1 TEXT, source_name_ch1 TEXT, extract_protocol_ch1 TEXT,"
        " description TEXT, data_processing TEXT, submission_date TEXT,"
        " supplementary_file TEXT)"
    )
    cur.execute("CREATE TABLE gse (gse TEXT, pubmed_id TEXT)")

    cur.executemany(
        "INSERT INTO gpl VALUES (?, ?)",
        [
            ("GPL1", "high-throughput sequencing"),
            ("GPL2", "array"),  # filtered out by the technology query
            ("GPL3", "high-throughput sequencing"),
        ],
    )
    cur.executemany(
        "INSERT INTO gsm VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                "GSM1", "GSE1", "GPL1", "Homo sapiens", "t1",
                "c1", "s1", "e1", "d1",
                "Genome build: hg19", "2015-01-01",
                "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1_bed.bed.gz",
            ),
            (
                "GSM2", "GSE1", "GPL1", "Homo sapiens\u00e9", "t2",
                "c2", "s2", "e2", "d2",
                "genome without a release hint", "2015-01-02",
                "https://x/GSM2_rawfastq.gz,"
                "ftp://x/GSM2_bedGraph.gz",  # comma-clumped; only bedGraph kept
            ),
            (
                "GSM3", "GSE2", "GPL3", "Mus musculus", "t3",
                "c3", "s3", "e3", "d3",
                "no version here either", "2015-01-03",
                "ftp://x/GSM3_mm9_align.bed.gz",  # mm9 picked from the URL
            ),
            (
                "GSM4", "GSE9", "GPL2", "Homo sapiens", "t4",
                "c4", "s4", "e4", "d4",
                "Genome build: hg19", "2015-01-04",
                "ftp://x/GSM4_bed.bed.gz",  # GPL2 -> must be filtered out
            ),
        ],
    )
    cur.executemany(
        "INSERT INTO gse VALUES (?, ?)",
        [("GSE1", "123"), ("GSE2", "456"), ("GSE3", None)],
    )
    con.commit()
    con.close()
    return path


def _ucsc(tmp_path):
    path = tmp_path / "ucsc_release.csv"
    _write_ucsc_release(path)
    return path


def test_build_gsminfo_end_to_end(tmp_path, geometadb):
    out = tmp_path / "gsminfo.csv.gz"
    rows = build_gsminfo.build_gsminfo(
        geometadb, out, ucsc_release_path=_ucsc(tmp_path)
    )
    assert rows == 3  # GSM1 + GSM2(bedGraph) + GSM3; GSM4 on GPL2 filtered out

    with gzip.open(out, "rt", encoding="utf-8") as handle:
        data = list(csv.DictReader(handle))
    assert [r["gsm"] for r in data] == ["GSM1", "GSM2", "GSM3"]
    assert set(build_gsminfo.OUTPUT_COLUMNS) == set(data[0].keys())

    # genomeVersion: data_processing first, then the URL
    assert data[0]["genomeVersion"] == "hg19"
    assert data[1]["genomeVersion"] == ""  # no hint anywhere
    assert data[2]["genomeVersion"] == "mm9"

    # supplementary files: comma-clumped entries split and pattern-filtered
    assert data[0]["supplementary_file"].endswith("GSM1_bed.bed.gz")
    assert data[1]["supplementary_file"].endswith("GSM2_bedGraph.gz")

    # "_ch1" renames and plain (non-bundled) passthrough columns
    assert data[0]["organism"] == "Homo sapiens"
    assert data[0]["series_id"] == "GSE1"
    assert data[0]["gpl"] == "GPL1"

    # pubmed joined through the gse table; unknown series -> empty
    assert data[0]["pubmed_id"] == "123"
    assert data[2]["pubmed_id"] == "456"


def test_build_gsminfo_ascii_clean(tmp_path, geometadb):
    out = tmp_path / "gsminfo.csv"
    build_gsminfo.build_gsminfo(
        geometadb, out, ucsc_release_path=_ucsc(tmp_path)
    )
    with open(out, encoding="utf-8") as handle:
        data = list(csv.DictReader(handle))
    assert all("é" not in r["organism"] for r in data)  # é stripped like iconv


def test_build_gsminfo_keep_unicode(tmp_path, geometadb):
    out = tmp_path / "gsminfo_unicode.csv"
    build_gsminfo.build_gsminfo(
        geometadb, out, ucsc_release_path=_ucsc(tmp_path), keep_unicode=True
    )
    with open(out, encoding="utf-8") as handle:
        data = list(csv.DictReader(handle))
    assert any("é" in r["organism"] for r in data)


def test_build_gsminfo_cli(tmp_path, geometadb, capsys):
    out = tmp_path / "cli.csv.gz"
    rc = build_gsminfo.main(
        [str(geometadb), "-o", str(out), "--ucsc-release", str(_ucsc(tmp_path)), "--quiet"]
    )
    assert rc == 0
    assert out.exists()

    # missing file -> argparse error (exit 2)
    with pytest.raises(SystemExit):
        build_gsminfo.main([str(tmp_path / "missing.sqlite"), "-o", str(out)])


def test_ucsc_release_mapping_requires_species_match(tmp_path, geometadb):
    # an organism absent from ucsc_release yields empty genomeVersion
    con = sqlite3.connect(geometadb)
    con.execute(
        "UPDATE gsm SET organism_ch1='Gallus gallus' WHERE gsm='GSM1'"
    )
    con.commit()
    con.close()
    out = tmp_path / "gsminfo2.csv"
    build_gsminfo.build_gsminfo(geometadb, out, ucsc_release_path=_ucsc(tmp_path))
    with open(out, encoding="utf-8") as handle:
        data = list(csv.DictReader(handle))
    assert data[0]["gsm"] == "GSM1"
    assert data[0]["genomeVersion"] == ""
