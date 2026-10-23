"""Rebuild a ChIPseeker-style ``gsminfo`` snapshot *offline* from GEOmetadb.

This is the Python equivalent of ChIPseeker's internal ``prepareGSMInfo()``
(R, see ``R/GEO.R``).  Upstream runs it as a maintainer batch job: it
downloads the huge ``GEOmetadb.sqlite``, then re-fetches every GPL/GSM soft
file with GEOquery to find the supplementary BED files.  Here the same
result is produced **offline and with zero dependencies**: the user downloads
the GEOmetadb SQLite once (https://gbnci-abcc.ncifcrf.gov/geo/ or the URL in
``GEOmetadb::getSQLiteFile()``), and this script queries it directly — the
``gsm`` table already carries ``supplementary_file``, so no per-sample HTTP
requests are needed.

The output CSV has exactly the columns of the bundled
``pychipseeker/data/gsminfo.csv.gz`` and can be handed to every GEO helper
through ``gsminfo=...``::

    python scripts/build_gsminfo.py GEOmetadb.sqlite -o gsminfo.csv.gz
    python -c "import pychipseeker as pcs; \
               print(pcs.getGEOInfo('hg19', gsminfo='gsminfo.csv.gz').head())"

Usage::

    python scripts/build_gsminfo.py GEOmetadb.sqlite [-o gsminfo.csv.gz]
        [--ucsc-release ucsc_release.csv.gz] [--keep-unicode] [--quiet]

Only files whose URL matches ChIPseeker's ``getGSMsuppFile`` patterns
(``bed.gz``, ``Peak.gz``, ``bedGraph.gz`` — R ``grep`` regex semantics) are
kept.  ``genomeVersion`` is resolved with the upstream ``ucsc_release``
table (species match, then search in ``data_processing``, then in the file
URL, first match wins; empty string when nothing matches) exactly like R's
``getGenomicVersion``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sqlite3
import sys
from pathlib import Path

# Column order mirrors pychipseeker/data/gsminfo.csv.gz (upstream gsminfo.rda).
OUTPUT_COLUMNS = [
    "series_id",
    "gsm",
    "gpl",
    "organism",
    "title",
    "characteristics",
    "source_name",
    "extract_protocol",
    "description",
    "data_processing",
    "submission_date",
    "supplementary_file",
    "genomeVersion",
    "pubmed_id",
]

# R getGSMsuppFile: grep("bed.gz"), grep("Peak.gz"), grep("bedGraph.gz")
SUPP_PATTERNS = ("bed.gz", "Peak.gz", "bedGraph.gz")

DEFAULT_UCSC_RELEASE = Path(__file__).with_name("ucsc_release.csv.gz")


def _read_ucsc_release(path: Path | None = None) -> dict[str, list[str]]:
    """Load the upstream named-vector ``ucsc_release`` table.

    Returns ``{species: [ucsc_version, ...]}`` in file order (order matters:
    R keeps the first match).
    """

    path = Path(path) if path is not None else DEFAULT_UCSC_RELEASE
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found - pass --ucsc-release <ucsc_release.csv[.gz]>"
        )
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for needed in ("species", "ucsc_version"):
            if needed not in (reader.fieldnames or []):
                raise ValueError(f"{path} is missing the '{needed}' column")
        versions: dict[str, list[str]] = {}
        for row in reader:
            species = (row.get("species") or "").strip()
            version = (row.get("ucsc_version") or "").strip()
            if species and version:
                versions.setdefault(species, []).append(version)
    return versions


def _get_genomic_version(
    ucsc_release: dict[str, list[str]],
    data_processing: str,
    organism: str,
    supplementary_file: str,
) -> str:
    """ChIPseeker ``getGenomicVersion``: grep in data_processing, then URL."""

    versions = ucsc_release.get(organism, [])
    if not versions:
        return ""
    for haystack in (data_processing, supplementary_file):
        matches = [v for v in versions if re.search(v, haystack)]
        if matches:
            return matches[0]
    return ""


def _iter_supp_files(supplementary_file: str) -> list[str]:
    """Split GEOmetadb's (possibly comma-clumped) supplementary_file column."""

    if not supplementary_file:
        return []
    return [u.strip() for u in re.split(r",\s*", supplementary_file) if u.strip()]


def _matches_supp_patterns(url: str) -> bool:
    return any(re.search(p, url) for p in SUPP_PATTERNS)


def _ascii_clean(value: str) -> str:
    """R iconv(latin1 -> ASCII, sub=\"\"): drop non-ASCII characters."""

    if not value:
        return value
    return value.encode("ascii", "ignore").decode("ascii")


def _query_high_throughput_gpl(con: sqlite3.Connection) -> list[tuple[str, str]]:
    cur = con.execute("SELECT gpl, technology FROM gpl")
    gpls: list[tuple[str, str]] = []
    for gpl_id, technology in cur:
        if technology == "high-throughput sequencing" and gpl_id:
            gpls.append((str(gpl_id), str(gpl_id)))
    return gpls


def _query_gsm_rows(con, keys: list[tuple[str, str]]):
    """One batched query for every GSM of all high-throughput GPLs.

    Equivalent to upstream's per-GSM GEOquery lookups: gsm.organism_ch1 ->
    organism, characteristics_ch1 -> characteristics, source_name_ch1 ->
    source_name (the same "_ch1" rename R applies).
    """

    placeholders = ",".join("?" for _ in keys)
    cur = con.execute(
        "SELECT gsm, series_id, gpl, organism_ch1, title, characteristics_ch1,"
        " source_name_ch1, extract_protocol_ch1, description, data_processing,"
        " submission_date, supplementary_file"
        f" FROM gsm WHERE gpl IN ({placeholders})",
        [gpl_id for gpl_id, _ in keys],
    )
    return cur


def build_gsminfo(
    sqlite_path: Path,
    output: Path,
    ucsc_release_path: Path | None = None,
    keep_unicode: bool = False,
    quiet: bool = False,
) -> int:
    """Build a gsminfo CSV from a local GEOmetadb SQLite. Returns row count."""

    ucsc_release = _read_ucsc_release(ucsc_release_path)
    con = sqlite3.connect(f"file:{sqlite_path}", uri=True)
    try:
        gpls = _query_high_throughput_gpl(con)
        if not quiet:
            print(f"{len(gpls)} high-throughput sequencing GPLs", flush=True)
        pubmed: dict[str, str] = {}
        series_seen: list[str] = []
        rows: list[dict[str, str]] = []
        for raw in _query_gsm_rows(con, gpls):
            (
                gsm, series_id, gpl, organism, title, characteristics,
                source_name, extract_protocol, description, data_processing,
                submission_date, supplementary_file,
            ) = ["" if v is None else str(v) for v in raw]
            series_seen.append(series_id)
            for url in _iter_supp_files(supplementary_file):
                if not _matches_supp_patterns(url):
                    continue
                rows.append(
                    {
                        "series_id": series_id,
                        "gsm": gsm,
                        "gpl": gpl,
                        "organism": organism,
                        "title": title,
                        "characteristics": characteristics,
                        "source_name": source_name,
                        "extract_protocol": extract_protocol,
                        "description": description,
                        "data_processing": data_processing,
                        "submission_date": submission_date,
                        "supplementary_file": url,
                        "genomeVersion": _get_genomic_version(
                            ucsc_release, data_processing, organism, url
                        ),
                        "pubmed_id": "",
                    }
                )
        # pubmed_id join, mirroring the gse-table lookup in upstream code.
        unique_series = list(dict.fromkeys(series_seen))
        for i in range(0, len(unique_series), 500):
            chunk = unique_series[i : i + 500]
            placeholders = ",".join("?" for _ in chunk)
            cur = con.execute(
                f"SELECT gse, pubmed_id FROM gse WHERE gse IN ({placeholders})",
                chunk,
            )
            for gse, pid in cur:
                if gse and (pid is not None) and str(pid) and gse not in pubmed:
                    pubmed[str(gse)] = str(pid)
    finally:
        con.close()

    for row in rows:
        row["pubmed_id"] = pubmed.get(row["series_id"], "")
        if not keep_unicode:
            for key, value in row.items():
                row[key] = _ascii_clean(value)

    output = Path(output)
    write = gzip.open if output.suffix == ".gz" else open
    output.parent.mkdir(parents=True, exist_ok=True)
    with write(output, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    if not quiet:
        print(f"{len(rows)} supplementary files written to {output}", flush=True)
    return len(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("sqlite", help="path to GEOmetadb.sqlite")
    parser.add_argument(
        "-o", "--output", default="gsminfo.csv.gz",
        help="output CSV (.csv or .csv.gz, default: gsminfo.csv.gz)",
    )
    parser.add_argument(
        "--ucsc-release", default=None,
        help="ucsc_release CSV (default: scripts/ucsc_release.csv.gz bundled "
             "in this folder, or pass a copy elsewhere)",
    )
    parser.add_argument(
        "--keep-unicode", action="store_true",
        help="keep non-ASCII text (upstream strips it with iconv)",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    sqlite_path = Path(args.sqlite)
    if not sqlite_path.exists():
        parser.error(f"{sqlite_path} not found")
    try:
        count = build_gsminfo(
            sqlite_path,
            Path(args.output),
            ucsc_release_path=(
                Path(args.ucsc_release) if args.ucsc_release else None
            ),
            keep_unicode=args.keep_unicode,
            quiet=args.quiet,
        )
    except sqlite3.OperationalError as exc:
        print(f"error: {exc}\nIs {sqlite_path} really a GEOmetadb.sqlite?",
              file=sys.stderr)
        return 2
    if count == 0:
        print("warning: no BED-like supplementary files matched", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
