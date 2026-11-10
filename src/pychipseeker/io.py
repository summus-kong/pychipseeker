"""Readers for GTF/GFF and peak tables."""

from __future__ import annotations

import bz2
import gzip
import lzma
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd

GTF_COLUMNS = [
    "chrom",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes",
]

_FEATURE_ALIASES = {
    "five_prime_utr": "five_prime_utr",
    "5utr": "five_prime_utr",
    "5_utr": "five_prime_utr",
    "5'utr": "five_prime_utr",
    "5'-utr": "five_prime_utr",
    "5'_utr": "five_prime_utr",
    "5' utr": "five_prime_utr",
    "three_prime_utr": "three_prime_utr",
    "3utr": "three_prime_utr",
    "3_utr": "three_prime_utr",
    "3'utr": "three_prime_utr",
    "3'-utr": "three_prime_utr",
    "3'_utr": "three_prime_utr",
    "3' utr": "three_prime_utr",
    "utr": "utr",
    "cds": "cds",
    "exon": "exon",
    "gene": "gene",
}


def _extract_gtf(series: pd.Series, key: str) -> pd.Series:
    return series.str.extract(
        rf"(?:^|;)\s*{re.escape(key)}\s+\"([^\"]*)\"",
        expand=False,
    )


def _extract_gff(series: pd.Series, key: str) -> pd.Series:
    return series.str.extract(
        rf"(?:^|;)\s*{re.escape(key)}=([^;]*)",
        expand=False,
    )


def _strip_feature_prefix(values: pd.Series) -> pd.Series:
    return values.astype("string").str.replace(
        r"^(gene|transcript|mRNA|mrna):",
        "",
        regex=True,
    )


class _GtfTextReader:
    """File-like text reader that stops at a GFF3 ``##FASTA`` section.

    pandas does not accept arbitrary line iterators, so this wrapper exposes
    ``read``/``readline`` while streaming from an already-opened text handle.
    Keeping the parser streaming avoids loading multi-gigabyte GENCODE GTFs
    into memory just to skip embedded genome sequence.
    """

    def __init__(self, handle) -> None:
        self._handle = handle
        self._buffer = ""
        self._eof = False

    def _read_line(self) -> None:
        line = self._handle.readline()
        if not line:
            self._eof = True
            return
        if line.startswith("##FASTA"):
            self._eof = True
            return
        self._buffer += line

    def _fill(self, size: int) -> None:
        while not self._eof and len(self._buffer) < size:
            self._read_line()

    def read(self, size: int = -1) -> str:
        if size is None or size < 0:
            while not self._eof:
                self._read_line()
            text, self._buffer = self._buffer, ""
            return text
        self._fill(int(size))
        text, self._buffer = self._buffer[: int(size)], self._buffer[int(size) :]
        return text

    def readline(self, size: int = -1) -> str:
        while "\n" not in self._buffer and not self._eof:
            self._read_line()
        index = self._buffer.find("\n")
        if index < 0:
            text, self._buffer = self._buffer, ""
            return text
        text, self._buffer = self._buffer[: index + 1], self._buffer[index + 1 :]
        return text

    def readable(self) -> bool:
        return True


def _normalise_feature_names(
    values: Union[str, Iterable[str]],
) -> set[str]:
    if isinstance(values, str):
        values = [values]
    names = {
        str(value).strip().lower().replace("-", "_")
        for value in values
        if str(value).strip()
    }
    return {_FEATURE_ALIASES.get(name, name) for name in names}


def read_gtf(
    path: Union[str, Path],
    *,
    features: Optional[Union[str, Iterable[str]]] = None,
    chunksize: Optional[int] = None,
) -> pd.DataFrame:
    """Read a GTF/GFF file and extract common gene/transcript attributes.

    Coordinates are converted to 0-based, half-open intervals.  GTF is the
    primary supported format; GFF3 is parsed on a best-effort basis, including
    comma-separated ``Parent`` values and an embedded ``##FASTA`` section.

    ``features`` optionally restricts the returned rows to one or more feature
    types (for example ``"exon"``).  Gene/transcript rows are still parsed
    internally to recover parent identifiers.  ``chunksize`` streams large
    files through pandas without changing the result.
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if chunksize is not None and int(chunksize) <= 0:
        raise ValueError("chunksize must be positive")
    feature_filter = None if features is None else _normalise_feature_names(features)

    dtype_map = {
        "chrom": "string",
        "source": "string",
        "feature": "string",
        "start": "string",
        "end": "string",
        "score": "string",
        "strand": "string",
        "frame": "string",
        "attributes": "string",
    }
    frames: list[pd.DataFrame] = []
    with _open_text(path) as handle:
        reader = _GtfTextReader(handle)
        parsed = pd.read_csv(
            reader,
            sep="\t",
            comment="#",
            header=None,
            names=GTF_COLUMNS,
            dtype=dtype_map,
            low_memory=False,
            chunksize=None if chunksize is None else int(chunksize),
        )
        chunks = [parsed] if isinstance(parsed, pd.DataFrame) else parsed
        for chunk in chunks:
            chunk = chunk.copy()
            chunk["start"] = pd.to_numeric(chunk["start"], errors="coerce")
            chunk["end"] = pd.to_numeric(chunk["end"], errors="coerce")
            chunk = chunk.dropna(
                subset=["chrom", "feature", "start", "end"]
            ).copy()
            if chunk.empty:
                continue
            chunk["start"] = chunk["start"].astype(np.int64) - 1
            chunk["end"] = chunk["end"].astype(np.int64)
            chunk["feature"] = (
                chunk["feature"]
                .astype("string")
                .str.strip()
                .str.lower()
                .str.replace("-", "_", regex=False)
            )
            chunk["feature"] = chunk["feature"].replace(_FEATURE_ALIASES)
            if feature_filter is not None:
                keep = feature_filter | {"gene", "transcript", "mrna"}
                chunk = chunk[chunk["feature"].isin(keep)].copy()
            if not chunk.empty:
                frames.append(chunk)
    if not frames:
        raise ValueError(f"no records found in GTF: {path}")
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    attrs = df["attributes"].fillna("")

    # GTF-style attributes.
    df["gene_id"] = _extract_gtf(attrs, "gene_id")
    df["transcript_id"] = _extract_gtf(attrs, "transcript_id")
    df["gene_name"] = _extract_gtf(attrs, "gene_name")
    df["gene_type"] = _extract_gtf(attrs, "gene_type")
    df["transcript_name"] = _extract_gtf(attrs, "transcript_name")
    df["transcript_type"] = _extract_gtf(attrs, "transcript_type")
    df["exon_number"] = _extract_gtf(attrs, "exon_number")
    df["exon_id"] = _extract_gtf(attrs, "exon_id")

    # GFF3-style fallbacks.  In GFF3, gene rows use ID=gene:..., transcript
    # rows use ID=transcript:...;Parent=gene:..., and exon/UTR/CDS rows use
    # Parent=transcript:... .  Do not blindly use ID as gene_id on child rows:
    # that would create one "gene" per exon.
    tx_features = {"mrna", "transcript"}
    child_features = {"exon", "five_prime_utr", "three_prime_utr", "utr", "cds"}

    mask = df["gene_id"].isna() & (df["feature"] == "gene")
    if mask.any():
        df.loc[mask, "gene_id"] = _extract_gff(attrs[mask], "ID")

    tx_mask = df["feature"].isin(tx_features)
    mask = df["gene_id"].isna() & tx_mask
    if mask.any():
        df.loc[mask, "gene_id"] = _extract_gff(attrs[mask], "Parent")

    mask = df["transcript_id"].isna() & tx_mask
    if mask.any():
        df.loc[mask, "transcript_id"] = _extract_gff(attrs[mask], "ID")

    mask = df["transcript_name"].isna() & tx_mask
    if mask.any():
        df.loc[mask, "transcript_name"] = _extract_gff(attrs[mask], "Name")

    mask = df["exon_number"].isna()
    if mask.any():
        df.loc[mask, "exon_number"] = _extract_gff(attrs[mask], "rank")
    mask = df["exon_number"].isna()
    if mask.any():
        df.loc[mask, "exon_number"] = _extract_gff(attrs[mask], "exon_number")

    child_mask = df["feature"].isin(child_features)
    mask = df["transcript_id"].isna() & child_mask
    if mask.any():
        df.loc[mask, "transcript_id"] = _extract_gff(attrs[mask], "Parent")

    # GFF3 allows a feature to belong to several parents:
    # ``Parent=transcript:tx1,transcript:tx2``.  Duplicate the row once per
    # parent so each transcript gets its own exon/UTR/CDS record.
    for column, feature_set in (
        ("transcript_id", child_features),
        ("gene_id", tx_features),
    ):
        row_mask = df["feature"].isin(feature_set)
        values = df[column].astype("string")
        multiple = row_mask & values.str.contains(",", regex=False, na=False)
        if multiple.any():
            df[column] = df[column].astype(object)
            df.loc[multiple, column] = df.loc[multiple, column].str.split(",")
            df = df.explode(column, ignore_index=True)
            df[column] = df[column].astype("string").str.strip()
            attrs = df["attributes"].fillna("")

    mask = df["gene_name"].isna() & (df["feature"] == "gene")
    if mask.any():
        df.loc[mask, "gene_name"] = _extract_gff(attrs[mask], "Name")

    mask = df["gene_name"].isna()
    if mask.any():
        df.loc[mask, "gene_name"] = _extract_gff(attrs[mask], "gene")

    mask = df["gene_type"].isna()
    if mask.any():
        vals = _extract_gtf(attrs[mask], "gene_biotype")
        df.loc[mask, "gene_type"] = vals

    mask = df["transcript_type"].isna()
    if mask.any():
        vals = _extract_gtf(attrs[mask], "transcript_biotype")
        df.loc[mask, "transcript_type"] = vals

    mask = df["gene_type"].isna()
    if mask.any():
        df.loc[mask, "gene_type"] = _extract_gff(attrs[mask], "gene_biotype")

    mask = df["transcript_type"].isna()
    if mask.any():
        df.loc[mask, "transcript_type"] = _extract_gff(attrs[mask], "transcript_biotype")

    for col in ("gene_id", "transcript_id"):
        df[col] = _strip_feature_prefix(df[col])

    # If child rows only carry a transcript ID (common in GFF3), recover the
    # gene ID from transcript/mRNA rows.
    tx2gene = (
        df.dropna(subset=["transcript_id", "gene_id"])
        .drop_duplicates("transcript_id")
        .set_index("transcript_id")["gene_id"]
    )
    if not tx2gene.empty:
        mask = df["gene_id"].isna() & df["transcript_id"].notna()
        if mask.any():
            df.loc[mask, "gene_id"] = df.loc[mask, "transcript_id"].map(tx2gene)

    df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
    df["strand"] = df["strand"].fillna(".")

    if feature_filter is not None:
        df = df[df["feature"].isin(feature_filter)].reset_index(drop=True)
        if df.empty:
            raise ValueError(
                f"no records matching features={sorted(feature_filter)} in {path}"
            )
    return df


def _normalise_peak_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename = {}
    reserved = set(df.columns)
    for col in df.columns:
        key = str(col).strip().lower()
        target = None
        if key in {"chr", "chrom", "chromosome", "seqnames"}:
            target = "chrom"
        elif key in {"start", "chromstart"}:
            target = "start"
        elif key in {"end", "chromend", "stop"}:
            target = "end"
        elif key in {"name", "peak"}:
            target = "name"
        elif key in {"peak_id", "peakid"}:
            target = "peak_id"
        elif key == "score":
            target = "score"
        elif key == "strand":
            target = "strand"
        elif key in {"signalvalue", "signal_value", "signal"}:
            target = "signalValue"
        elif key in {"pvalue", "p_value", "pval"}:
            target = "pValue"
        elif key in {"qvalue", "q_value", "qval", "fdr"}:
            target = "qValue"
        if target is None:
            continue
        if target in reserved and target != col:
            continue
        rename[col] = target
        reserved.add(target)
    return df.rename(columns=rename)


_PEAK_VARIANTS = {
    "narrowpeak": ["signalValue", "pValue", "qValue", "peak"],
    "broadpeak": ["signalValue", "pValue", "qValue"],
    "gappedpeak": [
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
        "signalValue",
        "pValue",
        "qValue",
    ],
}


def _peak_variant(name: str) -> Optional[str]:
    """Return the BED variant implied by a filename, if any."""

    lower = name.lower()
    for suffix in (".gz", ".bz2", ".xz"):
        if lower.endswith(suffix):
            lower = lower[: -len(suffix)]
            break
    for variant in _PEAK_VARIANTS:
        if lower.endswith("." + variant) or lower.endswith("." + variant + "peak"):
            return variant
    if lower.endswith((".bed", ".bed6", ".bed12")):
        return "bed"
    return None


def _is_bed_header(values: Sequence[object]) -> bool:
    """True if the first BED row looks like a column header."""

    if len(values) < 3:
        return False
    first = str(values[0]).strip().lower()
    second = str(values[1]).strip().lower()
    third = str(values[2]).strip().lower()
    return first in {"chrom", "chr", "chromosome", "#chrom"} and second in {
        "start",
        "chromstart",
    } and third in {"end", "chromend", "stop"}


def _open_text(path: Path):
    lower = path.name.lower()
    if lower.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if lower.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8", errors="replace")
    if lower.endswith(".xz"):
        return lzma.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _bed_header_info(path: Path) -> tuple[Optional[list[str]], int]:
    """Return ``(header, skiprows)`` for a BED-like file.

    The scan stops at the first non-comment/non-track line.  A header such as
    ``chrom start end`` (optionally ``#chrom``) is recognised even when it has
    fewer columns than the data rows, which pandas' CSV sniffer would
    otherwise reject.
    """

    with _open_text(path) as fh:
        for index, line in enumerate(fh):
            text = line.strip()
            if not text:
                continue
            lower = text.lower()
            if lower.startswith(("track", "browser")):
                continue
            fields = text.split("\t")
            if _is_bed_header(fields):
                return [str(field).strip() for field in fields], index + 1
            if text.startswith("#"):
                continue
            return None, index
    return None, 0


def _coerce_peaks(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise column names and coerce coordinates for a peak frame."""

    df = _normalise_peak_columns(df)
    required = {"chrom", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        # Last resort for BED-like files with unnamed first three columns.
        if len(df.columns) >= 3 and not required.intersection(df.columns):
            first = list(df.columns[:3])
            df = df.rename(
                columns={first[0]: "chrom", first[1]: "start", first[2]: "end"}
            )
        else:
            raise ValueError(f"peak table is missing columns: {sorted(missing)}")

    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    if df["start"].isna().any() or df["end"].isna().any():
        raise ValueError("peak start/end contain non-numeric values")

    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    return df.reset_index(drop=True)


def _bed_frame(rows: Sequence[Sequence[object]]) -> pd.DataFrame:
    """Build a BED-like frame from ``(chrom, start, end, *fields)`` rows."""

    if not rows:
        raise ValueError("no peak records were found")
    n_cols = max(len(row) for row in rows)
    standard = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    columns = standard[: min(n_cols, len(standard))]
    columns += [f"extra_{i + 1}" for i in range(n_cols - len(columns))]
    padded = [list(row) + [None] * (n_cols - len(row)) for row in rows]
    return pd.DataFrame(padded, columns=columns)


def _read_bigbed(path: Path) -> pd.DataFrame:
    """Read a bigBed file via ``bigBedToBed`` or the optional pybigtools."""

    converter = shutil.which("bigBedToBed") or shutil.which("bigbedtobed")
    if converter is not None:
        tmpdir = tempfile.mkdtemp(prefix="pychipseeker_bigbed_")
        try:
            bed_path = Path(tmpdir) / "peaks.bed"
            subprocess.run(
                [converter, str(path), str(bed_path)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            return read_peaks(bed_path, fmt="bed")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    try:
        import pybigtools  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "reading .bb/.bigBed files requires the UCSC bigBedToBed utility "
            "on PATH or the optional pybigtools package "
            "(`pip install pychipseeker[bigbed]`)"
        ) from exc

    handle = pybigtools.open(str(path))
    rows: list[list[object]] = []
    try:
        chroms = handle.chroms()
        for chrom, size in chroms.items():
            records = handle.records(str(chrom), 0, int(size))
            if records is None:
                continue
            for record in records:
                fields = list(record)
                if len(fields) >= 3 and isinstance(fields[2], str):
                    rest = [value for value in fields[2].split("\t") if value != ""]
                else:
                    rest = fields[2:]
                rows.append([str(chrom), int(fields[0]), int(fields[1]), *rest])
    finally:
        close = getattr(handle, "close", None)
        if close is not None:
            close()

    return _coerce_peaks(_bed_frame(rows))


def _read_bam(path: Path) -> pd.DataFrame:
    """Read aligned reads from a BAM/CRAM file with pysam."""

    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "reading .bam/.cram files requires the optional pysam package "
            "(`pip install pychipseeker[bam]`; on Windows use "
            "`conda install -c bioconda pysam`)"
        ) from exc

    rows: list[dict[str, object]] = []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_start is None:
                continue
            start = int(read.reference_start)
            end = read.reference_end
            if end is None:
                end = start + max(1, int(read.query_length or 1))
            rows.append(
                {
                    "chrom": str(read.reference_name),
                    "start": start,
                    "end": int(end),
                    "name": str(read.query_name or "."),
                    "score": int(read.mapping_quality),
                    "strand": "-" if read.is_reverse else "+",
                }
            )
    if not rows:
        raise ValueError(f"no aligned reads were found in {path}")
    return _coerce_peaks(pd.DataFrame.from_records(rows))


def read_peaks(
    path: Union[str, Path],
    fmt: Optional[str] = None,
) -> pd.DataFrame:
    """Read peaks from BED, narrowPeak/broadPeak/gappedPeak, bigBed, BAM,
    TSV or CSV.

    BED files are interpreted as 0-based, half-open.  The ENCODE peak
    variants get their standard column names (``signalValue``, ``pValue``,
    ``qValue``, ``peak``/``blockCount``/...); plain BED files keep extra
    columns as ``extra_1`` ... ``extra_n``.  A header row such as
    ``chrom start end`` is skipped, and ``track``/``browser`` lines are
    ignored.

    ``.bb``/``.bigBed`` files need either the UCSC ``bigBedToBed`` utility or
    the optional ``pybigtools`` package; ``.bam``/``.cram`` files need the
    optional ``pysam`` package.  For tabular files, the columns
    ``chrom``/``start``/``end`` (case-insensitive) are required.  ``fmt`` may
    be ``bed``, ``narrowpeak``, ``broadpeak``, ``gappedpeak``, ``bigbed``,
    ``bam``, ``tsv`` or ``csv`` to override the suffix.
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.name.lower()
    if fmt is not None and str(fmt).lower() in {"bam", "cram"}:
        return _read_bam(path)
    if fmt is not None and str(fmt).lower() in {"bigbed", "bb"}:
        return _read_bigbed(path)
    if fmt is None:
        if suffix.endswith((".bam", ".cram")):
            return _read_bam(path)
        if suffix.endswith((".bb", ".bigbed")):
            return _read_bigbed(path)

    explicit = None if fmt is None else str(fmt).lower()
    explicit_variant: Optional[str] = None
    if explicit in _PEAK_VARIANTS:
        explicit_variant = explicit
        explicit = "bed"
    elif explicit in {"bed", "csv", "tsv"}:
        explicit_variant = None
    elif explicit is not None:
        raise ValueError(
            "fmt must be one of bed, narrowpeak, broadpeak, gappedpeak, "
            "bigbed, bam, tsv or csv"
        )

    variant = explicit_variant or _peak_variant(suffix)
    if fmt is None:
        if variant is not None:
            fmt = "bed"
        elif suffix.endswith((".csv", ".csv.gz")):
            fmt = "csv"
        else:
            fmt = "tsv"
    elif explicit in {"tsv", "csv"}:
        fmt = explicit
    elif explicit == "bed":
        fmt = "bed"

    if fmt == "bed":
        header, skiprows = _bed_header_info(path)
        raw = pd.read_csv(
            path,
            sep="\t",
            header=None,
            skiprows=skiprows,
            comment="#",
            compression="infer",
            dtype={0: "string"},
        )
        if raw.shape[1] >= 1:
            first_col = raw.iloc[:, 0].astype(str).str.strip().str.lower()
            raw = raw[~first_col.str.startswith(("track", "browser"))].copy()
        if raw.shape[1] < 3:
            raise ValueError(f"BED file must have at least 3 columns: {path}")

        n_cols = raw.shape[1]
        base = ["chrom", "start", "end", "name", "score", "strand"]
        tail = _PEAK_VARIANTS.get(variant or "", [])
        standard = base + tail
        if header is None:
            columns = standard[: min(n_cols, len(standard))]
            n_extra = max(0, n_cols - len(columns))
            columns += [f"extra_{i + 1}" for i in range(n_extra)]
        else:
            columns = []
            used: set[str] = set()
            for column in header[:n_cols]:
                key = column.strip().lower().lstrip("#")
                if key in {"chr", "chrom", "chromosome", "seqnames"}:
                    name = "chrom"
                elif key in {"start", "chromstart"}:
                    name = "start"
                elif key in {"end", "chromend", "stop"}:
                    name = "end"
                elif key in {"name", "peak_id", "peakid"}:
                    name = "name"
                elif key == "score":
                    name = "score"
                elif key == "strand":
                    name = "strand"
                elif key in {"signalvalue", "signal_value", "signal"}:
                    name = "signalValue"
                elif key in {"pvalue", "p_value", "pval"}:
                    name = "pValue"
                elif key in {"qvalue", "q_value", "qval", "fdr"}:
                    name = "qValue"
                elif key in {"peak", "summit"}:
                    name = "peak"
                elif key:
                    name = column
                else:
                    name = f"extra_{len(used) + 1}"
                if name in used:
                    name = f"extra_{len(used) + 1}"
                used.add(name)
                columns.append(name)
            if len(columns) < n_cols:
                columns += standard[len(columns) : n_cols]
            columns += [
                f"extra_{i + 1}" for i in range(n_cols - len(columns))
            ]

        raw.columns = columns
        raw = raw.reset_index(drop=True)
        df = raw
    else:
        sep = "," if fmt == "csv" else "\t"
        df = pd.read_csv(path, sep=sep, header=0, compression="infer")
        df = _normalise_peak_columns(df)
        required = {"chrom", "start", "end"}
        if not required.issubset(df.columns):
            # Try a headerless file as a convenience.
            df = pd.read_csv(path, sep=sep, header=None, compression="infer")
            if df.shape[1] >= 3:
                first = list(df.columns[:3])
                df = df.rename(
                    columns={first[0]: "chrom", first[1]: "start", first[2]: "end"}
                )
            df = _normalise_peak_columns(df)
            if not required.issubset(df.columns):
                raise ValueError(
                    f"peak table must contain chrom/start/end columns: {path}"
                )

    return _coerce_peaks(df)


def write_table(df: pd.DataFrame, path: Union[str, Path]) -> None:
    """Write a DataFrame to TSV/CSV based on the filename suffix."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df.to_csv(path, sep=sep, index=False)
