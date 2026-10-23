"""Readers for GTF/GFF and peak tables."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

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


def read_gtf(path: Union[str, Path]) -> pd.DataFrame:
    """Read a GTF/GFF file and extract common gene/transcript attributes.

    Coordinates are converted to 0-based, half-open intervals.  GTF is the
    primary supported format; GFF3 is parsed on a best-effort basis.
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=GTF_COLUMNS,
        dtype={
            "chrom": "string",
            "source": "string",
            "feature": "string",
            "start": "string",
            "end": "string",
            "score": "string",
            "strand": "string",
            "frame": "string",
            "attributes": "string",
        },
        compression="infer",
        low_memory=False,
    )

    # GFF3 files may contain a ##FASTA section; sequence rows do not have
    # genomic coordinates.  Convert coordinates defensively and drop rows that
    # are not valid feature records.
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "feature", "start", "end"]).copy()

    if df.empty:
        raise ValueError(f"no records found in GTF: {path}")

    # GTF is 1-based, closed.  Convert to 0-based, half-open.
    df["start"] = df["start"].astype(np.int64) - 1
    df["end"] = df["end"].astype(np.int64)

    attrs = df["attributes"].fillna("")

    # Normalise feature names early so the GFF3 fallbacks can be
    # feature-aware.
    df["feature"] = (
        df["feature"]
        .astype("string")
        .str.strip()
        .str.lower()
        .str.replace("-", "_", regex=False)
    )
    df["feature"] = df["feature"].replace(_FEATURE_ALIASES)

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
        if target is None:
            continue
        if target in reserved and target != col:
            continue
        rename[col] = target
        reserved.add(target)
    return df.rename(columns=rename)


def read_peaks(
    path: Union[str, Path],
    fmt: Optional[str] = None,
) -> pd.DataFrame:
    """Read peaks from BED, TSV or CSV.

    BED files are interpreted as 0-based, half-open.  For tabular files, the
    columns ``chrom``/``start``/``end`` (case-insensitive) are required.
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.name.lower()
    if fmt is None:
        if suffix.endswith((".bed", ".bed.gz")):
            fmt = "bed"
        elif suffix.endswith((".csv", ".csv.gz")):
            fmt = "csv"
        else:
            fmt = "tsv"

    if fmt == "bed":
        raw = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            compression="infer",
            dtype={0: "string"},
        )
        if raw.shape[1] >= 1:
            first_col = raw.iloc[:, 0].astype(str).str.lower()
            raw = raw[~first_col.str.startswith(("track", "browser"))].copy()
        if raw.shape[1] < 3:
            raise ValueError(f"BED file must have at least 3 columns: {path}")

        n_extra = max(0, raw.shape[1] - 6)
        extra = raw.iloc[:, 6:].copy() if n_extra else None
        cols = ["chrom", "start", "end", "name", "score", "strand"]
        if raw.shape[1] < 6:
            cols = cols[: raw.shape[1]]
        raw = raw.iloc[:, :6].copy()
        raw.columns = cols
        if "name" not in raw:
            raw["name"] = "."
        if "score" not in raw:
            raw["score"] = "."
        if "strand" not in raw:
            raw["strand"] = "."
        if extra is not None:
            # Preserve extra BED columns rather than silently dropping them.
            extra.columns = [f"extra_{i + 1}" for i in range(n_extra)]
            raw = pd.concat(
                [raw.reset_index(drop=True), extra.reset_index(drop=True)],
                axis=1,
            )

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

    df = _normalise_peak_columns(df)
    required = {"chrom", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        # Last resort for BED-like files with unnamed first three columns.
        if len(df.columns) >= 3 and not required.intersection(df.columns):
            first = list(df.columns[:3])
            df = df.rename(columns={first[0]: "chrom", first[1]: "start", first[2]: "end"})
        else:
            raise ValueError(f"peak table is missing columns: {sorted(missing)}")

    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    if df["start"].isna().any() or df["end"].isna().any():
        raise ValueError("peak start/end contain non-numeric values")

    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    return df.reset_index(drop=True)


def write_table(df: pd.DataFrame, path: Union[str, Path]) -> None:
    """Write a DataFrame to TSV/CSV based on the filename suffix."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df.to_csv(path, sep=sep, index=False)
