"""ChIPseeker-compatible peak annotation.

This module implements the core annotation semantics of R/Bioconductor
ChIPseeker using only numpy/pandas:

* transcript-level (default) and gene-level nearest features;
* TSS-point versus full-feature overlap;
* strand-aware, signed distance to the TSS calculated from both peak ends;
* ChIPseeker's annotation priority, promoter bins, downstream bins and
  "Distal Intergenic" category;
* annotation/exon/intron alignment so that the reported feature is the one
  that supplied the genomic annotation;
* many-to-many peak-to-gene mapping (``seq2gene``) and flanking-gene output.

Coordinates in the public Python API are 0-based, half-open (BED convention).
ChIPseeker/GRanges coordinates are 1-based, closed; conversion is done
internally only where needed to reproduce R's distance formulas.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd

from .intervals import BinnedIntervalIndex, SpanIndex, TSSIndex
from .io import read_gtf, read_peaks

Number = Union[int, float]

_CORE_INTERVAL_COLUMNS = [
    "gene_id",
    "gene_name",
    "transcript_id",
    "transcript_name",
    "chrom",
    "start",
    "end",
    "strand",
]

_DEFAULT_PRIORITY = (
    "Promoter",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Downstream",
    "Intergenic",
)

_PRIORITY_ALIASES = {
    "promoter": "Promoter",
    "5utr": "5UTR",
    "5' utr": "5UTR",
    "5'utr": "5UTR",
    "five_prime_utr": "5UTR",
    "five prime utr": "5UTR",
    "3utr": "3UTR",
    "3' utr": "3UTR",
    "3'utr": "3UTR",
    "three_prime_utr": "3UTR",
    "three prime utr": "3UTR",
    "exon": "Exon",
    "intron": "Intron",
    "downstream": "Downstream",
    "intergenic": "Intergenic",
    "distal intergenic": "Intergenic",
}

_DETAIL_KEYS = {
    "Promoter": "Promoter",
    "5UTR": "fiveUTR",
    "3UTR": "threeUTR",
    "Exon": "Exon",
    "Intron": "Intron",
}


def _normalise_strand(value: object) -> str:
    text = str(value)
    return text if text in {"+", "-"} else "*"


def _format_kb(value: int) -> str:
    value = int(value)
    if value % 1000 == 0:
        return f"{value // 1000}kb"
    return f"{value / 1000:g}kb"


def _signed_distance(pos: int, tss: int, strand: str) -> int:
    """Signed distance from a 0-based position to a 0-based TSS."""

    delta = int(pos) - int(tss)
    return -delta if _normalise_strand(strand) == "-" else delta


def _distance_to_tss(peak_start: int, peak_end: int, tss: int, strand: str) -> int:
    """ChIPseeker distance using the closer of the two peak ends.

    ``peak_end`` is the half-open end.  Internally we use the inclusive end
    (``peak_end - 1``) because GRanges coordinates are closed.
    """

    d_start = _signed_distance(int(peak_start), int(tss), strand)
    d_end = _signed_distance(int(peak_end) - 1, int(tss), strand)
    return int(d_start if abs(d_start) <= abs(d_end) else d_end)


def _promoter_bins_from_region(tss_region: Sequence[int]) -> tuple[int, ...]:
    pm = max(abs(int(tss_region[0])), abs(int(tss_region[1])))
    if pm / 1000 < 2:
        return ()
    n = int(np.ceil(pm / 1000))
    return tuple(i * 1000 for i in range(1, n + 1))


def promoter_label(distance: float, bins: Sequence[int]) -> str:
    """Return a ChIPseeker promoter subcategory label."""

    if not bins:
        return "Promoter"
    abs_d = abs(float(distance))
    b0, b1, b2 = (list(bins) + [None, None, None])[:3]
    if abs_d <= b0:
        return f"Promoter (<={_format_kb(int(b0))})"
    if b1 is not None and abs_d <= b1:
        return f"Promoter ({int(b0) // 1000}-{int(b1) // 1000}kb)"
    if b2 is not None and abs_d <= b2:
        return f"Promoter ({int(b1) // 1000}-{int(b2) // 1000}kb)"
    # ChIPseeker only labels bins inside the tssRegion.  If a custom bin list
    # is shorter than the window, keep the generic label.
    return "Promoter"


def _promoter_sort_key(label: str) -> tuple[int, float]:
    """Order promoter labels like ChIPseeker's R locale sort.

    R places ``Promoter (<=1kb)`` before ``Promoter (1-2kb)``; a plain Python
    string sort would put the digit forms first.
    """

    import re

    text = str(label)
    if text == "Promoter":
        return (0, 0.0)
    match = re.search(r"<=([0-9.]+)kb", text)
    if match:
        return (1, float(match.group(1)))
    match = re.search(r"([0-9.]+)-([0-9.]+)kb", text)
    if match:
        return (1, float(match.group(2)))
    return (2, 0.0)


def group_annotation(annotation: str) -> str:
    """Collapse detailed labels into major ChIPseeker classes."""

    label = str(annotation)
    if label.startswith("Promoter"):
        return "Promoter"
    if label.startswith("Exon"):
        return "Exon"
    if label.startswith("Intron"):
        return "Intron"
    if label.startswith("Downstream"):
        return "Downstream"
    if label in {"5' UTR", "3' UTR", "UTR"}:
        return label
    if "Intergenic" in label:
        return "Intergenic"
    return label


def _normalise_priority(priority: Sequence[str]) -> tuple[str, ...]:
    result: list[str] = []
    for item in priority:
        key = str(item).strip()
        canonical = _PRIORITY_ALIASES.get(key.lower())
        if canonical is None:
            raise ValueError(
                "genomic_annotation_priority must contain any order of "
                "Promoter, 5UTR, 3UTR, Exon, Intron, Downstream, Intergenic"
            )
        result.append(canonical)
    if sorted(result) != sorted(_DEFAULT_PRIORITY):
        raise ValueError(
            "genomic_annotation_priority must be a permutation of "
            "Promoter, 5UTR, 3UTR, Exon, Intron, Downstream, Intergenic"
        )
    return tuple(result)


def _ensure_columns(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    df = df.copy()
    for col in columns:
        if col not in df.columns:
            df[col] = None
    return df


def _prepare_utr_table(df: pd.DataFrame, feature: str) -> pd.DataFrame:
    sub = df[df["feature"] == feature].copy()
    if sub.empty:
        return pd.DataFrame(columns=_CORE_INTERVAL_COLUMNS)
    return _ensure_columns(sub, _CORE_INTERVAL_COLUMNS)


def _split_generic_utr(
    df: pd.DataFrame,
    utr5: pd.DataFrame,
    utr3: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split ``UTR`` records into 5' and 3' UTRs using CDS bounds."""

    generic = df[df["feature"] == "utr"].copy()
    if generic.empty:
        return utr5, utr3
    cds = df[df["feature"] == "cds"].copy()
    if cds.empty or "transcript_id" not in generic:
        return utr5, utr3
    cds = cds.dropna(subset=["transcript_id", "start", "end"])
    if cds.empty:
        return utr5, utr3

    bounds = (
        cds.groupby("transcript_id", as_index=False)
        .agg(cds_start=("start", "min"), cds_end=("end", "max"))
    )
    generic = generic.merge(bounds, on="transcript_id", how="left")
    plus = generic["strand"].astype(str) != "-"
    valid = generic["cds_start"].notna() & generic["cds_end"].notna()

    five = valid & (
        (plus & (generic["end"] <= generic["cds_start"]))
        | (~plus & (generic["start"] >= generic["cds_end"]))
    )
    three = valid & (
        (plus & (generic["start"] >= generic["cds_end"]))
        | (~plus & (generic["end"] <= generic["cds_start"]))
    )

    if five.any():
        utr5 = pd.concat(
            [utr5, _ensure_columns(generic[five], _CORE_INTERVAL_COLUMNS)],
            ignore_index=True,
        )
    if three.any():
        utr3 = pd.concat(
            [utr3, _ensure_columns(generic[three], _CORE_INTERVAL_COLUMNS)],
            ignore_index=True,
        )
    return utr5, utr3


def _add_tss_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df["strand"] = df["strand"].map(_normalise_strand)
    df["tss"] = np.where(df["strand"] == "-", df["end"] - 1, df["start"]).astype(
        np.int64
    )
    df["length"] = (df["end"] - df["start"]).astype(np.int64)
    return df


def _sort_by_transcript_order(
    df: pd.DataFrame,
    transcript_order: Dict[str, int],
) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["_tx_order"] = (
        df["transcript_id"].astype(str).map(transcript_order).fillna(10**12)
    )
    return (
        df.sort_values(["_tx_order", "start", "end"], kind="stable")
        .drop(columns="_tx_order")
        .reset_index(drop=True)
    )


class GeneAnnotation:
    """Gene/transcript models and interval indexes used for annotation."""

    def __init__(
        self,
        genes: pd.DataFrame,
        transcripts: pd.DataFrame,
        exons: pd.DataFrame,
        introns: pd.DataFrame,
        utr5: pd.DataFrame,
        utr3: pd.DataFrame,
    ) -> None:
        self.genes = genes.reset_index(drop=True).copy()
        self.transcripts = transcripts.reset_index(drop=True).copy()
        self.exons = exons.reset_index(drop=True).copy()
        self.introns = introns.reset_index(drop=True).copy()
        self.utr5 = utr5.reset_index(drop=True).copy()
        self.utr3 = utr3.reset_index(drop=True).copy()

        self._features: Dict[str, pd.DataFrame] = {
            "gene": self.genes,
            "transcript": self.transcripts,
        }
        self._tss_indexes: Dict[str, TSSIndex] = {}
        self._span_indexes: Dict[str, SpanIndex] = {}
        for level, table in self._features.items():
            self._tss_indexes[level] = TSSIndex(table)
            self._span_indexes[level] = SpanIndex(table)

        self.gene_index = BinnedIntervalIndex(1_000_000).build(self.genes)
        self.transcript_index = BinnedIntervalIndex(1_000_000).build(self.transcripts)
        self.exon_index = BinnedIntervalIndex(100_000).build(self.exons)
        self.intron_index = BinnedIntervalIndex(100_000).build(self.introns)
        self.utr5_index = BinnedIntervalIndex(100_000).build(self.utr5)
        self.utr3_index = BinnedIntervalIndex(100_000).build(self.utr3)

        # Backwards-compatible gene-level attributes.
        self.tss_index = self._tss_indexes["gene"]
        self.span_index = self._span_indexes["gene"]

        self._positions: Dict[str, Dict[str, int]] = {}
        for level, table in self._features.items():
            id_col = "gene_id" if level == "gene" else "transcript_id"
            self._positions[level] = {
                str(v): i for i, v in enumerate(table[id_col].astype(str))
            }

        self._category_tables = {
            "5UTR": self.utr5,
            "3UTR": self.utr3,
            "Exon": self.exons,
            "Intron": self.introns,
        }
        self._category_indexes = {
            "5UTR": self.utr5_index,
            "3UTR": self.utr3_index,
            "Exon": self.exon_index,
            "Intron": self.intron_index,
        }

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def from_gtf(cls, gtf: Union[str, Path, pd.DataFrame]) -> "GeneAnnotation":
        if isinstance(gtf, (str, Path)):
            df = read_gtf(gtf)
        elif isinstance(gtf, pd.DataFrame):
            df = gtf.copy()
        else:
            raise TypeError("gtf must be a path or a pandas.DataFrame")

        if "feature" not in df.columns:
            raise ValueError("annotation table must contain a 'feature' column")

        df = df.copy()
        df["feature"] = (
            df["feature"]
            .astype(str)
            .str.strip()
            .str.lower()
            .str.replace("-", "_", regex=False)
        )
        for col in (
            "gene_id",
            "gene_name",
            "gene_type",
            "transcript_id",
            "transcript_name",
            "transcript_type",
            "exon_number",
            "strand",
            "chrom",
        ):
            if col not in df.columns:
                df[col] = None

        df["chrom"] = df["chrom"].astype(str)
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df = df.dropna(subset=["chrom", "start", "end"]).copy()
        df["start"] = df["start"].astype(np.int64)
        df["end"] = df["end"].astype(np.int64)
        df = df[df["end"] > df["start"]].copy()
        df["strand"] = df["strand"].map(_normalise_strand)
        if df.empty:
            raise ValueError("no valid feature records were found in the annotation")

        gene_rows = df[df["feature"] == "gene"].copy()
        tx_rows = df[df["feature"].isin({"transcript", "mrna"})].copy()
        exon_rows = df[df["feature"] == "exon"].copy()

        # A transcript_id is essential for the transcript-level model.  If
        # exons only carry a gene_id, use the gene as a one-transcript model.
        if not exon_rows.empty and exon_rows["transcript_id"].isna().any():
            missing_tx = exon_rows["transcript_id"].isna()
            if exon_rows.loc[missing_tx, "gene_id"].isna().any():
                warnings.warn(
                    "exon records without transcript_id/gene_id will be ignored",
                    RuntimeWarning,
                )
            exon_rows.loc[missing_tx, "transcript_id"] = exon_rows.loc[
                missing_tx, "gene_id"
            ]

        if tx_rows.empty:
            if not exon_rows.empty and exon_rows["transcript_id"].notna().any():
                tx_rows = (
                    exon_rows.dropna(subset=["transcript_id"])
                    .groupby("transcript_id", as_index=False)
                    .agg(
                        chrom=("chrom", "first"),
                        start=("start", "min"),
                        end=("end", "max"),
                        strand=("strand", "first"),
                        gene_id=("gene_id", "first"),
                        gene_name=("gene_name", "first"),
                        transcript_name=("transcript_name", "first"),
                        transcript_type=("transcript_type", "first"),
                        gene_type=("gene_type", "first"),
                    )
                )
            elif not gene_rows.empty:
                tx_rows = gene_rows.copy()
                tx_rows["transcript_id"] = tx_rows["gene_id"]
                tx_rows["transcript_name"] = tx_rows["transcript_id"]
                tx_rows["transcript_type"] = tx_rows.get("gene_type", None)
            else:
                raise ValueError(
                    "no transcript or exon records were found in the annotation"
                )

        # Recover gene IDs for transcript rows that only have a transcript ID.
        tx2gene = (
            tx_rows.dropna(subset=["transcript_id", "gene_id"])
            .drop_duplicates("transcript_id")
            .set_index("transcript_id")["gene_id"]
        )
        if exon_rows["gene_id"].isna().any() and not tx2gene.empty:
            missing_gene = exon_rows["gene_id"].isna()
            exon_rows.loc[missing_gene, "gene_id"] = exon_rows.loc[
                missing_gene, "transcript_id"
            ].map(tx2gene)

        if gene_rows.empty:
            source = tx_rows if not tx_rows.empty else exon_rows
            if source.empty:
                raise ValueError("no gene or transcript records were found")
            gene_rows = (
                source.dropna(subset=["gene_id"])
                .groupby("gene_id", as_index=False)
                .agg(
                    chrom=("chrom", "first"),
                    start=("start", "min"),
                    end=("end", "max"),
                    strand=("strand", "first"),
                    gene_name=("gene_name", "first"),
                    gene_type=("gene_type", "first"),
                )
            )

        gene_rows = gene_rows.dropna(subset=["gene_id"]).copy()
        gene_rows["gene_id"] = gene_rows["gene_id"].astype(str)
        gene_rows["gene_name"] = gene_rows["gene_name"].fillna(gene_rows["gene_id"])
        gene_rows = gene_rows.drop_duplicates("gene_id").reset_index(drop=True)
        gene_rows = _add_tss_columns(gene_rows)
        gene_name_map = (
            gene_rows.set_index("gene_id")["gene_name"].astype(str).to_dict()
        )
        gene_type_map = (
            gene_rows.set_index("gene_id")["gene_type"].to_dict()
            if "gene_type" in gene_rows.columns
            else {}
        )

        def _fill_gene_metadata(table: pd.DataFrame) -> pd.DataFrame:
            """Recover gene symbols/types on child records.

            ``read_gtf`` falls back to ``gene_id`` when a child record has no
            ``gene_name``; replace that placeholder with the gene-level symbol.
            """

            if table.empty:
                return table
            table = table.copy()
            ids = table["gene_id"].astype(str)
            names = table["gene_name"].astype(object)
            bad_name = names.isna() | names.astype(str).isin({"None", "nan", ""})
            bad_name = bad_name | (names.astype(str) == ids)
            mapped_names = ids.map(gene_name_map)
            table["gene_name"] = names.where(~bad_name, mapped_names).fillna(ids)

            if "gene_type" in table.columns and gene_type_map:
                types = table["gene_type"].astype(object)
                bad_type = types.isna() | types.astype(str).isin({"None", "nan", ""})
                mapped_types = ids.map(gene_type_map)
                table["gene_type"] = types.where(~bad_type, mapped_types)
            return table

        tx_rows = tx_rows.dropna(subset=["transcript_id"]).copy()
        tx_rows["transcript_id"] = tx_rows["transcript_id"].astype(str)
        tx_rows["gene_id"] = tx_rows["gene_id"].fillna(
            tx_rows["transcript_id"].map(tx2gene)
        )
        # Some GTF/GFF files place the gene on the exon rows only.
        if tx_rows["gene_id"].isna().any() and not exon_rows.empty:
            exon_gene = (
                exon_rows.dropna(subset=["transcript_id", "gene_id"])
                .drop_duplicates("transcript_id")
                .set_index("transcript_id")["gene_id"]
            )
            missing = tx_rows["gene_id"].isna()
            tx_rows.loc[missing, "gene_id"] = tx_rows.loc[missing, "transcript_id"].map(
                exon_gene
            )
        tx_rows = tx_rows.dropna(subset=["gene_id"]).copy()
        tx_rows["gene_id"] = tx_rows["gene_id"].astype(str)
        tx_rows = _fill_gene_metadata(tx_rows)
        tx_rows["transcript_name"] = tx_rows["transcript_name"].fillna(
            tx_rows["transcript_id"]
        )
        tx_rows = tx_rows.drop_duplicates("transcript_id").reset_index(drop=True)
        tx_rows = _add_tss_columns(tx_rows)
        transcripts = tx_rows[
            [
                "transcript_id",
                "transcript_name",
                "transcript_type",
                "gene_id",
                "gene_name",
                "gene_type",
                "chrom",
                "start",
                "end",
                "strand",
                "tss",
                "length",
            ]
        ].copy()

        genes = gene_rows[
            [
                "gene_id",
                "gene_name",
                "gene_type",
                "chrom",
                "start",
                "end",
                "strand",
                "tss",
                "length",
            ]
        ].copy()

        transcript_order = {
            str(row.transcript_id): i
            for i, row in enumerate(transcripts.itertuples(index=False))
        }

        # ---- exons -----------------------------------------------------
        exons = exon_rows.copy()
        exons = exons.dropna(subset=["transcript_id", "gene_id"]).copy()
        exons["transcript_id"] = exons["transcript_id"].astype(str)
        exons["gene_id"] = exons["gene_id"].astype(str)
        exons = exons.drop_duplicates(
            ["transcript_id", "chrom", "start", "end", "strand"]
        )
        if exons.empty:
            exons = pd.DataFrame(
                columns=[
                    "transcript_id",
                    "gene_id",
                    "gene_name",
                    "transcript_name",
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "exon_number",
                    "exon_rank",
                    "exon_total",
                    "annotation_id",
                ]
            )
        else:
            exons = _fill_gene_metadata(exons)
            exons["transcript_name"] = exons["transcript_name"].fillna(
                exons["transcript_id"]
            )
            exons = _sort_by_transcript_order(exons, transcript_order)
            exons["exon_total"] = exons.groupby("transcript_id")[
                "transcript_id"
            ].transform("size")
            exons["_grank"] = exons.groupby("transcript_id").cumcount() + 1
            calc_rank = np.where(
                exons["strand"] == "-",
                exons["exon_total"] - exons["_grank"] + 1,
                exons["_grank"],
            )
            given = pd.to_numeric(exons["exon_number"], errors="coerce")
            exons["exon_rank"] = (
                given.where(given.notna() & (given > 0), pd.Series(calc_rank, index=exons.index))
                .astype(int)
            )
            exons["annotation_id"] = (
                exons["transcript_id"] + "/" + exons["gene_id"]
            )
            exons = exons.drop(columns="_grank")

        # ---- introns ---------------------------------------------------
        intron_records: list[dict[str, object]] = []
        if not exons.empty:
            for tx_id, sub in exons.groupby("transcript_id", sort=False):
                sub = sub.sort_values("start", kind="stable")
                gaps: list[tuple[int, int]] = []
                starts = sub["start"].to_numpy(dtype=np.int64)
                ends = sub["end"].to_numpy(dtype=np.int64)
                for i in range(len(sub) - 1):
                    gap_start = int(ends[i])
                    gap_end = int(starts[i + 1])
                    if gap_end > gap_start:
                        gaps.append((gap_start, gap_end))
                total = len(gaps)
                if total == 0:
                    continue
                first = sub.iloc[0]
                for i, (gap_start, gap_end) in enumerate(gaps, start=1):
                    rank = total - i + 1 if first["strand"] == "-" else i
                    intron_records.append(
                        {
                            "transcript_id": tx_id,
                            "gene_id": first["gene_id"],
                            "gene_name": first["gene_name"],
                            "transcript_name": first["transcript_name"],
                            "chrom": first["chrom"],
                            "start": gap_start,
                            "end": gap_end,
                            "strand": first["strand"],
                            "intron_rank": rank,
                            "intron_total": total,
                            "annotation_id": first["annotation_id"],
                        }
                    )
        introns = pd.DataFrame(intron_records)
        if introns.empty:
            introns = pd.DataFrame(
                columns=[
                    "transcript_id",
                    "gene_id",
                    "gene_name",
                    "transcript_name",
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "intron_rank",
                    "intron_total",
                    "annotation_id",
                ]
            )
        else:
            introns = _sort_by_transcript_order(introns, transcript_order)

        # ---- UTRs ------------------------------------------------------
        utr5 = _prepare_utr_table(df, "five_prime_utr")
        utr3 = _prepare_utr_table(df, "three_prime_utr")
        utr5, utr3 = _split_generic_utr(df, utr5, utr3)

        def _prepare_utrs(table: pd.DataFrame) -> pd.DataFrame:
            if table.empty:
                return pd.DataFrame(
                    columns=[
                        "transcript_id",
                        "gene_id",
                        "gene_name",
                        "transcript_name",
                        "chrom",
                        "start",
                        "end",
                        "strand",
                        "annotation_id",
                    ]
                )
            table = table.dropna(subset=["transcript_id", "gene_id"]).copy()
            table["transcript_id"] = table["transcript_id"].astype(str)
            table["gene_id"] = table["gene_id"].astype(str)
            table = _fill_gene_metadata(table)
            table["transcript_name"] = table["transcript_name"].fillna(
                table["transcript_id"]
            )
            table = table.drop_duplicates(
                ["transcript_id", "chrom", "start", "end", "strand"]
            )
            table["annotation_id"] = table["transcript_id"] + "/" + table["gene_id"]
            return _sort_by_transcript_order(table, transcript_order)

        utr5 = _prepare_utrs(utr5)
        utr3 = _prepare_utrs(utr3)

        return cls(genes, transcripts, exons, introns, utr5, utr3)

    # ------------------------------------------------------------------
    # feature access
    # ------------------------------------------------------------------
    def features(self, level: str = "transcript") -> pd.DataFrame:
        if level not in self._features:
            raise ValueError("level must be 'transcript' or 'gene'")
        return self._features[level]

    def tss_index_for(self, level: str) -> TSSIndex:
        return self._tss_indexes[level]

    def span_index_for(self, level: str) -> SpanIndex:
        return self._span_indexes[level]

    def feature_by_id(self, level: str, feature_id: str) -> Optional[pd.Series]:
        pos = self._positions[level].get(str(feature_id))
        if pos is None:
            return None
        return self.features(level).iloc[pos]

    def feature_row(self, level: str, row_index: int) -> pd.Series:
        return self.features(level).iloc[int(row_index)]

    def category_table(self, category: str) -> pd.DataFrame:
        return self._category_tables[category]

    def category_index(self, category: str) -> BinnedIntervalIndex:
        return self._category_indexes[category]

    # ------------------------------------------------------------------
    # backwards-compatible gene helpers
    # ------------------------------------------------------------------
    def gene_row(self, gene_id: str) -> Optional[pd.Series]:
        pos = self._positions["gene"].get(str(gene_id))
        if pos is None:
            return None
        return self.genes.iloc[pos]

    def nearest_gene(self, chrom: str, pos: int) -> tuple[Optional[pd.Series], float]:
        result = self.tss_index.nearest(chrom, pos)
        if result is None:
            return None, float("nan")
        gene_idx, raw_delta = result
        row = self.genes.iloc[int(gene_idx)]
        distance = -raw_delta if _normalise_strand(row["strand"]) == "-" else raw_delta
        return row, float(distance)

    def choose_gene(
        self, gene_ids: Iterable[str], pos: int
    ) -> tuple[Optional[pd.Series], float]:
        best_row: Optional[pd.Series] = None
        best_key: Optional[tuple[float, float, str]] = None
        for gene_id in set(map(str, gene_ids)):
            row = self.gene_row(gene_id)
            if row is None:
                continue
            distance = _signed_distance(pos, int(row["tss"]), str(row["strand"]))
            key = (abs(distance), distance, str(gene_id))
            if best_key is None or key < best_key:
                best_key = key
                best_row = row
        if best_row is None:
            return None, float("nan")
        return best_row, float(best_key[1])  # type: ignore[index]

    def promoter_gene(
        self,
        chrom: str,
        pos: int,
        upstream: int,
        downstream: int,
    ) -> tuple[Optional[pd.Series], float]:
        """Legacy helper: nearest gene whose promoter contains ``pos``."""

        max_window = max(int(upstream), int(downstream))
        gene_idx, _ = self.tss_index.within(
            chrom,
            int(pos) - max_window,
            int(pos) + max_window + 1,
        )
        best_row: Optional[pd.Series] = None
        best_key: Optional[tuple[float, float, str]] = None
        best_distance = float("nan")
        for gi in gene_idx:
            row = self.genes.iloc[int(gi)]
            distance = _signed_distance(int(pos), int(row["tss"]), str(row["strand"]))
            if not (-int(upstream) <= distance <= int(downstream)):
                continue
            key = (abs(distance), distance, str(row["gene_id"]))
            if best_key is None or key < best_key:
                best_key = key
                best_row = row
                best_distance = float(distance)
        return best_row, best_distance

    def promoter_label(self, distance: float, bins: Sequence[int]) -> str:
        return promoter_label(distance, bins)


@dataclass
class _CategoryHit:
    annotation: str
    feature_id: str
    gene_id: str
    detail_key: str


@dataclass
class PeakAnnotation:
    """Result of :func:`annotate_peak`."""

    peaks: pd.DataFrame
    annotation: pd.DataFrame
    gene_annotation: GeneAnnotation
    tss_region: tuple[int, int] = (-3000, 3000)
    level: str = "transcript"
    promoter_bins: tuple[int, ...] = ()
    genomic_annotation_priority: tuple[str, ...] = _DEFAULT_PRIORITY
    downstream_distance: int = 300
    detail_genomic_annotation: pd.DataFrame = field(default_factory=pd.DataFrame)
    has_genomic_annotation: bool = True
    add_flank_gene_info: bool = False
    flank_distance: int = 5000

    @property
    def upstream(self) -> int:
        return abs(int(self.tss_region[0]))

    @property
    def downstream(self) -> int:
        return int(self.tss_region[1])

    @property
    def detail(self) -> pd.DataFrame:
        """Alias for the ChIPseeker ``detailGenomicAnnotation`` slot."""

        return self.detail_genomic_annotation

    @property
    def annoStat(self) -> pd.DataFrame:
        """Alias for the ChIPseeker ``annoStat`` slot."""

        return self.anno_stat()

    @property
    def peakNum(self) -> int:
        return len(self.peaks)

    @property
    def hasGenomicAnnotation(self) -> bool:
        return bool(self.has_genomic_annotation)

    def __repr__(self) -> str:
        return (
            f"PeakAnnotation(n_peaks={len(self.annotation)}, level={self.level!r}, "
            f"tssRegion=({self.tss_region[0]}, {self.tss_region[1]}))"
        )

    def to_dataframe(self, include_peak_metadata: bool = False) -> pd.DataFrame:
        df = self.annotation.copy()
        if include_peak_metadata:
            base = {"peak_id", "chrom", "start", "end", "name", "score", "strand"}
            extras = [
                c
                for c in self.peaks.columns
                if c not in base and c not in self.annotation.columns
            ]
            if extras:
                df = df.merge(
                    self.peaks[["peak_id", *extras]],
                    on="peak_id",
                    how="left",
                )
        return df

    def to_chipseeker(self) -> pd.DataFrame:
        """Return a table with ChIPseeker-style column names.

        Coordinates remain 0-based, half-open in the Python table.  ChIPseeker
        users will normally convert ``peakStart``/``geneStart`` with ``+1`` and
        ``peakEnd``/``geneEnd`` as-is when moving between the two conventions.
        """

        mapping = {
            "peak_id": "peakId",
            "chrom": "peakChr",
            "start": "peakStart",
            "end": "peakEnd",
            "width": "peakWidth",
            "center": "peakCenter",
            "strand": "peakStrand",
            "gene_id": "geneId",
            "gene_name": "geneName",
            "gene_chrom": "geneChr",
            "gene_start": "geneStart",
            "gene_end": "geneEnd",
            "gene_strand": "geneStrand",
            "gene_length": "geneLength",
            "gene_type": "geneBiotype",
            "transcript_id": "transcriptId",
            "transcript_name": "transcriptName",
            "transcript_type": "transcriptBiotype",
            "distance_to_tss": "distanceToTSS",
            "nearest_gene_id": "nearestGeneId",
            "nearest_gene_name": "nearestGeneName",
            "nearest_transcript_id": "nearestTranscriptId",
            "nearest_distance": "nearestDistance",
            "flank_gene_ids": "flank_geneIds",
            "flank_gene_distances": "flank_gene_distances",
            "flank_tx_ids": "flank_txIds",
        }
        return self.annotation.rename(columns=mapping).copy()

    def as_data_frame(self) -> pd.DataFrame:
        """ChIPseeker ``as.data.frame(csAnno)`` analogue."""

        return self.to_chipseeker()

    asDataFrame = as_data_frame

    def __len__(self) -> int:
        return len(self.annotation)

    def as_granges(self) -> pd.DataFrame:
        """Return a GRanges-like table (``seqnames/start/end/strand``).

        Python does not have Bioconductor ``GRanges``; this table keeps the
        ChIPseeker metadata columns and uses the conventional GRanges names.
        Coordinates remain 0-based, half-open.
        """

        df = self.to_chipseeker()
        rename = {
            "peakChr": "seqnames",
            "peakStart": "start",
            "peakEnd": "end",
            "peakStrand": "strand",
        }
        out = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
        leading = [c for c in ("seqnames", "start", "end", "strand") if c in out.columns]
        rest = [c for c in out.columns if c not in leading]
        return out[leading + rest]

    # ChIPseeker spelling.
    asGRanges = as_granges

    def get_anno_stat(self) -> pd.DataFrame:
        return self.anno_stat()

    getAnnoStat = get_anno_stat

    def drop_anno(
        self,
        distance_to_tss_cutoff: float = 10000,
    ) -> "PeakAnnotation":
        """Keep peaks within ``abs(distanceToTSS) < cutoff``.

        Mirrors ``ChIPseeker::dropAnno`` (the function filters by distance; it
        does not remove annotation columns).
        """

        cutoff = float(distance_to_tss_cutoff)
        if cutoff < 0:
            raise ValueError("distance_to_tss_cutoff must be non-negative")
        if "distance_to_tss" not in self.annotation.columns:
            raise ValueError("annotation has no distance_to_tss column")
        mask = self.annotation["distance_to_tss"].abs() < cutoff
        annotation = self.annotation.loc[mask].reset_index(drop=True)
        detail = (
            self.detail_genomic_annotation.loc[mask].reset_index(drop=True)
            if len(self.detail_genomic_annotation) == len(mask)
            else self.detail_genomic_annotation
        )
        peaks = self.peaks[
            self.peaks["peak_id"].isin(annotation["peak_id"])
        ].reset_index(drop=True)
        return PeakAnnotation(
            peaks=peaks,
            annotation=annotation,
            gene_annotation=self.gene_annotation,
            tss_region=self.tss_region,
            level=self.level,
            promoter_bins=self.promoter_bins,
            genomic_annotation_priority=self.genomic_annotation_priority,
            downstream_distance=self.downstream_distance,
            detail_genomic_annotation=detail,
            has_genomic_annotation=self.has_genomic_annotation,
            add_flank_gene_info=self.add_flank_gene_info,
            flank_distance=self.flank_distance,
        )

    dropAnno = drop_anno

    def anno_stat(
        self,
        ignore_promoter_subcategory: bool = False,
        ignore_1st_exon: bool = False,
        ignore_1st_intron: bool = False,
        ignore_downstream: bool = False,
    ) -> pd.DataFrame:
        return get_genomic_anno_stat(
            self.annotation["annotation"],
            downstream_distance=self.downstream_distance,
            ignore_promoter_subcategory=ignore_promoter_subcategory,
            ignore_1st_exon=ignore_1st_exon,
            ignore_1st_intron=ignore_1st_intron,
            ignore_downstream=ignore_downstream,
        )

    def summary(
        self,
        groupby: Optional[str] = None,
        normalize: bool = False,
    ) -> pd.DataFrame:
        """Convenience summary.

        With ``groupby=None`` this returns the major annotation classes.  Use
        :meth:`anno_stat` for the detailed ChIPseeker ``annoStat`` table.
        """

        df = self.to_dataframe(include_peak_metadata=True)
        if "annotation" not in df.columns:
            raise ValueError("annotation was not assigned")
        df["annotation_group"] = df["annotation"].map(group_annotation)

        if groupby is None:
            counts = df["annotation_group"].value_counts()
            out = counts.rename_axis("annotation").reset_index(name="count")
            out["percent"] = out["count"] / out["count"].sum() * 100
            return out

        if groupby not in df.columns:
            raise KeyError(f"groupby column not found: {groupby}")
        out = (
            df.groupby([groupby, "annotation_group"], dropna=False)
            .size()
            .unstack(fill_value=0)
        )
        if normalize:
            out = out.div(out.sum(axis=1).replace(0, np.nan), axis=0) * 100
        return out

    def seq2gene(
        self,
        tss_region: Optional[Sequence[int]] = None,
        flank_distance: int = 5000,
        same_strand: bool = False,
    ) -> list[str]:
        return seq2gene(
            self.peaks,
            self.gene_annotation,
            tss_region=tuple(self.tss_region) if tss_region is None else tuple(tss_region),
            flank_distance=flank_distance,
            same_strand=same_strand,
        )

    def seq2gene_links(
        self,
        tss_region: Optional[Sequence[int]] = None,
        flank_distance: int = 5000,
        same_strand: bool = False,
    ) -> pd.DataFrame:
        return seq2gene_links(
            self.peaks,
            self.gene_annotation,
            tss_region=tuple(self.tss_region) if tss_region is None else tuple(tss_region),
            flank_distance=flank_distance,
            same_strand=same_strand,
        )

    # Convenience plotting methods (imported lazily to keep core light).
    def plot_anno_pie(self, **kwargs):
        from .plotting import plot_anno_pie

        return plot_anno_pie(self, **kwargs)

    def plot_anno_bar(self, **kwargs):
        from .plotting import plot_anno_bar

        return plot_anno_bar(self, **kwargs)

    def plot_dist_to_tss(self, **kwargs):
        from .plotting import plot_dist_to_tss

        return plot_dist_to_tss(self, **kwargs)

    plotAnnoPie = plot_anno_pie
    plotAnnoBar = plot_anno_bar
    plotDistToTSS = plot_dist_to_tss

    def vennpie(self, **kwargs):
        from .compat import vennpie

        return vennpie(self, **kwargs)

    def upsetplot(self, **kwargs):
        from .compat import upsetplot

        return upsetplot(self, **kwargs)

    def get_tag_matrix(
        self,
        peak=None,
        upstream: int = 3000,
        downstream: int = 3000,
        type: str = "start_site",
        by: Optional[str] = None,
        **kwargs,
    ):
        """Build a ChIPseeker-style tag matrix using this annotation."""

        from .tagmatrix import get_tag_matrix

        if peak is None:
            peak = self.peaks
        return get_tag_matrix(
            peak=peak,
            upstream=upstream,
            downstream=downstream,
            type=type,
            by=by or self.level,
            annotation=self.gene_annotation,
            **kwargs,
        )

    getTagMatrix = get_tag_matrix

    def plot_peak_prof(self, **kwargs):
        from .tagmatrix import plot_peak_prof

        if "TxDb" not in kwargs and "annotation" not in kwargs and "gtf" not in kwargs:
            kwargs["annotation"] = self.gene_annotation
        return plot_peak_prof(**kwargs)

    plotPeakProf = plot_peak_prof

    def peak_heatmap(self, peak=None, **kwargs):
        from .tagmatrix import peak_heatmap

        if "TxDb" not in kwargs and "annotation" not in kwargs and "gtf" not in kwargs:
            kwargs["annotation"] = self.gene_annotation
        if peak is None:
            peak = self.peaks
        return peak_heatmap(peak, **kwargs)

    peakHeatmap = peak_heatmap

    def _tss_matrix(
        self,
        bigwig: Optional[Union[str, Path]] = None,
        *,
        matrix: Optional[np.ndarray] = None,
        signal: Optional[Union[str, Path, pd.DataFrame]] = None,
        window: int = 2000,
        n_bins: int = 100,
        min_value: Optional[float] = 0.0,
    ) -> np.ndarray:
        """Resolve a precomputed matrix, bedGraph signal, or bigWig."""

        provided = sum(x is not None for x in (matrix, signal, bigwig))
        if provided != 1:
            raise ValueError(
                "provide exactly one of matrix=, signal=, or a bigWig path"
            )

        if matrix is not None:
            arr = np.asarray(matrix, dtype=float)
            if arr.ndim != 2:
                raise ValueError("matrix must be 2-dimensional (features x bins)")
            return arr

        genes = self.gene_annotation.features(self.level)
        if signal is not None:
            from .plotting import extract_tss_matrix_from_signal

            return extract_tss_matrix_from_signal(
                signal,
                genes,
                window=window,
                n_bins=n_bins,
                min_value=min_value,
            )

        from .plotting import extract_tss_matrix

        return extract_tss_matrix(
            bigwig,
            genes,
            window=window,
            n_bins=n_bins,
            min_value=0.0 if min_value is None else float(min_value),
        )

    def plot_tss_profile(
        self,
        bigwig: Optional[Union[str, Path]] = None,
        *,
        matrix: Optional[np.ndarray] = None,
        signal: Optional[Union[str, Path, pd.DataFrame]] = None,
        **kwargs,
    ):
        from .plotting import plot_tss_profile

        window = int(kwargs.pop("window", 2000))
        n_bins = int(kwargs.pop("n_bins", 100))
        min_value = kwargs.pop("min_value", 0.0)
        matrix = self._tss_matrix(
            bigwig,
            matrix=matrix,
            signal=signal,
            window=window,
            n_bins=n_bins,
            min_value=min_value,
        )
        return plot_tss_profile(matrix, window=window, **kwargs)

    def plot_tss_heatmap(
        self,
        bigwig: Optional[Union[str, Path]] = None,
        *,
        matrix: Optional[np.ndarray] = None,
        signal: Optional[Union[str, Path, pd.DataFrame]] = None,
        **kwargs,
    ):
        from .plotting import plot_tss_heatmap

        window = int(kwargs.pop("window", 2000))
        n_bins = int(kwargs.pop("n_bins", 100))
        min_value = kwargs.pop("min_value", 0.0)
        matrix = self._tss_matrix(
            bigwig,
            matrix=matrix,
            signal=signal,
            window=window,
            n_bins=n_bins,
            min_value=min_value,
        )
        return plot_tss_heatmap(matrix, window=window, **kwargs)


# ----------------------------------------------------------------------
# peak table normalisation
# ----------------------------------------------------------------------
def _normalise_peaks(peaks: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    if isinstance(peaks, (str, Path)):
        df = read_peaks(peaks)
    elif isinstance(peaks, pd.DataFrame):
        df = peaks.copy()
    else:
        raise TypeError("peaks must be a path or a pandas.DataFrame")

    required = {"chrom", "start", "end"}
    if not required.issubset(df.columns) and len(df.columns) >= 3:
        first = list(df.columns[:3])
        df = df.rename(columns={first[0]: "chrom", first[1]: "start", first[2]: "end"})

    original_len = len(df)

    rename = {}
    reserved_targets = set(df.columns)
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
        if target in reserved_targets and target != col:
            continue
        rename[col] = target
        reserved_targets.add(target)

    df = df.rename(columns=rename)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"peak table is missing columns: {sorted(missing)}")

    if "name" not in df.columns and "peak_id" in df.columns:
        df["name"] = df["peak_id"]

    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    if df["start"].isna().any() or df["end"].isna().any():
        raise ValueError("peak start/end contain non-numeric values")

    df = df[df["end"] > df["start"]].copy()
    if len(df) == 0 and original_len > 0:
        raise ValueError("no valid peak intervals (end must be greater than start)")
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)

    index_ids = None
    if df.index.name is not None or not isinstance(df.index, pd.RangeIndex):
        index_ids = df.index.astype(str).tolist()

    if "name" in df.columns:
        raw_names = df["name"].tolist()
        ids = []
        for i, value in enumerate(raw_names):
            if pd.isna(value):
                ids.append(index_ids[i] if index_ids is not None else f"peak_{i + 1}")
                continue
            pid = str(value).strip()
            if pid in {"", ".", "nan", "None"}:
                ids.append(index_ids[i] if index_ids is not None else f"peak_{i + 1}")
            else:
                ids.append(pid)
    elif "peak_id" in df.columns and df["peak_id"].notna().any():
        ids = df["peak_id"].astype(str).tolist()
    elif index_ids is not None:
        ids = index_ids
    else:
        ids = [f"peak_{i + 1}" for i in range(len(df))]

    seen: dict[str, int] = {}
    unique_ids: list[str] = []
    for pid in ids:
        if pid in seen:
            seen[pid] += 1
            unique_ids.append(f"{pid}_{seen[pid]}")
        else:
            seen[pid] = 1
            unique_ids.append(pid)
    df["peak_id"] = unique_ids

    if "name" not in df.columns:
        df["name"] = df["peak_id"]
    else:
        invalid_names = df["name"].isna() | df["name"].astype(str).str.strip().isin(
            {"", ".", "nan", "None"}
        )
        if invalid_names.any():
            names = df["name"].astype(object).copy()
            names.loc[invalid_names] = df.loc[invalid_names, "peak_id"]
            df["name"] = names

    if "score" not in df.columns:
        df["score"] = "."
    if "strand" not in df.columns:
        df["strand"] = "*"
    df["strand"] = df["strand"].map(_normalise_strand)

    base = ["peak_id", "chrom", "start", "end", "name", "score", "strand"]
    extras = [c for c in df.columns if c not in base]
    return df[base + extras].reset_index(drop=True)


# ----------------------------------------------------------------------
# nearest-feature / annotation primitives
# ----------------------------------------------------------------------
def _nearest_feature(
    model: GeneAnnotation,
    level: str,
    chrom: str,
    start: int,
    end: int,
    peak_strand: str = "*",
    same_strand: bool = False,
    ignore_overlap: bool = False,
    ignore_upstream: bool = False,
    ignore_downstream: bool = False,
    overlap: str = "TSS",
) -> tuple[Optional[int], float]:
    """Return ``(feature_row_index, distance)`` for one peak.

    ``follow`` is the upstream neighbour and ``precede`` the downstream
    neighbour, matching ``GenomicRanges``.  The nearest of the two is chosen by
    absolute distance to the peak start/end, as in ChIPseeker.
    """

    features = model.features(level)
    if features.empty:
        return None, float("nan")

    tss_index = model.tss_index_for(level)
    span_index = model.span_index_for(level)
    follow_idx, precede_idx = tss_index.follow_precede(
        chrom,
        start,
        end,
        strand=peak_strand,
        same_strand=same_strand,
    )

    def _dist(feature_idx: Optional[int], pos: int) -> float:
        if feature_idx is None:
            return float("inf")
        row = features.iloc[int(feature_idx)]
        return float(_signed_distance(pos, int(row["tss"]), str(row["strand"])))

    ps_d = _dist(follow_idx, start)
    pe_d = _dist(precede_idx, end - 1)

    if ignore_upstream:
        chosen = 2 if precede_idx is not None else None
    elif ignore_downstream:
        chosen = 1 if follow_idx is not None else None
    elif ps_d == float("inf") and pe_d == float("inf"):
        chosen = None
    else:
        chosen = 1 if abs(ps_d) <= abs(pe_d) else 2

    if chosen == 1:
        feature_idx, distance = follow_idx, ps_d
    elif chosen == 2:
        feature_idx, distance = precede_idx, pe_d
    else:
        feature_idx, distance = None, float("nan")

    if feature_idx is None:
        return None, float("nan")

    if not ignore_overlap:
        if overlap == "all":
            hits = span_index.overlapping(
                chrom,
                start,
                end,
                strand=peak_strand,
                same_strand=same_strand,
            )
            if len(hits):
                feature_idx = int(hits[0])
                row = features.iloc[feature_idx]
                distance = float(
                    _distance_to_tss(start, end, int(row["tss"]), str(row["strand"]))
                )
        else:
            hits = tss_index.overlapping(
                chrom,
                start,
                end,
                strand=peak_strand,
                same_strand=same_strand,
            )
            if len(hits):
                feature_idx = int(hits[0])
                distance = 0.0

    return int(feature_idx), float(distance)


def _category_hit(
    model: GeneAnnotation,
    category: str,
    chrom: str,
    start: int,
    end: int,
    peak_strand: str,
    same_strand: bool,
) -> Optional[_CategoryHit]:
    table = model.category_table(category)
    if table.empty:
        return None
    index = model.category_index(category)
    rows = index.query(
        chrom,
        start,
        end,
        strand=peak_strand,
        same_strand=same_strand,
    )
    if len(rows) == 0:
        return None
    row = table.iloc[int(rows[0])]

    if category == "Exon":
        label = (
            f"Exon ({row['annotation_id']}, exon {int(row['exon_rank'])} "
            f"of {int(row['exon_total'])})"
        )
        detail_key = "Exon"
    elif category == "Intron":
        label = (
            f"Intron ({row['annotation_id']}, intron {int(row['intron_rank'])} "
            f"of {int(row['intron_total'])})"
        )
        detail_key = "Intron"
    elif category == "5UTR":
        label = "5' UTR"
        detail_key = "fiveUTR"
    else:
        label = "3' UTR"
        detail_key = "threeUTR"

    return _CategoryHit(
        annotation=label,
        feature_id=str(row["transcript_id"]),
        gene_id=str(row["gene_id"]),
        detail_key=detail_key,
    )


def _downstream_annotation(
    model: GeneAnnotation,
    level: str,
    chrom: str,
    start: int,
    end: int,
    peak_strand: str,
    same_strand: bool,
    downstream_distance: int,
) -> tuple[Optional[str], Optional[int]]:
    """Return a ChIPseeker ``Downstream`` label and the upstream-gene delta."""

    follow_idx, _ = model.span_index_for(level).follow_precede(
        chrom,
        start,
        end,
        strand=peak_strand,
        same_strand=same_strand,
    )
    if follow_idx is None:
        return None, None

    feature = model.features(level).iloc[int(follow_idx)]
    if str(feature["strand"]) == "+":
        dd = int(start) - (int(feature["end"]) - 1)
    else:
        dd = (int(end) - 1) - int(feature["start"])

    dsd = int(downstream_distance)
    if dd == 0 or abs(dd) > dsd:
        return None, int(dd)

    if dsd <= 1000:
        return f"Downstream (<={dsd}bp)", int(dd)

    for i in range(1, int(dsd / 1000) + 1):
        if abs(dd) <= i * 1000:
            if i == 1:
                return "Downstream (<1kb)", int(dd)
            return f"Downstream ({i - 1}-{i}kb)", int(dd)

    return f"Downstream ({dsd / 1000:g}kb-{dsd}bp)", int(dd)


# ----------------------------------------------------------------------
# main annotation API
# ----------------------------------------------------------------------
def annotate_peak(
    peaks: Optional[Union[str, Path, pd.DataFrame]] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame, GeneAnnotation]] = None,
    tss_region: Sequence[int] = (-3000, 3000),
    level: str = "transcript",
    assign_genomic_annotation: bool = True,
    genomic_annotation_priority: Sequence[str] = _DEFAULT_PRIORITY,
    add_flank_gene_info: bool = False,
    flank_distance: int = 5000,
    same_strand: bool = False,
    ignore_overlap: bool = False,
    ignore_upstream: bool = False,
    ignore_downstream: bool = False,
    overlap: str = "TSS",
    verbose: bool = False,
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
    promoter_bins: Optional[Sequence[int]] = None,
    downstream_distance: int = 300,
    drop_unmapped: bool = True,
    align_annotation_feature: bool = True,
    TxDb: Optional[Union[str, Path, pd.DataFrame, GeneAnnotation]] = None,
    **kwargs,
) -> PeakAnnotation:
    """Annotate peaks with ChIPSeeker-compatible genomic features.

    Parameters mirror ``ChIPseeker::annotatePeak``.  ``upstream``/``downstream``
    are kept as convenience aliases for an asymmetric ``tss_region``.
    Coordinates are 0-based, half-open (BED).

    ``annotatePeak`` is a bare alias of this function, so the R-style keyword
    spellings (``tssRegion``, ``genomicAnnotationPriority``, ``peak``,
    ``annotation``, ...) and the ``annoDb``/``columns`` gene-annotation
    shortcut are handled here directly -- there is no separate wrapper.

    ``align_annotation_feature`` follows the ChIPseeker development fix for
    issue #252: when a peak overlaps an exon/intron/UTR from one isoform but the
    nearest TSS belongs to another isoform, report the gene/transcript and
    distance of the isoform that supplied the genomic annotation.  Set it to
    ``False`` to reproduce ChIPseeker 1.48 and earlier, where the annotation
    text and transcript columns can refer to different isoforms.
    """

    # ChIPseeker-compatible keyword spellings.  ``annotatePeak`` is a bare
    # alias of this function, so R-style names and the ``peak=`` /
    # ``annotation=`` spellings are accepted here rather than in a wrapper.
    if peaks is None:
        peaks = kwargs.pop("peak", None)
    if peaks is None:
        raise ValueError("provide peaks=... (a BED-like path or peak DataFrame)")
    annotation_source = kwargs.pop("annotation", None)
    if annotation_source is not None:
        if gtf is not None or TxDb is not None:
            raise ValueError("provide only one of gtf=, TxDb= or annotation=")
        gtf = annotation_source
    tss_region = kwargs.pop("tssRegion", tss_region)
    assign_genomic_annotation = kwargs.pop(
        "assignGenomicAnnotation", assign_genomic_annotation
    )
    genomic_annotation_priority = kwargs.pop(
        "genomicAnnotationPriority", genomic_annotation_priority
    )
    add_flank_gene_info = kwargs.pop("addFlankGeneInfo", add_flank_gene_info)
    flank_distance = kwargs.pop("flankDistance", flank_distance)
    same_strand = kwargs.pop("sameStrand", same_strand)
    ignore_overlap = kwargs.pop("ignoreOverlap", ignore_overlap)
    ignore_upstream = kwargs.pop("ignoreUpstream", ignore_upstream)
    ignore_downstream = kwargs.pop("ignoreDownstream", ignore_downstream)
    anno_db = kwargs.pop("annoDb", None)
    anno_columns = kwargs.pop("columns", None)
    if kwargs:
        unexpected = ", ".join(sorted(kwargs))
        raise TypeError(
            f"annotate_peak() got unexpected keyword argument(s): {unexpected}"
        )
    if genomic_annotation_priority is None:
        genomic_annotation_priority = _DEFAULT_PRIORITY

    if level not in {"transcript", "gene"}:
        raise ValueError("level must be 'transcript' or 'gene'")
    if gtf is None:
        gtf = TxDb
    if gtf is None:
        raise ValueError(
            "provide gtf=... or TxDb=..., or a bundled assembly name such as "
            f"{', '.join(supported_genomes())}"
        )
    from .genomes import is_assembly_name, load_genome, supported_genomes

    if is_assembly_name(gtf):
        gtf = load_genome(gtf)
    if str(overlap) not in {"TSS", "all"}:
        raise ValueError("overlap must be 'TSS' or 'all'")

    if upstream is not None or downstream is not None:
        up = 3000 if upstream is None else int(upstream)
        down = 3000 if downstream is None else int(downstream)
        tss_region = (-up, down)
    tss_region = (int(tss_region[0]), int(tss_region[1]))
    if tss_region[1] < tss_region[0]:
        raise ValueError("tss_region must be (upstream-negative, downstream-positive)")

    if promoter_bins is None:
        bins = _promoter_bins_from_region(tss_region)
    else:
        bins = tuple(int(x) for x in promoter_bins)
        if any(b <= 0 for b in bins) or list(bins) != sorted(bins):
            raise ValueError("promoter_bins must be increasing positive integers")

    priority = _normalise_priority(genomic_annotation_priority)
    downstream_distance = int(downstream_distance)
    flank_distance = int(flank_distance)
    if downstream_distance < 0 or flank_distance < 0:
        raise ValueError("downstream_distance and flank_distance must be non-negative")

    model = gtf if isinstance(gtf, GeneAnnotation) else GeneAnnotation.from_gtf(gtf)
    peaks_df = _normalise_peaks(peaks)
    if peaks_df.empty:
        empty = pd.DataFrame(columns=_annotation_columns(level, assign_genomic_annotation))
        detail = pd.DataFrame(columns=list(_DETAIL_KEYS.values()) + [
            "genic", "Intergenic", "downstream", "distal_intergenic"
        ])
        return _finish_annotation(
            PeakAnnotation(
                peaks=peaks_df,
                annotation=empty,
                gene_annotation=model,
                tss_region=tss_region,
                level=level,
                promoter_bins=bins,
                genomic_annotation_priority=priority,
                downstream_distance=downstream_distance,
                detail_genomic_annotation=detail,
                has_genomic_annotation=assign_genomic_annotation,
                add_flank_gene_info=add_flank_gene_info,
                flank_distance=flank_distance,
            ),
            anno_db,
            anno_columns,
        )

    records: list[dict[str, object]] = []
    details: list[dict[str, bool]] = []
    dropped = 0
    total = len(peaks_df)

    for i, row in enumerate(peaks_df.itertuples(index=False), start=1):
        chrom = str(row.chrom)
        start = int(row.start)
        end = int(row.end)
        peak_strand = _normalise_strand(row.strand)

        nearest_idx, nearest_distance = _nearest_feature(
            model,
            level,
            chrom,
            start,
            end,
            peak_strand=peak_strand,
            same_strand=same_strand,
            ignore_overlap=ignore_overlap,
            ignore_upstream=ignore_upstream,
            ignore_downstream=ignore_downstream,
            overlap=overlap,
        )
        if nearest_idx is None:
            dropped += 1
            continue

        nearest_row = model.features(level).iloc[nearest_idx]
        promoter_hit = (
            assign_genomic_annotation
            and tss_region[0] <= nearest_distance <= tss_region[1]
        )

        hits: dict[str, _CategoryHit] = {}
        detail = {
            "genic": False,
            "Intergenic": False,
            "Promoter": bool(promoter_hit),
            "fiveUTR": False,
            "threeUTR": False,
            "Exon": False,
            "Intron": False,
            "downstream": False,
            "distal_intergenic": False,
        }
        if assign_genomic_annotation:
            for category in ("5UTR", "3UTR", "Exon", "Intron"):
                hit = _category_hit(
                    model,
                    category,
                    chrom,
                    start,
                    end,
                    peak_strand,
                    same_strand,
                )
                if hit is not None:
                    hits[category] = hit
                    detail[hit.detail_key] = True

        annotation: Optional[str] = None
        annotation_feature_id: Optional[str] = None
        annotation_gene_id: Optional[str] = None
        if assign_genomic_annotation:
            for category in reversed(priority):
                if category == "Promoter":
                    if promoter_hit:
                        annotation = promoter_label(nearest_distance, bins)
                        annotation_feature_id = None
                        annotation_gene_id = None
                elif category in {"Downstream", "Intergenic"}:
                    if annotation is None:
                        annotation = "Intergenic"
                else:
                    hit = hits.get(category)
                    if hit is not None:
                        annotation = hit.annotation
                        annotation_feature_id = hit.feature_id
                        annotation_gene_id = hit.gene_id
        else:
            annotation = None

        downstream_dd: Optional[int] = None
        if assign_genomic_annotation and annotation == "Intergenic":
            downstream_label, downstream_dd = _downstream_annotation(
                model,
                level,
                chrom,
                start,
                end,
                peak_strand,
                same_strand,
                downstream_distance,
            )
            if downstream_label is not None:
                annotation = downstream_label
            if downstream_dd is not None:
                detail["downstream"] = bool(0 < downstream_dd < downstream_distance)
            if annotation == "Intergenic":
                annotation = "Distal Intergenic"
                detail["distal_intergenic"] = True
        elif assign_genomic_annotation and annotation is not None:
            if annotation.startswith("Distal"):
                detail["distal_intergenic"] = True

        if assign_genomic_annotation and annotation is not None:
            if annotation == "Intergenic":
                annotation = "Distal Intergenic"
            if str(annotation).startswith("Distal"):
                detail["distal_intergenic"] = True

        # Align the reported feature with the feature that supplied the
        # exon/intron/UTR annotation (ChIPseeker issue #252).
        feature_idx = nearest_idx
        distance = nearest_distance
        alignment_id: Optional[str] = None
        if assign_genomic_annotation and align_annotation_feature:
            if level == "transcript":
                alignment_id = annotation_feature_id
            else:
                alignment_id = annotation_gene_id
        if alignment_id is not None:
            aligned_row = model.feature_by_id(level, alignment_id)
            if aligned_row is not None:
                id_col = "transcript_id" if level == "transcript" else "gene_id"
                pos = model._positions[level].get(str(aligned_row[id_col]))
                if pos is not None:
                    feature_idx = int(pos)
                    distance = float(
                        _distance_to_tss(
                            start,
                            end,
                            int(aligned_row["tss"]),
                            str(aligned_row["strand"]),
                        )
                    )

        feature = model.features(level).iloc[feature_idx]
        nearest_feature = model.features(level).iloc[nearest_idx]
        gene_id = str(feature["gene_id"])
        gene_name = str(feature["gene_name"])
        gene_chrom = str(feature["chrom"])
        gene_start = int(feature["start"])
        gene_end = int(feature["end"])
        gene_strand = str(feature["strand"])
        gene_length = int(feature["length"])
        gene_type = None if pd.isna(feature.get("gene_type")) else str(feature.get("gene_type"))
        transcript_id = None
        transcript_name = None
        transcript_type = None
        if level == "transcript":
            transcript_id = str(feature["transcript_id"])
            transcript_name = str(feature["transcript_name"])
            transcript_type = (
                None
                if pd.isna(feature.get("transcript_type"))
                else str(feature.get("transcript_type"))
            )

        detail["genic"] = bool(detail["Exon"] or detail["Intron"])
        detail["Intergenic"] = not detail["genic"]

        records.append(
            {
                "peak_id": str(row.peak_id),
                "chrom": chrom,
                "start": start,
                "end": end,
                "width": end - start,
                "center": (start + end) // 2,
                "strand": peak_strand,
                "annotation": annotation,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_chrom": gene_chrom,
                "gene_start": gene_start,
                "gene_end": gene_end,
                "gene_strand": gene_strand,
                "gene_length": gene_length,
                "gene_type": gene_type,
                "transcript_id": transcript_id,
                "transcript_name": transcript_name,
                "transcript_type": transcript_type,
                "distance_to_tss": distance,
                "nearest_gene_id": str(nearest_feature["gene_id"]),
                "nearest_gene_name": str(nearest_feature["gene_name"]),
                "nearest_transcript_id": (
                    str(nearest_feature["transcript_id"])
                    if level == "transcript"
                    else None
                ),
                "nearest_distance": nearest_distance,
            }
        )
        details.append(detail)

        if verbose and total > 0 and i % 10000 == 0:
            print(f"annotated {i:,}/{total:,} peaks")

    if dropped:
        warnings.warn(
            f"{dropped} of {total} peaks were dropped because no feature was "
            "found on their chromosome; check chromosome naming/seqlevels",
            RuntimeWarning,
        )
    if not records:
        raise ValueError(
            f"all {total} peaks were dropped: none of the peak chromosomes "
            "matches the annotation seqlevels"
        )

    annotation_df = pd.DataFrame.from_records(
        records, columns=_annotation_columns(level, assign_genomic_annotation)
    )
    detail_df = pd.DataFrame.from_records(details)

    if add_flank_gene_info:
        flank = get_flanking_genes(
            model,
            peaks_df[peaks_df["peak_id"].isin(annotation_df["peak_id"])],
            level=level,
            distance=flank_distance,
        )
        for col in ("flank_gene_ids", "flank_gene_distances"):
            if col not in annotation_df.columns:
                annotation_df[col] = None
        if level == "transcript" and "flank_tx_ids" not in annotation_df.columns:
            annotation_df["flank_tx_ids"] = None
        mapping = flank.set_index("peak_id")
        for col in ("flank_gene_ids", "flank_gene_distances", "flank_tx_ids"):
            if col in mapping.columns:
                annotation_df[col] = annotation_df["peak_id"].map(mapping[col])

    return _finish_annotation(
        PeakAnnotation(
            peaks=peaks_df,
            annotation=annotation_df.reset_index(drop=True),
            gene_annotation=model,
            tss_region=tss_region,
            level=level,
            promoter_bins=bins,
            genomic_annotation_priority=priority,
            downstream_distance=downstream_distance,
            detail_genomic_annotation=detail_df.reset_index(drop=True),
            has_genomic_annotation=assign_genomic_annotation,
            add_flank_gene_info=add_flank_gene_info,
            flank_distance=flank_distance,
        ),
        anno_db,
        anno_columns,
    )


def _finish_annotation(
    annotation: PeakAnnotation,
    anno_db,
    columns: Optional[Sequence[str]] = None,
) -> PeakAnnotation:
    """Apply the ChIPseeker ``annoDb``/``columns`` shortcut, if requested.

    Imported lazily because :mod:`pychipseeker.compat` imports this module;
    ``annotate_db`` lives there to keep the OrgDb-like helpers together.
    """

    if anno_db is None:
        return annotation
    from .compat import annotate_db

    return annotate_db(annotation, anno_db, columns=columns)


def _annotation_columns(level: str, assign: bool) -> list[str]:
    cols = [
        "peak_id",
        "chrom",
        "start",
        "end",
        "width",
        "center",
        "strand",
    ]
    if assign:
        cols.append("annotation")
    cols += [
        "gene_id",
        "gene_name",
        "gene_chrom",
        "gene_start",
        "gene_end",
        "gene_strand",
        "gene_length",
        "gene_type",
    ]
    if level == "transcript":
        cols += ["transcript_id", "transcript_name", "transcript_type"]
    cols += [
        "distance_to_tss",
        "nearest_gene_id",
        "nearest_gene_name",
        "nearest_transcript_id",
        "nearest_distance",
    ]
    return cols


# ----------------------------------------------------------------------
# annotation statistics
# ----------------------------------------------------------------------
def get_genomic_anno_stat(
    annotations: Sequence[str],
    downstream_distance: int = 300,
    ignore_promoter_subcategory: bool = False,
    ignore_1st_exon: bool = False,
    ignore_1st_intron: bool = False,
    ignore_downstream: bool = False,
) -> pd.DataFrame:
    """Replicate ``ChIPseeker:::getGenomicAnnoStat``.

    Returns a ``Feature``/``Frequency`` table where ``Frequency`` is a
    percentage, as in ChIPseeker's ``annoStat`` slot.
    """

    anno = pd.Series(annotations).astype(object).copy()
    if ignore_1st_exon:
        first_exon = None
        other_exon = "Exon"
    else:
        first_exon = "1st Exon"
        other_exon = "Other Exon"
    if ignore_1st_intron:
        first_intron = None
        other_intron = "Intron"
    else:
        first_intron = "1st Intron"
        other_intron = "Other Intron"

    if not ignore_1st_exon:
        mask = anno.astype(str).str.contains(r"exon 1 of", regex=True, na=False)
        anno[mask] = first_exon
    anno[anno.astype(str).str.match(r"Exon \(", na=False)] = other_exon

    if not ignore_1st_intron:
        mask = anno.astype(str).str.contains(r"intron 1 of", regex=True, na=False)
        anno[mask] = first_intron
    anno[anno.astype(str).str.match(r"Intron \(", na=False)] = other_intron

    if ignore_downstream:
        downstream_label = None
        intergenic_label = "Intergenic"
    else:
        dsd = int(downstream_distance)
        if dsd > 1000:
            dsd_label = f"{round(dsd / 1000, 1):g}kb"
        else:
            dsd_label = str(dsd)
        downstream_label = f"Downstream (<={dsd_label})"
        intergenic_label = "Distal Intergenic"

    if downstream_label is not None:
        anno[anno.astype(str).str.startswith("Downstream", na=False)] = downstream_label
    else:
        anno[anno.astype(str).str.startswith("Downstream", na=False)] = intergenic_label
    anno[anno.astype(str).str.startswith("Distal", na=False)] = intergenic_label

    if ignore_promoter_subcategory:
        promoter_labels = ["Promoter"] if (anno.astype(str).str.startswith("Promoter")).any() else []
        anno[anno.astype(str).str.startswith("Promoter", na=False)] = "Promoter"
    else:
        promoter_labels = sorted(
            {
                str(v)
                for v in anno
                if str(v).startswith("Promoter")
            },
            key=_promoter_sort_key,
        )

    levels = [
        *promoter_labels,
        "5' UTR",
        "3' UTR",
        first_exon,
        other_exon,
        first_intron,
        other_intron,
        downstream_label,
        intergenic_label,
    ]
    levels = [v for v in levels if v is not None]
    counts = anno.value_counts()
    counts = counts.reindex([v for v in levels if v in counts.index])
    extra = [v for v in counts.index if v not in levels]
    if extra:
        counts = pd.concat([counts, pd.Series(0, index=extra)]).loc[counts.index]
    total = int(counts.sum())
    frequency = (counts / total * 100) if total else counts.astype(float)
    out = pd.DataFrame(
        {
            "Feature": counts.index.astype(str),
            "Frequency": frequency.to_numpy(dtype=float),
            "Count": counts.to_numpy(dtype=int),
        }
    )
    return out


# ----------------------------------------------------------------------
# many-to-many peak -> gene mapping
# ----------------------------------------------------------------------
def seq2gene(
    peaks: Union[str, Path, pd.DataFrame],
    gene_annotation: Optional[GeneAnnotation] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    tss_region: Sequence[int] = (-3000, 3000),
    flank_distance: int = 5000,
    same_strand: bool = False,
    level: str = "gene",
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
) -> list[str]:
    """Return unique gene IDs, following ``ChIPseeker::seq2gene``.

    Host genes (exon/intron overlap) and genes within ``flank_distance`` of the
    nearest TSS are combined.  Use :func:`seq2gene_links` for a peak-gene table.
    """

    links = seq2gene_links(
        peaks,
        gene_annotation=gene_annotation,
        gtf=gtf,
        tss_region=tss_region,
        flank_distance=flank_distance,
        same_strand=same_strand,
        level=level,
        upstream=upstream,
        downstream=downstream,
    )
    if links.empty:
        return []
    return links["gene_id"].drop_duplicates().astype(str).tolist()


def seq2gene_links(
    peaks: Union[str, Path, pd.DataFrame],
    gene_annotation: Optional[GeneAnnotation] = None,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    tss_region: Sequence[int] = (-3000, 3000),
    flank_distance: int = 5000,
    same_strand: bool = False,
    level: str = "gene",
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
) -> pd.DataFrame:
    """Many-to-many peak-to-gene links used for enrichment analysis."""

    if gene_annotation is None:
        if gtf is None:
            raise ValueError("provide gene_annotation or gtf")
        gene_annotation = GeneAnnotation.from_gtf(gtf)
    if level not in {"gene", "transcript"}:
        raise ValueError("level must be 'gene' or 'transcript'")

    if upstream is not None or downstream is not None:
        up = 3000 if upstream is None else int(upstream)
        down = 3000 if downstream is None else int(downstream)
        tss_region = (-up, down)
    tss_region = (int(tss_region[0]), int(tss_region[1]))
    flank_distance = int(flank_distance)

    peaks_df = _normalise_peaks(peaks)
    links: dict[tuple[str, str], dict[str, object]] = {}

    def add_link(
        peak_row,
        gene_row: pd.Series,
        link_type: str,
        distance: float,
    ) -> None:
        key = (str(peak_row.peak_id), str(gene_row["gene_id"]))
        record = {
            "peak_id": str(peak_row.peak_id),
            "chrom": str(peak_row.chrom),
            "start": int(peak_row.start),
            "end": int(peak_row.end),
            "gene_id": str(gene_row["gene_id"]),
            "gene_name": str(gene_row["gene_name"]),
            "distance_to_tss": float(distance),
            "link_type": link_type,
        }
        priority = {"host": 3, "promoter": 2, "flanking": 1}
        existing = links.get(key)
        if existing is None or priority[link_type] > priority[str(existing["link_type"])]:
            links[key] = record

    genes = gene_annotation.genes
    for peak_row in peaks_df.itertuples(index=False):
        chrom = str(peak_row.chrom)
        start = int(peak_row.start)
        end = int(peak_row.end)
        peak_strand = _normalise_strand(peak_row.strand)

        host_genes: set[str] = set()
        for category in ("Exon", "Intron"):
            rows = gene_annotation.category_index(category).query(
                chrom,
                start,
                end,
                strand=peak_strand,
                same_strand=same_strand,
            )
            if len(rows):
                table = gene_annotation.category_table(category)
                host_genes.add(str(table.iloc[int(rows[0])]["gene_id"]))

        nearest_idx, distance = _nearest_feature(
            gene_annotation,
            "gene",
            chrom,
            start,
            end,
            peak_strand=peak_strand,
            same_strand=same_strand,
            overlap="TSS",
        )
        nearest_gene_id = None
        if nearest_idx is not None:
            nearest_row = genes.iloc[int(nearest_idx)]
            nearest_gene_id = str(nearest_row["gene_id"])
            within_flank = abs(distance) < flank_distance
            within_tss = tss_region[0] < distance < tss_region[1]
            if within_flank or within_tss:
                link_type = "promoter" if within_tss else "flanking"
                add_link(peak_row, nearest_row, link_type, distance)

        for gene_id in host_genes:
            gene_row = gene_annotation.gene_row(gene_id)
            if gene_row is None:
                continue
            if nearest_gene_id == gene_id:
                host_distance = float(distance)
            else:
                host_distance = float(
                    _distance_to_tss(
                        start,
                        end,
                        int(gene_row["tss"]),
                        str(gene_row["strand"]),
                    )
                )
            add_link(peak_row, gene_row, "host", host_distance)

    columns = [
        "peak_id",
        "chrom",
        "start",
        "end",
        "gene_id",
        "gene_name",
        "distance_to_tss",
        "link_type",
    ]
    if not links:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(list(links.values())).reindex(columns=columns)


# ----------------------------------------------------------------------
# flanking genes
# ----------------------------------------------------------------------
def get_flanking_genes(
    model: GeneAnnotation,
    peaks: Union[str, Path, pd.DataFrame],
    level: str = "transcript",
    distance: int = 5000,
) -> pd.DataFrame:
    """Return ChIPseeker-style semicolon-delimited flanking feature info."""

    if level not in {"transcript", "gene"}:
        raise ValueError("level must be 'transcript' or 'gene'")
    peaks_df = _normalise_peaks(peaks)
    features = model.features(level)
    records: list[dict[str, object]] = []
    distance = int(distance)

    for peak_row in peaks_df.itertuples(index=False):
        chrom = str(peak_row.chrom)
        start = int(peak_row.start)
        end = int(peak_row.end)
        ext_start = max(0, start - distance)
        ext_end = end + distance
        rows = model.span_index_for(level).overlapping(chrom, ext_start, ext_end)
        if len(rows) == 0:
            continue

        gene_ids: list[str] = []
        tx_ids: list[str] = []
        distances: list[str] = []
        for row_idx in rows:
            feature = features.iloc[int(row_idx)]
            f_start = int(feature["start"])
            f_end_incl = int(feature["end"]) - 1
            peak_end_incl = end - 1
            if f_start < end and f_end_incl >= start:
                d = 0
            elif str(feature["strand"]) == "+":
                ps_d = start - f_start
                pe_d = peak_end_incl - f_start
                d = pe_d if abs(ps_d) > abs(pe_d) else ps_d
            else:
                ps_d = f_end_incl - peak_end_incl
                pe_d = f_end_incl - start
                d = pe_d if abs(ps_d) > abs(pe_d) else ps_d

            gene_ids.append(str(feature["gene_id"]))
            if level == "transcript":
                tx_ids.append(str(feature["transcript_id"]))
            distances.append(str(int(d)))

        record = {
            "peak_id": str(peak_row.peak_id),
            "flank_gene_ids": ";".join(gene_ids),
            "flank_gene_distances": ";".join(distances),
        }
        if level == "transcript":
            record["flank_tx_ids"] = ";".join(tx_ids)
        records.append(record)

    columns = ["peak_id", "flank_gene_ids", "flank_gene_distances"]
    if level == "transcript":
        columns.append("flank_tx_ids")
    if not records:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame.from_records(records, columns=columns)
