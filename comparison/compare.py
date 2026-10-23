#!/usr/bin/env python
"""Compare pychipseeker output with ChIPseeker output on the bundled fixture.

The R side is produced by ``run_chipseeker.R``; the Python side by
``run_pychipseeker.py``.  ``compare.py`` can run both automatically.

Usage:
    python compare.py                 # run R + Python for all configs
    python compare.py --no-r          # compare existing output files
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

HERE = Path(__file__).resolve().parent

CONFIGS = {
    "default": dict(level="transcript", tss="-3000,3000", overlap="TSS"),
    "small": dict(level="transcript", tss="-100,100", overlap="TSS"),
    "all": dict(level="transcript", tss="-3000,3000", overlap="all"),
    "gene": dict(level="gene", tss="-3000,3000", overlap="TSS"),
    # ChIPseeker <=1.48 does not align the reported transcript to the isoform
    # that supplied an exon/intron hit.  This mode reproduces that behavior.
    "release": dict(
        level="transcript", tss="-3000,3000", overlap="TSS", align=False
    ),
}


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype={"peakId": str})


def _decode_r_gene_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Undo ChIPseeker 1.48.0's factor-to-integer conversion in as.data.frame.

    ``geneChr`` and ``geneStrand`` are factors in the GRanges metadata; the
    release version writes integer codes.  We can recover chromosome names from
    ``seqnames`` for this fixture (all nearest genes are on the peak chromosome)
    and strand codes from the known levels (1 = '+', 2 = '-').
    """

    df = df.copy()
    if "geneChr" in df.columns:
        codes = pd.to_numeric(df["geneChr"], errors="coerce")
        if codes.notna().all() and codes.isin([1, 2]).all():
            # This fixture has chr1/chr2 only.  Prefer seqnames when the gene
            # lies on the same chromosome, which is true for every test peak.
            df["geneChr"] = df["seqnames"].astype(str)
    if "geneStrand" in df.columns:
        codes = pd.to_numeric(df["geneStrand"], errors="coerce")
        if codes.notna().all() and codes.isin([1, 2]).all():
            df["geneStrand"] = codes.map({1: "+", 2: "-"})
    return df


def _compare_values(
    py: pd.DataFrame,
    r: pd.DataFrame,
    py_col: str,
    r_col: str,
    label: str,
    failures: list[str],
    ignore_peak_ids: set[str] | None = None,
) -> None:
    if ignore_peak_ids:
        py = py[~py["peak_id"].isin(ignore_peak_ids)]
        r = r[~r["peakId"].isin(ignore_peak_ids)]
    left = py[["peak_id", py_col]].rename(columns={py_col: "py"})
    right = r[["peakId", r_col]].rename(columns={"peakId": "peak_id", r_col: "r"})
    merged = left.merge(right, on="peak_id", how="outer", indicator=True)
    mismatched_rows = merged[merged["_merge"] != "both"]
    if not mismatched_rows.empty:
        failures.append(
            f"{label}: peak ID mismatch ({len(mismatched_rows)} row(s)): "
            f"{mismatched_rows['peak_id'].tolist()}"
        )
        return
    if pd.api.types.is_numeric_dtype(merged["py"]) and pd.api.types.is_numeric_dtype(
        merged["r"]
    ):
        bad = ~(
            (merged["py"].isna() & merged["r"].isna())
            | (
                (merged["py"].astype(float) - merged["r"].astype(float)).abs()
                < 1e-9
            )
        )
    else:
        bad = merged["py"].astype(str) != merged["r"].astype(str)
    if bad.any():
        examples = merged.loc[
            bad, ["peak_id", "py", "r"]
        ].head(10).to_dict("records")
        failures.append(
            f"{label}: {int(bad.sum())} mismatch(es); examples={examples}"
        )


def compare_one(
    py_prefix: Path, r_prefix: Path, *, align_feature: bool = True
) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    notes: list[str] = []
    py = _read_table(py_prefix.with_suffix(".tsv"))
    r = _read_table(r_prefix.with_suffix(".tsv"))
    r = _decode_r_gene_metadata(r)

    # ChIPseeker <=1.48 reports the nearest-TSS isoform even when the genomic
    # annotation came from a different isoform.  Our default follows the
    # development fix (issue #252); validate that intended difference
    # explicitly, then compare all other peaks and fields exactly.
    known: set[str] = set()
    if align_feature and "transcript_id" in py.columns:
        known = {"p_shared_exon"}
        py_row = py.set_index("peak_id").loc["p_shared_exon"]
        r_row = r.set_index("peakId").loc["p_shared_exon"]
        checks = {
            "annotation": str(py_row["annotation"]) == str(r_row["annotation"]),
            "py transcript is annotation isoform": py_row["transcript_id"] == "T4",
            "R transcript is nearest isoform": r_row["transcriptId"] == "T5",
            "py distance is aligned": float(py_row["distance_to_tss"])
            == float(r_row["distanceToTSS"]) + 1000.0,
            "py geneStart is aligned": int(py_row["gene_start"])
            == int(r_row["geneStart0"]) - 1000,
            "py geneLength is aligned": int(py_row["gene_length"])
            == int(r_row["geneLength"]) + 1000,
            "nearest transcript retained": py_row["nearest_transcript_id"] == "T5",
            "nearest distance retained": float(py_row["nearest_distance"])
            == float(r_row["distanceToTSS"]),
            "gene id retained": py_row["gene_id"] == r_row["geneId"],
        }
        bad = [name for name, ok in checks.items() if not ok]
        if bad:
            failures.append(f"known isoform-alignment check failed: {bad}")
        else:
            notes.append(
                "p_shared_exon: intentional devel fix (#252); ChIPseeker 1.48 "
                "would report T5/distance 3600, we report aligned T4/distance 4600"
            )

    _compare_values(
        py, r, "annotation", "annotation", "annotation", failures, known
    )
    _compare_values(
        py,
        r,
        "distance_to_tss",
        "distanceToTSS",
        "distanceToTSS",
        failures,
        known,
    )
    _compare_values(py, r, "gene_id", "geneId", "geneId", failures, known)
    _compare_values(
        py, r, "gene_start", "geneStart0", "geneStart0", failures, known
    )
    _compare_values(py, r, "gene_end", "geneEnd", "geneEnd", failures, known)
    _compare_values(
        py, r, "gene_length", "geneLength", "geneLength", failures, known
    )
    if "transcript_id" in py.columns and "transcriptId" in r.columns:
        _compare_values(
            py,
            r,
            "transcript_id",
            "transcriptId",
            "transcriptId",
            failures,
            known,
        )

    py_stat = _read_table(py_prefix.with_suffix(".annoStat.tsv"))
    r_stat = _read_table(r_prefix.with_suffix(".annoStat.tsv"))
    py_stat = py_stat.set_index("Feature")
    r_stat = r_stat.set_index("Feature")
    if list(py_stat.index) != list(r_stat.index):
        failures.append(
            "annoStat feature order differs:\n"
            f"  Python: {list(py_stat.index)}\n"
            f"  R:      {list(r_stat.index)}"
        )
    common = py_stat.index.intersection(r_stat.index)
    delta = (
        py_stat.loc[common, "Frequency"].astype(float)
        - r_stat.loc[common, "Frequency"].astype(float)
    ).abs()
    if (delta > 1e-9).any():
        examples = delta[delta > 1e-9].to_dict()
        failures.append(f"annoStat frequency differs: {examples}")

    py_detail = _read_table(py_prefix.with_suffix(".detail.tsv"))
    r_detail = _read_table(r_prefix.with_suffix(".detail.tsv"))
    bool_cols = [
        "genic",
        "Intergenic",
        "Promoter",
        "fiveUTR",
        "threeUTR",
        "Exon",
        "Intron",
        "downstream",
        "distal_intergenic",
    ]
    py_detail = py_detail.set_index("peakId")
    r_detail = r_detail.set_index("peakId")
    if set(py_detail.index) != set(r_detail.index):
        failures.append("detail peak IDs differ")
    else:
        r_detail = r_detail.loc[py_detail.index]
        for col in bool_cols:
            if col not in py_detail.columns or col not in r_detail.columns:
                failures.append(f"detail column missing: {col}")
                continue
            py_bool = py_detail[col].astype(str).str.lower().isin({"true", "1"})
            r_bool = r_detail[col].astype(str).str.lower().isin({"true", "1"})
            bad = py_bool != r_bool
            if bad.any():
                failures.append(
                    f"detail {col}: {int(bad.sum())} mismatch(es); "
                    f"peaks={py_detail.index[bad].tolist()}"
                )
    return failures, notes


def _find_rscript() -> Path:
    import os

    value = os.environ.get("RSCRIPT")
    if value:
        return Path(value)
    found = shutil.which("Rscript")
    if found:
        return Path(found)

    # Common Windows install locations.  Do not assume a specific R version.
    roots = [
        Path(r"C:\Program Files\R"),
        Path.home() / "AppData" / "Local" / "R",
        Path.home() / "scoop" / "apps" / "r",
    ]
    candidates: list[Path] = []
    for root in roots:
        if root.is_dir():
            candidates.extend(root.glob("R-*/bin/Rscript.exe"))
            candidates.extend(root.glob("*/bin/Rscript.exe"))
    if candidates:
        return sorted(set(candidates), key=lambda p: str(p), reverse=True)[0]

    raise FileNotFoundError(
        "Rscript not found; set RSCRIPT to your Rscript.exe, or use --no-r"
    )


def _r_has_comparison_packages(rscript: Path) -> bool:
    """Return True only when the optional R comparison dependencies exist."""

    code = (
        "quit(status=!(requireNamespace('ChIPseeker', quietly=TRUE) && "
        "requireNamespace('txdbmaker', quietly=TRUE)))"
    )
    try:
        result = subprocess.run(
            [str(rscript), "-e", code],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return False
    return result.returncode == 0


def run_pipeline(configs: Iterable[str], *, no_r: bool = False) -> int:
    rscript = None if no_r else _find_rscript()
    if rscript is not None and not _r_has_comparison_packages(rscript):
        print(
            f"R comparison skipped: {rscript} does not have both "
            "ChIPseeker and txdbmaker installed.\n"
            "The Python implementation and test suite do not need R. "
            "Install the R packages yourself if you want to run this comparison."
        )
        return 2
    failed = False
    for name in configs:
        cfg = CONFIGS[name]
        py_prefix = HERE / f"py_{name}"
        r_prefix = HERE / f"r_{name}"
        if not no_r:
            subprocess.run(
                [
                    str(rscript),
                    str(HERE / "run_chipseeker.R"),
                    str(HERE / "example.gtf"),
                    str(HERE / "example_peaks.bed"),
                    str(r_prefix),
                    cfg["level"],
                    cfg["tss"],
                    cfg["overlap"],
                ],
                check=True,
            )
        subprocess.run(
            [
                sys.executable,
                str(HERE / "run_pychipseeker.py"),
                str(HERE / "example.gtf"),
                str(HERE / "example_peaks.bed"),
                str(py_prefix),
                cfg["level"],
                cfg["tss"],
                cfg["overlap"],
                "true" if cfg.get("align", True) else "false",
            ],
            check=True,
            env={**__import__("os").environ, "PYTHONPATH": str(HERE.parent / "src")},
        )

        failures, notes = compare_one(
            py_prefix, r_prefix, align_feature=cfg.get("align", True)
        )
        if failures:
            failed = True
            print(f"[FAIL] {name}")
            for failure in failures:
                print("  -", failure)
        else:
            print(f"[OK]   {name}: annotation, distance, IDs, annoStat and detail match")
            for note in notes:
                print(f"       note: {note}")
    return 1 if failed else 0


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        choices=sorted(CONFIGS),
        action="append",
        dest="configs",
        help="configuration(s) to run; default: all",
    )
    parser.add_argument(
        "--no-r",
        action="store_true",
        help="only compare existing r_*/py_* output files",
    )
    args = parser.parse_args(argv)
    configs = args.configs or list(CONFIGS)
    try:
        return run_pipeline(configs, no_r=args.no_r)
    except FileNotFoundError as exc:
        print(f"R comparison skipped: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
