#!/usr/bin/env python
"""Run pychipseeker on a GTF/BED pair and write comparison-friendly tables.

Usage:
    python run_pychipseeker.py <gtf> <bed> <out_prefix> [level] [tssRegion] [overlap]
"""

from __future__ import annotations

import sys
from pathlib import Path

import pychipseeker as pcs


def main(argv: list[str]) -> int:
    if len(argv) < 4:
        raise SystemExit(
            "usage: run_pychipseeker.py <gtf> <bed> <out_prefix> "
            "[level] [tssRegion] [overlap]"
        )
    gtf, bed, out_prefix = argv[1], argv[2], argv[3]
    level = argv[4] if len(argv) >= 5 else "transcript"
    tss_text = argv[5] if len(argv) >= 6 else "-3000,3000"
    overlap = argv[6] if len(argv) >= 7 else "TSS"
    align_text = argv[7] if len(argv) >= 8 else "true"
    align_annotation_feature = str(align_text).lower() not in {
        "0",
        "false",
        "no",
        "off",
    }

    tss_region = tuple(int(x) for x in tss_text.split(","))
    if len(tss_region) != 2:
        raise SystemExit("tssRegion must be two comma-separated integers")

    ann = pcs.annotate_peak(
        peaks=bed,
        gtf=gtf,
        tss_region=tss_region,
        level=level,
        overlap=overlap,
        align_annotation_feature=align_annotation_feature,
        verbose=False,
    )

    out = Path(out_prefix)
    df = ann.to_dataframe(include_peak_metadata=True)
    df.to_csv(str(out) + ".tsv", sep="\t", index=False)

    ann.anno_stat().to_csv(
        str(out) + ".annoStat.tsv", sep="\t", index=False
    )

    detail = ann.detail_genomic_annotation.copy()
    detail.insert(0, "peakId", df["peak_id"].astype(str).to_numpy())
    detail.to_csv(str(out) + ".detail.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
