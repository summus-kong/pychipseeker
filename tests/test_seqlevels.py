"""Chromosome alias harmonisation and assembly detection."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pychipseeker import annotate_peak
from pychipseeker.genomes import detect_assembly
from pychipseeker.seqlevels import (
    alias_map,
    canonical_chrom,
    chrom_aliases,
    chromosome_style,
    harmonize_frame,
)


def test_canonical_chrom_and_aliases() -> None:
    assert canonical_chrom("chr1") == "1"
    assert canonical_chrom("ChrX") == "x"
    assert canonical_chrom("MT") == "m"
    assert canonical_chrom("chrM") == "m"
    aliases = chrom_aliases("chr1")
    assert aliases[0] == "chr1"
    assert "1" in aliases
    assert "chrM" in chrom_aliases("MT")


def test_alias_map_and_harmonize_frame() -> None:
    mapping = alias_map(["chr1", "chrM", "scaffold_1"], ["1", "MT"])
    assert mapping == {"chr1": "1", "chrM": "MT"}

    frame = pd.DataFrame(
        {"chrom": ["chr1", "chrM", "scaffold_1"], "start": [0, 0, 0], "end": [1, 1, 1]}
    )
    out, applied = harmonize_frame(frame, ["1", "MT"])
    assert out["chrom"].tolist() == ["1", "MT", "scaffold_1"]
    assert applied["chr1"] == "1"
    assert chromosome_style(["chr1", "chr2"]) == "ucsc"
    assert chromosome_style(["1", "2"]) == "ensembl"
    assert chromosome_style(["chr1", "2"]) == "mixed"


def test_detect_assembly_from_sizes_and_names() -> None:
    assert (
        detect_assembly(
            chrom_sizes={"chr1": 248956422, "chr2": 242193529}
        )
        == "hg38"
    )
    assert (
        detect_assembly(chrom_sizes={"1": 195471971, "2": 182113224}) == "mm10"
    )
    assert detect_assembly(peaks=pd.DataFrame({"chrom": ["I", "II"]})) == "sacCer3"
    with pytest.raises(ValueError, match="hg38/hg19"):
        detect_assembly(peaks=pd.DataFrame({"chrom": ["chr1", "chr2"]}))


def test_annotate_peak_harmonises_chromosome_names(tmp_path: Path) -> None:
    source = Path(__file__).resolve().parents[1] / "examples" / "example.gtf"
    gtf = tmp_path / "ensembl.gtf"
    gtf.write_text(
        "\n".join(
            line
            for line in source.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ).replace("chr", "")
        + "\n",
        encoding="utf-8",
    )
    peaks = pd.DataFrame(
        {"chrom": ["chr1"], "start": [1000], "end": [1100], "name": ["p"]}
    )
    annotated = annotate_peak(peaks, gtf, upstream=500, downstream=500)
    assert len(annotated.annotation) == 1
    assert annotated.annotation.iloc[0]["gene_id"] == "G1"
