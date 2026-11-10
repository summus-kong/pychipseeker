from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture
def examples_dir() -> Path:
    return EXAMPLES


@pytest.fixture
def small_bam(tmp_path) -> Path:
    """A tiny indexed BAM used by the real-pysam regression tests.

    Layout (0-based, 1000 bp chromosomes)::

        chr1    0-50    normal
        chr1  100-150   normal  <- peak 90-200
        chr1  120-170   normal  <- peak 90-200
        chr1  400-450   normal
        chr1  130-180   duplicate
        chr1  300-350   secondary
        chr1  500-550   normal, low MAPQ
        chr2    0-50    normal
        unmapped

    Six reads pass the default filters, two of them overlap the peak.
    """

    pysam = pytest.importorskip("pysam")
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}],
    }
    path = tmp_path / "small.bam"
    reads = [
        ("chr1", 0, 50, 60, 0),
        ("chr1", 100, 50, 60, 0),
        ("chr1", 120, 50, 60, 0),
        ("chr1", 400, 50, 60, 0),
        ("chr1", 130, 50, 60, 1024),
        ("chr1", 300, 50, 60, 256),
        ("chr1", 500, 50, 5, 0),
        ("chr2", 0, 50, 60, 0),
    ]
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for index, (chrom, start, length, mapq, flag) in enumerate(reads):
            read = pysam.AlignedSegment()
            read.query_name = f"read{index}"
            read.query_sequence = "A" * length
            read.flag = flag
            read.reference_id = 0 if chrom == "chr1" else 1
            read.reference_start = start
            read.mapping_quality = mapq
            read.cigarstring = f"{length}M"
            handle.write(read)

        unmapped = pysam.AlignedSegment()
        unmapped.query_name = "unmapped"
        unmapped.query_sequence = "A" * 50
        unmapped.flag = 4
        unmapped.reference_id = -1
        unmapped.reference_start = -1
        unmapped.mapping_quality = 0
        handle.write(unmapped)

    pysam.index(str(path))
    return path
