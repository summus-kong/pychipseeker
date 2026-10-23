import math

import numpy as np
import pandas as pd
import pytest

from pychipseeker import annotate_peak, seq2gene, seq2gene_links
from pychipseeker.annotation import GeneAnnotation, get_genomic_anno_stat
from pychipseeker.io import read_gtf


def _example(examples_dir, **kwargs):
    params = dict(upstream=10, downstream=10)
    params.update(kwargs)
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        **params,
    )


def test_gene_and_transcript_models(examples_dir):
    ann = _example(examples_dir)
    genes = ann.gene_annotation.genes.set_index("gene_id")
    assert genes.loc["G1", "tss"] == 1000
    assert genes.loc["G2", "tss"] == 6999
    assert genes.loc["G1", "strand"] == "+"
    assert genes.loc["G2", "strand"] == "-"
    assert genes.loc["G1", "gene_name"] == "GeneA"

    transcripts = ann.gene_annotation.transcripts.set_index("transcript_id")
    assert transcripts.loc["T1", "tss"] == 1000
    assert transcripts.loc["T2", "tss"] == 6999
    assert transcripts.loc["T1", "gene_name"] == "GeneA"

    exons = ann.gene_annotation.exons
    t2 = exons[exons["transcript_id"] == "T2"].sort_values("exon_rank")
    assert list(t2["exon_rank"]) == [1, 2, 3]
    assert list(t2["start"]) == [6800, 6000, 5000]


def test_annotation_categories(examples_dir):
    ann = _example(examples_dir)
    df = ann.to_dataframe().set_index("peak_id")
    expected = {
        "P_promoter_plus": "Promoter",
        "P_5utr_plus": "5' UTR",
        "P_exon_plus": "Exon (T1/G1, exon 1 of 3)",
        "P_intron_plus": "Intron (T1/G1, intron 1 of 2)",
        "P_3utr_plus": "3' UTR",
        "P_promoter_minus": "Promoter",
        # The peak end is exactly the TSS, so this is inside the -10..10
        # promoter window used by the test.
        "P_5utr_minus": "Promoter",
        "P_exon_minus": "Exon (T2/G2, exon 2 of 3)",
        "P_intron_minus": "Intron (T2/G2, intron 2 of 2)",
        "P_intergenic": "Distal Intergenic",
    }
    for peak_id, label in expected.items():
        assert df.loc[peak_id, "annotation"] == label, peak_id
    assert "P_unknown_chr" not in df.index


def test_unknown_chromosome_is_dropped_with_warning(examples_dir):
    with pytest.warns(RuntimeWarning, match="dropped"):
        ann = _example(examples_dir)
    assert len(ann.annotation) == 10


def test_strand_aware_distance_uses_peak_ends(examples_dir):
    ann = _example(examples_dir)
    df = ann.to_dataframe().set_index("peak_id")
    assert df.loc["P_promoter_plus", "distance_to_tss"] == 0
    assert df.loc["P_promoter_minus", "distance_to_tss"] == 0
    assert df.loc["P_5utr_plus", "distance_to_tss"] == 20
    assert df.loc["P_5utr_minus", "distance_to_tss"] == 10
    assert df.loc["P_exon_plus", "distance_to_tss"] == 100
    assert df.loc["P_3utr_plus", "distance_to_tss"] == 1920
    assert df.loc["P_exon_minus", "distance_to_tss"] == 850
    assert df.loc["P_intron_minus", "distance_to_tss"] == 1400


def test_default_promoter_bins(examples_dir):
    ann = _example(examples_dir, upstream=3000, downstream=3000)
    df = ann.to_dataframe().set_index("peak_id")
    assert df.loc["P_promoter_plus", "annotation"] == "Promoter (<=1kb)"
    assert df.loc["P_5utr_plus", "annotation"] == "Promoter (<=1kb)"
    assert df.loc["P_exon_plus", "annotation"] == "Promoter (<=1kb)"
    assert df.loc["P_intron_plus", "annotation"] == "Promoter (<=1kb)"
    assert df.loc["P_3utr_plus", "annotation"] == "Promoter (1-2kb)"
    assert df.loc["P_intron_minus", "annotation"] == "Promoter (1-2kb)"
    assert df.loc["P_intergenic", "annotation"] == "Distal Intergenic"


def test_previous_center_priority_example(examples_dir):
    ann = _example(examples_dir, upstream=3000, downstream=3000)
    df = ann.to_dataframe().set_index("peak_id")
    # Exon 1-3kb from TSS is still a promoter because promoter has priority.
    assert df.loc["P_exon_plus", "annotation"].startswith("Promoter")


def test_gene_assignment_and_nearest(examples_dir):
    ann = _example(examples_dir)
    df = ann.to_dataframe().set_index("peak_id")
    assert df.loc["P_intron_plus", "gene_name"] == "GeneA"
    assert df.loc["P_intron_plus", "transcript_id"] == "T1"
    assert df.loc["P_intron_minus", "gene_name"] == "GeneB"
    assert df.loc["P_exon_minus", "transcript_id"] == "T2"
    assert df.loc["P_intergenic", "gene_name"] == "GeneB"
    assert df.loc["P_intergenic", "distance_to_tss"] < 0


def test_summary_and_anno_stat(examples_dir):
    ann = _example(examples_dir)
    summary = ann.summary().set_index("annotation")
    assert summary.loc["Promoter", "count"] == 3
    assert summary.loc["5' UTR", "count"] == 1
    assert summary.loc["3' UTR", "count"] == 1
    assert summary.loc["Exon", "count"] == 2
    assert summary.loc["Intron", "count"] == 2
    assert summary.loc["Intergenic", "count"] == 1
    assert math.isclose(summary["percent"].sum(), 100.0)

    stat = ann.anno_stat().set_index("Feature")
    assert stat.loc["Promoter", "Count"] == 3
    assert stat.loc["1st Exon", "Count"] == 1
    assert stat.loc["Other Exon", "Count"] == 1
    assert stat.loc["Distal Intergenic", "Count"] == 1
    assert math.isclose(stat["Frequency"].sum(), 100.0)


def test_chipseeker_column_names(examples_dir):
    ann = _example(examples_dir)
    df = ann.to_chipseeker()
    for col in [
        "peakChr",
        "peakStart",
        "peakEnd",
        "geneId",
        "geneName",
        "transcriptId",
        "distanceToTSS",
    ]:
        assert col in df.columns


def test_seq2gene(examples_dir):
    ann = _example(examples_dir)
    genes = ann.seq2gene()
    assert set(genes) == {"G1", "G2"}
    links = ann.seq2gene_links()
    assert not links.empty
    assert set(links.loc[links["peak_id"] == "P_promoter_plus", "gene_id"]) == {"G1"}
    assert set(links.loc[links["peak_id"] == "P_intron_plus", "gene_id"]) == {"G1"}
    assert "host" in set(links["link_type"])


def test_seq2gene_standalone(examples_dir):
    gene_ann = GeneAnnotation.from_gtf(read_gtf(examples_dir / "example.gtf"))
    genes = seq2gene(
        examples_dir / "example_peaks.bed",
        gene_ann,
        upstream=10,
        downstream=10,
    )
    assert set(genes) == {"G1", "G2"}
    links = seq2gene_links(
        examples_dir / "example_peaks.bed",
        gene_ann,
        upstream=10,
        downstream=10,
    )
    assert not links.empty


def test_gene_level(examples_dir):
    ann = _example(examples_dir, level="gene")
    df = ann.to_dataframe().set_index("peak_id")
    assert df.loc["P_exon_plus", "gene_id"] == "G1"
    assert "transcript_id" not in df.columns
    assert df.loc["P_exon_plus", "distance_to_tss"] == 100


def test_invalid_peak_table():
    with pytest.raises(ValueError):
        annotate_peak(
            pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [10]}),
            pd.DataFrame(
                {
                    "feature": ["gene"],
                    "gene_id": ["G1"],
                    "gene_name": ["A"],
                    "chrom": ["chr1"],
                    "start": [0],
                    "end": [100],
                    "strand": ["+"],
                }
            ),
        )


def test_overlap_and_same_strand():
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "gene"],
            "gene_id": ["NEG", "POS"],
            "gene_name": ["Neg", "Pos"],
            "chrom": ["chr1", "chr1"],
            "start": [150, 1000],
            "end": [300, 1200],
            "strand": ["-", "+"],
        }
    )
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 290],
            "end": [200, 310],
            "name": ["p_star", "p_plus"],
            "strand": ["*", "+"],
        }
    )
    ann_all = annotate_peak(
        peaks.iloc[[0]],
        gtf,
        tss_region=(-100, 100),
        overlap="all",
        same_strand=True,
    )
    row = ann_all.to_dataframe().iloc[0]
    assert row["gene_id"] == "NEG"
    assert row["distance_to_tss"] == 100

    ann_plus = annotate_peak(
        peaks.iloc[[1]],
        gtf,
        tss_region=(-100, 100),
        overlap="all",
        same_strand=True,
    )
    row = ann_plus.to_dataframe().iloc[0]
    assert row["gene_id"] == "POS"
    assert row["distance_to_tss"] == -691

    ann_ignore = annotate_peak(
        peaks.iloc[[1]],
        gtf,
        tss_region=(-100, 100),
        overlap="TSS",
        ignore_overlap=True,
    )
    row = ann_ignore.to_dataframe().iloc[0]
    assert row["gene_id"] == "POS"
    assert row["distance_to_tss"] == -691


def test_downstream_annotation():
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "gene"],
            "gene_id": ["G1", "G2"],
            "gene_name": ["A", "B"],
            "chrom": ["chr1", "chr1"],
            "start": [1000, 5000],
            "end": [2000, 6000],
            "strand": ["+", "+"],
        }
    )
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [2100, 10000],
            "end": [2200, 10100],
            "name": ["down", "distal"],
        }
    )
    ann = annotate_peak(
        peaks,
        gtf,
        tss_region=(-100, 100),
        downstream_distance=3000,
    )
    df = ann.to_dataframe().set_index("peak_id")
    assert df.loc["down", "annotation"] == "Downstream (<1kb)"
    assert df.loc["distal", "annotation"] == "Distal Intergenic"
    assert ann.detail_genomic_annotation.loc[0, "downstream"]
    assert ann.detail_genomic_annotation.loc[1, "distal_intergenic"]


def test_flanking_genes(examples_dir):
    ann = _example(examples_dir, upstream=10, downstream=10)
    flank = ann.gene_annotation
    from pychipseeker import get_flanking_genes

    table = get_flanking_genes(
        flank,
        pd.DataFrame(
            {"chrom": ["chr1"], "start": [2900], "end": [2950], "name": ["p1"]}
        ),
        level="transcript",
        distance=1000,
    )
    assert table.loc[0, "flank_gene_ids"] == "G1"
    assert table.loc[0, "flank_gene_distances"] == "0"


def test_get_genomic_anno_stat():
    stat = get_genomic_anno_stat(
        [
            "Promoter (<=1kb)",
            "Promoter (1-2kb)",
            "Exon (T1/G1, exon 1 of 3)",
            "Exon (T1/G1, exon 2 of 3)",
            "Intron (T1/G1, intron 1 of 2)",
            "Downstream (<=300bp)",
            "Distal Intergenic",
        ]
    ).set_index("Feature")
    assert stat.loc["1st Exon", "Count"] == 1
    assert stat.loc["Other Exon", "Count"] == 1
    assert stat.loc["1st Intron", "Count"] == 1
    assert stat.loc["Downstream (<=300)", "Count"] == 1
    assert math.isclose(stat["Frequency"].sum(), 100.0)
