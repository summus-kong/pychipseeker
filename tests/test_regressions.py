import pandas as pd

from pychipseeker import annotate_peak
from pychipseeker.io import read_gtf


def test_gff3_transcript_to_gene_mapping(tmp_path):
    path = tmp_path / "example.gff3"
    path.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tID=gene:G1;Name=GeneA\n"
        "chr1\ttest\tmRNA\t1001\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\ttest\texon\t1001\t1200\t.\t+\t.\tParent=transcript:T1\n",
        encoding="utf-8",
    )
    df = read_gtf(path)
    exon = df[df["feature"] == "exon"].iloc[0]
    gene = df[df["feature"] == "gene"].iloc[0]
    assert exon["gene_id"] == "G1"
    assert exon["transcript_id"] == "T1"
    assert gene["gene_name"] == "GeneA"


def test_gff3_fasta_section_is_ignored(tmp_path):
    path = tmp_path / "example.gff3"
    path.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tID=gene:G1;Name=GeneA\n"
        "##FASTA\n"
        ">chr1\n"
        "ACGTACGTACGT\n",
        encoding="utf-8",
    )
    df = read_gtf(path)
    assert len(df) == 1
    assert df.iloc[0]["feature"] == "gene"


def test_promoter_uses_nearest_tss_only():
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "gene"],
            "gene_id": ["G1", "G2"],
            "gene_name": ["GeneA", "GeneB"],
            "chrom": ["chr1", "chr1"],
            "start": [5100, 6000],
            "end": [5500, 7000],
            "strand": ["+", "+"],
        }
    )
    peaks = pd.DataFrame(
        {"chrom": ["chr1"], "start": [5500], "end": [5600], "name": ["p1"]}
    )

    # ChIPseeker only tests the distance to the *nearest* TSS (G1 here).  The
    # peak is 400 bp downstream of G1, outside the asymmetric -1000/+100
    # promoter, even though G2 is 400 bp upstream and inside its upstream
    # window.  This is a deliberate difference from "any TSS in the window".
    ann = annotate_peak(peaks, gtf, upstream=1000, downstream=100)
    row = ann.to_dataframe().iloc[0]
    assert row["annotation"] == "Distal Intergenic"
    assert row["gene_name"] == "GeneA"
    assert row["nearest_gene_name"] == "GeneA"


def test_groupby_uses_peak_metadata(examples_dir):
    peaks = pd.read_csv(
        examples_dir / "example_peaks.bed",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    peaks["cluster"] = ["A"] * 5 + ["B"] * 6
    ann = annotate_peak(
        peaks,
        examples_dir / "example.gtf",
        upstream=10,
        downstream=10,
    )
    df = ann.to_dataframe(include_peak_metadata=True)
    assert "cluster" in df.columns
    summary = ann.summary(groupby="cluster")
    assert summary.loc["A", "Promoter"] == 1
    assert summary.loc["A", "5' UTR"] == 1
    assert summary.loc["B", "Promoter"] == 2
    assert summary.loc["B", "Intergenic"] == 1


def test_bed_without_name_gets_generated_ids(tmp_path, examples_dir):
    path = tmp_path / "peaks.bed"
    path.write_text("chr1\t950\t1050\n", encoding="utf-8")
    ann = annotate_peak(
        path,
        examples_dir / "example.gtf",
        upstream=10,
        downstream=10,
    )
    assert ann.to_dataframe().loc[0, "peak_id"] == "peak_1"


def test_dataframe_index_is_used_as_peak_id(examples_dir):
    peaks = pd.DataFrame(
        {"chrom": ["chr1"], "start": [950], "end": [1050]},
        index=pd.Index(["peakA"], name="peak"),
    )
    ann = annotate_peak(
        peaks,
        examples_dir / "example.gtf",
        upstream=10,
        downstream=10,
    )
    assert ann.to_dataframe().loc[0, "peak_id"] == "peakA"


