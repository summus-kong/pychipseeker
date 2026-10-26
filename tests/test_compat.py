import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pychipseeker import annotate_peak
from pychipseeker.compat import (
    as_granges,
    combine_cs_anno,
    covplot,
    drop_anno,
    enrich_anno_overlap,
    enrich_peak_overlap,
    get_anno_stat,
    overlap,
    read_peak_file,
    shuffle,
    vennpie,
    vennplot,
    vennplot_peakfile,
)


@pytest.fixture
def annotation(examples_dir):
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        upstream=1000,
        downstream=1000,
    )


def test_read_peak_file_and_as_granges(examples_dir, annotation):
    from pychipseeker import readPeakFile

    assert readPeakFile is read_peak_file
    df = read_peak_file(examples_dir / "example_peaks.bed")
    assert {"chrom", "start", "end"}.issubset(df.columns)
    gr = read_peak_file(examples_dir / "example_peaks.bed", as_type="GRanges")
    assert {"seqnames", "start", "end"}.issubset(gr.columns)
    gr_alias = readPeakFile(examples_dir / "example_peaks.bed", **{"as": "GRanges"})
    assert {"seqnames", "start", "end"}.issubset(gr_alias.columns)
    ann_gr = as_granges(annotation)
    assert {"seqnames", "start", "end", "distanceToTSS"}.issubset(ann_gr.columns)


def test_annotation_statistics_and_drop(annotation):
    stat = get_anno_stat(annotation)
    assert {"Feature", "Frequency", "Count"}.issubset(stat.columns)
    filtered = drop_anno(annotation, 1000)
    assert filtered.peakNum == len(filtered.annotation)
    assert filtered.peakNum < annotation.peakNum
    assert filtered.annotation["distance_to_tss"].abs().max() < 1000
    # dropAnno filters peaks; it does not remove annotation columns.
    assert "annotation" in filtered.annotation.columns


def test_chipseeker_annotatepeak_wrapper(examples_dir):
    from pychipseeker import annotatePeak, annotate_peak

    # A bare alias: the R-style keyword spellings live in the engine itself.
    assert annotatePeak is annotate_peak
    ann = annotatePeak(
        examples_dir / "example_peaks.bed",
        TxDb=examples_dir / "example.gtf",
        tssRegion=(-1000, 1000),
    )
    assert ann.tss_region == (-1000, 1000)
    assert ann.peakNum == 11
    ann2 = annotatePeak(
        examples_dir / "example_peaks.bed",
        gtf=examples_dir / "example.gtf",
        annoDb=pd.DataFrame(
            {"gene_id": ["G1", "G2", "G3"], "SYMBOL": ["A", "B", "C"]}
        ),
        columns=["SYMBOL"],
    )
    assert "SYMBOL" in ann2.annotation.columns
    assert set(ann2.annotation["SYMBOL"].dropna()).issubset({"A", "B", "C"})
    assert "A" in set(ann2.annotation["SYMBOL"].dropna())

    # R spellings of the peak/source arguments and the default priority.
    ann3 = annotatePeak(
        peak=examples_dir / "example_peaks.bed",
        annotation=examples_dir / "example.gtf",
        tssRegion=(-1000, 1000),
        genomicAnnotationPriority=None,
    )
    assert ann3.peakNum == ann.peakNum
    assert ann3.genomic_annotation_priority == ann.genomic_annotation_priority

    with pytest.raises(TypeError, match="unexpected keyword"):
        annotatePeak(
            examples_dir / "example_peaks.bed",
            gtf=examples_dir / "example.gtf",
            bogus=1,
        )


def test_combine_annotations(annotation):
    left = annotate_peak(
        annotation.peaks.iloc[:5],
        annotation.gene_annotation,
        tss_region=annotation.tss_region,
    )
    right = annotate_peak(
        annotation.peaks.iloc[5:],
        annotation.gene_annotation,
        tss_region=annotation.tss_region,
    )
    combined = combine_cs_anno(left, right)
    assert combined.peakNum == annotation.peakNum
    assert len(combined.annotation) == len(annotation.annotation)
    assert len(combined.detail_genomic_annotation) == len(annotation.annotation)


def test_overlap_and_venn_plots(annotation):
    table = overlap([{1, 2, 3}, {2, 3, 4}])
    rows = table.set_index(["Set1", "Set2"])["Weight"].to_dict()
    assert rows[(1, 1)] == 2
    assert rows[(1, 0)] == 1
    assert rows[(0, 1)] == 1

    fig, ax = vennplot([{1, 2, 3}, {2, 3, 4}], labels=["A", "B"])
    assert ax is not None
    plt.close(fig)

    fig, ax = vennplot([{1, 2}, {2, 3}, {3, 4}], labels=["A", "B", "C"])
    assert ax is not None
    plt.close(fig)

    fig, ax = vennpie(annotation)
    assert "Genomic annotation" == ax.get_title()
    plt.close(fig)


def test_vennplot_peakfile(examples_dir):
    fig, ax = vennplot_peakfile(
        [examples_dir / "example_peaks.bed", examples_dir / "example_peaks.bed"],
        labels=["A", "B"],
    )
    assert ax is not None
    plt.close(fig)


def test_upsetplot(annotation):
    from pychipseeker import upsetplot

    fig, ax = upsetplot(annotation)
    assert "Annotation intersections" == ax.get_title()
    plt.close(fig)


def test_peak_heatmap_multiple_sets(annotation):
    from pychipseeker import peakHeatmap_multiple_Sets

    fig, axes = peakHeatmap_multiple_Sets(
        [annotation.peaks, annotation.peaks.iloc[:5]],
        labels=["all", "subset"],
        annotation=annotation.gene_annotation,
        upstream=100,
        downstream=100,
        nbin=20,
    )
    assert len(axes) == 2
    plt.close(fig)


def test_plot_peak_prof_multi_windows(annotation):
    from pychipseeker import getPromoters, plotPeakProf_MultiWindows

    windows1 = getPromoters(
        annotation=annotation.gene_annotation, upstream=100, downstream=100
    )
    windows2 = getPromoters(
        annotation=annotation.gene_annotation, upstream=200, downstream=200
    )
    fig, ax = plotPeakProf_MultiWindows(
        annotation.peaks,
        [windows1, windows2],
        windows_name=["small", "large"],
        nbin=20,
    )
    assert ax is not None
    plt.close(fig)


def test_covplot(annotation):
    fig, axes = covplot(annotation.peaks)
    assert fig is not None
    assert axes is not None
    plt.close(fig)


def test_shuffle_preserves_width_and_chromosome(annotation):
    shuffled = shuffle(annotation.peaks, seed=1)
    assert len(shuffled) == len(annotation.peaks)
    assert set(shuffled["chrom"]) == set(annotation.peaks["chrom"])
    np.testing.assert_array_equal(
        shuffled["end"].to_numpy() - shuffled["start"].to_numpy(),
        annotation.peaks["end"].to_numpy() - annotation.peaks["start"].to_numpy(),
    )


def test_overlap_enrichment(annotation):
    result = enrich_peak_overlap(
        annotation.peaks,
        [annotation.peaks, annotation.peaks.iloc[:3]],
        nShuffle=3,
        seed=2,
    )
    assert {"qSample", "tSample", "qLen", "tLen", "N_OL", "pvalue", "p.adjust"}.issubset(
        result.columns
    )
    assert result["N_OL"].iloc[0] == len(annotation.peaks)
    assert (result["pvalue"] >= 0).all()
    assert (result["pvalue"] <= 1).all()


def test_anno_enrichment(annotation):
    result = enrich_anno_overlap(
        annotation.peaks,
        [annotation.peaks, annotation.peaks.iloc[:3]],
        annotation=annotation.gene_annotation,
    )
    assert len(result) == 2
    assert (result["N_OL"] >= 0).all()
    assert (result["pvalue"] >= 0).all()
    assert (result["pvalue"] <= 1).all()


# ----------------------------------------------------------------------
# ChIPseeker GEO helper tests
# ----------------------------------------------------------------------
def _geo_snapshot() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "series_id": ["GSE1", "GSE1", "GSE2", "GSE3"],
            "gsm": ["GSM1", "GSM2", "GSM3", "GSM4"],
            "gpl": ["GPL1", "GPL1", "GPL2", "GPL2"],
            "organism": ["Homo sapiens", "Homo sapiens", "Mus musculus", "Homo sapiens"],
            "title": ["t1", "t2", "t3", "t4"],
            "supplementary_file": [
                "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1_bed.bed.gz",
                "ftp://example.org/GSM2_Peaks.narrowPeak.gz",
                None,
                "https://example.org/GSM4_x.bedGraph.gz",
            ],
            "genomeVersion": ["hg19", "hg19", "mm9", None],
            "pubmed_id": [123.0, 123.0, None, 456.0],
        }
    )


def test_get_geo_species_snapshot():
    from pychipseeker.compat import getGEOspecies

    species = getGEOspecies(_geo_snapshot())
    assert list(species.columns) == ["species", "Freq"]
    # R's as.data.frame(table(species)) is alphabetically sorted.
    assert list(species["species"]) == ["Homo sapiens", "Mus musculus"]
    assert list(species["Freq"]) == [3, 1]


def test_get_geo_genome_version_snapshot():
    from pychipseeker.compat import getGEOgenomeVersion

    gv = getGEOgenomeVersion(_geo_snapshot())
    assert list(gv.columns) == ["organism", "genomeVersion", "Freq"]
    genome_versions = list(gv["genomeVersion"])
    assert "hg19" in genome_versions and "mm9" in genome_versions
    freq = dict(zip(gv["genomeVersion"], gv["Freq"]))
    assert freq["hg19"] == 2
    assert freq["mm9"] == 1


def test_get_geo_info_simplify():
    from pychipseeker.compat import getGEOInfo

    hg19 = getGEOInfo("hg19", gsminfo=_geo_snapshot())
    assert list(hg19["gsm"]) == ["GSM1", "GSM2"]
    assert [
        "series_id",
        "gsm",
        "organism",
        "title",
        "supplementary_file",
        "genomeVersion",
        "pubmed_id",
    ] == list(hg19.columns)

    full = getGEOInfo("hg19", simplify=False, gsminfo=_geo_snapshot())
    assert {"gpl", "title", "genomeVersion"}.issubset(full.columns)
    assert len(full) == 2
    full = getGEOInfo("mm9", simplify=False, gsminfo=_geo_snapshot())
    assert list(full["gsm"]) == ["GSM3"]
    assert {"gpl", "title", "genomeVersion"}.issubset(full.columns)


def test_get_geo_info_missing_genome():
    from pychipseeker.compat import getGEOInfo

    empty = getGEOInfo("hg38", gsminfo=_geo_snapshot())
    assert len(empty) == 0


def test_get_geo_info_bundled_snapshot():
    from pychipseeker.compat import getGEOInfo, getGEOgenomeVersion, getGEOspecies

    species = getGEOspecies()
    assert {"species", "Freq"}.issubset(species.columns)
    assert len(species) > 0
    genome_versions = getGEOgenomeVersion()
    assert {"organism", "genomeVersion", "Freq"}.issubset(genome_versions.columns)
    hg19 = getGEOInfo("hg19")
    assert len(hg19) > 0
    assert (hg19["genomeVersion"] == "hg19").all()


def test_download_gsm_bed_files(monkeypatch, tmp_path):
    import urllib.request

    from pychipseeker import compat

    class _Response:
        def __init__(self, payload):
            self._payload = payload

        def read(self):
            return self._payload

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    requested = []

    def fake_urlopen(url, timeout=None):
        requested.append(url)
        return _Response(b"peak-data")

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

    snapshot = _geo_snapshot()
    out = compat.downloadGSMbedFiles(["GSM1", "GSM3"], destDir=tmp_path, gsminfo=snapshot)
    # GSM3 has no supplementary file; GSM1 is downloaded and ftp->https
    assert len(out) == 1
    dest = out[0]
    assert dest.parent == tmp_path and dest.name == "GSM1_bed.bed.gz"
    assert requested[0].startswith("https://")
    assert dest.read_bytes() == b"peak-data"

    # existing destination files are not re-downloaded
    again = compat.downloadGSMbedFiles("GSM1", destDir=tmp_path, gsminfo=snapshot)
    assert again == []
    assert len(requested) == 1


def test_download_geo_bed_files_skips_failures(monkeypatch, tmp_path):
    import urllib.request

    from pychipseeker import compat

    def failing_urlopen(url, timeout=None):
        raise OSError("boom")

    monkeypatch.setattr(urllib.request, "urlopen", failing_urlopen)
    snapshot = _geo_snapshot()
    out = compat.downloadGEObedFiles("hg19", destDir=tmp_path, gsminfo=snapshot)
    assert out == []


def test_get_sample_files_labels_and_paths():
    from pychipseeker import read_peak_file
    from pychipseeker.compat import getSampleFiles

    files = getSampleFiles()
    assert set(files) == {"ARmo_0M", "ARmo_1nM", "ARmo_100nM", "CBX6_BF", "CBX7_BF"}
    for path in files.values():
        assert path.exists()
    bed = read_peak_file(files["CBX7_BF"])
    assert {"chrom", "start", "end"}.issubset(bed.columns)
    assert len(bed) > 0

    # custom directory support preserves the R naming rules
    custom = getSampleFiles(files["CBX7_BF"].parent)
    assert set(custom) == set(files)
