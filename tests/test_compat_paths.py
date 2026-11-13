"""Branches of the ChIPseeker-compatibility helpers not covered elsewhere."""

from __future__ import annotations

import copy
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pychipseeker import annotate_peak, compat
from pychipseeker.compat import (
    _annotation_gene_ids,
    _as_annotation_df,
    _chrom_sizes_from_peaks,
    _count_query_overlaps,
    _hypergeom_sf,
    _load_gsminfo,
    _log_comb,
    _p_adjust,
    _set_elements,
    annotate_db,
    as_granges,
    combine_cs_anno,
    covplot,
    detail_genomic_annotation,
    drop_anno,
    enrich_peak_overlap,
    get_anno_stat,
    get_sample_files,
    overlap,
    read_peak_file,
    shuffle,
    upsetplot,
    vennpie,
    vennplot,
)


@pytest.fixture
def annotation(examples_dir):
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        upstream=1000,
        downstream=1000,
    )


def test_as_annotation_df_accepts_both_inputs(annotation) -> None:
    from_annotation = _as_annotation_df(annotation)
    assert {"peakChr", "peakStart", "peakEnd"}.issubset(from_annotation.columns)

    frame = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    copied = _as_annotation_df(frame)
    assert copied is not frame
    assert copied.equals(frame)

    with pytest.raises(TypeError, match="PeakAnnotation or pandas.DataFrame"):
        _as_annotation_df([1, 2, 3])


def test_as_granges_and_detail_accept_dataframes() -> None:
    frame = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    granges = as_granges(frame)
    assert {"seqnames", "start", "end"}.issubset(granges.columns)

    detail = detail_genomic_annotation(frame)
    assert detail is not frame
    assert detail.equals(frame)

    with pytest.raises(TypeError, match="PeakAnnotation or pandas.DataFrame"):
        as_granges("peaks.bed")
    with pytest.raises(TypeError, match="PeakAnnotation or pandas.DataFrame"):
        detail_genomic_annotation("peaks.bed")


def test_read_peak_file_accepts_a_dataframe() -> None:
    frame = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    loaded = read_peak_file(frame)
    assert loaded is not frame
    assert loaded.equals(frame)


def test_get_anno_stat_dataframe_variants() -> None:
    stat = pd.DataFrame({"Feature": ["Promoter"], "Frequency": [1]})
    assert get_anno_stat(stat).equals(stat)

    labelled = pd.DataFrame(
        {
            "annotation": [
                "Promoter (<=1kb)",
                "Exon (1 of 2)",
                "Promoter (<=1kb)",
            ]
        }
    )
    counts = get_anno_stat(labelled).set_index("Feature")["Frequency"]
    assert counts["Promoter (<=1kb)"] == pytest.approx(200 / 3)
    assert counts.sum() == pytest.approx(100.0)

    # Column-name matching is case-insensitive.
    upper = pd.DataFrame({"Annotation": ["Exon (1 of 2)", "Intron (1 of 2)"]})
    assert "Other Exon" in set(get_anno_stat(upper)["Feature"])

    with pytest.raises(TypeError, match="annotation labels"):
        get_anno_stat(pd.DataFrame({"chrom": ["chr1"]}))


def test_drop_anno_dataframe_distance_spellings() -> None:
    snake = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [0, 10, 20],
            "end": [5, 15, 25],
            "distance_to_tss": [10, 5000, -2000],
        }
    )
    kept = drop_anno(snake, 3000)
    assert kept["distance_to_tss"].tolist() == [10, -2000]
    assert list(kept.index) == [0, 1]

    camel = snake.rename(columns={"distance_to_tss": "distanceToTSS"})
    assert len(drop_anno(camel, 3000)) == 2

    with pytest.raises(ValueError, match="distance_to_tss"):
        drop_anno(snake.drop(columns=["distance_to_tss"]))
    with pytest.raises(TypeError, match="PeakAnnotation or pandas.DataFrame"):
        drop_anno(object())


def test_annotate_db_dataframe_and_mapping(annotation) -> None:
    table = pd.DataFrame({"geneId": ["G1", "G2"], "SYMBOL": ["A", "B"]})
    result = annotate_db(copy.deepcopy(annotation), table, columns=["SYMBOL"])
    assert set(result.annotation["SYMBOL"].dropna()) == {"A", "B"}

    mapping = {"G1": {"SYMBOL": "A", "ENTREZID": "1"}, "G2": {"SYMBOL": "B"}}
    mapped = annotate_db(copy.deepcopy(annotation), mapping, columns=["SYMBOL"])
    assert "ENTREZID" not in mapped.annotation.columns
    assert set(mapped.annotation["SYMBOL"].dropna()) == {"A", "B"}

    scalars = annotate_db(copy.deepcopy(annotation), {"G1": "A"})
    assert set(scalars.annotation["SYMBOL"].dropna()) == {"A"}

    with pytest.raises(TypeError, match="PeakAnnotation"):
        annotate_db(annotation.annotation, table)
    with pytest.raises(TypeError, match="DataFrame or a mapping"):
        annotate_db(copy.deepcopy(annotation), [("G1", "A")])

    without_id = copy.deepcopy(annotation)
    without_id.annotation = without_id.annotation.drop(columns=["gene_id"])
    with pytest.raises(ValueError, match="gene_id"):
        annotate_db(without_id, table)


def test_combine_cs_anno_guards(annotation, examples_dir) -> None:
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
    assert combine_cs_anno([left, right]).peakNum == annotation.peakNum

    with pytest.raises(ValueError, match="at least one"):
        combine_cs_anno()
    with pytest.raises(TypeError, match="PeakAnnotation"):
        combine_cs_anno(left, pd.DataFrame())

    other_tss = annotate_peak(
        annotation.peaks,
        annotation.gene_annotation,
        tss_region=(-500, 500),
    )
    with pytest.raises(ValueError, match="tss_region"):
        combine_cs_anno(left, other_tss)

    other_level = annotate_peak(
        annotation.peaks,
        annotation.gene_annotation,
        tss_region=annotation.tss_region,
        level="gene",
    )
    with pytest.raises(ValueError, match="same level"):
        combine_cs_anno(left, other_level)

    no_genomic = annotate_peak(
        annotation.peaks,
        annotation.gene_annotation,
        tss_region=annotation.tss_region,
        assign_genomic_annotation=False,
    )
    with pytest.raises(ValueError, match="has_genomic_annotation"):
        combine_cs_anno(left, no_genomic)


def test_set_elements_variants(annotation) -> None:
    assert _set_elements(annotation) == set(annotation.annotation["peak_id"].astype(str))

    no_ids = copy.deepcopy(annotation)
    no_ids.annotation = no_ids.annotation.drop(columns=["peak_id"])
    assert _set_elements(no_ids) == set(no_ids.annotation.index.astype(str))

    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    assert _set_elements(peaks) == {"chr1:0-10"}
    assert _set_elements({"a": 1, "b": 2}) == {"a", "b"}
    assert _set_elements(["a", "b"]) == {"a", "b"}


def test_overlap_and_vennplot_guards() -> None:
    with pytest.raises(ValueError, match="empty"):
        overlap([])

    with pytest.raises(ValueError, match="empty"):
        vennplot([])
    with pytest.raises(ValueError, match="labels length"):
        vennplot([{1}, {2}], labels=["only"])


def test_vennplot_single_set_and_large_overlap() -> None:
    fig, ax = vennplot([{"a", "b", "c"}], labels=["one"])
    assert ax is not None
    plt.close(fig)

    fig, ax = vennplot({"A": {1, 2}, "B": {2, 3}, "C": {3, 4}, "D": {4, 5}})
    assert ax is not None
    plt.close(fig)

    fig, ax = plt.subplots()
    returned_fig, returned_ax = vennplot([{1}, {2}], ax=ax)
    assert returned_fig is fig and returned_ax is ax
    plt.close(fig)


def test_vennpie_guards_and_existing_axes(annotation) -> None:
    with pytest.raises(TypeError, match="PeakAnnotation"):
        vennpie(annotation.annotation)

    empty = copy.deepcopy(annotation)
    empty.detail_genomic_annotation = empty.detail_genomic_annotation.iloc[0:0]
    with pytest.raises(ValueError, match="detailGenomicAnnotation"):
        vennpie(empty)

    fig, ax = plt.subplots()
    returned_fig, returned_ax = vennpie(annotation, ax=ax)
    assert returned_fig is fig and returned_ax is ax
    plt.close(fig)


def test_covplot_list_chrs_xlim_and_lower(annotation) -> None:
    peaks = annotation.peaks
    chroms = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chrX"],
            "start": [0, 0, 0],
            "end": [100, 200, 50],
        }
    )
    fig, axes = covplot([peaks, chroms], title="combined")
    assert len(axes) == 4
    plt.close(fig)

    fig, axes = covplot(
        chroms,
        chrs=["chr1", "chr2", "chrX"],
        xlim=(0, 150),
        lower=0.5,
        weightCol="end",
        figsize=(5, 4),
    )
    assert len(axes) == 4
    assert all(len(ax.lines) >= 1 for ax in axes[:3])
    assert not axes[3].axison
    plt.close(fig)

    fig, ax = covplot(chroms, chrs=["chr2"])
    assert ax.get_title() == "chr2"
    plt.close(fig)

    with pytest.raises(ValueError, match="no chromosomes"):
        covplot(peaks, chrs=["chrZ"])


def test_chrom_sizes_from_peaks_empty() -> None:
    empty = pd.DataFrame({"chrom": [], "start": [], "end": []})
    assert _chrom_sizes_from_peaks(empty) == {}


def test_shuffle_n_and_size_sources(
    annotation, examples_dir, monkeypatch: pytest.MonkeyPatch
) -> None:
    peaks = annotation.peaks
    sampled = shuffle(peaks, n=3, seed=1)
    assert len(sampled) == 3

    expanded = shuffle(peaks, n=len(peaks) + 2, seed=1)
    assert len(expanded) == len(peaks) + 2

    from_gene_annotation = shuffle(peaks, TxDb=annotation.gene_annotation, seed=2)
    assert len(from_gene_annotation) == len(peaks)

    from_gtf = shuffle(peaks, TxDb=examples_dir / "example.gtf", seed=2)
    assert len(from_gtf) == len(peaks)

    monkeypatch.setattr(
        "pychipseeker.genomes.load_genome", lambda name: annotation.gene_annotation
    )
    from_assembly = shuffle(peaks, TxDb="hg38", seed=2)
    assert len(from_assembly) == len(peaks)

    monkeypatch.setattr(
        "pychipseeker.genomes.load_annotation", lambda spec: annotation.gene_annotation
    )
    from_spec = shuffle(peaks, TxDb="ensembl:mouse@110", seed=2)
    assert len(from_spec) == len(peaks)

    sizes = pd.DataFrame({"chrom": ["chr1", "chrX"], "size": [5000, 5000]})
    from_frame = shuffle(peaks, chrom_sizes=sizes, seed=3)
    assert set(from_frame["chrom"]) == {"chr1", "chrX"}


def test_count_query_overlaps_empty_input() -> None:
    empty = pd.DataFrame({"chrom": [], "start": [], "end": []})
    target = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [10]})
    assert _count_query_overlaps(empty, target) == 0
    assert _count_query_overlaps(target, empty) == 0


def test_p_adjust_methods() -> None:
    np.testing.assert_allclose(_p_adjust([0.01, 0.2], "none"), [0.01, 0.2])
    np.testing.assert_allclose(_p_adjust([0.01, 0.2], "nop"), [0.01, 0.2])
    assert _p_adjust([], "bonferroni").size == 0
    np.testing.assert_allclose(
        _p_adjust([0.01, 0.2], "bonferroni"), [0.02, 0.4]
    )
    with pytest.raises(ValueError, match="pAdjustMethod"):
        _p_adjust([0.1], "holm")


def test_enrich_peak_overlap_dataframe_sizes_and_no_shuffle(annotation) -> None:
    sizes = pd.DataFrame({"chrom": ["chr1", "chrX"], "end": [5000, 5000]})
    result = enrich_peak_overlap(
        annotation.peaks,
        annotation.peaks,
        chrom_sizes=sizes,
        nShuffle=2,
        seed=1,
    )
    assert result["N_OL"].iloc[0] == len(annotation.peaks)
    assert 0 < result["pvalue"].iloc[0] <= 1

    static = enrich_peak_overlap(
        annotation.peaks, annotation.peaks, nShuffle=0, seed=1
    )
    assert np.isnan(static["pvalue"].iloc[0])
    assert static["p.adjust"].iloc[0] == 1.0


def test_hypergeom_helpers_bounds() -> None:
    assert _log_comb(5, -1) == float("-inf")
    assert _log_comb(5, 6) == float("-inf")
    assert _hypergeom_sf(0, 100, 10, 20) == 1.0
    assert _hypergeom_sf(11, 100, 10, 20) == 0.0


def test_annotation_gene_ids_errors(annotation, monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(ValueError, match="provide TxDb"):
        _annotation_gene_ids(annotation.peaks)

    stub = SimpleNamespace(
        annotation=pd.DataFrame({"peak_id": [1]}),
        gene_annotation=SimpleNamespace(genes=pd.DataFrame({"gene_id": ["G1"]})),
    )
    monkeypatch.setattr(compat, "annotate_peak", lambda *a, **k: stub)
    with pytest.raises(ValueError, match="gene_id"):
        _annotation_gene_ids(annotation.peaks, annotation=annotation.gene_annotation)


def test_upsetplot_guards_and_options(annotation) -> None:
    with pytest.raises(TypeError, match="PeakAnnotation"):
        upsetplot(annotation.annotation)

    empty = copy.deepcopy(annotation)
    empty.detail_genomic_annotation = empty.detail_genomic_annotation.iloc[0:0]
    with pytest.raises(ValueError, match="detailGenomicAnnotation"):
        upsetplot(empty)

    no_flags = copy.deepcopy(annotation)
    no_flags.detail_genomic_annotation = pd.DataFrame({"label": ["a", "b"]})
    with pytest.raises(ValueError, match="boolean annotation columns"):
        upsetplot(no_flags)

    fig, ax = upsetplot(annotation, order_by="degree")
    assert "Annotation intersections" == ax.get_title()
    plt.close(fig)

    fig, ax = plt.subplots()
    returned_fig, returned_ax = upsetplot(annotation, ax=ax)
    assert returned_fig is fig and returned_ax is not None
    plt.close(fig)


def test_load_gsminfo_dataframe_path_and_errors(tmp_path: Path) -> None:
    frame = pd.DataFrame({"organism": ["human"], "genomeVersion": ["hg38"]})
    assert _load_gsminfo(frame).equals(frame)
    assert _load_gsminfo(frame) is not frame

    csv = tmp_path / "gsminfo.csv"
    frame.to_csv(csv, index=False)
    assert _load_gsminfo(csv)["organism"].tolist() == ["human"]

    gz = tmp_path / "gsminfo.csv.gz"
    frame.to_csv(gz, index=False, compression="gzip")
    assert _load_gsminfo(str(gz))["organism"].tolist() == ["human"]

    with pytest.raises(FileNotFoundError):
        _load_gsminfo(tmp_path / "missing.csv")
    with pytest.raises(TypeError, match="gsminfo must be"):
        _load_gsminfo([["human"]])


def test_download_geo_default_destination(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import urllib.request

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self):
            return b"data"

    monkeypatch.setattr(urllib.request, "urlopen", lambda url, timeout=None: Response())
    monkeypatch.setattr(Path, "cwd", classmethod(lambda cls: tmp_path))
    info = pd.DataFrame({"supplementary_file": ["https://example.org/sample.bed.gz"]})
    written = compat._download_geo_internal(info, None)
    assert written == [tmp_path / "sample.bed.gz"]
    assert written[0].read_bytes() == b"data"


def test_get_sample_files_missing_directory(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="sample data directory"):
        get_sample_files(tmp_path / "missing")
