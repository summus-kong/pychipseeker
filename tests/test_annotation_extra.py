"""Branch coverage for :mod:`pychipseeker.annotation`.

The main annotation semantics are exercised in ``test_annotation.py``; this
module covers the smaller helpers, GTF normalisation fallbacks, legacy
accessors and the parallel-annotation plumbing.
"""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pychipseeker.annotation as ann_mod
from pychipseeker import annotate_peak, get_flanking_genes, seq2gene, seq2gene_links
from pychipseeker.annotation import (
    GeneAnnotation,
    PeakAnnotation,
    _annotate_peaks_frame,
    _annotation_worker,
    _annotation_worker_init,
    _annotation_worker_init_from_source,
    _downstream_annotation,
    _ensure_columns,
    _format_kb,
    _load_annotation_model,
    _nearest_feature,
    _normalise_peaks,
    _normalise_priority,
    _promoter_sort_key,
    _sort_by_transcript_order,
    _split_generic_utr,
    get_genomic_anno_stat,
    group_annotation,
    promoter_label,
)


def _example(examples_dir: Path, **kwargs) -> PeakAnnotation:
    params = dict(upstream=10, downstream=10)
    params.update(kwargs)
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        **params,
    )


def _config(**overrides) -> dict[str, object]:
    config: dict[str, object] = {
        "level": "transcript",
        "tss_region": (-100, 100),
        "bins": (1000,),
        "priority": ann_mod._DEFAULT_PRIORITY,
        "assign_genomic_annotation": True,
        "same_strand": False,
        "ignore_overlap": False,
        "ignore_upstream": False,
        "ignore_downstream": False,
        "overlap": "TSS",
        "align_annotation_feature": True,
        "downstream_distance": 300,
        "verbose": False,
        "total": 0,
    }
    config.update(overrides)
    return config


# ----------------------------------------------------------------------
# label / priority helpers
# ----------------------------------------------------------------------
def test_format_kb_and_promoter_labels() -> None:
    assert _format_kb(1000) == "1kb"
    assert _format_kb(1500) == "1.5kb"

    assert promoter_label(10, ()) == "Promoter"
    assert promoter_label(10, (None,)) == "Promoter"  # type: ignore[arg-type]
    assert promoter_label(1500, (1000,)) == "Promoter"
    assert promoter_label(5000, (1000, 2000)) == "Promoter"
    assert promoter_label(1500, (1000, 2000, 3000)) == "Promoter (1-2kb)"
    assert promoter_label(2500, (1000, 2000, 3000)) == "Promoter (2-3kb)"

    assert _promoter_sort_key("Promoter") == (0, 0.0)
    assert _promoter_sort_key("Promoter (<=1kb)") == (1, 1.0)
    assert _promoter_sort_key("Promoter (1-2kb)") == (1, 2.0)
    assert _promoter_sort_key("Something else") == (2, 0.0)


def test_group_annotation_and_priority_validation() -> None:
    assert group_annotation("Downstream (<=300bp)") == "Downstream"
    assert group_annotation("Exon (T1/G1, exon 1 of 3)") == "Exon"
    assert group_annotation("Intron (T1/G1, intron 1 of 2)") == "Intron"
    assert group_annotation("5' UTR") == "5' UTR"
    assert group_annotation("Unclassified") == "Unclassified"

    assert _normalise_priority(list(reversed(ann_mod._DEFAULT_PRIORITY))) == tuple(
        reversed(ann_mod._DEFAULT_PRIORITY)
    )
    with pytest.raises(ValueError, match="must contain any order"):
        _normalise_priority(["Promoter", "Nope"])
    with pytest.raises(ValueError, match="permutation"):
        _normalise_priority(["Promoter", "Promoter"])


def test_ensure_columns_and_empty_transcript_sort() -> None:
    table = pd.DataFrame({"start": [1]})
    filled = _ensure_columns(table, ["start", "end"])
    assert list(filled.columns) == ["start", "end"]
    assert filled["end"].isna().all()
    assert "end" not in table.columns

    empty = pd.DataFrame(columns=["transcript_id", "start", "end"])
    assert _sort_by_transcript_order(empty, {"T1": 0}).empty


# ----------------------------------------------------------------------
# generic UTR splitting
# ----------------------------------------------------------------------
def test_generic_utr_split_uses_cds_bounds() -> None:
    gtf = pd.DataFrame(
        {
            "feature": [
                "gene",
                "gene",
                "transcript",
                "transcript",
                "cds",
                "cds",
                "utr",
                "utr",
                "utr",
                "utr",
            ],
            "gene_id": ["GP", "GM", "GP", "GM", "GP", "GM", "GP", "GP", "GM", "GM"],
            "gene_name": ["P", "M", "P", "M", "P", "M", "P", "P", "M", "M"],
            "transcript_id": [
                None,
                None,
                "TP",
                "TM",
                "TP",
                "TM",
                "TP",
                "TP",
                "TM",
                "TM",
            ],
            "chrom": ["chr1"] * 10,
            "start": [100, 1200, 100, 1200, 200, 1300, 100, 400, 1500, 1200],
            "end": [500, 1600, 500, 1600, 400, 1500, 200, 500, 1600, 1300],
            "strand": ["+", "-", "+", "-", "+", "-", "+", "+", "-", "-"],
        }
    )
    model = GeneAnnotation.from_gtf(gtf)

    five = model.utr5.set_index("gene_id")
    three = model.utr3.set_index("gene_id")
    assert five.loc["GP", "start"] == 100
    assert three.loc["GP", "start"] == 400
    assert five.loc["GM", "start"] == 1500
    assert three.loc["GM", "start"] == 1200


def test_generic_utr_without_cds_is_ignored() -> None:
    base = pd.DataFrame(
        {
            "feature": ["gene", "gene"],
            "gene_id": ["G1", "G2"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [500, 600],
            "strand": ["+", "+"],
        }
    )
    utr_only = pd.concat(
        [
            base,
            pd.DataFrame(
                {
                    "feature": ["utr"],
                    "gene_id": ["G1"],
                    "transcript_id": ["T1"],
                    "chrom": ["chr1"],
                    "start": [100],
                    "end": [150],
                    "strand": ["+"],
                }
            ),
        ],
        ignore_index=True,
    )
    model = GeneAnnotation.from_gtf(utr_only)
    assert model.utr5.empty and model.utr3.empty

    utr_no_tx = pd.concat(
        [
            base,
            pd.DataFrame(
                {
                    "feature": ["cds", "utr"],
                    "gene_id": ["G1", "G1"],
                    "chrom": ["chr1", "chr1"],
                    "start": [200, 100],
                    "end": [400, 150],
                    "strand": ["+", "+"],
                }
            ),
        ],
        ignore_index=True,
    )
    model = GeneAnnotation.from_gtf(utr_no_tx)
    assert model.utr5.empty and model.utr3.empty


def test_split_generic_utr_returns_inputs_when_empty() -> None:
    empty = pd.DataFrame(columns=["feature", "transcript_id"])
    utr5 = pd.DataFrame({"gene_id": ["G1"]})
    utr3 = pd.DataFrame({"gene_id": ["G2"]})
    assert _split_generic_utr(empty, utr5, utr3) == (utr5, utr3)


# ----------------------------------------------------------------------
# GTF normalisation fallbacks
# ----------------------------------------------------------------------
def test_from_gtf_rejects_invalid_inputs() -> None:
    with pytest.raises(TypeError, match="path or a pandas"):
        GeneAnnotation.from_gtf(123)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="'feature' column"):
        GeneAnnotation.from_gtf(pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [1]}))
    with pytest.raises(ValueError, match="no valid feature records"):
        GeneAnnotation.from_gtf(
            pd.DataFrame(
                {
                    "feature": ["gene"],
                    "gene_id": ["G1"],
                    "chrom": ["chr1"],
                    "start": [10],
                    "end": [10],
                    "strand": ["+"],
                }
            )
        )
    with pytest.raises(ValueError, match="no transcript or exon records"):
        GeneAnnotation.from_gtf(
            pd.DataFrame(
                {
                    "feature": ["cds"],
                    "gene_id": ["G1"],
                    "chrom": ["chr1"],
                    "start": [0],
                    "end": [10],
                    "strand": ["+"],
                }
            )
        )


def test_from_gtf_reconstructs_transcripts_from_exons() -> None:
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "gene", "exon", "exon"],
            "gene_id": ["G1", None, "G1", None],
            "gene_name": ["A", None, "A", None],
            "chrom": ["chr1"] * 4,
            "start": [100, 200, 120, 300],
            "end": [500, 600, 180, 350],
            "strand": ["+"] * 4,
        }
    )
    with pytest.warns(RuntimeWarning, match="will be ignored"):
        model = GeneAnnotation.from_gtf(gtf)
    assert list(model.genes["gene_id"]) == ["G1"]
    assert list(model.transcripts["transcript_id"]) == ["G1"]
    assert list(model.exons["transcript_id"]) == ["G1"]


def test_from_gtf_gene_only_and_gene_recovery() -> None:
    gene_only = pd.DataFrame(
        {
            "feature": ["gene", "gene"],
            "gene_id": ["G1", "G2"],
            "gene_name": ["A", "B"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [500, 600],
            "strand": ["+", "-"],
        }
    )
    model = GeneAnnotation.from_gtf(gene_only)
    assert set(model.transcripts["transcript_id"]) == {"G1", "G2"}
    assert model.transcripts["transcript_type"].isna().all()

    exons_only = pd.DataFrame(
        {
            "feature": ["transcript", "exon"],
            "gene_id": ["G1", None],
            "gene_name": ["A", None],
            "transcript_id": ["T1", "T1"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 120],
            "end": [500, 180],
            "strand": ["+", "+"],
        }
    )
    model = GeneAnnotation.from_gtf(exons_only)
    assert list(model.genes["gene_id"]) == ["G1"]
    assert list(model.exons["gene_id"]) == ["G1"]
    assert list(model.transcripts["transcript_id"]) == ["T1"]


def test_from_gtf_drops_transcripts_without_gene_id() -> None:
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "transcript"],
            "gene_id": ["G1", None],
            "gene_name": ["A", None],
            "transcript_id": [None, "T1"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 120],
            "end": [500, 180],
            "strand": ["+", "+"],
        }
    )
    model = GeneAnnotation.from_gtf(gtf)
    assert list(model.genes["gene_id"]) == ["G1"]
    assert model.transcripts.empty


def test_from_gtf_recovers_transcript_gene_from_exons() -> None:
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "transcript", "exon"],
            "gene_id": ["G1", None, "G1"],
            "gene_name": ["A", None, "A"],
            "transcript_id": [None, "T1", "T1"],
            "chrom": ["chr1"] * 3,
            "start": [100, 100, 120],
            "end": [500, 500, 180],
            "strand": ["+", "+", "+"],
        }
    )
    model = GeneAnnotation.from_gtf(gtf)
    transcripts = model.transcripts.set_index("transcript_id")
    assert transcripts.loc["T1", "gene_id"] == "G1"
    assert transcripts.loc["T1", "gene_name"] == "A"


# ----------------------------------------------------------------------
# accessors
# ----------------------------------------------------------------------
def test_gene_annotation_accessors(examples_dir: Path) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    assert model.features("gene") is model.genes
    with pytest.raises(ValueError, match="level must be"):
        model.features("isoform")

    assert model.feature_by_id("transcript", "T1")["gene_id"] == "G1"
    assert model.feature_by_id("transcript", "missing") is None
    assert model.feature_row("transcript", 0)["transcript_id"] == "T1"
    assert model.gene_row("G2")["gene_name"] == "GeneB"
    assert model.gene_row("missing") is None

    row, distance = model.nearest_gene("chr1", 1000)
    assert row is not None and row["gene_id"] == "G1"
    assert distance == 0.0
    row, distance = model.nearest_gene("chrZ", 1000)
    assert row is None and math.isnan(distance)

    row, distance = model.choose_gene(["G1", "G2", "missing"], 1000)
    assert row is not None and row["gene_id"] == "G1"
    assert distance == 0.0
    row, distance = model.choose_gene([], 1000)
    assert row is None and math.isnan(distance)

    row, distance = model.promoter_gene("chr1", 1000, 2000, 2000)
    assert row is not None and row["gene_id"] == "G1"
    assert distance == 0.0
    row, distance = model.promoter_gene("chr1", 50_000, 10, 10)
    assert row is None and math.isnan(distance)
    row, distance = model.promoter_gene("chr1", 1100, 100, 0)
    assert row is None and math.isnan(distance)

    assert model.promoter_label(500, (1000,)) == "Promoter (<=1kb)"


def test_peak_annotation_aliases_and_properties(examples_dir: Path) -> None:
    ann = _example(examples_dir)
    assert ann.upstream == 10
    assert ann.downstream == 10
    assert ann.detail is ann.detail_genomic_annotation
    assert ann.annoStat.equals(ann.anno_stat())
    assert ann.peakNum == len(ann.peaks)
    assert ann.hasGenomicAnnotation is True
    assert "PeakAnnotation(n_peaks=" in repr(ann)

    assert ann.get_anno_stat().equals(ann.anno_stat())
    pd.testing.assert_frame_equal(ann.as_data_frame(), ann.to_chipseeker())
    granges = ann.as_granges()
    assert list(granges.columns[:4]) == ["seqnames", "start", "end", "strand"]


def test_drop_anno_filters_and_validates(examples_dir: Path) -> None:
    ann = _example(examples_dir, upstream=3000, downstream=3000)
    dropped = ann.drop_anno(1000)
    assert len(dropped) < len(ann)
    assert (dropped.annotation["distance_to_tss"].abs() < 1000).all()
    assert len(dropped.peaks) == len(dropped.annotation)

    with pytest.raises(ValueError, match="non-negative"):
        ann.drop_anno(-1)
    bare = replace(ann, annotation=ann.annotation.drop(columns=["distance_to_tss"]))
    with pytest.raises(ValueError, match="distance_to_tss"):
        bare.drop_anno()


def test_summary_groupby_and_errors(examples_dir: Path) -> None:
    ann = _example(examples_dir, upstream=3000, downstream=3000)
    plain = ann.summary()
    assert math.isclose(plain["percent"].sum(), 100.0)
    grouped = ann.summary(groupby="strand")
    assert "Promoter" in grouped.columns
    normalized = ann.summary(groupby="strand", normalize=True)
    assert np.allclose(normalized.sum(axis=1).dropna(), 100.0)

    with pytest.raises(KeyError, match="groupby"):
        ann.summary(groupby="missing_column")
    bare = replace(ann, annotation=ann.annotation.drop(columns=["annotation"]))
    with pytest.raises(ValueError, match="not assigned"):
        bare.summary()


def test_peak_annotation_delegates(examples_dir: Path) -> None:
    ann = _example(examples_dir, upstream=3000, downstream=3000)
    fig, ax = ann.vennpie()
    assert ax is not None
    ann.upsetplot()
    ann.plot_anno_pie()
    ann.plot_anno_bar()
    ann.plot_dist_to_tss()
    tag_matrix = ann.get_tag_matrix()
    ann.plot_peak_prof(tagMatrix=tag_matrix)
    ann.peak_heatmap()
    ann.plot_tss_profile(matrix=np.zeros((3, 4)), window=100, n_bins=4)
    ann.plot_tss_heatmap(matrix=np.zeros((3, 4)), window=100, n_bins=4)
    plt.close("all")


def test_tss_matrix_sources() -> None:
    model = GeneAnnotation.from_gtf(
        pd.DataFrame(
            {
                "feature": ["gene"],
                "gene_id": ["G1"],
                "gene_name": ["A"],
                "chrom": ["1"],
                "start": [10_000],
                "end": [20_000],
                "strand": ["+"],
            }
        )
    )
    ann = PeakAnnotation(
        peaks=pd.DataFrame(columns=["peak_id"]),
        annotation=pd.DataFrame(columns=["peak_id"]),
        gene_annotation=model,
    )
    matrix = np.zeros((2, 3))
    np.testing.assert_array_equal(ann._tss_matrix(matrix=matrix), matrix)

    with pytest.raises(ValueError, match="exactly one"):
        ann._tss_matrix()
    with pytest.raises(ValueError, match="exactly one"):
        ann._tss_matrix(matrix=matrix, signal="track.bedgraph")
    with pytest.raises(ValueError, match="2-dimensional"):
        ann._tss_matrix(matrix=np.zeros(3))

    pytest.importorskip("pyBigWig")
    bigwig = Path(__file__).parent / "data" / "test.bw"
    if not bigwig.exists():  # pragma: no cover - packaging guard
        pytest.skip("tests/data/test.bw is missing")
    tss = ann._tss_matrix(bigwig=bigwig, window=100, n_bins=4)
    assert tss.shape == (1, 4)

    signal = pd.DataFrame(
        {"chrom": ["1"], "start": [9_000], "end": [21_000], "value": [1.0]}
    )
    from_signal = ann._tss_matrix(signal=signal, window=100, n_bins=4)
    assert from_signal.shape == (1, 4)


# ----------------------------------------------------------------------
# peak table normalisation
# ----------------------------------------------------------------------
def test_normalise_peaks_ids_and_renames() -> None:
    with pytest.raises(TypeError, match="path or a pandas"):
        _normalise_peaks(123)  # type: ignore[arg-type]

    positional = _normalise_peaks(
        pd.DataFrame({"a": ["chr1"], "b": [0], "c": [10], "foo": ["bar"]})
    )
    assert positional.loc[0, "peak_id"] == "peak_1"
    assert positional.loc[0, "foo"] == "bar"

    both_chrom = _normalise_peaks(
        pd.DataFrame({"chr": ["chr1"], "chrom": ["chr1"], "start": [0], "end": [10]})
    )
    assert both_chrom.loc[0, "chrom"] == "chr1"

    with pytest.raises(ValueError, match="missing columns"):
        _normalise_peaks(pd.DataFrame({"chrom": ["chr1"]}))
    with pytest.raises(ValueError, match="non-numeric"):
        _normalise_peaks(
            pd.DataFrame({"chrom": ["chr1"], "start": ["x"], "end": [10]})
        )

    peak_id_only = _normalise_peaks(
        pd.DataFrame(
            {"chrom": ["chr1"], "start": [0], "end": [10], "peak_id": ["p1"]}
        )
    )
    assert peak_id_only.loc[0, "name"] == "p1"

    missing_names = _normalise_peaks(
        pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [0, 10], "end": [10, 20], "name": [np.nan, "."]}
        )
    )
    assert list(missing_names["peak_id"]) == ["peak_1", "peak_2"]

    invalid_names = _normalise_peaks(
        pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [0, 10], "end": [10, 20], "name": ["", "p2"]}
        )
    )
    assert list(invalid_names["name"]) == ["peak_1", "p2"]

    duplicate_names = _normalise_peaks(
        pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [0, 10], "end": [10, 20], "name": ["p", "p"]}
        )
    )
    assert list(duplicate_names["peak_id"]) == ["p", "p_2"]

    custom_index = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [10], "name": [np.nan]},
        index=pd.Index(["row-a"]),
    )
    from_index = _normalise_peaks(custom_index)
    assert from_index.loc[0, "peak_id"] == "row-a"

    no_name_index = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [10]},
        index=pd.Index(["anonymous"]),
    )
    assert _normalise_peaks(no_name_index).loc[0, "peak_id"] == "anonymous"


# ----------------------------------------------------------------------
# annotate_peak validation / short options
# ----------------------------------------------------------------------
def test_annotate_peak_argument_validation(examples_dir: Path) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [1000], "end": [1100]})

    with pytest.raises(ValueError, match="provide peaks"):
        annotate_peak(gtf=model)
    with pytest.raises(ValueError, match="only one of"):
        annotate_peak(peaks, gtf=model, annotation=model)
    with pytest.raises(TypeError, match="unexpected keyword"):
        annotate_peak(peaks, gtf=model, bogus=1)
    with pytest.raises(ValueError, match="level must be"):
        annotate_peak(peaks, gtf=model, level="isoform")
    with pytest.raises(ValueError, match="provide gtf"):
        annotate_peak(peaks)
    with pytest.raises(ValueError, match="overlap must be"):
        annotate_peak(peaks, gtf=model, overlap="bad")
    with pytest.raises(ValueError, match="tss_region must be"):
        annotate_peak(peaks, gtf=model, tss_region=(10, -10))
    with pytest.raises(ValueError, match="promoter_bins"):
        annotate_peak(peaks, gtf=model, promoter_bins=[2000, 1000])
    with pytest.raises(ValueError, match="non-negative"):
        annotate_peak(peaks, gtf=model, downstream_distance=-1)


def test_annotate_peak_resolves_annotation_shortcuts(
    examples_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pychipseeker import genomes

    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [1000], "end": [1100]})

    monkeypatch.setattr(genomes, "detect_assembly", lambda peaks=None: model)
    assert len(annotate_peak(peaks, gtf="auto").annotation) == 1

    monkeypatch.setattr(genomes, "is_assembly_name", lambda name: name == "toy")
    monkeypatch.setattr(genomes, "load_genome", lambda name: model)
    assert len(annotate_peak(peaks, gtf="toy").annotation) == 1

    monkeypatch.setattr(genomes, "is_genome_spec", lambda spec: spec == "toy.gtf")
    monkeypatch.setattr(genomes, "load_annotation", lambda spec: model)
    assert len(annotate_peak(peaks, gtf="toy.gtf").annotation) == 1


def test_annotate_peak_shortcut_options(examples_dir: Path) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    peaks = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [1000],
            "end": [1100],
            "extra": ["kept"],
        }
    )
    ann = annotate_peak(
        peaks,
        annotation=model,
        genomic_annotation_priority=None,
    )
    assert ann.annotation.loc[0, "gene_id"] == "G1"
    with_metadata = ann.to_dataframe(include_peak_metadata=True)
    assert with_metadata.loc[0, "extra"] == "kept"

    unassigned = annotate_peak(peaks, gtf=model, assign_genomic_annotation=False)
    assert "annotation" not in unassigned.annotation.columns
    assert unassigned.annotation.loc[0, "gene_id"] == "G1"

    db = pd.DataFrame({"gene_id": ["G1"], "SYMBOL": ["A"]})
    enriched = annotate_peak(peaks, gtf=model, annoDb=db)
    assert enriched.annotation.loc[0, "SYMBOL"] == "A"


def test_annotate_peak_ignore_neighbour_flags(examples_dir: Path) -> None:
    upstream = _example(examples_dir, ignore_upstream=True)
    downstream = _example(examples_dir, ignore_downstream=True)
    assert 0 < len(upstream.annotation) <= 10
    assert 0 < len(downstream.annotation) <= 10


def test_annotate_peak_empty_and_all_dropped(examples_dir: Path) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    empty = annotate_peak(
        pd.DataFrame(columns=["chrom", "start", "end"]), gtf=model
    )
    assert empty.annotation.empty
    assert empty.detail_genomic_annotation.empty

    unknown = pd.DataFrame({"chrom": ["chrZ"], "start": [10], "end": [20]})
    with pytest.warns(RuntimeWarning, match="dropped"):
        with pytest.raises(ValueError, match="all .* peaks were dropped"):
            annotate_peak(unknown, gtf=model)


def test_annotate_peak_adds_flanking_info(examples_dir: Path) -> None:
    ann = _example(examples_dir, add_flank_gene_info=True, flank_distance=1000)
    assert ann.add_flank_gene_info
    assert "flank_gene_ids" in ann.annotation.columns
    assert ann.annotation["flank_gene_ids"].notna().any()

    gene_level = _example(
        examples_dir, level="gene", add_flank_gene_info=True, flank_distance=1000
    )
    assert "flank_tx_ids" not in gene_level.annotation.columns


# ----------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------
def test_genomic_anno_stat_ignore_flags() -> None:
    labels = [
        "Promoter (<=1kb)",
        "Promoter (1-2kb)",
        "Exon (T1/G1, exon 1 of 3)",
        "Exon (T1/G1, exon 2 of 3)",
        "Intron (T1/G1, intron 1 of 2)",
        "Intron (T1/G1, intron 2 of 2)",
        "Downstream (<=300bp)",
        "Distal Intergenic",
    ]
    stat = get_genomic_anno_stat(
        labels,
        downstream_distance=3000,
        ignore_1st_exon=True,
        ignore_1st_intron=True,
        ignore_downstream=True,
        ignore_promoter_subcategory=True,
    ).set_index("Feature")
    assert stat.loc["Promoter", "Count"] == 2
    assert stat.loc["Exon", "Count"] == 2
    assert stat.loc["Intron", "Count"] == 2
    assert stat.loc["Intergenic", "Count"] == 2
    assert "1st Exon" not in stat.index
    assert "1st Intron" not in stat.index
    assert not any(index.startswith("Downstream") for index in stat.index)

    detailed = get_genomic_anno_stat(labels, downstream_distance=3000).set_index(
        "Feature"
    )
    assert detailed.loc["Downstream (<=3kb)", "Count"] == 1
    assert detailed.loc["Distal Intergenic", "Count"] == 1


# ----------------------------------------------------------------------
# seq2gene helpers
# ----------------------------------------------------------------------
def test_seq2gene_empty_and_validation(examples_dir: Path) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    far = pd.DataFrame({"chrom": ["chr1"], "start": [50_000], "end": [50_100]})
    assert seq2gene(far, model, tss_region=(0, 0), flank_distance=1) == []
    links = seq2gene_links(far, model, tss_region=(0, 0), flank_distance=1)
    assert links.empty
    assert list(links.columns) == [
        "peak_id",
        "chrom",
        "start",
        "end",
        "gene_id",
        "gene_name",
        "distance_to_tss",
        "link_type",
    ]

    with pytest.raises(ValueError, match="provide gene_annotation"):
        seq2gene_links(far)
    with pytest.raises(ValueError, match="level must be"):
        seq2gene_links(far, model, level="isoform")

    from_gtf = seq2gene_links(
        far,
        gtf=examples_dir / "example.gtf",
        upstream=10,
        downstream=10,
    )
    assert from_gtf.empty


def test_seq2gene_links_host_replacements(monkeypatch: pytest.MonkeyPatch) -> None:
    gtf = pd.DataFrame(
        {
            "feature": ["gene", "gene", "exon", "exon"],
            "gene_id": ["G1", "G2", "G1", "G2"],
            "gene_name": ["A", "B", "A", "B"],
            "transcript_id": [None, None, "T1", "T2"],
            "chrom": ["chr1"] * 4,
            "start": [1000, 3000, 1000, 3000],
            "end": [2000, 4000, 2000, 4000],
            "strand": ["+", "+", "+", "+"],
        }
    )
    model = GeneAnnotation.from_gtf(gtf)
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [1900], "end": [3100]})

    links = seq2gene_links(peaks, model)
    host = links[links["link_type"] == "host"].set_index("gene_id")
    assert set(host.index) == {"G1", "G2"}
    assert host.loc["G1", "distance_to_tss"] == 900.0
    assert host.loc["G2", "distance_to_tss"] == 0.0

    missing_gene = lambda gene_id: None  # noqa: E731 - tiny test double
    monkeypatch.setattr(model, "gene_row", missing_gene)
    skipped = seq2gene_links(peaks, model)
    assert "host" not in set(skipped["link_type"])


def test_nearest_feature_with_empty_model() -> None:
    columns = [
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
    genes = pd.DataFrame(columns=columns)
    transcripts = pd.DataFrame(
        columns=[
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
    )
    model = GeneAnnotation(
        genes, transcripts, genes.copy(), genes.copy(), genes.copy(), genes.copy()
    )
    feature_idx, distance = _nearest_feature(model, "gene", "chr1", 0, 10)
    assert feature_idx is None
    assert math.isnan(distance)


def test_downstream_annotation_labels() -> None:
    model = GeneAnnotation.from_gtf(
        pd.DataFrame(
            {
                "feature": ["gene"],
                "gene_id": ["G1"],
                "chrom": ["chr1"],
                "start": [1000],
                "end": [2000],
                "strand": ["+"],
            }
        )
    )
    candidate = lambda start: _downstream_annotation(  # noqa: E731
        model, "transcript", "chr1", start, start + 100, "*", False, 300
    )
    assert candidate(2099) == ("Downstream (<=300bp)", 100)
    assert candidate(10_000) == (None, 8001)

    long_window = lambda start, dsd: _downstream_annotation(  # noqa: E731
        model, "transcript", "chr1", start, start + 100, "*", False, dsd
    )
    assert long_window(3499, 3000) == ("Downstream (1-2kb)", 1500)
    assert long_window(2500, 3000) == ("Downstream (<1kb)", 501)
    assert long_window(3399, 1500) == ("Downstream (1.5kb-1500bp)", 1400)


def test_get_flanking_genes_branches() -> None:
    model = GeneAnnotation.from_gtf(
        pd.DataFrame(
            {
                "feature": ["gene"],
                "gene_id": ["G1"],
                "gene_name": ["A"],
                "chrom": ["chr1"],
                "start": [1000],
                "end": [2000],
                "strand": ["-"],
            }
        )
    )
    with pytest.raises(ValueError, match="level must be"):
        get_flanking_genes(model, pd.DataFrame(), level="isoform")

    far = pd.DataFrame({"chrom": ["chr1"], "start": [5000], "end": [5100]})
    empty = get_flanking_genes(model, far, level="gene", distance=1)
    assert empty.empty
    assert list(empty.columns) == [
        "peak_id",
        "flank_gene_ids",
        "flank_gene_distances",
    ]

    nearby = pd.DataFrame({"chrom": ["chr1"], "start": [2050], "end": [2100]})
    table = get_flanking_genes(model, nearby, level="gene", distance=200)
    assert table.loc[0, "flank_gene_ids"] == "G1"
    assert table.loc[0, "flank_gene_distances"] == "-51"

    plus_model = GeneAnnotation.from_gtf(
        pd.DataFrame(
            {
                "feature": ["gene"],
                "gene_id": ["G1"],
                "gene_name": ["A"],
                "chrom": ["chr1"],
                "start": [1000],
                "end": [2000],
                "strand": ["+"],
            }
        )
    )
    plus_table = get_flanking_genes(plus_model, nearby, level="gene", distance=200)
    assert plus_table.loc[0, "flank_gene_distances"] == "1050"


# ----------------------------------------------------------------------
# parallel annotation plumbing
# ----------------------------------------------------------------------
def test_worker_helpers(examples_dir: Path, tmp_path: Path) -> None:
    from pychipseeker.genomes import save_bundle

    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    config = _config()
    _annotation_worker_init(model, config)
    assert ann_mod._WORKER_MODEL is model
    assert _annotation_worker(pd.DataFrame()) == ([], [], 0)

    source = ("gtf", str(examples_dir / "example.gtf"))
    _annotation_worker_init_from_source(source, config)
    assert ann_mod._WORKER_MODEL is not None
    assert _annotation_worker(pd.DataFrame()) == ([], [], 0)

    loaded = _load_annotation_model(source)
    assert {"G1", "G2"}.issubset(set(loaded.genes["gene_id"]))
    with pytest.raises(ValueError, match="unknown annotation source"):
        _load_annotation_model(("bogus", "x"))

    bundle = tmp_path / "bundle.npz"
    save_bundle(model, bundle)
    npz_model = _load_annotation_model(("npz", str(bundle)))
    assert npz_model.source == ("npz", str(bundle))


def test_annotate_peaks_frame_parallel_paths(
    examples_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    rng = np.random.default_rng(3)
    n_peaks = 2100
    starts = rng.integers(500, 9000, n_peaks)
    peaks = _normalise_peaks(
        pd.DataFrame(
            {
                "chrom": np.repeat("chr1", n_peaks),
                "start": starts,
                "end": starts + rng.integers(50, 300, n_peaks),
            }
        )
    )

    captured: dict[str, int] = {}

    def fake_parallel_map(func, items, *, n_jobs, initializer=None, initargs=()):
        captured["workers"] = n_jobs
        captured["chunks"] = len(items)
        if initializer is not None:
            initializer(*initargs)
        return [func(item) for item in items]

    monkeypatch.setattr(ann_mod, "parallel_map", fake_parallel_map)
    monkeypatch.setattr(
        ann_mod,
        "effective_jobs",
        lambda n_jobs, tasks=None: 4 if n_jobs == -1 else 1,
    )

    model.source = ("gtf", str(examples_dir / "example.gtf"))
    config = _config(total=n_peaks, verbose=True)
    records, details, dropped = _annotate_peaks_frame(
        model, peaks, config, n_jobs=-1
    )
    assert captured["workers"] == 4
    assert captured["chunks"] > 1
    assert len(records) + dropped == n_peaks
    assert len(details) == len(records)

    model.source = None
    records, details, dropped = _annotate_peaks_frame(
        model, peaks, _config(total=n_peaks), n_jobs=-1
    )
    assert len(records) + dropped == n_peaks


def test_annotate_rows_defaults_total_and_verbose(
    examples_dir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    rng = np.random.default_rng(5)
    n_peaks = 10_000
    starts = rng.integers(500, 9000, n_peaks)
    peaks = _normalise_peaks(
        pd.DataFrame(
            {
                "chrom": np.repeat("chr1", n_peaks),
                "start": starts,
                "end": starts + rng.integers(50, 300, n_peaks),
            }
        )
    )
    config = _config(verbose=True)
    config.pop("total")
    records, details, dropped = ann_mod._annotate_rows(model, peaks, **config)
    assert len(records) + dropped == n_peaks
    assert len(details) == len(records)
    assert "annotated 10,000/10,000 peaks" in capsys.readouterr().out
