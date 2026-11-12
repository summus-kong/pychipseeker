"""Tests for the high-level analysis summaries."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pychipseeker import (
    annotate_peak,
    annotation_summary,
    distance_summary,
    gene_type_summary,
    peaks_per_gene,
    signal_groups,
    top_genes,
    tss_enrichment_by_quartile,
)


@pytest.fixture
def annotation(examples_dir):
    return annotate_peak(
        examples_dir / "example_peaks.bed",
        examples_dir / "example.gtf",
        upstream=1000,
        downstream=1000,
    )


def test_annotation_summary_counts_and_percent(annotation):
    summary = annotation_summary(annotation)
    assert list(summary.columns)[:4] == [
        "annotation",
        "count",
        "percent",
        "median_width_bp",
    ]
    assert int(summary["count"].sum()) == len(annotation)
    assert np.isclose(summary["percent"].sum(), 100.0)
    assert (summary["count"].to_numpy()[:-1] >= summary["count"].to_numpy()[1:]).all()


def test_peaks_per_gene_and_top_genes(annotation):
    per_gene = peaks_per_gene(annotation)
    assert int(per_gene["n_peaks"].sum()) == len(annotation)
    assert per_gene.iloc[0]["gene_id"] == "G1"
    assert int(per_gene.iloc[0]["n_peaks"]) == 5

    top = top_genes(annotation, n=2)
    assert list(top["rank"]) == [1, 2]
    assert len(top) == 2
    assert top.iloc[0]["gene_id"] == "G1"


def test_gene_type_summary_folds_tail(annotation):
    summary = gene_type_summary(annotation, top=1)
    assert int(summary["count"].sum()) == len(annotation)
    assert np.isclose(summary["percent"].sum(), 100.0)
    assert summary.iloc[0]["gene_type"] == "protein_coding"

    frame = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start": [0, 100, 200, 300],
            "end": [50, 150, 250, 350],
            "annotation": ["Promoter"] * 4,
            "gene_type": ["protein_coding", "protein_coding", "lincRNA", "miRNA"],
        }
    )
    folded = gene_type_summary(frame, top=1)
    assert list(folded["gene_type"]) == ["protein_coding", "Other"]
    assert int(folded["count"].sum()) == 4
    assert int(folded.iloc[-1]["count"]) == 2


def test_distance_summary_metrics(annotation):
    summary = distance_summary(annotation, thresholds=(1000,))
    assert list(summary.columns) == ["Metric", "Value"]
    metrics = dict(zip(summary["Metric"], summary["Value"]))
    assert metrics["Peaks"] == len(annotation)
    assert 0.0 <= metrics["Within 1 kb"] <= 100.0
    assert metrics["Median distance (bp)"] >= 0.0


def _synthetic_tag(n_rows: int = 8, n_bins: int = 20):
    from pychipseeker import TagMatrix

    matrix = np.zeros((n_rows, n_bins), dtype=float)
    matrix[:, 9:11] = 10.0
    matrix[:, :2] = 1.0
    matrix[:, -2:] = 1.0
    rows = pd.DataFrame({"gene_id": [f"G{i}" for i in range(n_rows)]})
    frame = pd.DataFrame(
        {
            "chrom": ["chr1"] * n_rows,
            "start": np.arange(n_rows) * 100,
            "end": np.arange(n_rows) * 100 + 50,
            "annotation": ["Promoter"] * n_rows,
            "gene_id": [f"G{i}" for i in range(n_rows)],
            "gene_name": [f"Gene{i}" for i in range(n_rows)],
            "gene_type": ["protein_coding"] * n_rows,
            "score": np.arange(1, n_rows + 1, dtype=float),
        }
    )
    tag = TagMatrix(matrix, upstream=1000, downstream=1000, rows=rows)
    return tag, frame


def test_signal_groups_are_balanced():
    tag, frame = _synthetic_tag(8, 20)
    groups = signal_groups(tag, frame, n_groups=4)
    assert len(groups) == np.asarray(tag).shape[0]
    assert set(groups["group"]) == {0, 1, 2, 3}
    counts = groups["group"].value_counts().to_numpy()
    assert counts.max() - counts.min() <= 1
    assert groups.loc[groups["group"] == 0, "label"].iloc[0].endswith("lowest")
    assert groups.loc[groups["group"] == 3, "label"].iloc[0].endswith("highest")


def test_signal_groups_caps_at_row_count():
    tag, frame = _synthetic_tag(3, 20)
    groups = signal_groups(tag, frame, n_groups=4)
    assert set(groups["group"]) == {0, 1, 2}


def test_tss_enrichment_by_quartile():
    tag, frame = _synthetic_tag(8, 20)
    summary = tss_enrichment_by_quartile(tag, frame, n_groups=4)
    assert len(summary) == 4
    assert int(summary["n"].sum()) == np.asarray(tag).shape[0]
    assert np.isfinite(summary["tss_enrichment"]).all()
    assert np.allclose(summary["tss_enrichment"], 10.0)
    assert np.allclose(summary["max_coverage"], 10.0)


def test_analysis_argument_validation(annotation):
    with pytest.raises(TypeError):
        annotation_summary("not-an-annotation")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        top_genes(annotation, n=0)
    with pytest.raises(ValueError):
        gene_type_summary(annotation, top=0)


def test_signal_groups_falls_back_to_coverage():
    matrix = np.arange(12, dtype=float).reshape(4, 3)
    frame = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start": [0, 100, 200, 300],
            "end": [50, 150, 250, 350],
            "annotation": ["Promoter"] * 4,
        }
    )
    groups = signal_groups(matrix, frame, n_groups=2)
    assert len(groups) == 4
    assert set(groups["group"]) == {0, 1}
