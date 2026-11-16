"""Branch coverage for :mod:`pychipseeker.analysis` edge cases."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pychipseeker import TagMatrix
from pychipseeker.analysis import (
    _human_bp,
    annotation_summary,
    distance_summary,
    gene_type_summary,
    peaks_per_gene,
    signal_groups,
    tss_enrichment_by_quartile,
)


def _frame(**columns) -> pd.DataFrame:
    base = {
        "chrom": ["chr1"] * 4,
        "start": [0, 100, 200, 300],
        "end": [50, 150, 250, 350],
        "annotation": ["Promoter"] * 4,
    }
    base.update(columns)
    return pd.DataFrame(base)


def test_annotation_frame_validation() -> None:
    with pytest.raises(ValueError, match="missing columns"):
        annotation_summary(pd.DataFrame({"chrom": ["chr1"], "start": [0]}))

    with pytest.raises(ValueError, match="no valid peaks"):
        annotation_summary(
            pd.DataFrame(
                {
                    "chrom": ["chr1"],
                    "start": [100],
                    "end": [100],
                    "annotation": ["Promoter"],
                }
            )
        )


def test_peaks_per_gene_gene_name_and_errors() -> None:
    per_gene = peaks_per_gene(_frame(gene_name=["A", "A", "B", "B"]))
    assert set(per_gene["gene_id"]) == {"A", "B"}
    assert list(per_gene["gene_name"]) == ["A", "B"]

    from_id = peaks_per_gene(_frame(gene_id=["G1", "G1", "G2", "G2"]))
    assert set(from_id["gene_name"]) == {"G1", "G2"}

    with pytest.raises(ValueError, match="neither gene_id nor gene_name"):
        peaks_per_gene(_frame())


def test_gene_type_summary_requires_gene_type() -> None:
    with pytest.raises(ValueError, match="no gene_type column"):
        gene_type_summary(_frame())


def test_distance_summary_fallbacks_and_labels() -> None:
    nearest = distance_summary(
        _frame(nearest_distance=[0, 1_000_000, 2_500, np.nan]),
        thresholds=(1_000_000, 2_500),
    )
    metrics = dict(zip(nearest["Metric"], nearest["Value"]))
    assert metrics["Within 1 Mb"] >= 0.0
    assert metrics["Within 2500 bp"] >= 0.0

    with pytest.raises(ValueError, match="no distance_to_tss"):
        distance_summary(_frame())

    with pytest.raises(ValueError, match="no finite distance"):
        distance_summary(_frame(distance_to_tss=[np.nan, np.nan, np.nan, np.nan]))

    assert _human_bp(2_500) == "2500 bp"
    assert _human_bp(2_000_000) == "2 Mb"


def _tag(matrix: np.ndarray, rows: pd.DataFrame | None = None) -> TagMatrix:
    return TagMatrix(
        matrix,
        upstream=1000 if rows is not None else None,
        downstream=1000 if rows is not None else None,
        rows=rows,
    )


def test_signal_groups_validation() -> None:
    frame = _frame()
    with pytest.raises(ValueError, match="at least 2"):
        signal_groups(np.zeros((4, 3)), frame, n_groups=1)

    with pytest.raises(ValueError, match="2-dimensional"):
        signal_groups(np.zeros(4), frame)

    with pytest.raises(ValueError, match="no rows"):
        signal_groups(np.zeros((0, 3)), frame)

    with pytest.raises(ValueError, match="no finite signal"):
        signal_groups(np.full((3, 4), np.nan), frame)


def test_signal_groups_falls_back_when_rows_do_not_match() -> None:
    matrix = np.ones((4, 6), dtype=float)
    matrix[:, 2:4] = 5.0
    rows = pd.DataFrame({"gene_id": ["missing"] * 4})
    groups = signal_groups(
        _tag(matrix, rows),
        _frame(gene_id=["G1"] * 4, score=[1.0, 2.0, 3.0, 4.0]),
        n_groups=2,
    )
    assert set(groups["group"]) == {0, 1}
    assert np.isfinite(groups["score"]).all()


def test_tss_enrichment_validation_and_center_fallback() -> None:
    frame = _frame(gene_id=["G1"] * 4)
    matrix = np.ones((4, 10), dtype=float)

    with pytest.raises(ValueError, match="center_fraction"):
        tss_enrichment_by_quartile(matrix, frame, center_fraction=0.0)
    with pytest.raises(ValueError, match="flank_fraction"):
        tss_enrichment_by_quartile(matrix, frame, flank_fraction=1.5)

    # A plain ndarray has no upstream/downstream, so the centre falls back to
    # the matrix midpoint.
    summary = tss_enrichment_by_quartile(matrix, frame, n_groups=2)
    assert len(summary) == 2
    assert np.allclose(summary["tss_enrichment"], 1.0)
