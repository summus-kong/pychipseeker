"""Packaging metadata and public API surface consistency."""

from __future__ import annotations

import importlib.metadata
from pathlib import Path

import pytest

import pychipseeker as pcs

_REPO_ROOT = Path(__file__).resolve().parents[1]


def test_dunder_all_resolves_and_is_unique() -> None:
    missing = [name for name in pcs.__all__ if not hasattr(pcs, name)]
    assert missing == []
    assert len(pcs.__all__) == len(set(pcs.__all__))


def test_installed_metadata_matches_module_version() -> None:
    try:
        installed = importlib.metadata.version("pychipseeker")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        pytest.skip("pychipseeker is not installed in this environment")
    assert installed == pcs.__version__, (
        "installed metadata and pychipseeker.__version__ disagree; reinstall "
        "with `pip install -e .` to refresh the distribution metadata"
    )


def test_installed_metadata_uses_spdx_license() -> None:
    try:
        meta = importlib.metadata.metadata("pychipseeker")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        pytest.skip("pychipseeker is not installed in this environment")
    assert meta.get("License-Expression") == "MIT", (
        "installed metadata does not use the SPDX license expression; reinstall "
        "with `pip install -e .` to refresh the distribution metadata"
    )
    assert "License :: OSI Approved :: MIT License" not in meta.get_all(
        "Classifier", []
    )


def test_changelog_documents_current_version() -> None:
    text = (_REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{pcs.__version__}]" in text


@pytest.mark.parametrize(
    ("public", "alias"),
    [
        ("annotate_peak", "annotatePeak"),
        ("read_peak_file", "readPeakFile"),
        ("as_granges", "asGRanges"),
        ("combine_cs_anno", "combine_csAnno"),
        ("get_anno_stat", "getAnnoStat"),
        ("detail_genomic_annotation", "detailGenomicAnnotation"),
        ("get_promoters", "getPromoters"),
        ("get_bio_region", "getBioRegion"),
        ("get_tag_matrix", "getTagMatrix"),
        ("make_bio_region_from_granges", "makeBioRegionFromGranges"),
        ("plot_avg_prof", "plotAvgProf"),
        ("plot_peak_prof2", "plotAvgProf2"),
        ("plot_peak_prof", "plotPeakProf"),
        ("plot_peak_prof2", "plotPeakProf2"),
        ("plot_peak_prof_multi_windows", "plotPeakProf_MultiWindows"),
        ("plot_multi_prof", "plotMultiProf"),
        ("tag_heatmap", "tagHeatmap"),
        ("peak_heatmap", "peakHeatmap"),
        ("peak_heatmap_multiple_sets", "peakHeatmap_multiple_Sets"),
        ("peak_profile_heatmap", "peak_Profile_Heatmap"),
        ("plot_anno_pie", "plotAnnoPie"),
        ("plot_anno_bar", "plotAnnoBar"),
        ("plot_dist_to_tss", "plotDistToTSS"),
        ("plot_tss_profile", "plotProfile"),
        ("plot_tss_heatmap", "plotHeatmap"),
        ("get_geo_species", "getGEOspecies"),
        ("get_geo_genome_version", "getGEOgenomeVersion"),
        ("get_geo_info", "getGEOInfo"),
        ("download_geo_bed_files", "downloadGEObedFiles"),
        ("download_gsm_bed_files", "downloadGSMbedFiles"),
        ("get_sample_files", "getSampleFiles"),
        ("enrich_peak_overlap", "enrichPeakOverlap"),
        ("enrich_anno_overlap", "enrichAnnoOverlap"),
    ],
)
def test_chipseeker_aliases_share_objects(public: str, alias: str) -> None:
    assert getattr(pcs, public) is getattr(pcs, alias)
