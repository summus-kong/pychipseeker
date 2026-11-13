"""Offline enrichment helpers and plotting normalisation edge cases."""

from __future__ import annotations

import math
import sys
import types

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pychipseeker import enrichment
from pychipseeker.enrichment import (
    _as_gene_sets,
    _as_ranked_frame,
    _enrichment_plot_frame,
    _hypergeom_sf,
    _normal_quantile,
    _p_adjust,
    enrich_go,
    enrich_kegg,
    enrich_peak_genes,
    enrichr,
    gsea,
    gsea_peak_genes,
    ora,
    plot_enrichment,
    read_gmt,
)


def test_p_adjust_methods_and_empty() -> None:
    assert _p_adjust([]) == []
    assert _p_adjust([0.02, 0.5], "none") == [0.02, 0.5]
    assert _p_adjust([0.02, 0.5], "raw") == [0.02, 0.5]
    np.testing.assert_allclose(
        _p_adjust([0.02, 0.5], "bonferroni"), [0.04, 1.0]
    )
    with pytest.raises(ValueError, match="p_adjust"):
        _p_adjust([0.1], "holm")


def test_hypergeom_sf_bounds() -> None:
    assert _hypergeom_sf(0, 100, 10, 20) == 1.0
    assert _hypergeom_sf(11, 100, 10, 20) == 0.0
    assert 0.0 < _hypergeom_sf(3, 100, 10, 20) < 1.0


def test_read_gmt_skips_malformed_lines(tmp_path) -> None:
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "\n"
        "only_one_field\n"
        "\tdesc\tg1\n"
        "term A\tdesc\tg1\t\tg2\n"
        "empty\tdesc\n",
        encoding="utf-8",
    )
    assert read_gmt(gmt) == {"term A": ["g1", "g2"]}


def test_as_gene_sets_accepts_mapping_and_rejects_other(tmp_path) -> None:
    path = tmp_path / "sets.gmt"
    path.write_text("term\tdesc\ta\tb\n", encoding="utf-8")
    assert _as_gene_sets(path) == {"term": {"a", "b"}}
    assert _as_gene_sets({"t": ["a", None, "b"]}) == {"t": {"a", "b"}}
    with pytest.raises(TypeError, match="gene_sets must be"):
        _as_gene_sets([["a"]])  # type: ignore[arg-type]


def test_ora_defaults_background_and_rejects_no_overlap() -> None:
    sets = {"A": ["g1", "g2"], "B": ["g2", "g3"]}
    result = ora(["g1", "g2"], sets)
    assert result.iloc[0]["Term"] == "A"
    assert result.iloc[0]["BackgroundSize"] == 3

    with pytest.raises(ValueError, match="present in the background"):
        ora(["gX"], sets)
    with pytest.raises(ValueError, match="no gene set"):
        ora(["gX"], sets, background=["gX"])
    with pytest.raises(ValueError, match="p_adjust"):
        ora(["g1"], sets, p_adjust="holm")


def test_ora_haldane_correction_keeps_odds_finite() -> None:
    result = ora(
        ["a", "b"],
        {"all": ["a", "b"]},
        background=["a", "b", "c"],
    )
    row = result.iloc[0]
    assert math.isfinite(row["OddsRatio"])
    assert row["ORLower"] <= row["OddsRatio"] <= row["ORUpper"]


def test_normal_quantile_tails_and_validation() -> None:
    assert _normal_quantile(0.5) == pytest.approx(0.0, abs=1e-6)
    assert _normal_quantile(0.01) == pytest.approx(-2.326347874, rel=1e-5)
    assert _normal_quantile(0.99) == pytest.approx(2.326347874, rel=1e-5)
    with pytest.raises(ValueError, match="between 0 and 1"):
        _normal_quantile(0.0)
    with pytest.raises(ValueError, match="between 0 and 1"):
        _normal_quantile(1.0)


def _install_fake_gseapy(monkeypatch, results: pd.DataFrame, calls: dict | None = None):
    class EnrichrResult:
        def __init__(self, table):
            self.results = table

    def fake_enrichr(**kwargs):
        if calls is not None:
            calls.update(kwargs)
        return EnrichrResult(results.copy())

    fake = types.ModuleType("gseapy")
    fake.enrichr = fake_enrichr
    monkeypatch.setitem(sys.modules, "gseapy", fake)
    return fake


def test_enrichr_background_and_recomputed_adjustment(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    calls: dict = {}
    table = pd.DataFrame(
        {"Term": ["A", "B"], "P-value": [0.01, 0.2], "Adjusted P-value": [0.02, 0.2]}
    )
    _install_fake_gseapy(monkeypatch, table, calls)

    out = enrichr(
        ["g1", "g2"],
        background=["g1", "g2", None, "g1"],
        p_adjust="bonferroni",
        outdir=tmp_path,
        cutoff=0.5,
    )
    assert calls["background"] == ["g1", "g2"]
    assert calls["outdir"] == str(tmp_path)
    assert calls["cutoff"] == 0.5
    np.testing.assert_allclose(out["Adjusted P-value"], [0.02, 0.4])

    with pytest.raises(ValueError, match="background is empty"):
        enrichr(["g1"], background=[None, ""])


def test_enrichr_leaves_empty_results_untouched(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_gseapy(monkeypatch, pd.DataFrame(columns=["Term", "P-value"]))
    out = enrichr(["g1"], p_adjust="BH")
    assert out.empty
    assert "Adjusted P-value" not in out.columns


def test_enrich_go_and_kegg_forward_libraries(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []

    def fake_enrichr(gene_list, **kwargs):
        calls.append((gene_list, kwargs))
        return pd.DataFrame()

    monkeypatch.setattr(enrichment, "enrichr", fake_enrichr)
    enrich_go(["g1"], organism="human")
    enrich_kegg(["g2"], organism="mouse")
    assert calls[0][1]["gene_sets"] == "GO_Biological_Process_2023"
    assert calls[1][1]["gene_sets"] == "KEGG_2019_Mouse"
    enrich_kegg(["g3"], organism="human")
    assert calls[2][1]["gene_sets"] == "KEGG_2021_Human"


class FakeLinksAnnotation:
    def __init__(self, links: pd.DataFrame):
        self.links = links
        self.calls = []

    def seq2gene_links(self, tss_region=None):
        self.calls.append(tss_region)
        return self.links


def test_enrich_peak_genes_local_and_remote(monkeypatch: pytest.MonkeyPatch) -> None:
    links = pd.DataFrame({"gene_name": ["A", "B", "A"], "gene_id": ["g1", "g2", "g1"]})
    annotation = FakeLinksAnnotation(links)
    captured = {}
    monkeypatch.setattr(
        enrichment, "ora", lambda genes, **kwargs: captured.update(genes=genes, **kwargs)
    )
    enrich_peak_genes(annotation, upstream=1000, local=True, p_adjust="none")
    assert annotation.calls == [(-1000, 0)]
    assert captured["genes"] == ["A", "B"]

    monkeypatch.setattr(
        enrichment, "enrichr", lambda genes, **kwargs: captured.update(genes=genes)
    )
    enrich_peak_genes(annotation, downstream=500)
    assert annotation.calls[-1] == (0, 500)
    assert captured["genes"] == ["A", "B"]

    empty = FakeLinksAnnotation(pd.DataFrame())
    with pytest.raises(ValueError, match="no peak-gene links"):
        enrich_peak_genes(empty)
    with pytest.raises(KeyError, match="gene column"):
        enrich_peak_genes(annotation, gene_column="SYMBOL")


def test_as_ranked_frame_variants() -> None:
    frame = _as_ranked_frame(
        pd.DataFrame({"Gene": ["a", "b"], "Rank": [1.0, 2.0]})
    )
    assert frame.to_dict("list") == {"gene": ["b", "a"], "score": [2.0, 1.0]}

    unknown = _as_ranked_frame(pd.DataFrame({"id": ["a", "b"], "value": [1, 3]}))
    assert unknown["gene"].tolist() == ["b", "a"]

    with pytest.raises(ValueError, match="at least two columns"):
        _as_ranked_frame(pd.DataFrame({"gene": ["a"]}))
    with pytest.raises(ValueError, match="ranked gene list is empty"):
        _as_ranked_frame(pd.DataFrame({"gene": [""], "score": [1]}))

    ranked = _as_ranked_frame(["a", "b", "a", None])
    assert ranked["gene"].tolist() == ["a", "b"]
    assert ranked["score"].tolist() == [2.0, 1.0]

    mapped = _as_ranked_frame({"a": 1.0, "b": 2.0})
    assert mapped["gene"].tolist() == ["b", "a"]


def test_gsea_wrapper_argument_variants(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []

    class Result:
        res2d = pd.DataFrame({"Term": ["A"], "NES": [1.5]})

    def fake_prerank(**kwargs):
        calls.append(kwargs)
        return Result()

    fake = types.ModuleType("gseapy")
    fake.prerank = fake_prerank
    monkeypatch.setitem(sys.modules, "gseapy", fake)

    gsea({"a": 2.0, "b": 1.0}, "hallmark.gmt")
    assert calls[-1]["gene_sets"] == "hallmark.gmt"

    gsea(["a", "b"], {"A": ["a", None, "b"]})
    assert calls[-1]["gene_sets"] == {"A": ["a", "b"]}

    gsea(["a", "b"], ["A", "B"])
    assert calls[-1]["gene_sets"] == ["A", "B"]
    assert calls[-1]["min_size"] == 15
    assert calls[-1]["max_size"] == 500
    assert calls[-1]["verbose"] is False

    class ResultsOnly:
        results = pd.DataFrame({"Term": ["B"]})

    fake.prerank = lambda **kwargs: ResultsOnly()
    assert gsea(["a"], "x.gmt")["Term"].tolist() == ["B"]

    fake.prerank = lambda **kwargs: object()
    with pytest.raises(RuntimeError, match="no GSEA result"):
        gsea(["a"], "x.gmt")


def test_gsea_peak_genes_restricts_to_links(monkeypatch: pytest.MonkeyPatch) -> None:
    links = pd.DataFrame({"gene_name": ["A", "B"]})
    annotation = FakeLinksAnnotation(links)
    captured = {}

    def fake_gsea(frame, **kwargs):
        captured["frame"] = frame
        return pd.DataFrame({"Term": ["A"]})

    monkeypatch.setattr(enrichment, "gsea", fake_gsea)
    out = gsea_peak_genes(annotation, {"A": 2.0, "B": 1.0, "C": 3.0}, upstream=100)
    assert out["Term"].tolist() == ["A"]
    assert captured["frame"]["gene"].tolist() == ["A", "B"]
    assert annotation.calls == [(-100, 0)]

    with pytest.raises(ValueError, match="no peak-gene links"):
        gsea_peak_genes(FakeLinksAnnotation(pd.DataFrame()), {"A": 1.0})
    with pytest.raises(KeyError, match="gene column"):
        gsea_peak_genes(annotation, {"A": 1.0}, gene_column="SYMBOL")
    with pytest.raises(ValueError, match="no ranked genes"):
        gsea_peak_genes(annotation, {"C": 1.0})


def test_enrichment_plot_frame_normalises_columns() -> None:
    with pytest.raises(ValueError, match="result is empty"):
        _enrichment_plot_frame(pd.DataFrame())

    indexed = pd.DataFrame({"PValue": [0.01, 0.2]}, index=["A", "B"])
    indexed.index.name = "Term"
    frame = _enrichment_plot_frame(indexed)
    assert frame["Term"].tolist() == ["B", "A"]

    overlap = pd.DataFrame(
        {
            "Term": ["A", "B"],
            "Adjusted P-value": [0.01, 0.1],
            "Overlap": ["3/10", "1/20"],
        }
    )
    frame = _enrichment_plot_frame(overlap)
    assert frame.set_index("Term")["count"].to_dict() == {"A": 3, "B": 1}
    assert frame.set_index("Term")["ratio"].to_dict() == {"A": 0.3, "B": 0.05}

    ratio = pd.DataFrame(
        {
            "Term": ["A"],
            "padj": [0.05],
            "GeneRatio": ["2/8"],
            "GeneSetSize": [40],
        }
    )
    assert _enrichment_plot_frame(ratio).iloc[0]["ratio"] == pytest.approx(2 / 40)

    numeric = pd.DataFrame(
        {"Term": ["A"], "FDR q-val": [0.0], "Overlap": [4], "GeneSetSize": [8]}
    )
    frame = _enrichment_plot_frame(numeric)
    assert frame.iloc[0]["count"] == 4
    assert frame.iloc[0]["ratio"] == pytest.approx(0.5)
    assert np.isfinite(frame.iloc[0]["neglog10"])

    with pytest.raises(ValueError, match="no p-value"):
        _enrichment_plot_frame(pd.DataFrame({"Term": ["A"], "NES": [1.0]}))
    with pytest.raises(ValueError, match="no finite p-values"):
        _enrichment_plot_frame(
            pd.DataFrame({"Term": ["A"], "PValue": [np.nan]})
        )


def test_plot_enrichment_accepts_axes_and_rejects_kind() -> None:
    result = pd.DataFrame(
        {"Term": ["A", "B"], "PValue": [0.01, 0.2], "Overlap": ["1/10", "2/10"]}
    )
    fig, ax = plt.subplots()
    returned_fig, returned_ax = plot_enrichment(result, kind="dot", ax=ax, top=1)
    assert returned_fig is fig and returned_ax is ax
    plt.close(fig)

    with pytest.raises(ValueError, match="kind must be"):
        plot_enrichment(result, kind="volcano")
