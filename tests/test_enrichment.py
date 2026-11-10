import math
import sys
import types

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pychipseeker.enrichment import (
    enrichr,
    gsea,
    ora,
    plot_enrichment,
    read_gmt,
)


def test_enrichr_empty_gene_list():
    with pytest.raises(ValueError):
        enrichr([])


def test_enrichr_requires_gseapy_when_missing():
    try:
        import gseapy  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="gseapy"):
            enrichr(["TP53"])
    else:
        pytest.skip("gseapy is installed")


def test_ora_hypergeometric_and_bh():
    gene_sets = {
        "term A": [f"g{i}" for i in range(1, 11)],
        "term B": [f"g{i}" for i in range(11, 16)],
        "term C": ["g1", "g11", "g21", "g22"],
    }
    query = ["g1", "g2", "g3", "g11", "g12"]
    result = ora(
        query,
        gene_sets,
        background=[f"g{i}" for i in range(1, 101)],
    )
    assert list(result["Term"]) == ["term A", "term C", "term B"]
    first = result.iloc[0]
    assert first["Overlap"] == 3
    assert first["GeneSetSize"] == 10
    assert first["BackgroundSize"] == 100
    assert first["QuerySize"] == 5
    assert first["FoldEnrichment"] == pytest.approx(6.0)
    expected = sum(
        math.comb(10, k) * math.comb(90, 5 - k) / math.comb(100, 5)
        for k in range(3, 6)
    )
    assert first["PValue"] == pytest.approx(expected)
    # BH values are monotone in the original p-value order.
    adjusted = result.set_index("Term")["AdjustedPValue"]
    assert adjusted["term A"] <= adjusted["term C"] <= adjusted["term B"]
    assert adjusted["term A"] >= first["PValue"]
    assert first["OddsRatio"] == pytest.approx((3 * 88) / (2 * 7))
    assert first["ORLower"] < first["OddsRatio"] < first["ORUpper"]


def test_ora_gmt_and_min_overlap(tmp_path):
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "term A\tdesc\tg1\tg2\tg3\n"
        "term B\tdesc\tg1\tg4\tg5\n"
        "term C\tdesc\tg9\tg10\n",
        encoding="utf-8",
    )
    sets = read_gmt(gmt)
    assert sets["term A"] == ["g1", "g2", "g3"]

    result = ora(
        ["g1", "g2"],
        gmt,
        background=["g1", "g2", "g3", "g4", "g5"],
        min_overlap=2,
    )
    assert list(result["Term"]) == ["term A"]
    assert result.iloc[0]["Genes"] == "g1;g2"


def test_plot_enrichment_bar_and_dot():
    result = ora(
        ["g1", "g2", "g3"],
        {
            "term A": [f"g{i}" for i in range(1, 11)],
            "term B": ["g1", "g20", "g21"],
        },
        background=[f"g{i}" for i in range(1, 101)],
    )
    fig, ax = plot_enrichment(result, top=2)
    assert len(ax.patches) == 2
    assert ax.get_xlabel() == "-log10(adjusted p-value)"
    plt.close(fig)

    fig, ax = plot_enrichment(result, top=2, kind="dot")
    assert len(ax.collections) == 1
    assert len(ax.get_yticks()) == 2
    plt.close(fig)


def test_gsea_wrapper_uses_gseapy(monkeypatch):
    calls = {}

    class Result:
        res2d = pd.DataFrame(
            {
                "Term": ["term A"],
                "NES": [1.8],
                "NOM p-val": [0.01],
                "FDR q-val": [0.05],
            }
        )

    def fake_prerank(**kwargs):
        calls.update(kwargs)
        return Result()

    fake = types.ModuleType("gseapy")
    fake.prerank = fake_prerank
    monkeypatch.setitem(sys.modules, "gseapy", fake)

    table = gsea(
        {"g1": 3.0, "g2": 2.0, "g3": 1.0},
        {"term A": ["g1", "g2"]},
        permutation_num=10,
    )
    assert list(table["Term"]) == ["term A"]
    assert calls["rnk"]["gene"].tolist() == ["g1", "g2", "g3"]
    assert calls["min_size"] == 15
    assert calls["permutation_num"] == 10


def test_gsea_requires_gseapy_when_missing():
    try:
        import gseapy  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="gseapy"):
            gsea(["g1"], {"term": ["g1"]})
    else:
        pytest.skip("gseapy is installed")


def test_ora_rejects_empty_inputs():
    with pytest.raises(ValueError, match="gene_list"):
        ora([], {"term": ["g1"]})
    with pytest.raises(ValueError, match="gene_sets"):
        ora(["g1"], {})
    with pytest.raises(ValueError, match="background"):
        ora(["g1"], {"term": ["g1"]}, background=[])
    with pytest.raises(ValueError, match="no gene set"):
        ora(["gX"], {"term": ["g1"]}, background=["g1", "gX"])
