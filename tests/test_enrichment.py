import pytest

from pychipseeker.enrichment import enrichr


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
