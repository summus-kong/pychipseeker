"""Tests for the bundled ``TxDb="assembly"`` genome-annotation lookup."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import pychipseeker as pcs
from pychipseeker.annotation import GeneAnnotation
from pychipseeker.genomes import (
    TABLES,
    canonical_assembly,
    is_assembly_name,
    load_genome,
    save_bundle,
    supported_genomes,
)

_REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_GTF = _REPO_ROOT / "examples" / "example.gtf"


def _bundle_available(name: str) -> bool:
    from importlib.resources import as_file, files

    try:
        with as_file(
            files("pychipseeker").joinpath("data", f"genome_{name}.npz")
        ) as path:
            return Path(path).is_file()
    except (FileNotFoundError, ModuleNotFoundError):
        return False


def _stringize(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalise a frame for comparison: missing data becomes empty text."""

    out = frame.reset_index(drop=True).copy()
    for col in out.columns:
        if pd.api.types.is_integer_dtype(out[col]):
            out[col] = out[col].astype(np.int64)
        else:
            out[col] = (
                out[col].astype(object)
                .where(~pd.isna(out[col]), "")
                .astype(str)
            )
    return out


@pytest.fixture(scope="module")
def example_model() -> GeneAnnotation:
    return GeneAnnotation.from_gtf(EXAMPLE_GTF)


@pytest.fixture()
def fake_bundled_genome(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Install a tiny example-derived bundle under a patched table path."""

    model = GeneAnnotation.from_gtf(EXAMPLE_GTF)
    bundle = tmp_path / "genome_demo.npz"
    save_bundle(model, bundle, drop=())

    from pychipseeker import genomes as genomes_mod

    def patched_path(_assembly: str):
        return bundle

    monkeypatch.setattr(genomes_mod, "_bundle_path", patched_path)
    monkeypatch.setitem(genomes_mod.GENOME_ALIASES, "demo", "demo")
    return "demo"


def test_assembly_name_checks() -> None:
    assert is_assembly_name("hg38") and is_assembly_name("GRCh38")
    assert is_assembly_name(" mm10 ") and is_assembly_name("GRCm37")
    assert not is_assembly_name("sacCer3")
    assert not is_assembly_name(None) and not is_assembly_name(42)


def test_canonical_assembly_and_errors() -> None:
    assert canonical_assembly("GRCh38") == "hg38"
    assert canonical_assembly("mm9") == "mm9"
    assert supported_genomes() == ("hg38", "hg19", "mm10", "mm9")
    with pytest.raises(ValueError, match="hg38.*hg19.*mm10.*mm9"):
        canonical_assembly("sacCer3")
    with pytest.raises(ValueError, match="TxDb must be"):
        canonical_assembly(None)  # type: ignore[arg-type]


def test_load_genome_roundtrip(
    fake_bundled_genome: str, example_model: GeneAnnotation
) -> None:
    model = load_genome(fake_bundled_genome, refresh=True)
    assert isinstance(model, GeneAnnotation)
    for table in TABLES:
        got = _stringize(getattr(model, table))
        want = _stringize(getattr(example_model, table))
        pd.testing.assert_frame_equal(got, want)
    # cached afterwards
    assert load_genome(fake_bundled_genome) is model


def test_annotate_peak_with_assembly_name(fake_bundled_genome: str) -> None:
    gtf_frame = pd.read_csv(EXAMPLE_GTF, sep="\t", comment="#", header=None)
    chroms = gtf_frame[0].unique().tolist()
    peaks = pd.DataFrame(
        {
            "chrom": chroms,
            "start": [0] * len(chroms),
            "end": [100] * len(chroms),
        }
    )
    ann = pcs.annotate_peak(peaks=peaks, TxDb=fake_bundled_genome)
    table = ann.to_chipseeker()
    assert not table.empty
    assert table["annotation"].notna().any()


def test_get_tag_matrix_with_assembly_name(fake_bundled_genome: str) -> None:
    gtf_frame = pd.read_csv(EXAMPLE_GTF, sep="\t", comment="#", header=None)
    genes = gtf_frame[gtf_frame[2] == "gene"]
    first = genes.iloc[0]
    peaks = pd.DataFrame(
        {
            "chrom": [str(first[0])],
            "start": [int(first[3]) - 1],  # TSS of a known gene (0-based)
            "end": [int(first[3])],
        }
    )
    promoters = pcs.get_promoters(
        TxDb=fake_bundled_genome, upstream=1000, downstream=1000
    )
    tag_matrix = pcs.get_tag_matrix(peak=peaks, windows=promoters, nbin=10)
    # tag-matrix rows follow the peaks, columns follow the bins, as in R
    assert tag_matrix.shape == (len(peaks), 10)


def test_lookup_skips_download_when_installed_exists(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, example_model: GeneAnnotation
) -> None:
    from pychipseeker import genomes as genomes_mod

    fake = tmp_path / "genome_hg38.npz"
    save_bundle(example_model, fake)
    monkeypatch.setattr(genomes_mod, "_installed_bundle", lambda n: fake)
    downloads: list[str] = []
    monkeypatch.setattr(
        genomes_mod,
        "_download_bundle",
        lambda n: downloads.append(n) or fake,
    )
    load_genome("hg38", refresh=True)
    assert downloads == []


def test_lookup_downloads_when_not_installed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, example_model: GeneAnnotation
) -> None:
    from pychipseeker import genomes as genomes_mod

    fake = tmp_path / "genome_hg19.npz"
    save_bundle(example_model, fake)
    monkeypatch.setattr(genomes_mod, "_installed_bundle", lambda n: None)
    monkeypatch.setattr(genomes_mod, "_download_bundle", lambda n: fake)
    assert load_genome("hg19") is load_genome("hg19")


def _peak_df(x: int = 1_000_000) -> pd.DataFrame:
    return pd.DataFrame({"chrom": ["chr1"], "start": [x], "end": [x + 249]})


@pytest.mark.skipif(not _bundle_available("mm10"), reason="mm10 bundle not built yet")
def test_bundled_mm10_load() -> None:
    model = load_genome("mm10")
    assert isinstance(model, GeneAnnotation)
    for table in TABLES:
        assert len(getattr(model, table)) > 0 or table.endswith("3"), table
    for col in ("gene_id", "chrom", "start", "end", "strand", "tss"):
        assert col in model.genes.columns


@pytest.mark.skipif(not _bundle_available("mm10"), reason="mm10 bundle not built yet")
def test_bundled_mm10_annotate_and_tagmatrix() -> None:
    peaks = _peak_df()
    ann = pcs.annotate_peak(peaks=peaks, TxDb="mm10")
    table = ann.to_chipseeker()
    assert table["annotation"].notna().all()

    promoters = pcs.get_promoters(TxDb="mm10", upstream=1000, downstream=1000)
    assert len(promoters) > 1000
    tag = pcs.get_tag_matrix(peak=peaks, windows=promoters, nbin=5)
    assert tag.shape[1] == 5
