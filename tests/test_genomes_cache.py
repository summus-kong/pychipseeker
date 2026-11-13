"""Genome cache, spec parsing, assembly detection and provenance paths."""

from __future__ import annotations

import hashlib
import json
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

from pychipseeker import genomes
from pychipseeker.genomes import (
    _bundle_path,
    _cache_bundle,
    _download_bundle,
    _download_file,
    _installed_bundle,
    _installed_provenance,
    _npz_model,
    cached_genomes,
    clear_genome_cache,
    detect_assembly,
    download_genome,
    download_gtf,
    genome_provenance,
    load_annotation,
    load_genome,
    parse_genome_spec,
    prune_genome_cache,
    save_bundle,
)

EXAMPLE_GTF = Path(__file__).resolve().parents[1] / "examples" / "example.gtf"


@pytest.fixture()
def cache(tmp_path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setattr(genomes, "cache_dir", lambda: tmp_path)
    return tmp_path


def test_parse_genome_spec_aliases_and_validation() -> None:
    human = parse_genome_spec("ensembl:homo_sapiens")
    assert human["provider"] == "ensembl"
    assert human["release"] == genomes.DEFAULT_ENSEMBL_RELEASE
    assert human["url"].endswith("Homo_sapiens.GRCh38.110.gtf.gz")

    directory = parse_genome_spec(" ensembl:mus_musculus@100 ")
    assert directory["url"].endswith("Mus_musculus.GRCm39.100.gtf.gz")

    mouse = parse_genome_spec("gencode:GRCm39")
    assert mouse["release"] == "M33"
    assert mouse["url"].endswith("GRCm39.gencode.vM33.annotation.gtf.gz")

    human_default = parse_genome_spec("gencode:human")
    assert human_default["release"] == "44"

    ftp = parse_genome_spec("url:ftp://example.org/a.gtf.gz")
    assert ftp["url"] == "ftp://example.org/a.gtf.gz"

    with pytest.raises(ValueError, match="http"):
        parse_genome_spec("url:example.org/a.gtf.gz")
    with pytest.raises(ValueError, match="unknown GENCODE"):
        parse_genome_spec("gencode:unicorn")
    with pytest.raises(TypeError, match="must be a string"):
        parse_genome_spec(None)  # type: ignore[arg-type]


def test_detect_assembly_from_names_and_sizes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with pytest.raises(ValueError, match="without peaks"):
        detect_assembly()

    with pytest.raises(ValueError, match="chrom column"):
        detect_assembly(peaks=pd.DataFrame({"start": [0]}))

    yeast = pd.DataFrame(
        {"chrom": ["I", "II", "III"], "start": [0, 0, 0], "end": [10, 10, 10]}
    )
    assert detect_assembly(yeast) == "sacCer3"

    fly = pd.DataFrame(
        {"chrom": ["2L", "3L", "X"], "start": [0, 0, 0], "end": [10, 10, 10]}
    )
    assert detect_assembly(fly) == "dm6"

    human = pd.DataFrame(
        {"chrom": ["chr20"], "start": [0], "end": [10]}
    )
    with pytest.raises(ValueError, match="cannot distinguish hg38"):
        detect_assembly(human)
    human_short = pd.DataFrame(
        {"chrom": ["chr1"], "start": [0], "end": [10]}
    )
    with pytest.raises(ValueError, match="hg38/hg19"):
        detect_assembly(human_short)
    with pytest.raises(ValueError, match="could not detect"):
        detect_assembly(pd.DataFrame({"chrom": ["chrUn"], "start": [0], "end": [1]}))

    hg38 = detect_assembly(
        chrom_sizes={"1": 248956422, "2": 242193529}
    )
    assert hg38 == "hg38"

    monkeypatch.setattr(
        genomes, "_assembly_size_candidates", lambda sizes, tolerance: ["hg38", "hg19"]
    )
    with pytest.raises(ValueError, match="several assemblies"):
        detect_assembly(chrom_sizes={"1": 1})


def test_detect_assembly_reads_bigwig_header(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Handle:
        def chroms(self):
            return {
                "1": genomes.ASSEMBLY_CHROM_SIZES["mm10"]["1"],
                "2": genomes.ASSEMBLY_CHROM_SIZES["mm10"]["2"],
            }

        def close(self):
            return None

    fake = types.ModuleType("pyBigWig")
    fake.open = lambda path: Handle()
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)
    assert detect_assembly(bigwig="fake.bw") == "mm10"

    fake.open = lambda path: None
    with pytest.raises(ValueError, match="could not open bigWig"):
        detect_assembly(bigwig="missing.bw")


def test_cache_dir_platform_branches(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    windows_root = tmp_path / "win"
    monkeypatch.setattr(genomes.sys, "platform", "win32")
    monkeypatch.setenv("LOCALAPPDATA", str(windows_root))
    assert genomes.cache_dir() == windows_root / "pychipseeker" / "genomes"
    assert genomes.cache_dir().is_dir()

    monkeypatch.setattr(genomes.sys, "platform", "linux")
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("XDG_CACHE_HOME", str(xdg))
    assert genomes.cache_dir() == xdg / "pychipseeker" / "genomes"


def test_installed_and_cache_lookup(cache: Path) -> None:
    assert _installed_bundle("definitely_missing") is None
    assert _installed_provenance("definitely_missing") is None

    empty = cache / "genome_empty.npz"
    empty.write_bytes(b"")
    assert _cache_bundle("empty") is None
    good = cache / "genome_good.npz"
    good.write_bytes(b"data")
    assert _cache_bundle("good") == good


def test_download_file_retries_and_fails(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    attempts = []

    class Response:
        def __init__(self, payload):
            self._buffer = payload

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self, size=-1):
            data, self._buffer = self._buffer, b""
            return data

    def flaky(url, timeout=None):
        attempts.append(url)
        if len(attempts) == 1:
            raise OSError("first attempt")
        return Response(b"payload")

    monkeypatch.setattr(genomes.urllib.request, "urlopen", flaky)
    monkeypatch.setattr(genomes.time, "sleep", lambda seconds: None)
    dest = cache / "file.bin"
    _download_file("https://example.org/file.bin", dest, tries=2)
    assert dest.read_bytes() == b"payload"
    assert not (cache / "file.bin.part").exists()
    assert len(attempts) == 2

    def broken(url, timeout=None):
        raise OSError("down")

    monkeypatch.setattr(genomes.urllib.request, "urlopen", broken)
    with pytest.raises(RuntimeError, match="could not download"):
        _download_file("https://example.org/file.bin", cache / "other.bin", tries=1)
    assert not (cache / "other.bin").exists()


def test_download_bundle_checks_digest_and_falls_back(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"bundle-bytes"
    digest = hashlib.sha256(payload).hexdigest()

    def fake_download(url, dest, tries=4):
        dest.parent.mkdir(parents=True, exist_ok=True)
        if url.endswith(".json"):
            dest.write_text(json.dumps({"bundle_sha256": digest}), encoding="utf-8")
        else:
            dest.write_bytes(payload)

    monkeypatch.setattr(genomes, "_bundle_bases", lambda: ("https://base",))
    monkeypatch.setattr(genomes, "_download_file", fake_download)
    path = _download_bundle("demo")
    assert path.read_bytes() == payload

    def bad_download(url, dest, tries=4):
        dest.parent.mkdir(parents=True, exist_ok=True)
        if url.endswith(".json"):
            dest.write_text(json.dumps({"bundle_sha256": "0" * 64}), encoding="utf-8")
        else:
            dest.write_bytes(payload)

    monkeypatch.setattr(genomes, "_download_file", bad_download)
    with pytest.raises(RuntimeError, match="could not fetch genome bundle"):
        _download_bundle("bad")
    assert not (cache / "genome_bad.npz").exists()

    attempts = []

    def provenance_fails(url, dest, tries=4):
        attempts.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if url.endswith(".json") and "badbase" in url:
            raise OSError("no provenance")
        if url.endswith(".json"):
            raise OSError("no provenance")
        if "badbase" in url:
            raise OSError("no bundle")
        dest.write_bytes(payload)

    monkeypatch.setattr(
        genomes, "_bundle_bases", lambda: ("https://badbase", "https://goodbase")
    )
    monkeypatch.setattr(genomes, "_download_file", provenance_fails)
    recovered = _download_bundle("recovered")
    assert recovered.read_bytes() == payload
    assert any("goodbase" in url for url in attempts)


def test_bundle_path_prefers_installed_then_cache(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(genomes, "_installed_bundle", lambda name: "installed")
    assert _bundle_path("x") == "installed"

    monkeypatch.setattr(genomes, "_installed_bundle", lambda name: None)
    monkeypatch.setattr(genomes, "_cache_bundle", lambda name: "cached")
    assert _bundle_path("x") == "cached"

    monkeypatch.setattr(genomes, "_cache_bundle", lambda name: None)
    monkeypatch.setattr(genomes, "_download_bundle", lambda name: "downloaded")
    assert _bundle_path("x") == "downloaded"


def test_download_gtf_cache_hit_and_override(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = "url:https://example.org/genes.gtf.gz"
    cache_path, _, _ = genomes._gtf_cache_paths(spec)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(b"cached")
    assert download_gtf(spec) == cache_path

    override = cache / "manual.gtf"
    override.write_bytes(b"manual")
    assert download_gtf(spec, dest=override) == override

    fresh = "url:https://example.org/fresh.gtf"
    writes = []

    def fake_download(url, dest, tries=4):
        writes.append(url)
        dest.write_bytes(b"fresh")

    monkeypatch.setattr(genomes, "_download_file", fake_download)
    result = download_gtf(fresh)
    assert result.read_bytes() == b"fresh"
    assert writes == ["https://example.org/fresh.gtf"]


def test_load_annotation_refresh_reparses(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = "url:" + EXAMPLE_GTF.resolve().as_uri()
    monkeypatch.setattr(genomes, "_loaded_specs", {})

    first = load_annotation(spec)
    npz_path = next((cache / "gtf").glob("spec_*.npz"))
    mtime = npz_path.stat().st_mtime_ns

    def no_download(url, dest, tries=4):
        raise AssertionError("cached GTF should be reused")

    monkeypatch.setattr(genomes, "_download_file", no_download)
    second = load_annotation(spec, refresh=True)
    assert len(second.genes) == len(first.genes)
    assert npz_path.stat().st_mtime_ns >= mtime

    # A fresh process reads the npz without the GTF.
    monkeypatch.setattr(genomes, "_loaded_specs", {})
    third = load_annotation(spec)
    assert len(third.genes) == len(first.genes)


def test_build_manifest_bundle_and_load(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    name = "demogenome"
    monkeypatch.setitem(genomes.GENOME_ALIASES, name, name)
    monkeypatch.setitem(
        genomes.GENOME_SOURCES,
        name,
        {
            "provider": "test",
            "release": "1",
            "url": "https://example.org/demo.gtf",
            "description": "demo",
        },
    )

    def fake_download(url, dest, tries=4):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(EXAMPLE_GTF.read_bytes())

    monkeypatch.setattr(genomes, "_download_file", fake_download)
    npz_path = download_genome(name)
    assert npz_path.is_file()
    assert (cache / "genome_demogenome.json").is_file()
    assert (cache / "genome_demogenome.gtf").is_file()

    monkeypatch.setattr(genomes, "_loaded", {})
    model = load_genome(name)
    assert not model.genes.empty
    record = genome_provenance(name)
    assert record["provider"] == "test"
    assert record["description"] == "demo"


def test_genome_provenance_spec_and_manifest(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = "ensembl:mouse@110"
    record = genome_provenance(spec)
    assert record["downloaded"] is False
    assert record["provider"] == "ensembl"
    assert record["release"] == "110"

    manifest = genome_provenance("mm39")
    assert manifest["downloaded"] is False
    assert manifest["provider"] == "ensembl"


def test_cached_genomes_and_prune(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (cache / "genome_good.npz").write_bytes(b"data")
    (cache / "genome_good.json").write_text(
        json.dumps({"provider": "test", "release": "1"}), encoding="utf-8"
    )
    (cache / "genome_zero.npz").write_bytes(b"")
    (cache / "genome_orphan.json").write_text("{}", encoding="utf-8")
    (cache / "genome_part.npz.part").write_bytes(b"partial")
    gtf = cache / "gtf"
    gtf.mkdir()
    (gtf / "spec_x.npz").write_bytes(b"spec")
    (gtf / "spec_x.json").write_text(
        json.dumps({"assembly": "url:x", "provider": "url", "release": ""}),
        encoding="utf-8",
    )
    (gtf / "spec_bad.json").write_text("not json", encoding="utf-8")

    rows = {row["assembly"]: row for row in cached_genomes()}
    assert rows["good"]["source"] == "cache"
    assert rows["good"]["provider"] == "test"
    assert rows["url:x"]["source"] == "spec-cache"
    assert "spec_bad" not in rows

    removed = {path.name for path in prune_genome_cache()}
    assert {"genome_zero.npz", "genome_orphan.json", "genome_part.npz.part"} <= removed
    assert (cache / "genome_good.npz").exists()
    assert not (cache / "genome_zero.npz").exists()
    assert not (gtf / "spec_bad.json").exists()


def test_clear_genome_cache_specific_and_all(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(genomes.GENOME_ALIASES, "demo", "demo")
    (cache / "genome_demo.npz").write_bytes(b"data")
    (cache / "genome_demo.json").write_text("{}", encoding="utf-8")
    (cache / "genome_demo.gtf.gz").write_bytes(b"gtf")
    gtf = cache / "gtf"
    gtf.mkdir()
    spec = "url:https://example.org/genes.gtf"
    _, npz_path, _ = genomes._gtf_cache_paths(spec)
    npz_path.write_bytes(b"spec")
    npz_path.with_suffix(".json").write_text("{}", encoding="utf-8")

    removed = clear_genome_cache("demo")
    assert {path.name for path in removed} == {
        "genome_demo.npz",
        "genome_demo.json",
        "genome_demo.gtf.gz",
    }

    assert clear_genome_cache(spec)
    assert not npz_path.exists()
    assert clear_genome_cache() == []
    assert not gtf.exists()


def test_load_genome_bundle_roundtrip_and_errors(
    cache: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pychipseeker.annotation import GeneAnnotation

    monkeypatch.setattr(genomes, "_loaded", {})
    monkeypatch.setattr(genomes, "_installed_bundle", lambda name: None)
    monkeypatch.setattr(genomes, "_cache_bundle", lambda name: None)
    monkeypatch.setattr(genomes, "_download_bundle", lambda name: None)
    with pytest.raises(FileNotFoundError, match="no genome bundle"):
        load_genome("hg38", refresh=True)

    reference = GeneAnnotation.from_gtf(EXAMPLE_GTF)
    monkeypatch.setattr(genomes, "_loaded", {})
    monkeypatch.setattr(
        genomes, "_cache_bundle", lambda name: cache / f"genome_{name}.npz"
    )
    save_bundle(reference, cache / "genome_hg38.npz")
    loaded = load_genome("hg38")
    assert len(loaded.genes) == len(reference.genes)
    assert load_genome("hg38") is loaded

    # A custom alias outside BUNDLED_GENOMES/GENOME_SOURCES resolves through
    # the generic bundle lookup.
    monkeypatch.setitem(genomes.GENOME_ALIASES, "customassembly", "customassembly")
    monkeypatch.setattr(genomes, "_loaded", {})
    save_bundle(reference, cache / "genome_customassembly.npz")
    custom = load_genome("customassembly")
    assert len(custom.genes) == len(reference.genes)

    monkeypatch.setattr(genomes, "_cache_bundle", lambda name: None)
    monkeypatch.setattr(genomes, "_bundle_path", lambda name: None)
    monkeypatch.setattr(genomes, "_loaded", {})
    with pytest.raises(FileNotFoundError, match="no genome bundle"):
        load_genome("customassembly")


def test_npz_model_roundtrip_preserves_tables(cache: Path) -> None:
    from pychipseeker.annotation import GeneAnnotation

    reference = GeneAnnotation.from_gtf(EXAMPLE_GTF)
    path = cache / "genome_roundtrip.npz"
    save_bundle(reference, path)
    restored = _npz_model(path)
    assert restored.source == ("npz", str(path))
    assert len(restored.genes) == len(reference.genes)
    assert list(restored.genes.columns) == list(reference.genes.columns)
    assert restored.transcripts["chrom"].tolist() == reference.transcripts[
        "chrom"
    ].tolist()
