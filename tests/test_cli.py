
import json
from pathlib import Path

import pandas as pd
import pytest

from pychipseeker.cli import main
from pychipseeker.io import read_peaks


def _real_bigwig() -> Path:
    pytest.importorskip("pyBigWig")
    path = Path(__file__).parent / "data" / "test.bw"
    if not path.exists():  # pragma: no cover - packaging guard
        pytest.skip("tests/data/test.bw is missing")
    return path


def test_cli_callpeak_and_filter(tmp_path):
    output = tmp_path / "peaks.bed"
    rc = main(
        [
            "callpeak",
            "--bigwig",
            str(_real_bigwig()),
            "--output",
            str(output),
            "--threshold",
            "0.05",
            "--min-width",
            "1",
            "--max-gap",
            "0",
            "--chroms",
            "1",
        ]
    )
    assert rc == 0
    called = read_peaks(output)
    assert len(called) > 0

    filtered = tmp_path / "filtered.bed"
    rc = main(
        [
            "filter",
            "--peaks",
            str(output),
            "--output",
            str(filtered),
            "--min-width",
            "3",
        ]
    )
    assert rc == 0
    assert len(read_peaks(filtered)) <= len(called)


def test_cli_qc_and_cache_json(tmp_path, examples_dir, monkeypatch, capsys):
    qc_out = tmp_path / "qc.tsv"
    rc = main(
        [
            "qc",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--bigwig",
            str(_real_bigwig()),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--no-tss",
            "--output",
            str(qc_out),
        ]
    )
    assert rc == 0
    table = pd.read_csv(qc_out, sep="\t")
    assert "FRiP" in set(table["Metric"])
    capsys.readouterr()  # discard the QC table printed to stdout

    from pychipseeker import genomes as genomes_mod

    monkeypatch.setattr(genomes_mod, "cache_dir", lambda: tmp_path)
    (tmp_path / "genome_dm6.npz").write_bytes(b"npz")
    (tmp_path / "genome_dm6.json").write_text(
        json.dumps({"provider": "ensembl", "release": "110"}), encoding="utf-8"
    )
    rc = main(["cache", "list", "--json"])
    assert rc == 0
    rows = json.loads(capsys.readouterr().out)
    assert any(row["assembly"] == "dm6" for row in rows)

    rc = main(["cache", "info", "hg38"])
    assert rc == 0
    provenance = json.loads(capsys.readouterr().out)
    assert provenance["assembly"] == "hg38"

    (tmp_path / "genome_hg19.npz.part").write_bytes(b"partial")
    rc = main(["cache", "prune"])
    assert rc == 0
    assert not (tmp_path / "genome_hg19.npz.part").exists()


def test_cli_annotate_negative_tss_region(tmp_path, examples_dir):
    """Regression: --tss-region -3000,3000 must be accepted (space form)."""

    out = tmp_path / "annotated.tsv"
    rc = main(
        [
            "annotate",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(out),
            "--tss-region",
            "-3000,3000",
            "--summary",
            "--quiet",
        ]
    )
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 10
    assert "annotation" in df.columns


def test_cli_annotate_tss_region_equals_form(tmp_path, examples_dir):
    """The --tss-region=-3000,3000 (equals) form keeps working."""

    out = tmp_path / "annotated_eq.tsv"
    rc = main(
        [
            "annotate",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(out),
            "--tss-region=-3000,3000",
            "--quiet",
        ]
    )
    assert rc == 0


def test_cli_annotate(tmp_path, examples_dir):
    out = tmp_path / "annotated.tsv"
    plot_dir = tmp_path / "plots"
    rc = main(
        [
            "annotate",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(out),
            "--upstream",
            "10",
            "--downstream",
            "10",
            "--summary",
            "--plot-dir",
            str(plot_dir),
            "--quiet",
        ]
    )
    assert rc == 0
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 10
    assert "annotation" in df.columns
    assert (tmp_path / "annotated.annoStat.tsv").exists()
    assert (plot_dir / "annotation_pie.png").exists()
    assert (plot_dir / "annotation_bar.png").exists()
    assert (plot_dir / "distance_to_tss.png").exists()


def test_cli_plot(tmp_path, examples_dir):
    annotated = tmp_path / "annotated.tsv"
    main(
        [
            "annotate",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(annotated),
            "--upstream",
            "10",
            "--downstream",
            "10",
            "--quiet",
        ]
    )
    outdir = tmp_path / "plots2"
    rc = main(
        [
            "plot",
            "--annotated",
            str(annotated),
            "--outdir",
            str(outdir),
            "--prefix",
            "test",
        ]
    )
    assert rc == 0
    assert (outdir / "test_annotation_pie.png").exists()
    assert (outdir / "test_annotation_bar.png").exists()
    assert (outdir / "test_distance_to_tss.png").exists()


def test_cli_report(tmp_path, examples_dir):
    output = tmp_path / "report.pdf"
    rc = main(
        [
            "report",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(output),
            "--no-tag-matrix",
        ]
    )
    assert rc == 0
    assert output.read_bytes().startswith(b"%PDF")
    assert output.stat().st_size > 10_000


def test_cli_enrich_local(tmp_path, examples_dir):
    annotated = tmp_path / "annotated.tsv"
    main(
        [
            "annotate",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(annotated),
            "--quiet",
        ]
    )

    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "term A\tdescription\tGeneA\tGeneB\n"
        "term B\tdescription\tGeneA\tGeneC\n",
        encoding="utf-8",
    )
    background = tmp_path / "background.txt"
    background.write_text("GeneA\nGeneB\nGeneC\nGeneD\nGeneE\n", encoding="utf-8")
    outdir = tmp_path / "enrichment"
    rc = main(
        [
            "enrich",
            "--annotated",
            str(annotated),
            "--local",
            "--gmt",
            str(gmt),
            "--background",
            str(background),
            "--outdir",
            str(outdir),
        ]
    )
    assert rc == 0
    result = pd.read_csv(outdir / "enrichment.tsv", sep="\t")
    assert {"Term", "PValue", "AdjustedPValue"} <= set(result.columns)
    assert len(result) >= 1


def test_cli_cache_list_and_clear(tmp_path, monkeypatch, capsys):
    from pychipseeker import genomes as genomes_mod

    monkeypatch.setattr(genomes_mod, "cache_dir", lambda: tmp_path)
    (tmp_path / "genome_dm6.npz").write_bytes(b"npz")
    (tmp_path / "genome_dm6.json").write_text(
        '{"provider": "ensembl", "release": "110"}', encoding="utf-8"
    )

    assert main(["cache", "list"]) == 0
    listed = capsys.readouterr().out
    assert "dm6" in listed
    assert "cache" in listed

    assert main(["cache", "clear", "dm6"]) == 0
    assert not (tmp_path / "genome_dm6.npz").exists()


def test_cli_cache_download(tmp_path, monkeypatch, capsys):
    from pychipseeker import genomes as genomes_mod

    target = tmp_path / "genome_dm6.npz"
    target.write_bytes(b"npz")
    monkeypatch.setattr(genomes_mod, "download_genome", lambda _name: target)

    assert main(["cache", "download", "dm6"]) == 0
    assert "dm6" in capsys.readouterr().out
