"""Additional branch coverage for the ``pychipseeker`` command line interface."""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path

import pandas as pd
import pytest

from pychipseeker import cli
from pychipseeker.cli import _cmd_cache, _load_annotation_source, _read_gene_list, main


def test_parse_ints_rejects_wrong_arity(tmp_path: Path, examples_dir: Path) -> None:
    peaks = examples_dir / "example_peaks.bed"
    gtf = examples_dir / "example.gtf"

    with pytest.raises(SystemExit, match="two comma-separated"):
        main(
            [
                "annotate",
                "--peaks",
                str(peaks),
                "--gtf",
                str(gtf),
                "--output",
                str(tmp_path / "out.tsv"),
                "--tss-region",
                "1,2,3",
            ]
        )

    with pytest.raises(SystemExit, match="two comma-separated"):
        main(
            [
                "report",
                "--peaks",
                str(peaks),
                "--gtf",
                str(gtf),
                "--output",
                str(tmp_path / "out.pdf"),
                "--tss-region",
                "5",
            ]
        )


def test_read_gene_list_table_and_text(tmp_path: Path) -> None:
    table = tmp_path / "genes.tsv"
    pd.DataFrame({"gene_name": ["A", "B", "A", None], "score": [1, 2, 3, 4]}).to_csv(
        table, sep="\t", index=False
    )
    assert _read_gene_list(table) == ["A", "B"]

    gene_id_table = tmp_path / "genes.csv"
    pd.DataFrame({"gene_id": ["G1", "G2"]}).to_csv(gene_id_table, index=False)
    assert _read_gene_list(gene_id_table) == ["G1", "G2"]

    plain = tmp_path / "genes.txt"
    plain.write_text("A\n\nB\n", encoding="utf-8")
    assert _read_gene_list(plain) == ["A", "B"]


def test_enrich_validation_errors(tmp_path: Path, examples_dir: Path) -> None:
    missing_column = tmp_path / "no-gene.tsv"
    pd.DataFrame({"other": ["x"]}).to_csv(missing_column, sep="\t", index=False)
    with pytest.raises(SystemExit, match="gene column not found"):
        main(["enrich", "--annotated", str(missing_column), "--outdir", str(tmp_path)])

    empty_genes = tmp_path / "empty-genes.tsv"
    pd.DataFrame({"gene_name": [None, None]}).to_csv(empty_genes, sep="\t", index=False)
    with pytest.raises(SystemExit, match="no genes to enrich"):
        main(
            [
                "enrich",
                "--annotated",
                str(empty_genes),
                "--outdir",
                str(tmp_path),
            ]
        )

    annotated = tmp_path / "annotated.tsv"
    pd.DataFrame({"gene_name": ["GeneA"]}).to_csv(annotated, sep="\t", index=False)
    with pytest.raises(SystemExit, match="--local requires --gmt"):
        main(
            [
                "enrich",
                "--annotated",
                str(annotated),
                "--local",
                "--outdir",
                str(tmp_path),
            ]
        )


def test_enrich_remote_and_plot(
    tmp_path: Path,
    examples_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    annotated = tmp_path / "annotated.tsv"
    pd.DataFrame({"gene_name": ["GeneA", "GeneB"]}).to_csv(
        annotated, sep="\t", index=False
    )
    calls: dict[str, object] = {}

    def fake_enrichr(genes, **kwargs):
        calls["genes"] = genes
        calls.update(kwargs)
        return pd.DataFrame(
            {"Term": ["term A"], "PValue": [0.01], "AdjustedPValue": [0.05]}
        )

    monkeypatch.setattr(cli, "enrichr", fake_enrichr)
    rc = main(
        [
            "enrich",
            "--annotated",
            str(annotated),
            "--outdir",
            str(tmp_path / "remote"),
            "--plot",
        ]
    )
    assert rc == 0
    assert calls["genes"] == ["GeneA", "GeneB"]
    assert (tmp_path / "remote" / "enrichment.tsv").exists()
    assert (tmp_path / "remote" / "enrichment.png").exists()
    capsys.readouterr()

    def broken_enrichr(*args, **kwargs):
        raise ValueError("Enrichr is offline")

    monkeypatch.setattr(cli, "enrichr", broken_enrichr)
    with pytest.raises(SystemExit, match="offline"):
        main(
            [
                "enrich",
                "--annotated",
                str(annotated),
                "--outdir",
                str(tmp_path / "broken"),
            ]
        )


def test_filter_chromosome_options(tmp_path: Path, examples_dir: Path) -> None:
    output = tmp_path / "kept.bed"
    rc = main(
        [
            "filter",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--output",
            str(output),
            "--chroms",
            "chr1,chr2",
            "--exclude-chroms",
            "chr2",
        ]
    )
    assert rc == 0
    kept = pd.read_csv(output, sep="\t")
    assert set(kept["chrom"]) == {"chr1"}


def test_load_annotation_source_branches(
    examples_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pychipseeker import genomes

    sentinel = object()
    monkeypatch.setattr(genomes, "detect_assembly", lambda peaks=None: "toy")
    monkeypatch.setattr(genomes, "is_assembly_name", lambda name: name == "toy")
    monkeypatch.setattr(genomes, "load_genome", lambda name: sentinel)
    assert _load_annotation_source("auto") is sentinel

    monkeypatch.setattr(genomes, "is_assembly_name", lambda name: False)
    monkeypatch.setattr(genomes, "is_genome_spec", lambda spec: spec == "toy.gtf")
    monkeypatch.setattr(genomes, "load_annotation", lambda spec: sentinel)
    assert _load_annotation_source("toy.gtf") is sentinel

    monkeypatch.setattr(genomes, "is_genome_spec", lambda spec: False)
    model = _load_annotation_source(str(examples_dir / "example.gtf"))
    assert {"G1", "G2"}.issubset(set(model.genes["gene_id"]))


def test_qc_requires_signal_and_reports_source_errors(
    tmp_path: Path, examples_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(SystemExit, match="qc needs --bigwig or --bam"):
        main(["qc", "--peaks", str(examples_dir / "example_peaks.bed")])

    with pytest.raises(SystemExit):
        main(
            [
                "qc",
                "--peaks",
                str(examples_dir / "example_peaks.bed"),
                "--bigwig",
                str(examples_dir / "example_peaks.bed"),
                "--gtf",
                str(tmp_path / "missing.gtf"),
            ]
        )

    captured: dict[str, object] = {}

    def fake_qc(peaks, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame({"Metric": ["Peaks"], "Value": [1]})

    monkeypatch.setattr("pychipseeker.qc.qc_summary", fake_qc)
    output = tmp_path / "qc.tsv"
    rc = main(
        [
            "qc",
            "--peaks",
            str(examples_dir / "example_peaks.bed"),
            "--bigwig",
            str(examples_dir / "example_peaks.bed"),
            "--gtf",
            str(examples_dir / "example.gtf"),
            "--output",
            str(output),
        ]
    )
    assert rc == 0
    assert captured["genes"] is not None
    assert output.exists()


def test_cache_empty_and_invalid_actions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from pychipseeker import genomes

    monkeypatch.setattr(genomes, "cache_dir", lambda: tmp_path)
    monkeypatch.setattr(genomes, "cached_genomes", lambda: [])
    assert main(["cache", "list"]) == 0
    assert "No genome bundles found." in capsys.readouterr().out

    assert main(["cache", "prune"]) == 0
    assert "Nothing to prune." in capsys.readouterr().out

    assert main(["cache", "clear", "hg38"]) == 0
    assert "Nothing to remove." in capsys.readouterr().out

    with pytest.raises(SystemExit, match="cache clear needs"):
        main(["cache", "clear"])

    with pytest.raises(SystemExit, match="unknown cache action"):
        _cmd_cache(argparse.Namespace(cache_command="bogus"))


def test_report_accepts_custom_tss_region(tmp_path: Path, examples_dir: Path) -> None:
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
            "--tss-region=-1000,1000",
            "--no-tag-matrix",
        ]
    )
    assert rc == 0
    assert output.read_bytes().startswith(b"%PDF")


def test_module_entry_point(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["pychipseeker", "--version"])
    with pytest.raises(SystemExit) as excinfo:
        runpy.run_module("pychipseeker.cli", run_name="__main__")
    assert excinfo.value.code == 0
