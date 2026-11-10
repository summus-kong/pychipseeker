
import pandas as pd

from pychipseeker.cli import main


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
