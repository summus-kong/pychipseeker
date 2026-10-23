from pathlib import Path

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
