"""Focused coverage for the PDF report helpers in :mod:`pychipseeker.report`."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pychipseeker import annotate_peak
from pychipseeker import report as report_mod
from pychipseeker.annotation import GeneAnnotation, PeakAnnotation
from pychipseeker.report import _methods_page, _qc_page, _stats_page, _text_page


def _cell(cell) -> str:
    """Return the text rendered inside a matplotlib table cell."""

    return cell.get_text().get_text()


def test_text_page_writes_title_footer_and_lines() -> None:
    fig = _text_page("Report", [("Peaks", "10"), ("Level", "gene")], footer="footer")
    texts = [text.get_text() for text in fig.texts]
    assert "Report" in texts
    assert "Peaks" in texts and "10" in texts
    assert "footer" in texts

    bare = _text_page("Report", [])
    assert len(bare.texts) == 1
    plt.close("all")


def test_qc_page_formats_each_value_type() -> None:
    table = pd.DataFrame(
        {
            "Metric": [
                "Mitochondrial fraction",
                "FRiP",
                "Mean signal",
                "TSS enrichment",
                "Mean width (bp)",
                "Peaks",
                "Genome",
            ],
            "Value": [
                0.123456,
                0.5,
                2.34567,
                3.141592,
                125.55,
                12345,
                "mm10",
            ],
        }
    )
    fig = _qc_page(table)
    cells = fig.axes[0].tables[0]
    values = [_cell(cells[(row, 1)]) for row in range(1, len(table) + 1)]
    assert values == [
        "0.1235",
        "0.5000",
        "2.3457",
        "3.142",
        "125.5",
        "12,345",
        "mm10",
    ]
    plt.close(fig)


def test_methods_page_reports_signal_sources(examples_dir: Path) -> None:
    model = GeneAnnotation.from_gtf(examples_dir / "example.gtf")
    annotation = PeakAnnotation(
        peaks=pd.DataFrame(
            {
                "peak_id": ["p1"],
                "chrom": ["chr1"],
                "start": [100],
                "end": [200],
            }
        ),
        annotation=pd.DataFrame({"peak_id": ["p1"]}),
        gene_annotation=model,
    )
    full = _methods_page(
        annotation,
        tss_region=(-1000, 2000),
        bigwig=Path("signal.bw"),
        bam=Path("reads.bam"),
        n_signal=12,
        seed=7,
        qc=True,
    )
    body = full.texts[1].get_text()
    assert "bigWig: signal.bw" in body
    assert "12 transcripts sampled, seed 7" in body
    assert "BAM: reads.bam" in body
    assert "QC: peak widths" in body
    assert "annotation level: transcript" in body

    bare = _methods_page(
        annotation,
        tss_region=(-500, 500),
        bigwig=None,
        bam=None,
        n_signal=0,
        seed=1,
        qc=False,
    )
    bare_body = bare.texts[1].get_text()
    assert "no bigWig supplied" in bare_body
    assert "BAM:" not in bare_body
    assert "QC:" not in bare_body
    plt.close("all")


def test_stats_page_renders_features(examples_dir: Path) -> None:
    annotation = annotate_peak(
        examples_dir / "example_peaks.bed",
        gtf=examples_dir / "example.gtf",
    )
    fig = _stats_page(annotation)
    cells = fig.axes[0].tables[0]
    features = {
        _cell(cells[(row, 0)])
        for row in range(1, len(annotation.anno_stat()) + 1)
    }
    assert any(feature.startswith("Promoter") for feature in features)
    assert "Distal Intergenic" in features
    plt.close(fig)


def test_build_report_validates_inputs(tmp_path: Path, examples_dir: Path) -> None:
    peaks = examples_dir / "example_peaks.bed"
    gtf = examples_dir / "example.gtf"

    with pytest.raises(ValueError, match="provide gtf"):
        report_mod.build_report(peaks, output=tmp_path / "out.pdf")

    with pytest.raises(FileNotFoundError):
        report_mod.build_report(
            peaks,
            gtf=gtf,
            bigwig=tmp_path / "missing.bw",
            output=tmp_path / "out.pdf",
        )

    with pytest.raises(FileNotFoundError):
        report_mod.build_report(
            peaks,
            gtf=gtf,
            bam=tmp_path / "missing.bam",
            output=tmp_path / "out.pdf",
        )


def test_build_report_renders_tag_matrix_pages(
    tmp_path: Path, examples_dir: Path
) -> None:
    output = report_mod.build_report(
        examples_dir / "example_peaks.bed",
        gtf=examples_dir / "example.gtf",
        output=tmp_path / "tag.pdf",
        include_tag_matrix=True,
        sample=10,
    )
    assert output.exists()
    assert output.read_bytes().startswith(b"%PDF")
    plt.close("all")


def test_build_report_samples_and_plots_bigwig(
    tmp_path: Path, examples_dir: Path
) -> None:
    pytest.importorskip("pyBigWig")
    bigwig = Path(__file__).parent / "data" / "test.bw"
    if not bigwig.exists():  # pragma: no cover - packaging guard
        pytest.skip("tests/data/test.bw is missing")

    gtf = tmp_path / "signal.gtf"
    gtf.write_text(
        '1\ttest\tgene\t1001\t5000\t.\t+\t.\tgene_id "G1"; gene_name "A";\n'
        '1\ttest\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        '1\ttest\tgene\t10001\t15000\t.\t-\t.\tgene_id "G2"; gene_name "B";\n'
        '1\ttest\texon\t10001\t12000\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n',
        encoding="utf-8",
    )
    peaks = tmp_path / "signal_peaks.bed"
    peaks.write_text("1\t1100\t1200\n1\t10100\t10200\n", encoding="utf-8")

    output = report_mod.build_report(
        peaks,
        gtf=gtf,
        bigwig=bigwig,
        output=tmp_path / "signal.pdf",
        sample=1,
        upstream=200,
        downstream=200,
        include_tag_matrix=False,
    )
    assert output.exists()
    assert output.read_bytes().startswith(b"%PDF")

    # ``sample`` larger than the transcript table keeps every row.
    full_output = report_mod.build_report(
        peaks,
        gtf=gtf,
        bigwig=bigwig,
        output=tmp_path / "signal_full.pdf",
        sample=100,
        include_tag_matrix=False,
    )
    assert full_output.exists()
    plt.close("all")


def test_build_report_skips_empty_qc_reference_tables(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transcripts = pd.DataFrame(
        {
            "chrom": ["1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "tss": [100],
        }
    )

    class FakeModel:
        def features(self, level: str):
            if level == "gene":
                return None
            return transcripts

    class FakeAnnotation:
        peaks = pd.DataFrame(
            {
                "peak_id": ["p1"],
                "chrom": ["1"],
                "start": [100],
                "end": [200],
            }
        )
        level = "transcript"
        gene_annotation = FakeModel()

        def anno_stat(self) -> pd.DataFrame:
            return pd.DataFrame(
                {"Feature": ["Promoter"], "Frequency": [100.0], "Count": [1]}
            )

        def _figure(self):
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])
            return fig, ax

        plot_anno_pie = _figure
        plot_anno_bar = _figure
        plot_dist_to_tss = _figure

    captured: dict[str, object] = {}

    def fake_qc_summary(peaks, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame({"Metric": ["Peaks"], "Value": [1]})

    monkeypatch.setattr(report_mod, "annotate_peak", lambda **_: FakeAnnotation())
    monkeypatch.setattr("pychipseeker.qc.qc_summary", fake_qc_summary)

    output = report_mod.build_report(
        pd.DataFrame({"chrom": ["1"], "start": [100], "end": [200]}),
        gtf=pd.DataFrame({"chrom": ["1"]}),
        output=tmp_path / "qc.pdf",
        include_tag_matrix=False,
        qc=True,
    )
    assert output.exists()
    assert captured["chrom_sizes"] == {"1": 200}
    assert captured["genes"].equals(transcripts)
    plt.close("all")
