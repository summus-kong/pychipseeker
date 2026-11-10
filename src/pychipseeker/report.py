"""Multi-page PDF report for a peak annotation.

The report bundles the standard ChIPseeker-style figures (annotation
distribution, TSS distance, coverage, optional bigWig signal profiles and
tag-matrix heatmaps) plus the detailed ``annoStat`` table and a methods page.
It is the programmatic counterpart of ChIPseeker's ``plotAnnoPie`` /
``plotAnnoBar`` workflow and powers ``pychipseeker report``.
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Any, Optional, Sequence, Union

import numpy as np
import pandas as pd

from ._version import __version__
from .annotation import PeakAnnotation, annotate_peak

_PAGE = (11.0, 8.5)


def _pyplot():
    import matplotlib.pyplot as plt

    return plt


def _text_page(
    title: str,
    lines: Sequence[tuple[str, str]],
    footer: Optional[str] = None,
):
    plt = _pyplot()
    fig = plt.figure(figsize=_PAGE)
    fig.text(0.06, 0.90, title, fontsize=24, weight="bold")
    y = 0.80
    for label, value in lines:
        fig.text(0.06, y, str(label), fontsize=10, color="0.35")
        fig.text(0.34, y, str(value), fontsize=10)
        y -= 0.042
    if footer:
        fig.text(
            0.06,
            0.10,
            footer,
            fontsize=9,
            color="0.45",
            linespacing=1.5,
        )
    return fig


def _stats_page(annotation: PeakAnnotation):
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=_PAGE)
    ax.axis("off")
    ax.set_title("annoStat - detailed annotation categories", fontsize=14, loc="left")
    stat = annotation.anno_stat()
    columns = [c for c in ("Feature", "Frequency", "Count") if c in stat.columns]
    table = ax.table(
        cellText=stat[columns].round(3).astype(str).values,
        colLabels=columns,
        loc="center",
        cellLoc="center",
        colWidths=[0.42, 0.25, 0.18][: len(columns)],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.9)
    for cell in table.get_celld().values():
        cell.set_edgecolor("white")
    return fig


def _qc_page(table: pd.DataFrame):
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=_PAGE)
    ax.axis("off")
    ax.set_title("Quality control metrics", fontsize=14, loc="left")
    rows: list[tuple[str, str]] = []
    for row in table.itertuples(index=False):
        metric = str(row.Metric)
        value = row.Value
        if isinstance(value, (float, np.floating)):
            if "fraction" in metric.lower() or "FRiP" in metric:
                text = f"{float(value):.4f}"
            elif "signal" in metric.lower():
                text = f"{float(value):.4f}"
            elif "enrichment" in metric.lower():
                text = f"{float(value):.3f}"
            else:
                text = f"{float(value):,.1f}"
        elif isinstance(value, (int, np.integer)):
            text = f"{int(value):,}"
        else:
            text = str(value)
        rows.append((metric, text))
    table_obj = ax.table(
        cellText=rows,
        colLabels=["Metric", "Value"],
        loc="center",
        cellLoc="center",
        colWidths=[0.52, 0.28],
    )
    table_obj.auto_set_font_size(False)
    table_obj.set_fontsize(10)
    table_obj.scale(1.0, 1.8)
    for cell in table_obj.get_celld().values():
        cell.set_edgecolor("white")
    return fig


def _methods_page(
    annotation: PeakAnnotation,
    *,
    tss_region: tuple[int, int],
    bigwig: Optional[Path],
    bam: Optional[Path],
    n_signal: int,
    seed: int,
    qc: bool,
):
    signal_line = (
        f"   bigWig: {bigwig} ({n_signal:,} transcripts sampled, seed {seed})\n"
        if bigwig is not None
        else "   no bigWig supplied; signal panels omitted\n"
    )
    if bam is not None:
        signal_line += f"   BAM: {bam}\n"
    if qc:
        signal_line += "   QC: peak widths, blacklist, FRiP and TSS enrichment\n"
    body = (
        "1. Inputs\n"
        f"   peaks: {len(annotation.peaks):,}\n"
        f"   annotation level: {annotation.level}\n"
        f"   promoter window: {tss_region[0]:,} / +{tss_region[1]:,} bp\n\n"
        "2. Annotation\n"
        "   pychipseeker.annotate_peak with ChIPseeker-compatible categories\n"
        "   (Promoter, 5' UTR, 3' UTR, Exon, Intron, Downstream, Intergenic).\n\n"
        "3. Figures\n"
        "   plot_anno_pie / plot_anno_bar / plot_dist_to_tss / covplot /\n"
        "   plot_avg_prof / tag_heatmap / plot_tss_profile / plot_tss_heatmap.\n"
        f"{signal_line}\n"
        "4. Reproduce\n"
        "   pychipseeker report --peaks <peaks> --txdb <assembly> "
        "--output report.pdf\n\n"
        f"   pychipseeker {__version__}"
    )
    plt = _pyplot()
    fig = plt.figure(figsize=_PAGE)
    fig.text(0.06, 0.92, "Methods / reproducibility", fontsize=17, weight="bold")
    fig.text(
        0.06,
        0.86,
        body,
        fontsize=10,
        va="top",
        family="monospace",
        linespacing=1.45,
    )
    return fig


def build_report(
    peaks: Union[str, Path, pd.DataFrame],
    *,
    gtf: Optional[Union[str, Path, pd.DataFrame]] = None,
    TxDb: Optional[Union[str, Path, pd.DataFrame]] = None,
    bigwig: Optional[Union[str, Path]] = None,
    bam: Optional[Union[str, Path]] = None,
    blacklist: Optional[Union[str, Path, pd.DataFrame]] = None,
    output: Union[str, Path] = "pychipseeker_report.pdf",
    level: str = "transcript",
    tss_region: Sequence[int] = (-3000, 3000),
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
    sample: int = 3000,
    seed: int = 20261026,
    title: Optional[str] = None,
    include_tag_matrix: bool = True,
    qc: bool = False,
) -> Path:
    """Build a multi-page PDF report and return its path.

    Parameters mirror :func:`pychipseeker.annotate_peak`.  ``bigwig`` is
    optional; when supplied, TSS profile and heatmap pages are added using a
    reproducible subsample of ``sample`` transcripts (``seed``).  ``upstream``
    and ``downstream`` optionally define an asymmetric signal window.
    ``bam`` and ``blacklist`` feed the QC/FRiP pages enabled by ``qc=True``.
    """

    from .compat import covplot
    from .plotting import (
        extract_tss_matrix,
        plot_tss_heatmap,
        plot_tss_profile,
    )
    from .tagmatrix import (
        get_promoters,
        get_tag_matrix,
        plot_avg_prof,
        tag_heatmap,
    )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    source = gtf if gtf is not None else TxDb
    if source is None:
        raise ValueError("provide gtf=... or TxDb=...")
    annotation = annotate_peak(
        peaks=peaks,
        gtf=source,
        level=level,
        tss_region=tuple(int(x) for x in tss_region),
        verbose=False,
    )

    bigwig_path = None if bigwig is None else Path(bigwig)
    if bigwig_path is not None and not bigwig_path.exists():
        raise FileNotFoundError(bigwig_path)
    bam_path = None if bam is None else Path(bam)
    if bam_path is not None and not bam_path.exists():
        raise FileNotFoundError(bam_path)

    signal_up = (
        abs(int(tss_region[0])) if upstream is None else abs(int(upstream))
    )
    signal_down = (
        abs(int(tss_region[1])) if downstream is None else abs(int(downstream))
    )

    plt = _pyplot()
    pages: list[tuple[Any, str]] = []
    n_peaks = len(annotation.peaks)
    total = len(annotation.gene_annotation.features(level))
    n_signal = 0
    title_text = title or "ChIPseeker-style peak annotation report"
    stat = annotation.anno_stat()
    top = stat.iloc[0] if not stat.empty else None

    overview = [
        ("Peaks", f"{n_peaks:,}"),
        ("Annotation level", level),
        ("Promoter window", f"{int(tss_region[0]):,} / +{int(tss_region[1]):,} bp"),
        ("Features", f"{total:,}"),
        ("Top category", "" if top is None else f"{top['Feature']} ({top['Frequency']:.1f}%)"),
        ("Generated", _dt.date.today().isoformat()),
        ("pychipseeker", __version__),
    ]
    pages.append(
        (
            _text_page(
                title_text,
                overview,
                "Figures use the same matplotlib defaults as the ChIPseeker "
                "theme_grey translation. Signal panels may use a reproducible "
                "transcript subsample.",
            ),
            "Title",
        )
    )

    pages.append((annotation.plot_anno_pie()[0], "Annotation distribution (pie)"))
    pages.append((annotation.plot_anno_bar()[0], "Annotation distribution (bar)"))
    pages.append((annotation.plot_dist_to_tss()[0], "Distance to TSS"))
    pages.append((covplot(annotation.peaks, title="Peak coverage")[0], "Peak coverage"))

    if include_tag_matrix:
        try:
            promoters = get_promoters(
                TxDb=annotation.gene_annotation,
                upstream=abs(int(tss_region[0])),
                downstream=abs(int(tss_region[1])),
            )
            tag = get_tag_matrix(peak=annotation.peaks, windows=promoters, nbin=60)
            pages.append(
                (plot_avg_prof(tag, conf=None, title="Peak profile")[0], "Peak profile")
            )
            pages.append((tag_heatmap(tag, title="Peak heatmap")[0], "Peak heatmap"))
        except Exception:  # pragma: no cover - report best effort
            pass

    transcripts = annotation.gene_annotation.features("transcript")

    if qc:
        from .qc import qc_summary

        tables = [
            annotation.gene_annotation.features("gene"),
            annotation.gene_annotation.features("transcript"),
        ]
        chrom_sizes: dict[str, int] = {}
        for table in tables:
            if table is None or table.empty:
                continue
            for chrom, value in table.groupby("chrom", sort=False)["end"].max().items():
                chrom_sizes[str(chrom)] = max(
                    chrom_sizes.get(str(chrom), 0), int(value)
                )
        qc_table = qc_summary(
            annotation.peaks,
            bigwig=bigwig_path,
            bam=bam_path,
            genes=transcripts,
            blacklist=blacklist,
            chrom_sizes=chrom_sizes or None,
            tss_upstream=signal_up,
            tss_downstream=signal_down,
            tss_sample=sample,
            seed=int(seed),
        )
        pages.append((_qc_page(qc_table), "Quality control"))

    if bigwig_path is not None:
        if sample and int(sample) < len(transcripts):
            subset = transcripts.sample(n=int(sample), random_state=int(seed))
        else:
            subset = transcripts
        n_signal = len(subset)
        matrix = extract_tss_matrix(
            bigwig_path,
            subset,
            n_bins=100,
            upstream=signal_up,
            downstream=signal_down,
        )
        assert isinstance(matrix, np.ndarray)
        pages.append(
            (
                plot_tss_profile(
                    matrix,
                    upstream=signal_up,
                    downstream=signal_down,
                    title="Signal around TSS",
                    show_sem=False,
                )[0],
                f"Signal profile around TSS ({n_signal:,} sampled transcripts)",
            )
        )
        pages.append(
            (
                plot_tss_heatmap(
                    matrix,
                    upstream=signal_up,
                    downstream=signal_down,
                    cmap="Reds",
                    sort_rows=True,
                    title="Signal around TSS",
                )[0],
                f"Signal heatmap around TSS ({n_signal:,} sampled transcripts)",
            )
        )

    pages.append((_stats_page(annotation), "Annotation statistics"))
    pages.append(
        (
            _methods_page(
                annotation,
                tss_region=(int(tss_region[0]), int(tss_region[1])),
                bigwig=bigwig_path,
                bam=bam_path,
                n_signal=n_signal,
                seed=int(seed),
                qc=qc,
            ),
            "Methods",
        )
    )

    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(str(output)) as pdf:
        for fig, name in pages:
            fig.suptitle(name, fontsize=9, color="0.45", y=0.985)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    return output
