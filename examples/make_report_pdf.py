"""Build a complete multi-page PDF report from the bundled example data.

Run::

    python examples/make_report_pdf.py

Uses the virtual data shipped with pychipseeker —``examples/example_peaks.bed`` (11 hand-wired peaks covering every genomic
category) and ``examples/example.gtf`` — plus the bundled ENCODE H3K4me3
demo signal (``ENCSR000CDQ`` chr5 20 Mb subset, mm10), and writes every
visualization as one page of ``examples/output/pychipseeker_report.pdf``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import pychipseeker as pcs
from pychipseeker._version import __version__
from pychipseeker.tagmatrix import example_signal_path, example_signal_provenance

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"

#: pseudo-peak calling on the bundled demo track (fold-change threshold).
PEAK_THRESHOLD = 1.2
#: keep only the strongest peaks by mean fold change, for readable figures.
TOP_PEAKS = 120
BIN_BP = 25
GAP_BINS = 2  # merge nearby positive bins separated by <= 2 empty bins
MIN_BINS = 3  # discard putative peaks narrower than 3 bins

CHROM_SIZE = 151_834_684


def _title_page(stat: pd.DataFrame, df: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.93, "pychipseeker", ha="center", fontsize=28, weight="bold")
    fig.text(0.5, 0.895, "ChIPseeker-style peak annotation report", ha="center", fontsize=13)
    fig.text(
        0.5, 0.86,
        f"pychipseeker {__version__} · virtual example data + ENCODE demo signal",
        ha="center", fontsize=10, color="0.4",
    )

    lines = [
        "Input data (bundled, virtual):",
        f"  peaks : examples/example_peaks.bed  ({len(df)} synthetic peaks)",
        "ENCODE demo signal (bundled):",
        "  ENCSR000CDQ  H3K4me3 fold change over control",
        "  ENCF893ZGR (mm10) -> chr5 20 Mb subset, 166 KB",
    ]
    fig.text(0.12, 0.76, "\n".join(lines), fontsize=11, family="monospace")

    fig.text(0.12, 0.66, "Annotation statistics (ChIPseeker annoStat):", fontsize=12, weight="bold")
    ax = fig.add_axes((0.12, 0.30, 0.76, 0.30))
    ax.axis("off")
    table = ax.table(
        cellText=stat.round(2).astype(str).values,
        colLabels=[c.replace("_", " ").title() for c in stat.columns],
        loc="upper center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)

    genes = df.groupby("gene_name").size().sort_values(ascending=False)
    fig.text(
        0.12, 0.22,
        "Peaks assigned per gene:\n  "
        + ";  ".join(f"{g}: {n}" for g, n in genes.items()),
        fontsize=11, family="monospace",
    )

    kinds = df["annotation"].value_counts()
    fig.text(
        0.12, 0.10,
        "Categories hit by the demo peaks:\n  "
        + "\n  ".join(f"{k}: {n}" for k, n in kinds.items()),
        fontsize=11, family="monospace",
    )
    fig.text(
        0.5, 0.03,
        "promoter window = 10 bp (deliberately narrow so UTR/exon/intron categories appear)",
        ha="center", fontsize=9, style="italic", color="0.35",
    )
    return fig


def _final_page() -> plt.Figure:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.85, "Reproduce this report", ha="center", fontsize=16, weight="bold")
    fig.text(
        0.10, 0.58,
        "python examples/make_report_pdf.py\n\n"
        "from pychipseeker.report steps:\n"
        "  ann = pcs.annotate_peak(peaks=..., gtf=..., upstream=10, downstream=10)\n"
        "  ann.plot_anno_pie() / ann.plot_anno_bar() / ann.plot_dist_to_tss()\n"
        "  pcs.covplot(ann.peaks)\n"
        "  promoters = pcs.getPromoters(TxDb=ann.gene_annotation, upstream=1000, downstream=1000)\n"
        "  tag = pcs.getTagMatrix(peak=ann.peaks, windows=promoters, nbin=50)\n"
        "  pcs.plotAvgProf(tag, conf=None)  ·  pcs.tagHeatmap(tag)\n"
        "  ann.plot_tss_profile(signal=..., window=100, n_bins=50)\n"
        "  ann.plot_tss_heatmap(signal=..., window=100, n_bins=50)\n\n"
        "ENCODE demo signal (bundled bigWig):\n"
        "  bw_path = pcs.tagmatrix.example_signal_path()\n"
        "  prov    = pcs.tagmatrix.example_signal_provenance()\n"
        "  peaks   = _pseudo_peaks(bw_path, threshold=1.2)   # fold-change calling\n"
        "  cov     = pcs.covplot(peaks, weightCol='score')\n"
        "  tag     = pcs.get_tag_matrix_from_bigwig(bw_path, peaks +/- 1 kb, nbin=50)\n"
        "  pcs.plotAvgProf(tag, conf=0.95)  ·  pcs.tagHeatmap(tag)\n"
        "  mat = pcs.extract_tss_matrix(bw_path, summit_table, window=1000)",
        fontsize=10.5, family="monospace", va="top",
    )
    fig.text(0.5, 0.06, f"pychipseeker {__version__}", ha="center", fontsize=9, color="0.4")
    return fig


def _pseudo_peaks(bw_path: Path) -> pd.DataFrame:
    """Call pseudo-peaks on the bundled demo track.

    Adjacent 25 bp bins whose mean fold change is >= ``PEAK_THRESHOLD`` are
    merged across gaps of at most ``GAP_BINS`` empty bins; runs narrower than
    ``MIN_BINS`` bins are dropped.  ``score`` records the track mean.
    """

    import pyBigWig

    bw = pyBigWig.open(str(bw_path))
    try:
        intervals = bw.intervals("chr5") or []
    finally:
        bw.close()

    hot = [(start, end, float(v)) for start, end, v in intervals if v >= PEAK_THRESHOLD]
    if not hot:
        return pd.DataFrame(columns=["chrom", "start", "end", "score"])

    merged = []
    cur_start, cur_end, cur_vals = hot[0][0], hot[0][1], [hot[0][2]]
    for start, end, value in hot[1:]:
        if start - cur_end <= GAP_BINS * BIN_BP:
            cur_end = end
            cur_vals.append(value)
        else:
            merged.append((cur_start, cur_end, cur_vals))
            cur_start, cur_end, cur_vals = start, end, [value]
    merged.append((cur_start, cur_end, cur_vals))

    peaks = [
        {"chrom": "chr5", "start": start, "end": end, "score": float(np.mean(values))}
        for start, end, values in merged
        if (end - start) >= MIN_BINS * BIN_BP
    ]
    peaks = pd.DataFrame(peaks, columns=["chrom", "start", "end", "score"])
    if len(peaks) > TOP_PEAKS:
        peaks = (
            peaks.sort_values("score", ascending=False)
            .head(TOP_PEAKS)
            .sort_values("start")
            .reset_index(drop=True)
        )
    return peaks


def _summit_table(peaks: pd.DataFrame, bw_path: Path) -> pd.DataFrame:
    """Summit anchors (highest-signal bin of each pseudo-peak)."""

    if peaks.empty:
        return pd.DataFrame(columns=["chrom", "tss", "strand"])
    import pyBigWig

    bw = pyBigWig.open(str(bw_path))
    try:
        rows = []
        for chrom, start, end in peaks[["chrom", "start", "end"]].itertuples(index=False):
            bins = bw.intervals(chrom, int(start), int(end)) or []
            best_value = -np.inf
            for bin_start, _bin_end, value in bins:
                if value > best_value:
                    best_value, best_bin = value, bin_start
            if best_bin is not None:
                rows.append({"chrom": chrom, "tss": int(best_bin) + BIN_BP // 2, "strand": "+"})
    finally:
        bw.close()
    return pd.DataFrame(rows, columns=["chrom", "tss", "strand"])


def _encode_track_page(
    bw_path: Path,
    prov: dict,
    peaks: pd.DataFrame,
) -> plt.Figure:
    """Chromosome-signal line over the bundled 20 Mb chr5 window."""

    import pyBigWig

    subset = prov["subset"]
    start = int(subset["start"])
    end = int(subset["end"])
    bw = pyBigWig.open(str(bw_path))
    try:
        n_bins = 800
        stats = bw.stats("chr5", start, end, nBins=n_bins, type="mean")
    finally:
        bw.close()
    xs = np.linspace(start, end, n_bins)
    ys = np.asarray([np.nan if v is None else float(v) for v in stats], dtype=float)

    fig = plt.figure(figsize=(8.5, 6.2))
    ax = fig.add_axes((0.10, 0.30, 0.87, 0.60))
    ax.fill_between(xs, 0, ys, color="#3182bd", alpha=0.30, linewidth=0)
    ax.plot(xs, ys, color="#3182bd", linewidth=1.0)

    for _row in peaks.itertuples(index=False):
        ax.axvspan(_row.start, _row.end, color="#fb8072", alpha=0.35, linewidth=0)
    y_max = float(np.nanmax(ys))
    ax.axhline(PEAK_THRESHOLD, color="#fb8072", linewidth=1.0, linestyle="--")
    ax.text(
        start + (end - start) * 0.62,
        0.94 * y_max,
        f"pseudo-peak threshold = {PEAK_THRESHOLD}× fold change (dashed)",
        ha="left",
        va="top",
        fontsize=8.5,
        color="#b34a3a",
    )
    ax.set_xlim(start, end)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(f"chr5 position ({start:,} — {end:,}, mm10)")
    ax.set_ylabel("H3K4me3 fold change over control")
    ax.set_title("ENCODE H3K4me3 demo track — chr5 (20 Mb bundled subset)", loc="left")
    ax.set_facecolor("#EBEBEB")
    ax.grid(color="white", linewidth=1.0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0, labelsize=9, labelcolor="#4D4D4D")
    ax.xaxis.label.set_color("#333333")
    ax.yaxis.label.set_color("#333333")

    meta = (
        f"ENCODE experiment {prov['experiment']} · file {prov['source_file']} · "
        f"{prov['biosample']}\n"
        f"subset: {subset['total_bins']:,} × 25 bp bins re-binned to means; "
        f"{subset['retained_bins']:,} bins retained (fold change ≥ "
        f"{subset['min_value_threshold']}); bundle size "
        f"{prov['subset_size_bytes'] / 1024:.0f} KB\n"
        f"{len(peaks)} pseudo-peaks called (strongest of "
        f"{TOP_PEAKS} kept; shaded above); public ENCODE data, no usage restrictions"
    )
    fig.text(0.10, 0.20, meta, fontsize=9.5, va="top", color="0.25")
    return fig


def _encode_covplot_page(peaks: pd.DataFrame) -> plt.Figure:
    fig, _ax = pcs.covplot(
        peaks,
        weightCol="score",
        title="Pseudo-peak coverage over chr5 (shading weight = fold change)",
        xlab="chr5 position (bundled window)",
        ylab="Weighted coverage",
        figsize=(8.5, 3.4),
    )
    summary = (
        f"pseudo-peaks: {len(peaks)}\n"
        f"median width: {float(np.median(peaks['end'] - peaks['start'])):,.0f} bp\n"
        f"highest mean fold change in a peak: {peaks['score'].max():.1f}×"
    )
    fig.text(0.62, 0.80, summary, fontsize=9.5, family="monospace", va="top", color="0.25")
    return fig


def _summit_profile_page(matrix: np.ndarray) -> plt.Figure:
    fig = plt.figure(figsize=(8.5, 5.0))
    ax = fig.add_axes((0.12, 0.16, 0.85, 0.70))
    pcs.plot_tss_profile(
        matrix,
        window=1000,
        ax=ax,
        title="H3K4me3 around pseudo-peak summits (mean ± SEM)",
        color="#3182bd",
    )
    ax.set_xlabel("Distance to summit (bp)")
    return fig


def main() -> int:
    OUT.mkdir(exist_ok=True)
    pdf_path = OUT / "pychipseeker_report.pdf"

    ann = pcs.annotate_peak(
        peaks=HERE / "example_peaks.bed",
        gtf=HERE / "example.gtf",
        upstream=10,
        downstream=10,
    )
    stat = ann.summary()
    df = ann.to_dataframe()
    print(stat.to_string(index=False))

    # Virtual signal for the signal-based TSS profile/heatmap pages
    # (ChIPseeker accepts plain bedGraph-like DataFrames here).
    signal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [10000, 20000],
            "value": [1.0, 2.0],
        }
    )

    promoters = pcs.getPromoters(TxDb=ann.gene_annotation, upstream=1000, downstream=1000)
    tag = pcs.getTagMatrix(peak=ann.peaks, windows=promoters, nbin=50)

    # ------------------------------ ENCODE chapter ------------------------
    bw_path = example_signal_path()
    prov = example_signal_provenance()
    peaks = _pseudo_peaks(bw_path)
    print(f"pseudo-peaks: {len(peaks)} at fold change >= {PEAK_THRESHOLD}")
    summits = _summit_table(peaks, bw_path)

    flanked = peaks.copy()
    flanked["start"] = (peaks["start"] - 1000).clip(lower=0).astype(int)
    flanked["end"] = (peaks["end"] + 1000).clip(upper=CHROM_SIZE).astype(int)
    peak_tag = pcs.get_tag_matrix_from_bigwig(
        bw_path, flanked, nbin=50, value_type="mean"
    )

    summit_matrix = pcs.extract_tss_matrix(bw_path, summits, window=1000, n_bins=50)

    figures = [
        (_title_page(stat, df), "Title / overview"),
        (ann.plot_anno_pie()[0], "Annotation distribution (pie)"),
        (ann.plot_anno_bar()[0], "Annotation distribution (bar)"),
        (ann.plot_dist_to_tss()[0], "Distance to TSS"),
        (pcs.covplot(ann.peaks, title="Peak coverage over chromosomes")[0], "Peak coverage (covplot)"),
        (ann.plot_tss_profile(signal=signal, window=100, n_bins=50)[0], "TSS avg profile (signal)"),
        (ann.plot_tss_heatmap(signal=signal, window=100, n_bins=50)[0], "TSS heatmap (signal)"),
        (pcs.plot_avg_prof(tag, conf=None)[0], "Promoter avg profile (tagMatrix)"),
        (pcs.tag_heatmap(tag)[0], "Promoter tag heatmap (tagMatrix)"),
        (_encode_track_page(bw_path, prov, peaks), "ENCODE chapter — H3K4me3 chr5 track"),
        (_encode_covplot_page(peaks), "ENCODE chapter — pseudo-peak coverage"),
        (
            pcs.plot_avg_prof(
                peak_tag,
                conf=0.95,
                color="#3182bd",
                label="H3K4me3 (ENCODE)",
                ylab="Mean fold change in window",
                xlab="Pseudo-peak ± 1 kb (5'->3')",
            )[0],
            "ENCODE chapter — profile around pseudo-peaks",
        ),
        (
            pcs.tag_heatmap(
                peak_tag, palette="Reds", title="H3K4me3 around pseudo-peaks"
            )[0],
            "ENCODE chapter — tag heatmap",
        ),
        (
            _summit_profile_page(summit_matrix),
            "ENCODE chapter \u2014 summit-anchored profile",
        ),
        (
            pcs.plot_tss_heatmap(
                summit_matrix, window=1000, cmap="viridis"
            )[0],
            "ENCODE chapter \u2014 summit-anchored heatmap",
        ),
        (_final_page(), "How to reproduce"),
    ]

    with PdfPages(str(pdf_path)) as pdf:
        for fig, name in figures:
            fig.suptitle(name, fontsize=10, color="0.3", y=0.995)
            pdf.savefig(fig, bbox_inches="tight")
            print("page:", name)
            plt.close(fig)

    print(f"PDF written: {pdf_path} ({pdf_path.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
