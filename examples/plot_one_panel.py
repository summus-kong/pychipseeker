#!/usr/bin/env python3
"""Render one gallery panel for any peak file and save it as PDF (or PNG/SVG).

This script gives you the same 18 panel designs as the ENCODE gallery without
running the whole gallery.  Pick a letter from ``--list`` and point it at your
own peaks and genome annotation.

Examples
--------
List the panels::

    python examples/plot_one_panel.py --list

Render one panel as a vector PDF::

    python examples/plot_one_panel.py \
        --panel K --peaks my_peaks.bed --txdb mm10 \
        --outdir figures

Render the same panel as PDF and PNG::

    python examples/plot_one_panel.py \
        --panel A --peaks my_peaks.bed --txdb mm10 --formats pdf,png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt

import pychipseeker as pcs

try:  # running as a script: examples/ is already on sys.path
    import make_encode_plate as gallery
except ImportError:  # imported as a module from the repository root
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import make_encode_plate as gallery

REPO = Path(__file__).resolve().parents[1]
DEFAULT_PEAKS = REPO / "examples" / "encode_peaks.bed"


def _build_figure(letter, ann, tag, df, profiles, upstream, downstream):
    """Create the standalone matplotlib figure for one catalogue panel."""

    if letter == "A":
        return gallery._figure_pie(ann)
    if letter == "B":
        return gallery._figure_width(df)
    if letter == "C":
        return gallery._figure_distance(df)
    if letter == "D":
        return gallery._figure_profile(tag)
    if letter == "E":
        return gallery._figure_heatmap(tag)
    if letter == "F":
        return gallery._figure_table(ann, len(ann.peaks))
    if letter == "G":
        return gallery._figure_chrom(df)
    if letter == "H":
        return gallery._figure_width_box(df)
    if letter == "I":
        return gallery._figure_signal_violin(df)
    if letter == "J":
        return gallery._figure_nonpromoter_distance(df)
    if letter == "K":
        return gallery._figure_profile_quartiles(profiles, upstream, downstream)
    if letter == "L":
        return gallery._figure_width_signal(df)
    if letter == "M":
        return gallery._figure_width_ecdf(df)
    if letter == "N":
        return gallery._figure_top_genes(df)
    if letter == "O":
        return gallery._figure_peaks_per_gene(df)
    if letter == "P":
        return gallery._figure_tss_quartiles(tag, ann)
    if letter == "Q":
        return gallery._figure_width_distance(df)
    if letter == "R":
        return gallery._figure_gene_type(df)
    raise ValueError(f"unknown panel letter: {letter}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", default="K", help="panel letter A-R")
    parser.add_argument("--peaks", default=str(DEFAULT_PEAKS), help="peak file")
    parser.add_argument("--txdb", default="mm10", help="assembly name, GTF or TxDb spec")
    parser.add_argument(
        "--outdir",
        default=str(REPO / "examples" / "output" / "one-panel"),
        help="directory for the rendered figure",
    )
    parser.add_argument(
        "--formats",
        default="pdf",
        help="comma-separated formats: pdf, png, svg, eps, tif",
    )
    parser.add_argument("--upstream", type=int, default=3000)
    parser.add_argument("--downstream", type=int, default=3000)
    parser.add_argument("--list", action="store_true", help="print panels and exit")
    args = parser.parse_args()

    if args.list:
        print("letter  title")
        for letter, _stem, title, _blurb in gallery.PANEL_CATALOG:
            print(f"{letter:>6}  {title}")
        return 0

    letter = args.panel.strip().upper()
    if letter not in gallery.PANEL_BY_LETTER:
        by_stem = {item[1].lower(): item[0] for item in gallery.PANEL_CATALOG}
        letter = by_stem.get(args.panel.strip().lower(), "")
    if not letter:
        parser.error(f"unknown panel {args.panel!r}; use --list")

    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    allowed = {"pdf", "png", "svg", "eps", "tif"}
    unknown = sorted(set(formats) - allowed)
    if unknown:
        parser.error(f"unsupported format(s): {', '.join(unknown)}")
    if not formats:
        parser.error("at least one format is required")

    ann = pcs.annotate_peak(peaks=args.peaks, TxDb=args.txdb, verbose=False)
    df = ann.to_dataframe(include_peak_metadata=True).merge(
        ann.peaks[["peak_id", "score"]],
        on="peak_id",
        how="left",
    )

    tag = None
    profiles = []
    if letter in {"D", "E", "K", "P"}:
        promoters = pcs.get_promoters(
            TxDb=ann.gene_annotation,
            upstream=args.upstream,
            downstream=args.downstream,
        )
        tag = pcs.get_tag_matrix(
            peak=ann.peaks,
            windows=promoters,
            nbin=120,
        )
    if letter == "K":
        profiles = gallery._quartile_profiles(tag, df)

    fig = _build_figure(
        letter,
        ann,
        tag,
        df,
        profiles,
        args.upstream,
        args.downstream,
    )

    letter, stem, title, _blurb = gallery.PANEL_BY_LETTER[letter]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        path = outdir / f"{letter}_{stem}.{fmt}"
        fig.savefig(path, bbox_inches="tight")
        written.append(path)
    plt.close(fig)

    print(f"panel : {letter} - {title}")
    print(f"peaks : {len(ann.peaks):,}")
    for path in written:
        print(f"        {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
