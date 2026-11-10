"""Command-line interface for pychipseeker."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Optional, Sequence, Union

import pandas as pd

from ._version import __version__
from .annotation import GeneAnnotation, annotate_peak
from .enrichment import enrichr, ora
from .io import read_peaks, write_table
from .plotting import plot_anno_bar, plot_anno_pie, plot_dist_to_tss


def _read_annotation_table(path: Union[str, Path]) -> pd.DataFrame:
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def _save_plot(fig, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")


def _parse_ints(value: Optional[str]) -> Optional[tuple[int, ...]]:
    if value is None:
        return None
    values = tuple(int(x.strip()) for x in value.split(",") if x.strip())
    return values or None


def _cmd_annotate(args: argparse.Namespace) -> int:
    tss_region = None
    if args.tss_region:
        parsed = _parse_ints(args.tss_region)
        if parsed is None or len(parsed) != 2:
            raise SystemExit("--tss-region must be two comma-separated integers")
        tss_region = parsed

    priority = tuple(x.strip() for x in args.priority.split(",") if x.strip())
    promoter_bins = _parse_ints(args.promoter_bins)
    source = args.txdb if args.txdb is not None else args.gtf

    annotation = annotate_peak(
        peaks=args.peaks,
        gtf=source,
        tss_region=tss_region or (-3000, 3000),
        level=args.level,
        overlap=args.overlap,
        same_strand=args.same_strand,
        ignore_overlap=args.ignore_overlap,
        ignore_upstream=args.ignore_upstream,
        ignore_downstream=args.ignore_downstream,
        genomic_annotation_priority=priority,
        downstream_distance=args.downstream_distance,
        add_flank_gene_info=args.flank,
        flank_distance=args.flank_distance,
        upstream=args.upstream,
        downstream=args.downstream,
        promoter_bins=promoter_bins,
        drop_unmapped=not args.keep_unmapped,
        align_annotation_feature=not args.no_align_annotation_feature,
        verbose=not args.quiet,
        n_jobs=args.jobs,
    )
    df = annotation.to_chipseeker() if args.chipseeker_columns else annotation.to_dataframe()

    out = Path(args.output)
    write_table(df, out)
    print(f"Annotated {len(df):,} peaks -> {out}")

    if args.summary:
        summary = annotation.anno_stat()
        summary_path = out.with_name(out.stem + ".annoStat.tsv")
        summary.to_csv(summary_path, sep="\t", index=False)
        print(f"Annotation statistics -> {summary_path}")
        print(summary.to_string(index=False))

    if args.plot_dir:
        plot_dir = Path(args.plot_dir)
        plot_dir.mkdir(parents=True, exist_ok=True)
        fig, _ = annotation.plot_anno_pie()
        _save_plot(fig, plot_dir / "annotation_pie.png", args.dpi)
        fig.clf()
        fig, _ = annotation.plot_anno_bar()
        _save_plot(fig, plot_dir / "annotation_bar.png", args.dpi)
        fig.clf()
        fig, _ = annotation.plot_dist_to_tss()
        _save_plot(fig, plot_dir / "distance_to_tss.png", args.dpi)
        fig.clf()
        print(f"Plots -> {plot_dir}")

    return 0


def _cmd_plot(args: argparse.Namespace) -> int:
    df = _read_annotation_table(args.annotated)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix

    fig, _ = plot_anno_pie(df)
    _save_plot(fig, outdir / f"{prefix}_annotation_pie.png", args.dpi)
    fig.clf()

    fig, _ = plot_anno_bar(df)
    _save_plot(fig, outdir / f"{prefix}_annotation_bar.png", args.dpi)
    fig.clf()

    if "distance_to_tss" in df.columns or "distanceToTSS" in df.columns:
        fig, _ = plot_dist_to_tss(df)
        _save_plot(fig, outdir / f"{prefix}_distance_to_tss.png", args.dpi)
        fig.clf()

    print(f"Plots -> {outdir}")
    return 0


def _read_gene_list(path: Union[str, Path]) -> list[str]:
    path = Path(path)
    if path.suffix.lower() in {".tsv", ".csv"} and path.stat().st_size:
        df = _read_annotation_table(path)
        for column in ("gene", "gene_name", "gene_id", "Gene", "Genes"):
            if column in df.columns:
                return df[column].dropna().astype(str).unique().tolist()
    with open(path, encoding="utf-8", errors="replace") as fh:
        return [line.strip() for line in fh if line.strip()]


def _cmd_enrich(args: argparse.Namespace) -> int:
    df = _read_annotation_table(args.annotated)
    if args.gene_column not in df.columns:
        raise SystemExit(f"gene column not found: {args.gene_column}")
    genes = df[args.gene_column].dropna().astype(str).unique().tolist()
    if not genes:
        raise SystemExit("no genes to enrich")

    background = None if args.background is None else _read_gene_list(args.background)
    try:
        if args.local:
            if not args.gmt:
                raise SystemExit("--local requires --gmt with term/gene sets")
            result = ora(
                genes,
                gene_sets=args.gmt,
                background=background,
                p_adjust=args.p_adjust,
            )
        else:
            result = enrichr(
                genes,
                gene_sets=args.gene_sets,
                organism=args.organism,
                outdir=args.outdir,
                background=background,
                p_adjust=args.p_adjust,
            )
    except (ImportError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    out = Path(args.outdir) / "enrichment.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, sep="\t", index=False)
    print(f"Enrichment -> {out}")
    if args.plot:
        from .enrichment import plot_enrichment

        fig, _ = plot_enrichment(result, top=args.top)
        plot_path = out.parent / "enrichment.png"
        _save_plot(fig, plot_path, args.dpi)
        fig.clf()
        print(f"Enrichment plot -> {plot_path}")
    return 0


def _cmd_callpeak(args: argparse.Namespace) -> int:
    from .peaks import call_peaks, write_peaks

    chroms = None
    if args.chroms:
        chroms = tuple(value.strip() for value in args.chroms.split(",") if value.strip())
    peaks = call_peaks(
        args.bigwig,
        threshold=args.threshold,
        max_gap=args.max_gap,
        min_width=args.min_width,
        chroms=chroms,
        sample=args.sample,
        seed=args.seed,
        n_jobs=args.jobs,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_peaks(peaks, output)
    print(f"Called {len(peaks):,} peaks -> {output}")
    return 0


def _cmd_filter(args: argparse.Namespace) -> int:
    from .peaks import filter_peaks, write_peaks

    chroms = None
    if args.chroms:
        chroms = tuple(value.strip() for value in args.chroms.split(",") if value.strip())
    exclude = None
    if args.exclude_chroms:
        exclude = tuple(
            value.strip() for value in args.exclude_chroms.split(",") if value.strip()
        )
    before = len(read_peaks(args.peaks))
    peaks = filter_peaks(
        args.peaks,
        blacklist=args.blacklist,
        drop_mito=args.drop_mito,
        min_width=args.min_width,
        max_width=args.max_width,
        min_score=args.min_score,
        chroms=chroms,
        exclude_chroms=exclude,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_peaks(peaks, output)
    print(f"Kept {len(peaks):,}/{before:,} peaks -> {output}")
    return 0


def _load_annotation_source(
    source: str, peaks: Optional[Union[str, Path]] = None
) -> GeneAnnotation:
    from .genomes import (
        detect_assembly,
        is_assembly_name,
        is_genome_spec,
        load_annotation,
        load_genome,
    )

    if source.strip().lower() == "auto":
        source = detect_assembly(peaks=peaks)
    if is_assembly_name(source):
        return load_genome(source)
    if is_genome_spec(source):
        return load_annotation(source)
    return GeneAnnotation.from_gtf(source)


def _cmd_qc(args: argparse.Namespace) -> int:
    from .qc import qc_summary

    if args.bigwig is None and args.bam is None:
        raise SystemExit("qc needs --bigwig or --bam for FRiP and TSS metrics")
    genes = None
    source = args.txdb if args.txdb is not None else args.gtf
    if source is not None and not args.no_tss:
        try:
            annotation = _load_annotation_source(source, peaks=args.peaks)
        except (ValueError, FileNotFoundError) as exc:
            raise SystemExit(str(exc)) from exc
        genes = annotation.features("transcript")
    elif source is not None and args.no_tss:
        genes = None

    table = qc_summary(
        args.peaks,
        bigwig=args.bigwig,
        bam=args.bam,
        genes=genes,
        blacklist=args.blacklist,
        tss_upstream=args.tss_upstream,
        tss_downstream=args.tss_downstream,
        tss_sample=args.tss_sample,
        seed=args.seed,
        exclude_duplicates=not args.include_duplicates,
        n_jobs=args.jobs,
    )
    print(table.to_string(index=False))
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(output, sep="\t", index=False)
        print(f"QC -> {output}")
    return 0


def _cmd_report(args: argparse.Namespace) -> int:
    from .report import build_report

    tss_region = None
    if args.tss_region:
        parsed = _parse_ints(args.tss_region)
        if parsed is None or len(parsed) != 2:
            raise SystemExit("--tss-region must be two comma-separated integers")
        tss_region = parsed

    source = args.txdb if args.txdb is not None else args.gtf
    output = build_report(
        peaks=args.peaks,
        gtf=source,
        bigwig=args.bigwig,
        bam=args.bam,
        blacklist=args.blacklist,
        output=args.output,
        level=args.level,
        tss_region=tss_region or (-3000, 3000),
        upstream=args.upstream,
        downstream=args.downstream,
        sample=args.sample,
        seed=args.seed,
        title=args.title,
        include_tag_matrix=not args.no_tag_matrix,
        qc=args.qc,
        n_jobs=args.jobs,
    )
    print(f"Report -> {output} ({output.stat().st_size / 1024:.0f} KiB)")
    return 0


def _cmd_cache(args: argparse.Namespace) -> int:
    from .genomes import (
        cached_genomes,
        clear_genome_cache,
        download_genome,
        genome_provenance,
        prune_genome_cache,
    )

    action = args.cache_command
    if action == "list":
        rows = cached_genomes()
        if args.json:
            print(json.dumps(rows, indent=2, default=str))
            return 0
        if not rows:
            print("No genome bundles found.")
            return 0
        print(
            pd.DataFrame(rows)[["assembly", "source", "size_mb", "path"]].to_string(
                index=False
            )
        )
        return 0

    if action == "info":
        print(json.dumps(genome_provenance(args.assembly), indent=2, default=str))
        return 0

    if action == "prune":
        removed = prune_genome_cache()
        if not removed:
            print("Nothing to prune.")
        for path in removed:
            print(f"removed {path}")
        return 0

    if action == "download":
        for assembly in args.assemblies:
            path = download_genome(assembly)
            print(f"{assembly}: {path} ({path.stat().st_size / 1e6:.1f} MB)")
        return 0

    if action == "clear":
        assemblies = args.assemblies or ([None] if args.all else [])
        if not assemblies:
            raise SystemExit("cache clear needs ASSEMBLY... or --all")
        removed = []
        for assembly in assemblies:
            removed.extend(clear_genome_cache(assembly))
        if not removed:
            print("Nothing to remove.")
        for path in removed:
            print(f"removed {path}")
        return 0

    raise SystemExit(f"unknown cache action: {action}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pychipseeker",
        description="ChIPseeker-style peak annotation and visualization in Python.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command", required=True)

    p_ann = sub.add_parser("annotate", help="annotate peaks against a GTF or TxDb")
    p_ann.add_argument("--peaks", required=True, help="BED/TSV/CSV peak file")
    genome_source = p_ann.add_mutually_exclusive_group(required=True)
    genome_source.add_argument("--gtf", help="GTF/GFF annotation file")
    genome_source.add_argument(
        "--txdb",
        help=(
            "assembly name (hg38, hg19, mm10, mm9, dm6, danRer11, rn7, mm39, "
            "sacCer3), ensembl:/gencode:/url: spec, or alias"
        ),
    )
    p_ann.add_argument("--output", required=True, help="output annotation table")
    p_ann.add_argument("--plot-dir", default=None, help="directory for PNG plots")
    p_ann.add_argument("--summary", action="store_true", help="write ChIPseeker annoStat")
    p_ann.add_argument(
        "--level",
        choices=("transcript", "gene"),
        default="transcript",
        help="annotation level (ChIPseeker default: transcript)",
    )
    p_ann.add_argument(
        "--tss-region",
        default=None,
        help="promoter region, e.g. -3000,3000 (overrides --upstream/--downstream)",
    )
    p_ann.add_argument("--upstream", type=int, default=None, help="promoter upstream bp")
    p_ann.add_argument(
        "--downstream", type=int, default=None, help="promoter downstream bp"
    )
    p_ann.add_argument(
        "--promoter-bins",
        default=None,
        help="optional comma-separated promoter cutoffs; default follows tss-region",
    )
    p_ann.add_argument(
        "--overlap",
        choices=("TSS", "all"),
        default="TSS",
        help="overlap mode for nearest-feature selection",
    )
    p_ann.add_argument("--same-strand", action="store_true")
    p_ann.add_argument("--ignore-overlap", action="store_true")
    p_ann.add_argument("--ignore-upstream", action="store_true")
    p_ann.add_argument("--ignore-downstream", action="store_true")
    p_ann.add_argument(
        "--priority",
        default="Promoter,5UTR,3UTR,Exon,Intron,Downstream,Intergenic",
        help="genomic annotation priority",
    )
    p_ann.add_argument(
        "--downstream-distance",
        type=int,
        default=300,
        help="distance used for the Downstream category",
    )
    p_ann.add_argument("--flank", action="store_true", help="add flanking gene columns")
    p_ann.add_argument("--flank-distance", type=int, default=5000)
    p_ann.add_argument(
        "--keep-unmapped",
        action="store_true",
        help="keep peaks on chromosomes absent from the annotation",
    )
    p_ann.add_argument(
        "--chipseeker-columns",
        action="store_true",
        help="write ChIPseeker-style column names (geneId, distanceToTSS, ...)",
    )
    p_ann.add_argument(
        "--no-align-annotation-feature",
        action="store_true",
        help=(
            "reproduce ChIPseeker <=1.48 isoform metadata behavior instead of "
            "aligning gene/transcript fields to the exon/intron hit"
        ),
    )
    p_ann.add_argument("--dpi", type=int, default=300)
    p_ann.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="worker processes for large peak sets (-1 = all cores)",
    )
    p_ann.add_argument("--quiet", action="store_true")
    p_ann.set_defaults(func=_cmd_annotate)

    p_plot = sub.add_parser("plot", help="plot an existing annotation table")
    p_plot.add_argument("--annotated", required=True, help="annotation table from `annotate`")
    p_plot.add_argument("--outdir", required=True, help="output directory")
    p_plot.add_argument("--prefix", default="peak")
    p_plot.add_argument("--dpi", type=int, default=300)
    p_plot.set_defaults(func=_cmd_plot)

    p_enrich = sub.add_parser("enrich", help="Enrichr analysis on annotated genes")
    p_enrich.add_argument("--annotated", required=True)
    p_enrich.add_argument("--outdir", default="enrichment")
    p_enrich.add_argument("--gene-column", default="gene_name")
    p_enrich.add_argument("--gene-sets", default="GO_Biological_Process_2023")
    p_enrich.add_argument("--organism", default="human")
    p_enrich.add_argument(
        "--background",
        default=None,
        help="background gene list (one per line, or a gene table)",
    )
    p_enrich.add_argument(
        "--gmt",
        default=None,
        help="GMT file for offline over-representation analysis (--local)",
    )
    p_enrich.add_argument(
        "--local",
        action="store_true",
        help="run the built-in hypergeometric ORA instead of Enrichr",
    )
    p_enrich.add_argument(
        "--p-adjust",
        default="BH",
        help="p-value adjustment: BH (default), Bonferroni or none",
    )
    p_enrich.add_argument(
        "--plot",
        action="store_true",
        help="save a -log10(adjusted p) bar chart beside the table",
    )
    p_enrich.add_argument("--top", type=int, default=20, help="terms shown in the plot")
    p_enrich.add_argument("--dpi", type=int, default=300)
    p_enrich.set_defaults(func=_cmd_enrich)

    p_call = sub.add_parser(
        "callpeak", help="call broad peaks from a bigWig signal track"
    )
    p_call.add_argument("--bigwig", required=True, help="signal bigWig")
    p_call.add_argument("--output", required=True, help="output BED6+ file")
    p_call.add_argument("--threshold", type=float, default=3.0)
    p_call.add_argument("--max-gap", type=int, default=1000)
    p_call.add_argument("--min-width", type=int, default=500)
    p_call.add_argument(
        "--chroms", default=None, help="comma-separated chromosome subset"
    )
    p_call.add_argument("--sample", type=int, default=None)
    p_call.add_argument("--seed", type=int, default=None)
    p_call.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="worker processes, one task per chromosome (-1 = all cores)",
    )
    p_call.set_defaults(func=_cmd_callpeak)

    p_filter = sub.add_parser("filter", help="filter peaks by blacklist/width/score")
    p_filter.add_argument("--peaks", required=True)
    p_filter.add_argument("--output", required=True)
    p_filter.add_argument("--blacklist", default=None, help="BED blacklist intervals")
    p_filter.add_argument("--drop-mito", action="store_true")
    p_filter.add_argument("--min-width", type=int, default=None)
    p_filter.add_argument("--max-width", type=int, default=None)
    p_filter.add_argument("--min-score", type=float, default=None)
    p_filter.add_argument("--chroms", default=None, help="keep these chromosomes")
    p_filter.add_argument("--exclude-chroms", default=None)
    p_filter.set_defaults(func=_cmd_filter)

    p_qc = sub.add_parser("qc", help="peak, FRiP and TSS quality-control metrics")
    p_qc.add_argument("--peaks", required=True)
    p_qc.add_argument("--bigwig", default=None)
    p_qc.add_argument("--bam", default=None)
    qc_source = p_qc.add_mutually_exclusive_group(required=False)
    qc_source.add_argument("--gtf", default=None)
    qc_source.add_argument("--txdb", default=None)
    p_qc.add_argument("--blacklist", default=None)
    p_qc.add_argument("--output", default=None, help="optional TSV output")
    p_qc.add_argument("--tss-upstream", type=int, default=2000)
    p_qc.add_argument("--tss-downstream", type=int, default=2000)
    p_qc.add_argument("--tss-sample", type=int, default=2000)
    p_qc.add_argument("--seed", type=int, default=20261026)
    p_qc.add_argument("--no-tss", action="store_true")
    p_qc.add_argument(
        "--include-duplicates",
        action="store_true",
        help="count duplicate BAM reads in FRiP",
    )
    p_qc.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="worker processes for the BAM TSS matrix (-1 = all cores)",
    )
    p_qc.set_defaults(func=_cmd_qc)

    p_report = sub.add_parser(
        "report", help="build a multi-page PDF annotation report"
    )
    p_report.add_argument("--peaks", required=True, help="BED/TSV/CSV peak file")
    report_source = p_report.add_mutually_exclusive_group(required=True)
    report_source.add_argument("--gtf", help="GTF/GFF annotation file")
    report_source.add_argument("--txdb", help="bundled assembly name")
    p_report.add_argument("--bigwig", default=None, help="optional bigWig signal track")
    p_report.add_argument("--bam", default=None, help="optional BAM for FRiP/QC")
    p_report.add_argument("--blacklist", default=None, help="BED blacklist for QC")
    p_report.add_argument(
        "--qc", action="store_true", help="add a quality-control metrics page"
    )
    p_report.add_argument("--output", default="pychipseeker_report.pdf")
    p_report.add_argument(
        "--level", choices=("transcript", "gene"), default="transcript"
    )
    p_report.add_argument(
        "--tss-region", default=None, help="promoter region, e.g. -3000,3000"
    )
    p_report.add_argument(
        "--upstream", type=int, default=None, help="asymmetric signal window upstream"
    )
    p_report.add_argument(
        "--downstream", type=int, default=None, help="asymmetric signal window downstream"
    )
    p_report.add_argument(
        "--sample",
        type=int,
        default=3000,
        help="transcripts sampled for bigWig signal panels",
    )
    p_report.add_argument("--seed", type=int, default=20261026)
    p_report.add_argument("--title", default=None)
    p_report.add_argument(
        "--no-tag-matrix",
        action="store_true",
        help="skip the tagMatrix profile/heatmap pages",
    )
    p_report.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="worker processes for annotation/QC (-1 = all cores)",
    )
    p_report.set_defaults(func=_cmd_report)

    p_cache = sub.add_parser("cache", help="manage downloaded genome bundles")
    cache_sub = p_cache.add_subparsers(dest="cache_command", required=True)
    cache_list = cache_sub.add_parser("list", help="list installed and cached bundles")
    cache_list.add_argument("--json", action="store_true", help="machine-readable output")
    cache_info = cache_sub.add_parser(
        "info", help="show provenance for an assembly or download spec"
    )
    cache_info.add_argument("assembly")
    cache_sub.add_parser("prune", help="remove partial and orphaned cache files")
    cache_download = cache_sub.add_parser("download", help="download genome bundles")
    cache_download.add_argument("assemblies", nargs="+")
    cache_clear = cache_sub.add_parser("clear", help="remove cached bundles")
    cache_clear.add_argument("assemblies", nargs="*")
    cache_clear.add_argument(
        "--all", action="store_true", help="remove every cached bundle"
    )
    p_cache.set_defaults(func=_cmd_cache)

    _allow_negative_values(parser)
    _allow_negative_values(p_ann)
    _allow_negative_values(p_plot)
    _allow_negative_values(p_enrich)
    _allow_negative_values(p_call)
    _allow_negative_values(p_filter)
    _allow_negative_values(p_qc)
    _allow_negative_values(p_report)
    _allow_negative_values(p_cache)

    return parser


def _allow_negative_values(parser: argparse.ArgumentParser) -> None:
    """Teach argparse that comma-separated negative values are values.

    Stock argparse only treats a token like ``-3000`` as a value; the
    otherwise-valid ``--tss-region -3000,3000`` is misread as an unknown
    option ("expected one argument").  Extending the negative-number
    matcher to comma/dot forms fixes ``--tss-region -3000,3000``,
    ``-3000,-500`` and ``--tss-region=-3000.5,3000.5`` alike.
    """

    parser._negative_number_matcher = re.compile(
        r"^-\d+(?:\.\d+)?(?:,-?\d+(?:\.\d+)?)*$"
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
