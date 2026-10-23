"""Optional enrichment helpers.

The statistical work is delegated to gseapy/Enrichr; this module only provides
small convenience wrappers so that peak-to-gene results can be passed directly
to enrichment analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import pandas as pd


def _require_gseapy():
    try:
        import gseapy as gp
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "gseapy is required for enrichment analysis; install with "
            "`pip install pychipseeker[enrich]`"
        ) from exc
    return gp


def enrichr(
    gene_list: Iterable[str],
    gene_sets: Union[str, Sequence[str]] = "GO_Biological_Process_2023",
    organism: str = "human",
    outdir: Optional[Union[str, Path]] = None,
    **kwargs,
) -> pd.DataFrame:
    """Run Enrichr through gseapy and return the result table."""

    genes = [str(g) for g in gene_list if g is not None and str(g) != ""]
    # Preserve order while removing duplicates.
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("gene_list is empty")

    gp = _require_gseapy()
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism=organism,
        outdir=None if outdir is None else str(outdir),
        **kwargs,
    )
    return enr.results


def enrich_go(
    gene_list: Iterable[str],
    organism: str = "human",
    outdir: Optional[Union[str, Path]] = None,
    **kwargs,
) -> pd.DataFrame:
    """GO Biological Process enrichment via Enrichr."""

    return enrichr(
        gene_list,
        gene_sets="GO_Biological_Process_2023",
        organism=organism,
        outdir=outdir,
        **kwargs,
    )


def enrich_kegg(
    gene_list: Iterable[str],
    organism: str = "human",
    outdir: Optional[Union[str, Path]] = None,
    **kwargs,
) -> pd.DataFrame:
    """KEGG pathway enrichment via Enrichr."""

    gene_sets = "KEGG_2021_Human" if organism.lower() == "human" else "KEGG_2019_Mouse"
    return enrichr(
        gene_list,
        gene_sets=gene_sets,
        organism=organism,
        outdir=outdir,
        **kwargs,
    )


def enrich_peak_genes(
    annotation,
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
    include_gene_body: bool = True,
    gene_column: str = "gene_name",
    **kwargs,
) -> pd.DataFrame:
    """Map peaks to genes with ``seq2gene`` and run Enrichr.

    ``include_gene_body`` is retained for backwards compatibility; host-gene
    (exon/intron) links are always included by ``seq2gene_links``.
    """

    links = annotation.seq2gene_links(
        tss_region=(
            None
            if upstream is None and downstream is None
            else (-int(upstream or 0), int(downstream or 0))
        ),
    )
    if links.empty:
        raise ValueError("no peak-gene links were found")
    if gene_column not in links.columns:
        raise KeyError(f"gene column not found: {gene_column}")
    genes = links[gene_column].dropna().astype(str).unique().tolist()
    return enrichr(genes, **kwargs)
