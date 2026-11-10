"""Optional enrichment helpers.

The statistical work is delegated to gseapy/Enrichr when the optional
dependency is installed; this module also provides a dependency-free local
over-representation analysis (ORA) with hypergeometric tests and BH/Bonferroni
correction.  The local path is useful when no network access is available or
when a custom GMT file should be tested against an explicit background.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence, Union

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


def _clean_genes(gene_list: Iterable[str]) -> list[str]:
    genes = [str(g) for g in gene_list if g is not None and str(g) != ""]
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(genes))


def _p_adjust(pvalues: Sequence[float], method: str = "BH") -> list[float]:
    """Benjamini-Hochberg or Bonferroni adjustment, without SciPy."""

    values = [float(p) for p in pvalues]
    n = len(values)
    if n == 0:
        return []
    key = str(method).strip().lower()
    if key in {"none", "no", "raw"}:
        return values
    if key in {"bonferroni", "bonf"}:
        return [min(1.0, p * n) for p in values]
    if key not in {"bh", "fdr", "benjamini-hochberg"}:
        raise ValueError("p_adjust must be one of BH, Bonferroni or none")

    order = sorted(range(n), key=lambda i: values[i])
    adjusted = [1.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, values[index] * n / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """P(X >= k) for a hypergeometric distribution.

    ``M`` is the background size, ``n`` the query size and ``N`` the gene-set
    size.  The sum is evaluated in log space so that large backgrounds remain
    numerically stable.
    """

    if k <= 0:
        return 1.0
    upper = min(n, N)
    if k > upper:
        return 0.0
    terms = []
    for i in range(int(k), int(upper) + 1):
        terms.append(
            math.exp(
                _log_comb(n, i)
                + _log_comb(M - n, N - i)
                - _log_comb(M, N)
            )
        )
    return float(min(max(math.fsum(terms), 0.0), 1.0))


def read_gmt(path: Union[str, Path]) -> dict[str, list[str]]:
    """Read a GMT file into ``{term: [genes]}``."""

    gene_sets: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            term = fields[0].strip()
            if not term:
                continue
            genes = [gene.strip() for gene in fields[2:] if gene.strip()]
            if genes:
                gene_sets[term] = genes
    return gene_sets


def _as_gene_sets(
    gene_sets: Union[Mapping[str, Iterable[str]], str, Path],
) -> dict[str, set[str]]:
    if isinstance(gene_sets, (str, Path)):
        mapping: Mapping[str, Iterable[str]] = read_gmt(gene_sets)
    elif isinstance(gene_sets, Mapping):
        mapping = gene_sets
    else:
        raise TypeError("gene_sets must be a dict, GMT path or mapping")
    return {
        str(term): {str(gene) for gene in genes if gene is not None}
        for term, genes in mapping.items()
    }


def ora(
    gene_list: Iterable[str],
    gene_sets: Union[Mapping[str, Iterable[str]], str, Path],
    background: Optional[Iterable[str]] = None,
    p_adjust: str = "BH",
    min_overlap: int = 1,
) -> pd.DataFrame:
    """Local over-representation analysis with a hypergeometric test.

    Parameters
    ----------
    gene_list:
        Query genes (for example genes linked to peaks).
    gene_sets:
        Mapping ``{term: genes}`` or a GMT file path.
    background:
        Universe of genes.  Defaults to the union of all gene sets, which is
        the conservative choice when the assay only covers those genes.
    p_adjust:
        ``"BH"``/``"fdr"`` (default), ``"Bonferroni"`` or ``"none"``.
    min_overlap:
        Terms with fewer overlapping genes are dropped.

    Returns
    -------
    pandas.DataFrame
        Sorted by p-value with ``PValue``, ``AdjustedPValue``, effect size
        (``FoldEnrichment``) and the overlapping ``Genes``.
    """

    genes = _clean_genes(gene_list)
    if not genes:
        raise ValueError("gene_list is empty")
    sets = _as_gene_sets(gene_sets)
    if not sets:
        raise ValueError("gene_sets is empty")

    query = set(genes)
    if background is None:
        universe = set().union(*sets.values())
    else:
        universe = {str(gene) for gene in background if gene is not None}
        if not universe:
            raise ValueError("background is empty")
    query &= universe
    if not query:
        raise ValueError("no query genes are present in the background")

    M = len(universe)
    n = len(query)
    rows: list[dict[str, object]] = []
    for term, members in sets.items():
        term_genes = members & universe
        overlap = query & term_genes
        if len(overlap) < int(min_overlap):
            continue
        N = len(term_genes)
        k = len(overlap)
        expected = n * N / M if M else 0.0
        fold = (k / n) / (N / M) if N and M else float("nan")
        rows.append(
            {
                "Term": term,
                "Overlap": k,
                "GeneSetSize": N,
                "BackgroundSize": M,
                "QuerySize": n,
                "Expected": expected,
                "FoldEnrichment": fold,
                "PValue": _hypergeom_sf(k, M, n, N),
                "Genes": ";".join(sorted(overlap)),
            }
        )

    if not rows:
        raise ValueError("no gene set overlapped the query genes")
    table = pd.DataFrame.from_records(rows)
    table["AdjustedPValue"] = _p_adjust(table["PValue"].tolist(), p_adjust)
    return table.sort_values(
        ["PValue", "Term"], ascending=[True, True]
    ).reset_index(drop=True)


def enrichr(
    gene_list: Iterable[str],
    gene_sets: Union[str, Sequence[str]] = "GO_Biological_Process_2023",
    organism: str = "human",
    outdir: Optional[Union[str, Path]] = None,
    background: Optional[Iterable[str]] = None,
    p_adjust: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """Run Enrichr through gseapy and return the result table.

    ``background`` is forwarded to Enrichr (a custom background list is only
    supported by recent gseapy versions).  Enrichr already reports its own
    adjusted p-values; passing ``p_adjust="BH"`` recomputes them across the
    returned terms with the local implementation, which makes the output
    comparable to :func:`ora`.
    """

    genes = _clean_genes(gene_list)
    if not genes:
        raise ValueError("gene_list is empty")

    gp = _require_gseapy()
    call_kwargs = dict(kwargs)
    if background is not None:
        background_genes = _clean_genes(background)
        if not background_genes:
            raise ValueError("background is empty")
        call_kwargs["background"] = background_genes
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism=organism,
        outdir=None if outdir is None else str(outdir),
        **call_kwargs,
    )
    results = enr.results
    if p_adjust is not None and "P-value" in results.columns and not results.empty:
        results = results.copy()
        results["Adjusted P-value"] = _p_adjust(
            results["P-value"].tolist(), p_adjust
        )
    return results


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

    For offline analysis pass ``local=True`` to use :func:`ora` instead of
    Enrichr; all remaining keyword arguments are forwarded accordingly.
    """

    local = bool(kwargs.pop("local", False))
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
    if local:
        return ora(genes, **kwargs)
    return enrichr(genes, **kwargs)
