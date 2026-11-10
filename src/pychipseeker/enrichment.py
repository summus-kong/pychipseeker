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

import numpy as np
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
        odds_ratio, or_lower, or_upper = _odds_ratio_ci(k, n - k, N - k, M - n - N + k)
        rows.append(
            {
                "Term": term,
                "Overlap": k,
                "GeneSetSize": N,
                "BackgroundSize": M,
                "QuerySize": n,
                "Expected": expected,
                "FoldEnrichment": fold,
                "OddsRatio": odds_ratio,
                "ORLower": or_lower,
                "ORUpper": or_upper,
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


def _odds_ratio_ci(
    a: float,
    b: float,
    c: float,
    d: float,
    *,
    confidence: float = 0.95,
) -> tuple[float, float, float]:
    """Odds ratio and Woolfe log-normal confidence interval for a 2x2 table.

    ``a``/``b`` are query genes in/outside the term; ``c``/``d`` are the
    remaining background genes in/outside the term.  A Haldane-Anscombe
    correction of 0.5 is applied to zero cells so the estimate stays finite.
    """

    cells = [float(a), float(b), float(c), float(d)]
    if min(cells) < 0:
        return (float("nan"), float("nan"), float("nan"))
    if any(value == 0 for value in cells):
        cells = [value + 0.5 for value in cells]
    a, b, c, d = cells
    odds = (a * d) / (b * c) if b and c else float("nan")
    if not np.isfinite(odds) or odds <= 0:
        return (odds, float("nan"), float("nan"))
    z = float(_normal_quantile(0.5 + confidence / 2.0))
    se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    return (odds, math.exp(math.log(odds) - z * se), math.exp(math.log(odds) + z * se))


def _normal_quantile(probability: float) -> float:
    """Inverse standard-normal CDF (Acklam's rational approximation)."""

    p = float(probability)
    if not 0.0 < p < 1.0:
        raise ValueError("probability must be between 0 and 1")
    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    ]
    plow = 0.02425
    phigh = 1.0 - plow
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )
    q = p - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    )


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


def _as_ranked_frame(
    gene_list: Union[Mapping[str, float], pd.DataFrame, Sequence[str]],
) -> pd.DataFrame:
    """Coerce a ranked gene list to a two-column ``gene``/``score`` frame."""

    if isinstance(gene_list, pd.DataFrame):
        frame = gene_list.copy()
        if frame.shape[1] < 2:
            raise ValueError(
                "a ranked DataFrame needs at least two columns (gene, score)"
            )
        lower = {str(column).strip().lower(): column for column in frame.columns}
        gene_col = lower.get("gene") or lower.get("gene_name") or frame.columns[0]
        score_col = lower.get("score") or lower.get("rank") or frame.columns[1]
        frame = frame[[gene_col, score_col]].rename(
            columns={gene_col: "gene", score_col: "score"}
        )
    elif isinstance(gene_list, Mapping):
        frame = pd.DataFrame(
            {"gene": [str(gene) for gene in gene_list], "score": list(gene_list.values())}
        )
    else:
        genes = _clean_genes(gene_list)
        # First gene gets the highest rank, matching a descending sort.
        frame = pd.DataFrame(
            {"gene": genes, "score": np.arange(len(genes), 0, -1, dtype=float)}
        )

    frame["gene"] = frame["gene"].astype(str)
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame = frame.dropna(subset=["gene", "score"])
    frame = frame[frame["gene"] != ""]
    frame = frame.drop_duplicates("gene", keep="first")
    if frame.empty:
        raise ValueError("ranked gene list is empty")
    return frame.sort_values("score", ascending=False).reset_index(drop=True)


def gsea(
    gene_list: Union[Mapping[str, float], pd.DataFrame, Sequence[str]],
    gene_sets: Union[str, Sequence[str], Mapping[str, Iterable[str]]],
    *,
    min_size: int = 15,
    max_size: int = 500,
    permutation_num: int = 1000,
    seed: int = 20261026,
    outdir: Optional[Union[str, Path]] = None,
    **kwargs,
) -> pd.DataFrame:
    """Preranked Gene Set Enrichment Analysis through gseapy.

    ``gene_list`` may be a ``{gene: score}`` mapping, a two-column DataFrame,
    or a plain gene sequence (ranked by order, highest first).  ``gene_sets``
    accepts an Enrichr library name, a GMT path or a ``{term: genes}`` mapping.

    This is a thin wrapper: it requires the optional ``gseapy`` dependency and
    does not re-implement the permutation test.
    """

    frame = _as_ranked_frame(gene_list)
    gp = _require_gseapy()
    call_kwargs = dict(kwargs)
    call_kwargs.setdefault("min_size", int(min_size))
    call_kwargs.setdefault("max_size", int(max_size))
    call_kwargs.setdefault("permutation_num", int(permutation_num))
    call_kwargs.setdefault("seed", int(seed))
    call_kwargs.setdefault("verbose", False)
    if isinstance(gene_sets, (str, Path)):
        sets_arg: Union[str, list[str], dict[str, list[str]]] = str(gene_sets)
    elif isinstance(gene_sets, Mapping):
        sets_arg = {str(term): _clean_genes(genes) for term, genes in gene_sets.items()}
    else:
        sets_arg = list(gene_sets)

    result = gp.prerank(
        rnk=frame,
        gene_sets=sets_arg,
        outdir=None if outdir is None else str(outdir),
        **call_kwargs,
    )
    table = getattr(result, "res2d", None)
    if table is None:
        table = getattr(result, "results", None)
    if table is None:  # pragma: no cover - gseapy API drift
        raise RuntimeError("gseapy returned no GSEA result table")
    return pd.DataFrame(table).reset_index(drop=True)


def gsea_peak_genes(
    annotation,
    ranking: Union[Mapping[str, float], pd.DataFrame, Sequence[str]],
    *,
    gene_column: str = "gene_name",
    upstream: Optional[int] = None,
    downstream: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Restrict a ranked gene list to peak-linked genes and run GSEA."""

    frame = _as_ranked_frame(ranking)
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
    linked = set(links[gene_column].dropna().astype(str))
    subset = frame[frame["gene"].isin(linked)].reset_index(drop=True)
    if subset.empty:
        raise ValueError("no ranked genes are linked to the peaks")
    return gsea(subset, **kwargs)


def _enrichment_plot_frame(result: pd.DataFrame) -> pd.DataFrame:
    """Normalise ORA/Enrichr/GSEA tables to common plotting columns."""

    if not isinstance(result, pd.DataFrame):
        result = pd.DataFrame(result)
    if result.empty:
        raise ValueError("enrichment result is empty")

    columns = set(result.columns)
    if "Term" in columns:
        term = result["Term"].astype(str)
    elif result.index.name:
        term = result.index.astype(str)
    else:
        term = pd.Series([str(i) for i in result.index])

    for name in ("AdjustedPValue", "Adjusted P-value", "padj", "FDR q-val"):
        if name in columns:
            padj = pd.to_numeric(result[name], errors="coerce")
            break
    else:
        for name in ("PValue", "P-value", "pvalue", "NOM p-val"):
            if name in columns:
                padj = pd.to_numeric(result[name], errors="coerce")
                break
        else:
            raise ValueError("result has no p-value column")

    count: Optional[pd.Series] = None
    ratio: Optional[pd.Series] = None
    if "Overlap" in columns:
        raw = result["Overlap"]
        if pd.api.types.is_numeric_dtype(raw):
            count = pd.to_numeric(raw, errors="coerce")
        else:
            text = raw.astype(str)
            numerator, _, denominator = text.str.partition("/")
            count = pd.to_numeric(numerator, errors="coerce")
            denominator_values = pd.to_numeric(denominator, errors="coerce")
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio = count / denominator_values
    if "GeneRatio" in columns:
        text = result["GeneRatio"].astype(str)
        numerator, _, denominator = text.str.partition("/")
        count = pd.to_numeric(numerator, errors="coerce")
        denominator_values = pd.to_numeric(denominator, errors="coerce")
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = count / denominator_values

    if "GeneSetSize" in columns:
        denominator_values = pd.to_numeric(result["GeneSetSize"], errors="coerce")
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = count / denominator_values if count is not None else ratio

    frame = pd.DataFrame(
        {
            "Term": term.to_numpy(),
            "padj": padj.to_numpy(dtype=float),
            "count": (count if count is not None else pd.Series(np.ones(len(result)))).to_numpy(
                dtype=float
            ),
            "ratio": (
                ratio if ratio is not None else pd.Series(np.full(len(result), np.nan))
            ).to_numpy(dtype=float),
        }
    )
    frame = frame.dropna(subset=["padj"])
    frame = frame[np.isfinite(frame["padj"])]
    if frame.empty:
        raise ValueError("no finite p-values to plot")
    frame["padj"] = frame["padj"].clip(lower=np.finfo(float).tiny)
    frame["neglog10"] = -np.log10(frame["padj"])
    return frame.sort_values("neglog10", ascending=True).reset_index(drop=True)


def plot_enrichment(
    result: pd.DataFrame,
    top: int = 20,
    *,
    kind: str = "bar",
    ax=None,
    title: Optional[str] = None,
    figsize: tuple[float, float] = (7.0, 5.0),
    cmap: str = "Blues",
    **kwargs,
):
    """Plot an ORA/Enrichr/GSEA result table.

    ``kind="bar"`` draws a horizontal ``-log10(adjusted p)`` bar chart;
    ``kind="dot"`` draws the clusterProfiler-style dot plot with gene-ratio
    encoded by dot size and significance by colour.  The styling follows the
    same ggplot2 ``theme_grey`` defaults as the rest of pychipseeker.
    """

    if kind not in {"bar", "dot"}:
        raise ValueError("kind must be 'bar' or 'dot'")
    frame = _enrichment_plot_frame(result)
    frame = frame.tail(max(1, int(top))).reset_index(drop=True)

    from matplotlib import colormaps

    from .plotting import _GG_PANEL, _GG_TEXT, _style_axes, _title

    plt = __import__("matplotlib.pyplot", fromlist=["pyplot"])
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    y = np.arange(len(frame))
    cmap_obj = colormaps.get_cmap(cmap)
    values = frame["neglog10"].to_numpy(dtype=float)
    denom = float(values.max()) if values.size and values.max() > 0 else 1.0
    colors = [cmap_obj(0.25 + 0.75 * value / denom) for value in values]

    if kind == "bar":
        ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.4, **kwargs)
        ax.set_xlabel("-log10(adjusted p-value)")
        xmax = float(values.max()) if values.size else 1.0
        ax.set_xlim(0.0, xmax * 1.05 if xmax > 0 else 1.0)
    else:
        sizes = frame["count"].to_numpy(dtype=float)
        finite = sizes[np.isfinite(sizes)]
        if finite.size:
            sizes = np.where(np.isfinite(sizes), sizes, finite.min())
            scale_low, scale_high = float(finite.min()), float(finite.max())
            span = max(scale_high - scale_low, 1.0)
            point_sizes = 30.0 + 260.0 * (sizes - scale_low) / span
        else:  # pragma: no cover - defensive
            point_sizes = np.full(len(frame), 80.0)
        ax.scatter(
            values,
            y,
            s=point_sizes,
            c=colors,
            edgecolor="white",
            linewidth=0.5,
            **kwargs,
        )
        ax.set_xlabel("-log10(adjusted p-value)")
        xmax = float(values.max()) if values.size else 1.0
        ax.set_xlim(0.0, xmax * 1.1 if xmax > 0 else 1.0)

    ax.set_yticks(y)
    ax.set_yticklabels(frame["Term"].tolist(), fontsize=8)
    ax.tick_params(colors=_GG_TEXT, labelsize=8)
    ax.set_facecolor(_GG_PANEL)
    _style_axes(ax, grid_axis="both")
    _title(ax, title or ("Enrichment dot plot" if kind == "dot" else "Enrichment"))
    fig.tight_layout()
    return fig, ax
