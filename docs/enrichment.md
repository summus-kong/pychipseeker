# Enrichment

Three backends are available.

## Enrichr (online)

Requires the `enrich` extra (`gseapy`).

```python
import pychipseeker as pcs

genes = ann.to_chipseeker()["geneName"].dropna().unique().tolist()
res = pcs.enrichr(
    genes,
    gene_sets="GO_Biological_Process_2023",
    organism="mouse",
    background=background_genes,   # forwarded to recent gseapy versions
)
```

`enrich_peak_genes(ann, local=False, ...)` builds the peak-to-gene links with
`seq2gene_links` and then runs Enrichr.

## Local ORA (offline)

`ora()` implements the hypergeometric test and BH/Bonferroni correction without
SciPy or network access.

```python
res = pcs.ora(
    genes,
    gene_sets="gene_sets.gmt",     # or {term: [genes]}
    background=background_genes,   # defaults to the union of all sets
    p_adjust="BH",
    min_overlap=1,
)
```

Output columns: `Term`, `Overlap`, `GeneSetSize`, `BackgroundSize`, `QuerySize`,
`Expected`, `FoldEnrichment`, `OddsRatio`, `ORLower`, `ORUpper`, `PValue`,
`AdjustedPValue`, `Genes`.  The odds ratio uses a Woolfe log-normal 95%
confidence interval with a Haldane-Anscombe correction for zero cells.

Use `pcs.read_gmt(path)` to load a GMT file directly.

The local test is only as meaningful as the background supplied.  For
ChIP-seq work, pass the set of genes that could have been detected (for
example all genes with a valid TSS in the annotation), not the whole genome.

## Preranked GSEA

```python
ranks = {"G1": 3.2, "G2": 1.7, "G3": 0.4}
res = pcs.gsea(
    ranks,
    "GO_Biological_Process_2023",   # library, GMT path or {term: genes}
    permutation_num=1000,
    seed=20261026,
)
```

`gsea` accepts a `{gene: score}` mapping, a two-column DataFrame or a ranked
gene sequence (first gene highest).  `gsea_peak_genes(ann, ranks, ...)`
restricts the ranking to peak-linked genes first.  Both require the `enrich`
extra and return the `gseapy` result table.

## Plots

```python
fig, ax = pcs.plot_enrichment(res, top=20)              # -log10(padj) bars
fig, ax = pcs.plot_enrichment(res, top=20, kind="dot")  # clusterProfiler style
```

The plotter accepts ORA, Enrichr and GSEA tables and follows the same ggplot2
`theme_grey` defaults as the rest of pychipseeker.  From the command line,
`pychipseeker enrich ... --plot` writes `enrichment.png` beside the table.
