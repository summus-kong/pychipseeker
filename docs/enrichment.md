# Enrichment

Two backends are available.

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
`Expected`, `FoldEnrichment`, `PValue`, `AdjustedPValue`, `Genes`.

Use `pcs.read_gmt(path)` to load a GMT file directly.

The local test is only as meaningful as the background supplied.  For
ChIP-seq work, pass the set of genes that could have been detected (for
example all genes with a valid TSS in the annotation), not the whole genome.
