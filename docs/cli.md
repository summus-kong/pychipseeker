# Command line

```text
pychipseeker annotate --peaks PEAKS (--gtf GTF | --txdb TXDB) --output OUT
pychipseeker plot     --annotated TABLE --outdir DIR
pychipseeker enrich   --annotated TABLE [--gmt GMT --local]
pychipseeker report   --peaks PEAKS (--gtf GTF | --txdb TXDB) --output PDF
pychipseeker cache    list | download ASSEMBLY... | clear [ASSEMBLY...] [--all]
```

## annotate

```bash
pychipseeker annotate \
    --peaks peaks.narrowPeak \
    --txdb mm10 \
    --output peaks.annotated.tsv \
    --plot-dir plots \
    --summary
```

`--txdb` accepts a bundled assembly (`hg38`, `hg19`, `mm10`, `mm9`), an
on-demand assembly (`dm6`, `danRer11`, `rn7`, `mm39`, `sacCer3`), an
`ensembl:`/`gencode:`/`url:` spec, or an alias.  Use `--gtf` for a local
GTF/GFF file.

## report

```bash
pychipseeker report \
    --peaks peaks.narrowPeak \
    --txdb mm10 \
    --bigwig signal.bw \
    --sample 3000 \
    --output report.pdf
```

The PDF contains title/summary, pie, bar, TSS distance, coverage, tag-matrix
profile/heatmap, optional bigWig TSS profile/heatmap, `annoStat`, and methods
pages.  `--no-tag-matrix` skips the tag-matrix pages.

## cache

```bash
pychipseeker cache list
pychipseeker cache download dm6 sacCer3
pychipseeker cache clear dm6
pychipseeker cache clear --all
```

Installed wheel data is never removed; `clear` only touches the user cache.

## enrich

Online Enrichr:

```bash
pychipseeker enrich \
    --annotated peaks.annotated.tsv \
    --gene-column gene_name \
    --gene-sets GO_Biological_Process_2023 \
    --organism mouse \
    --background background_genes.txt \
    --outdir enrichment
```

Offline over-representation analysis with a GMT file:

```bash
pychipseeker enrich \
    --annotated peaks.annotated.tsv \
    --local \
    --gmt gene_sets.gmt \
    --background background_genes.txt \
    --p-adjust BH \
    --outdir enrichment
```

`--p-adjust` accepts `BH` (default), `Bonferroni` or `none`.
