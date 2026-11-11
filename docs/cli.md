# Command line

```text
pychipseeker annotate --peaks PEAKS (--gtf GTF | --txdb TXDB) --output OUT
pychipseeker plot     --annotated TABLE --outdir DIR
pychipseeker enrich   --annotated TABLE [--gmt GMT --local] [--plot]
pychipseeker callpeak --bigwig SIGNAL --output BED
pychipseeker filter   --peaks PEAKS --output BED
pychipseeker qc       --peaks PEAKS (--bigwig BW | --bam BAM)
pychipseeker report   --peaks PEAKS (--gtf GTF | --txdb TXDB) --output PDF
pychipseeker cache    list [--json] | info NAME | prune |
                      download ASSEMBLY... | clear [ASSEMBLY...] [--all]
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

`--txdb` accepts a prebuilt assembly (`hg38`, `hg19`, `mm10`, `mm9`), an
on-demand assembly (`dm6`, `danRer11`, `rn7`, `mm39`, `sacCer3`), an
`ensembl:`/`gencode:`/`url:` spec, `auto`, or an alias.  Use `--gtf` for a
local GTF/GFF file.

`--jobs N` annotates large peak sets in parallel (`-1` = all cores); the
default `1` keeps the run sequential.

## callpeak and filter

```bash
pychipseeker callpeak \
    --bigwig signal.bw \
    --output peaks.bed \
    --threshold 3 \
    --max-gap 1000 \
    --min-width 500

pychipseeker filter \
    --peaks peaks.bed \
    --blacklist ENCFF356LFX.bed.gz \
    --drop-mito \
    --min-width 200 \
    --max-width 10000 \
    --output peaks.clean.bed
```

`callpeak` is a dependency-free threshold caller for signal tracks (it needs
the `bigwig` extra).  `filter` also accepts `--min-score`, `--chroms` and
`--exclude-chroms`.  `callpeak --jobs -1` calls chromosomes in parallel; each
worker opens the bigWig read-only.

## qc

```bash
pychipseeker qc \
    --peaks peaks.clean.bed \
    --bam alignments.bam \
    --txdb mm10 \
    --blacklist ENCFF356LFX.bed.gz \
    --output qc.tsv
```

The table contains peak counts/widths, mitochondrial and blacklist fractions,
FRiP and TSS enrichment.  Supply `--bigwig` instead of `--bam` for signal
coverage, `--no-tss` to skip TSS enrichment, or `--include-duplicates` to count
duplicate BAM reads in FRiP.  `--jobs N` parallelises the BAM coverage matrix
over TSS features (`-1` = all cores); FRiP itself uses the C-level counter and
is already fast on multi-million-read files.

## report

```bash
pychipseeker report \
    --peaks peaks.narrowPeak \
    --txdb mm10 \
    --bigwig signal.bw \
    --blacklist ENCFF356LFX.bed.gz \
    --qc \
    --upstream 2000 \
    --downstream 4000 \
    --sample 3000 \
    --output report.pdf
```

The PDF contains title/summary, pie, bar, TSS distance, coverage, tag-matrix
profile/heatmap, optional bigWig TSS profile/heatmap, optional QC, `annoStat`,
and methods pages.  `--no-tag-matrix` skips the tag-matrix pages.
`--upstream`/`--downstream` define an asymmetric signal window and default to
the `--tss-region` values.  `--jobs N` forwards the worker count to annotation,
QC and the bigWig/tag-matrix extraction (`-1` = all cores).

## cache

```bash
pychipseeker cache list
pychipseeker cache list --json
pychipseeker cache info mm10
pychipseeker cache download dm6 sacCer3
pychipseeker cache prune
pychipseeker cache clear dm6
pychipseeker cache clear --all
```

`list --json` is machine-readable, `info` prints the provenance record, and
`prune` removes interrupted downloads (`.part`, zero-byte bundles) and
orphaned JSON files.  Installed package data is never removed; `clear` only
touches the user cache.

## enrich

Online Enrichr:

```bash
pychipseeker enrich \
    --annotated peaks.annotated.tsv \
    --gene-column gene_name \
    --gene-sets GO_Biological_Process_2023 \
    --organism mouse \
    --background background_genes.txt \
    --outdir enrichment \
    --plot
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

`--p-adjust` accepts `BH` (default), `Bonferroni` or `none`.  `--plot` writes
`enrichment.png` beside the table; `--top` controls how many terms are shown.
