# pychipseeker documentation

`pychipseeker` is a Python implementation of ChIPseeker-compatible peak
annotation, summarization and visualization.

- [Installation](installation.md)
- [Quick start](quickstart.md)
- [Command line](cli.md)
- [Genomes and TxDb](genomes.md)
- [Enrichment](enrichment.md)
- [Performance](performance.md)

The API mirrors the R package where practical: `annotatePeak`, `readPeakFile`,
`getPromoters`, `getTagMatrix`, `plotAnnoPie`, `plotAnnoBar`, `plotDistToTSS`,
`covplot`, `vennplot` and friends are all available alongside snake_case
aliases.  See the [README](https://github.com/summus-kong/pychipseeker#readme)
and [PARITY.md](https://github.com/summus-kong/pychipseeker/blob/main/PARITY.md)
for the full compatibility matrix.
