# Genomes and TxDb

`TxDb=` accepts four kinds of annotation sources:

1. a bundled assembly name (`hg38`, `hg19`, `mm10`, `mm9`);
2. a downloadable assembly name (`dm6`, `danRer11`, `rn7`, `mm39`,
   `sacCer3`);
3. a GTF/GFF path or pandas DataFrame;
4. an `ensembl:`, `gencode:` or `url:` download spec.

```python
import pychipseeker as pcs

ann = pcs.annotatePeak("peaks.bed", TxDb="mm10")          # bundled
ann = pcs.annotatePeak("peaks.bed", TxDb="dm6")           # Ensembl GTF, cached
ann = pcs.annotatePeak("peaks.bed", TxDb="ensembl:mouse@110")
ann = pcs.annotatePeak("peaks.bed", TxDb="gencode:human@44")
ann = pcs.annotatePeak("peaks.bed", TxDb="url:https://example.org/genes.gtf.gz")
```

## Cache locations

- Windows: `%LOCALAPPDATA%\pychipseeker\genomes`
- Linux/macOS: `$XDG_CACHE_HOME/pychipseeker/genomes` or
  `~/.cache/pychipseeker/genomes`

Downloaded assemblies are parsed once and stored as compact
`genome_<assembly>.npz` tables next to their GTF and provenance JSON.  Use:

```bash
pychipseeker cache list
pychipseeker cache download dm6
pychipseeker cache clear --all
```

or the Python API:

```python
pcs.supported_genomes()
pcs.download_genome("dm6")
pcs.cached_genomes()
pcs.clear_genome_cache()
pcs.genome_provenance("mm10")
```

## Air-gapped mirrors

Set `PYCHIPSEEKER_GENOME_BASE` to a directory or HTTP mirror containing the
bundle files; it is tried before GitHub Releases and raw GitHub:

```bash
export PYCHIPSEEKER_GENOME_BASE=https://mirror.example/pychipseeker
```

For on-demand assemblies, download the Ensembl GTF once and point `TxDb=` at
the local path (or use an internal `url:`/`file://` mirror).
