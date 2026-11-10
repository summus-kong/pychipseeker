# Genomes and TxDb

`TxDb=` accepts five kinds of annotation sources:

1. a prebuilt assembly name (`hg38`, `hg19`, `mm10`, `mm9`);
2. a downloadable assembly name (`dm6`, `danRer11`, `rn7`, `mm39`,
   `sacCer3`);
3. a GTF/GFF path or pandas DataFrame;
4. an `ensembl:`, `gencode:` or `url:` download spec;
5. `"auto"`, which runs `detect_assembly` against the peaks.  Peak-only
detection is deliberately conservative: names shared by human and mouse (for
example `chr1`, `chr2`) raise an error instead of guessing.  Use chromosome
sizes, a bigWig header or an explicit assembly when possible.

```python
import pychipseeker as pcs

ann = pcs.annotatePeak("peaks.bed", TxDb="mm10")          # prebuilt bundle
ann = pcs.annotatePeak("peaks.bed", TxDb="dm6")           # Ensembl GTF, cached
ann = pcs.annotatePeak("peaks.bed", TxDb="ensembl:mouse@110")
ann = pcs.annotatePeak("peaks.bed", TxDb="gencode:human@44")
ann = pcs.annotatePeak("peaks.bed", TxDb="url:https://example.org/genes.gtf.gz")
ann = pcs.annotatePeak("peaks.bed", TxDb="auto")
```

## Cache locations

- Windows: `%LOCALAPPDATA%\pychipseeker\genomes`
- Linux/macOS: `$XDG_CACHE_HOME/pychipseeker/genomes` or
  `~/.cache/pychipseeker/genomes`

The wheel does not embed the large `.npz` tables.  `hg19`, `hg38`, `mm9` and
`mm10` are downloaded once from the GitHub Release; other assemblies are built
from Ensembl GTFs.  Each bundle is stored as compact `genome_<assembly>.npz`
tables next to their source and provenance JSON.  Use:

```bash
pychipseeker cache list                 # human-readable table
pychipseeker cache list --json          # machine-readable
pychipseeker cache info mm10            # provenance record
pychipseeker cache download dm6
pychipseeker cache prune                # drop partial/orphaned files
pychipseeker cache clear --all
```

or the Python API:

```python
pcs.supported_genomes()
pcs.detect_assembly(bigwig="signal.bw")
pcs.download_genome("dm6")
pcs.cached_genomes()
pcs.genome_provenance("mm10")
pcs.prune_genome_cache()
pcs.clear_genome_cache()
```

## Air-gapped mirrors

Set `PYCHIPSEEKER_GENOME_BASE` to a directory or HTTP mirror containing the
bundle files; it is tried before GitHub Releases and raw GitHub:

```bash
export PYCHIPSEEKER_GENOME_BASE=https://mirror.example/pychipseeker
```

For on-demand assemblies, download the Ensembl GTF once and point `TxDb=` at
the local path (or use an internal `url:`/`file://` mirror).
