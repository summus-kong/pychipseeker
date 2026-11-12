# FAQ and troubleshooting

## Is this a port of R ChIPseeker?

It is an independent Python implementation that follows ChIPseeker's
annotation rules: promoter bins, feature priority, strand-aware signed TSS
distances, `annoStat` categories and the common plotting helpers. The
[PARITY.md](https://github.com/summus-kong/pychipseeker/blob/main/PARITY.md)
file lists what is covered and where behaviour differs.

You do **not** need R or Bioconductor. `TxDb=` accepts a GTF/GFF file, a
`GeneAnnotation`, a prebuilt assembly name, an Ensembl/GENCODE spec or a URL.

## Which input formats work?

BED, narrowPeak, broadPeak, gappedPeak (plain or gzip), TSV/CSV, bigBed
(optional `pybigtools`) and BAM/CRAM (optional `pysam`). For BED-like files,
the fourth column becomes `peak_id`; extra columns are preserved as peak
metadata.

## The genome download failed, or I am offline

The wheel does not embed genome tables. Prebuilt bundles are downloaded from
the GitHub release and cached locally. Options:

1. Retry with `pychipseeker cache download mm10`.
2. Point `PYCHIPSEEKER_GENOME_BASE` at an internal HTTP or directory mirror:

   ```bash
   export PYCHIPSEEKER_GENOME_BASE=https://mirror.example/pychipseeker
   ```

3. Use a local GTF/GFF and skip the bundle entirely:

   ```python
   ann = pcs.annotate_peak("peaks.bed", "genes.gtf")
   ```

4. Download the Ensembl GTF once, then reuse the local path or a `file://` /
   `url:` spec.

Cache locations and cache commands are documented in
[Genomes and TxDb](genomes.md).

## `TxDb="auto"` refuses to guess

Peak-only detection is deliberately conservative. Chromosome names such as
`chr1` are shared by human and mouse, so the package raises an error instead of
silently choosing an assembly. Pass an explicit assembly, chromosome sizes, or
a bigWig header.

## No peaks were annotated

Usually this is a chromosome naming mismatch (`chr1` vs `1`). pychipseeker
harmonises common UCSC/Ensembl styles, but a custom header may still fail.
Check both tables:

```python
print(sorted(ann.peaks["chrom"].unique())[:5])
print(sorted(ann.gene_annotation.features("gene")["chrom"].unique())[:5])
```

The annotation call also emits a warning when every peak on a chromosome is
dropped.

## `pyBigWig` or `pysam` will not install on Windows

Those projects do not publish wheels for every Windows/Python combination.
Install them with conda instead:

```bash
conda install -c conda-forge pybigwig
conda install -c bioconda pysam
```

Annotation, peak-coverage tag matrices, summaries, reports and genome tables
work without them. Tests that need the missing optional package skip rather
than fail.

## Multiprocessing hangs or restarts on Windows

Windows and macOS use the `spawn` start method. Guard script entry points:

```python
if __name__ == "__main__":
    ann = pcs.annotate_peak("peaks.bed", TxDb="mm10", n_jobs=-1)
```

Notebooks and interactive shells are fine for `n_jobs=1`. For small inputs,
parallelism costs more than it saves; the package already applies internal
thresholds. See [Performance](performance.md).

## PDF text is not editable, or fonts look wrong

PDF, SVG and EPS outputs are vector files. PDF uses Type 42 fonts so text stays
selectable and embeddable. If a journal requires outlined text or a specific
font, set matplotlib's font family before plotting:

```python
import matplotlib as mpl

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
```

Then save with `fig.savefig("figure.pdf", bbox_inches="tight")`.

## How do I export only PDF or only vector formats?

For your own figures, save explicitly:

```python
fig.savefig("figure.pdf", bbox_inches="tight")
fig.savefig("figure.svg", bbox_inches="tight")
```

For the ENCODE gallery:

```bash
python examples/make_encode_plate.py --formats pdf,svg
```

All plates and panels are always available as PDF; `--formats` only limits what
is written in a given run.

## How do I choose which panels to draw?

```bash
python examples/make_encode_plate.py --list-panels
python examples/make_encode_plate.py --panels A,C,K,P
```

`--panels` also accepts stems such as `width_ecdf` and the keyword `all`. The
result is `encode_custom_plate.*`. See [Figure gallery](gallery.md).

## Where are the example outputs?

`examples/output/` is generated and intentionally gitignored. Run
`python examples/make_encode_plate.py` to recreate the plates, the 18 panels,
`encode_figures.pdf` and `gallery.html`.

## Why do the signal quartiles put genes without a peak in Q1?

A gene without an assigned peak has no signal evidence. Treating it as zero
keeps it in the lowest group instead of filling the median, which would create
a large artificial block in the middle of the ranking. This makes the
centre-to-flank enrichment in panel P monotonic and easy to interpret.

## Can I use this in a publication?

Yes. The code is MIT licensed. For the example data, cite the ENCODE Consortium
and accession **ENCSR000CDQ** (mouse H3K4me3, mm10); the signal track used for
peak calling is **ENCFF893ZGR**. The package itself has no required citation —
link to the repository instead.

## Still stuck?

Open an issue with the command, the package version
(`python -c "import pychipseeker; print(pychipseeker.__version__)"`), your OS
and Python version, and a minimal reproducible example:
<https://github.com/summus-kong/pychipseeker/issues>.
