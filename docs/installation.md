# Installation

## Requirements

- Python 3.9 or newer
- NumPy, pandas and matplotlib (installed automatically)

The core has **no compiler requirement** and no SciPy dependency. Genome
bundles are downloaded on first use and cached locally.

## Install from PyPI

```bash
pip install pychipseeker
```

Verify the install:

```bash
python -c "import pychipseeker as pcs; print(pcs.__version__)"
pychipseeker --help
```

## Optional extras

Install only what you need:

| Extra | Adds | Used for |
|---|---|---|
| `bigwig` | `pyBigWig` (non-Windows) | `.bw` / `.bigWig` profiles, peak calling, FRiP |
| `bigbed` | `pybigtools` | reading `.bb` / `.bigBed` peak files |
| `bam` | `pysam` (non-Windows) | BAM/CRAM tag matrices, FRiP and TSS QC |
| `enrich` | `gseapy` | Enrichr and preranked GSEA |
| `all` | all of the above | full feature set |
| `dev` | pytest, pytest-cov, ruff, mypy | development |
| `docs` | MkDocs Material | building this documentation |

```bash
pip install "pychipseeker[all]"
```

!!! note "Windows and optional binary packages"

    `pyBigWig` and `pysam` do not publish wheels for every Windows/Python
    combination, so the `bigwig`, `bam` and `all` extras skip them on Windows.
    This keeps `pip install pychipseeker[all]` working there.

    If you need direct signal or BAM access on Windows, use conda:

    ```bash
    conda install -c conda-forge pybigwig
    conda install -c bioconda pysam
    ```

    Annotation, tag matrices, QC from existing matrices, summaries, reports and
    genome tables all work without them.

## Development install

```bash
git clone https://github.com/summus-kong/pychipseeker
cd pychipseeker
python -m pip install -e ".[dev]"
pytest -q
```

Add optional features as needed, for example `-e ".[dev,bigwig,enrich]"`.

## Build the documentation locally

```bash
python -m pip install -e ".[docs]"
mkdocs serve
```

Open <http://127.0.0.1:8000/pychipseeker/>. Use `mkdocs build --strict` before
opening a pull request; the same command runs in the Pages workflow.

## Genome data

The wheel does not embed the large genome tables. `hg19`, `hg38`, `mm9` and
`mm10` are downloaded once from the GitHub release; other assemblies are built
from Ensembl GTFs and cached. See [Genomes and TxDb](genomes.md) for cache
locations, offline mirrors and the CLI cache commands.

## Next step

Run the [Quick start](quickstart.md) on your own peaks, or browse the
[Figure gallery](gallery.md) to see what the plotting layer can produce.
