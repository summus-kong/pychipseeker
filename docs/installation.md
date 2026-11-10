# Installation

Core installation needs only NumPy, pandas and matplotlib:

```bash
pip install pychipseeker
```

Optional features are split into extras:

| Extra | Adds | Used for |
|---|---|---|
| `bigwig` | `pyBigWig` | direct `.bw`/`.bigWig` TSS profiles, peak calling, FRiP |
| `bigbed` | `pybigtools` | reading `.bb`/`.bigBed` peak files |
| `bam` | `pysam` | BAM/CRAM tag matrices, FRiP and TSS QC |
| `enrich` | `gseapy` | Enrichr and preranked GSEA |
| `all` | all of the above | full feature set |
| `dev` | pytest, pytest-cov, ruff, mypy | development |

```bash
pip install "pychipseeker[all]"
```

For development:

```bash
git clone https://github.com/summus-kong/pychipseeker
cd pychipseeker
pip install -e ".[dev,bigwig]"
pytest -q
```

`pyBigWig` has no official Windows wheel on every Python version; use
conda-forge (`conda install -c conda-forge pybigwig`) or install a C compiler
if pip must build it from source.  Genome tables do not need `pyBigWig`; only
direct signal-track access does.
