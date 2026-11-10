# Installation

Core installation needs only NumPy, pandas and matplotlib:

```bash
pip install pychipseeker
```

Optional features are split into extras:

| Extra | Adds | Used for |
|---|---|---|
| `bigwig` | `pyBigWig` (non-Windows) | direct `.bw`/`.bigWig` TSS profiles, peak calling, FRiP |
| `bigbed` | `pybigtools` | reading `.bb`/`.bigBed` peak files |
| `bam` | `pysam` (non-Windows) | BAM/CRAM tag matrices, FRiP and TSS QC |
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

`pyBigWig` and `pysam` have no PyPI wheels on every Windows/Python version,
so the `bigwig`, `bam` and `all` extras skip them on Windows.  This keeps
`pip install pychipseeker[all]` working there.  Install them separately with
conda-forge (`conda install -c conda-forge pybigwig`) or bioconda
(`conda install -c bioconda pysam`) if you need direct signal or BAM access;
otherwise all annotation, tag-matrix, QC and report features work without
them.  Genome tables do not need `pyBigWig`.
