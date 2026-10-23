# Third-party bundled content and licenses

pychipseeker's own code is **MIT** licensed (see `LICENSE`).  The following
data and derived files bundled in this repository come from third parties and
keep their original licenses; they are **not** covered by pychipseeker's MIT
license.

## From ChIPseeker (GPL-3)

Source: ChIPseeker for ChIP peak Annotation, Visualization and Comparison,
Copyright (C) 2015- Guangchuang YU, Guangzhou.  Distributed in R as an
R/Bioconductor package under the **GNU General Public License version 3**.
We include verbatim-evaluable transformations (CSV exports) of upstream data
files so that all bundled content is readable without R:

| Bundled file(s) | Upstream source |
|---|---|
| `src/pychipseeker/data/gsminfo.csv.gz` | ChIPseeker `data/gsminfo.rda` |
| `scripts/ucsc_release.csv.gz` | ChIPseeker `data/ucsc_release.rda` |
| `src/pychipseeker/data/GEO_sample_data/*.bed.gz` | ChIPseeker `inst/extdata/GEO_sample_data/` |

Upstream reference: https://github.com/GuangchuangYu/ChIPseeker
(ChIPseeker 1.48.0 / Bioconductor release 3.21).

The CSV exports were produced from the upstream `.rda` files without any
substantive modification of their content, using `pyreadr`.  They inherit the
upstream GPL-3 license: if you use or redistribute them, the GPL-3 terms
apply to those files, including a copy of the license.  See
https://www.gnu.org/licenses/gpl-3.0.html.

The GEO sample peak files are originally supplementary data from GEO
(NCBI Gene Expression Omnibus) samples GSM1174480-GSM1174482, GSM1295076 and
GSM1295077, as distributed in ChIPseeker's `extdata`.

## From pyBigWig (MIT)

| Bundled file(s) | Upstream source |
|---|---|
| `windows_build/io.c` | `libBigWig/io.c` as vendored in pyBigWig 0.3.26 |

pyBigWig and libBigWig are MIT-licensed, Copyright (c) 2015 Devon Ryan
(https://github.com/dpryan79/pyBigWig).  The bundled copy carries a patch
adding UTF-8 aware local file opening on Windows; the upstream MIT license
notice is preserved in the file header.

## Test data

| Bundled file(s) | Upstream source |
|---|---|
| `tests/data/test.bw` | pyBigWig test suite (`pyBigWigTest/test.bw`) |

MIT, same pyBigWig upstream as above.

## ENCODE demo signal

| Bundled file(s) | Upstream source |
|---|---|
| `src/pychipseeker/data/ENCODE_demo_signal.bw` | ENCODE experiment ENCSR000CDQ, file `ENCFF893ZGR` (H3K4me3 fold change over control, *Mus musculus* mm10), https://www.encodeproject.org/files/ENCFF893ZGR/ |
| `src/pychipseeker/data/ENCODE_demo_signal.json` | generated provenance record (accession, MD5/SHA-256 of source and bundle, subsetting window and binning parameters) |

The bundled bigWig is a small derived subset of `ENCFF893ZGR`: one 20 Mb
window of chr5 re-binned to 25 bp means, with near-background bins
(fold change < 0.05) left as gaps.  ENCODE portal data from released
experiments is distributed by ENCODE with no usage restrictions; the
subset was produced by `scripts/build_encode_subset.py` (which also
documents the exact provenance) and is regenerated deterministically.

## Recreating the bundled datasets

- `gsminfo.csv.gz` can be rebuilt from a locally downloaded `GEOmetadb.sqlite`
  with `python scripts/build_gsminfo.py` (offline, stdlib only).
- `ucsc_release.csv.gz` is additionally referenced by that script.
- `ENCODE_demo_signal.bw` is rebuilt from the downloaded ENCODE file with
  `python scripts/build_encode_subset.py --out-dir src/pychipseeker/data`;
  the writer used is `scripts/build_bigwig_py.py` (pure-Python bigWig
  writer, stdlib only, needed because the MinGW-built
  pyBigWig/libBigWig writer crashes on Windows).
