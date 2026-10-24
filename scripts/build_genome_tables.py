"""Build the bundled genome-annotation tables (``TxDb="assembly"`` bundles).

Downloads a GENCODE GTF, converts it through the regular
``GeneAnnotation.from_gtf`` pipeline and serialises the six derived frames
(genes / transcripts / exons / introns / 5'UTR / 3'UTR) as a compact
npz bundle plus a provenance JSON record, written next to the other
package data in ``src/pychipseeker/data/``.

Usage:

    python scripts/build_genome_tables.py --assembly hg38
    python scripts/build_genome_tables.py --assembly hg19 \
        --gtf /path/to/gencode.v19.annotation.gtf.gz

The bundle reproduces the same downstream semantics as ``from_gtf``, so
``TxDb="hg38"`` behaves exactly like ``GeneAnnotation.from_gtf(<bundled
GTF>)``.

GENCODE data are free for academic and commercial research use; see
``NOTICE.md`` and the provenance JSON written next to each bundle.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from pychipseeker.annotation import GeneAnnotation  # noqa: E402
from pychipseeker.genomes import TABLES, SUPPORTED_GENOMES, save_bundle  # noqa: E402

#: GENCODE release used per assembly (kept explicit for reproducibility).
DEFAULT_RELEASES: Dict[str, str] = {
    "hg38": "release_47",
    "hg19": "release_19",
    "mm10": "release_M25",
    "mm9": "release_M1",
}

#: filename per assembly (basic subsets where GENCODE provides them: closest
#: to the per-gene models ChIPseeker's TxDb packages carry and much smaller)
DEFAULT_FILES: Dict[str, str] = {
    "hg38": "gencode.v47.basic.annotation.gtf.gz",
    "hg19": "gencode.v19.annotation.gtf.gz",
    "mm10": "gencode.vM25.basic.annotation.gtf.gz",
    "mm9": "gencode.vM1.annotation.gtf.gz",
}

#: columns dropped from bundles: raw GTF decoration unused by the annotation
#: pipeline (annotation_id/rank/total/name columns are kept -- plotDistToTSS
#: style detail strings depend on them)
BUNDLE_DROP = (
    "attributes",
    "exon_id",
    "source",
    "score",
    "frame",
    "feature",
)

GENCODE_BASE = "https://ftp.ebi.ac.uk/pub/databases/gencode"

ALIASES: Dict[str, list] = {
    "hg38": ["hg38", "GRCh38"],
    "hg19": ["hg19", "GRCh37"],
    "mm10": ["mm10", "GRCm38"],
    "mm9": ["mm9", "GRCm37"],
}

LICENSE_NOTE = (
    "GENCODE data: free for academic and commercial research use; see "
    "https://www.gencodegenes.org/pages/license.html and NOTICE.md"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _gencode_species(assembly: str) -> str:
    return "mouse" if assembly in {"mm10", "mm9"} else "human"


def download_gtf(assembly: str, release: str, filename: str,
                 cache_dir: Path, url: Optional[str] = None) -> tuple[Path, str]:
    source_url = url or (
        f"{GENCODE_BASE}/Gencode_{_gencode_species(assembly)}/{release}/{filename}"
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / filename
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[cache] {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest, source_url
    print(f"[download] {source_url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Optional[Exception] = None
    for attempt in range(1, 8):
        try:
            with urllib.request.urlopen(source_url, timeout=120) as resp, open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh)
            break
        except Exception as exc:  # noqa: BLE001 - retry ANY network error
            last_error = exc
            print(f"  attempt {attempt} failed: {exc}")
            if attempt < 7:
                import time

                time.sleep(20 * attempt)
    else:
        raise RuntimeError(f"could not download {source_url}") from last_error
    tmp.replace(dest)
    return dest, source_url


def frame_counts(model: GeneAnnotation) -> Dict[str, int]:
    return {table: int(len(getattr(model, table))) for table in TABLES}


def build_bundle(
    assembly: str,
    gtf_path: Path,
    source_url: str,
    release: str,
    out_dir: Path,
    license_note: str = LICENSE_NOTE,
) -> Path:
    print(f"[convert] {gtf_path} -> GeneAnnotation ...")
    model = GeneAnnotation.from_gtf(gtf_path)
    counts = frame_counts(model)
    print(f"[frames] {counts}")

    npz_path = out_dir / f"genome_{assembly}.npz"
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    save_bundle(model, npz_path, drop=BUNDLE_DROP)
    provenance = {
        "assembly": assembly,
        "aliases": ALIASES[assembly],
        "gencode_release": release,
        "source_url": source_url,
        "source_sha256": sha256_file(gtf_path),
        "source_size_bytes": int(gtf_path.stat().st_size),
        "bundle": npz_path.name,
        "bundle_sha256": sha256_file(npz_path),
        "bundle_size_bytes": int(npz_path.stat().st_size),
        "built_at": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d"),
        "license": license_note,
        "counts": counts,
        "pipeline": "GENCODE GTF -> GeneAnnotation.from_gtf -> dict-encoded npz",
    }
    json_path = out_dir / f"genome_{assembly}.json"
    json_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"[built] {npz_path.name}: {npz_path.stat().st_size / 1e6:.2f} MB")
    print(f"[built] {json_path.name}")
    return npz_path


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True, choices=SUPPORTED_GENOMES)
    parser.add_argument("--gtf", type=Path, default=None,
                        help="local GENCODE GTF to convert (skips download)")
    parser.add_argument("--release", default=None,
                        help="GENCODE release name (default: known release)")
    parser.add_argument("--url", default=None,
                        help="explicit download URL (default: known file)")
    parser.add_argument("--out", type=Path,
                        default=REPO_ROOT / "src" / "pychipseeker" / "data")
    parser.add_argument("--cache-dir", type=Path,
                        default=REPO_ROOT / "build" / "gtf_cache")
    parser.add_argument("--license-note", default=None)
    args = parser.parse_args(argv)

    release = args.release or DEFAULT_RELEASES[args.assembly]
    if args.gtf is not None:
        gtf_path = args.gtf
        source_url = args.url or "(local file)"
    else:
        filename = (
            args.url.rsplit("/", 1)[-1]
            if args.url
            else DEFAULT_FILES[args.assembly]
        )
        gtf_path, source_url = download_gtf(
            args.assembly, release, filename, args.cache_dir, args.url
        )
    build_bundle(
        args.assembly,
        gtf_path,
        source_url=source_url,
        release=release,
        out_dir=args.out,
        license_note=args.license_note or LICENSE_NOTE,
    )


if __name__ == "__main__":
    main()
