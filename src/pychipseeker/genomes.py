"""Bundled genome annotation lookup — ChIPseeker's ``TxDb=`` analogue.

In R, ``ChIPseeker`` functions take a ``TxDb`` argument pointing at one of the
pre-built Bioconductor TxDb packages (``TxDb.Hsapiens.UCSC.hg19.knownGene``,
``TxDb.Mmusculus.UCSC.mm10.knownGene``, ...).  ``pychipseeker`` ships compact
bundled gene-annotation tables for the assemblies ChIPseeker users reach for
most often, so ``TxDb="hg38"`` (or an alias such as ``"GRCh38"``) works out of
the box, without any download or external dependency.

Supported assemblies and accepted aliases:

=========  =====================
bundle     accepted names
=========  =====================
``hg38``   hg38, GRCh38
``hg19``   hg19, GRCh37
``mm10``   mm10, GRCm38
``mm9``    mm9, GRCm37
=========  =====================

Matching is case-insensitive; see :func:`load_genome` and
:func:`is_assembly_name`.  A bundle stores the six frames produced by
``GeneAnnotation.from_gtf`` (genes / transcripts / exons / introns / 5'UTR /
3'UTR) dict-encoded inside ``data/genome_<assembly>.npz``; provenance for
each bundle lives next to it in ``data/genome_<assembly>.json`` and can be
read with :func:`genome_provenance`.
"""

from __future__ import annotations

import json
import hashlib
import os
import shutil
import sys
import time
import urllib.request
from typing import Any, Dict, Optional, Union
from pathlib import Path

import numpy as np
import pandas as pd

from .annotation import GeneAnnotation

#: assembly name -> bundled-table stem (without extension)
GENOME_ALIASES: Dict[str, str] = {
    "hg38": "hg38",
    "grch38": "hg38",
    "hg19": "hg19",
    "grch37": "hg19",
    "mm10": "mm10",
    "grcm38": "mm10",
    "mm9": "mm9",
    "grcm37": "mm9",
}

#: bundled assemblies, in the canonical order
SUPPORTED_GENOMES: tuple[str, ...] = ("hg38", "hg19", "mm10", "mm9")

#: the six derived frames carried by every bundle
TABLES: tuple[str, ...] = (
    "genes",
    "transcripts",
    "exons",
    "introns",
    "utr5",
    "utr3",
)

_loaded: Dict[str, GeneAnnotation] = {}


def is_assembly_name(value: Any) -> bool:
    """True when ``value`` is one of the supported bundled assembly names."""

    return isinstance(value, str) and value.strip().lower() in GENOME_ALIASES


def canonical_assembly(value: Any) -> str:
    """Map an assembly name or alias to the bundled assembly stem."""

    if not isinstance(value, str):
        raise ValueError(
            "TxDb must be a GeneAnnotation, a GTF path/DataFrame, or a bundled "
            f"assembly name; supported names: {', '.join(SUPPORTED_GENOMES)}"
        )
    name = value.strip().lower()
    if name not in GENOME_ALIASES:
        raise ValueError(
            f"unsupported TxDb assembly {value!r}; supported names: "
            f"{', '.join(SUPPORTED_GENOMES)}"
        )
    return GENOME_ALIASES[name]


# ---------------------------------------------------------------------------
# bundle serialisation (dict-coded frames shared with the build script)
# ---------------------------------------------------------------------------


def _encode_column(table: str, values: pd.Series) -> Dict[str, np.ndarray]:
    """Dict-encode one frame column into npz payload arrays.

    The pandas dtype tag travels in ``<table>__cols`` so the loader can
    rebuild the frames with the exact dtypes produced by ``from_gtf`` --
    including ``pd.NA`` missing data in extensible string columns.
    """

    name = str(values.name)
    prefix = f"{table}__{name}"
    tag = str(values.dtype)
    if tag == "int64":
        raw = values.to_numpy(dtype=object)
        missing = pd.isna(raw)
        if missing.any():
            arr = np.array(
                [0 if m else int(v) for v, m in zip(raw, missing)], dtype=np.int64
            )
        else:
            arr = raw.astype(np.int64)
        if int(np.max(arr, initial=0)) < np.iinfo(np.int32).max:
            arr = arr.astype(np.int32)
        if missing.any():
            return {
                prefix: arr,
                f"{prefix}__na": np.ascontiguousarray(missing, dtype=bool),
            }
        return {prefix: arr}
    flat = values.to_numpy(dtype=object)
    missing = pd.isna(flat)
    strings = np.array(
        ["" if m else str(v) for v, m in zip(flat, missing)], dtype=object
    )
    uniq, inverse = np.unique(strings, return_inverse=True)
    payload: Dict[str, np.ndarray] = {
        f"{prefix}__dict": uniq.astype(str),
        f"{prefix}__idx": np.ascontiguousarray(inverse, dtype=np.int32),
    }
    if missing.any():
        payload[f"{prefix}__na"] = np.ascontiguousarray(missing, dtype=bool)
    return payload


def _write_npz(dest, payload: Dict[str, np.ndarray]) -> None:
    """npz writer with DEFLATE level 9 (np.savez_compressed uses level 6)."""

    import zipfile

    with zipfile.ZipFile(
        dest, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9, allowZip64=True
    ) as zf:
        for key, arr in payload.items():
            with zf.open(key, "w", force_zip64=True) as fp:
                from numpy.lib import format as _npyformat

                _npyformat.write_array(
                    fp, np.ascontiguousarray(arr), allow_pickle=False
                )


def _encode_frame(frame: pd.DataFrame, table: str) -> Dict[str, np.ndarray]:
    specs = [[str(col), str(frame[col].dtype)] for col in frame.columns]
    out: Dict[str, np.ndarray] = {f"{table}__cols": json.dumps(specs).encode()}
    for col in frame.columns:
        out.update(_encode_column(table, frame[col]))
    return out


def save_bundle(model: GeneAnnotation, dest, drop: tuple[str, ...] = (
    "attributes",
)) -> None:
    """Serialise a :class:`GeneAnnotation` as a bundled npz table.

    Shared by the build script and the tests: the bundle reproduces the exact
    frames produced by ``GeneAnnotation.from_gtf`` (values, missing data and
    dtype tags); textual values are stored dict-encoded.  Columns named in
    ``drop`` (raw GTF ``attributes`` text by default -- bulky and unused by
    the annotation pipeline) are left out of the bundle.
    """

    payload: Dict[str, np.ndarray] = {"__tables": json.dumps(list(TABLES)).encode()}
    for table in TABLES:
        frame = getattr(model, table)
        if drop:
            dropcols = [c for c in frame.columns if str(c) in drop]
            if dropcols:
                frame = frame.drop(columns=dropcols)
        payload.update(_encode_frame(frame, table))
    _write_npz(dest, payload)


#: Public source of on-demand bundle downloads (raw files of this repository's
#: committed ``src/pychipseeker/data`` directory; override with the
#: ``PYCHIPSEEKER_GENOME_BASE`` environment variable, e.g. for a mirror).
DEFAULT_BUNDLE_BASE = (
    "https://raw.githubusercontent.com/summus-kong/pychipseeker/main/"
    "src/pychipseeker/data"
)


def cache_dir() -> Path:
    """Directory where downloaded genome bundles are stored.

    ``%LOCALAPPDATA%\\pychipseeker\\genomes`` on Windows and
    ``~/.cache/pychipseeker/genomes`` (honouring ``XDG_CACHE_HOME``)
    elsewhere.  Created on first call.
    """

    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home())
        root = Path(base) / "pychipseeker" / "genomes"
    else:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
        root = Path(base) / "pychipseeker" / "genomes"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _installed_bundle(name: str) -> Optional[Path]:
    """Bundle npz shipped inside the installed package, if any."""

    from importlib.resources import as_file, files

    try:
        with as_file(
            files("pychipseeker").joinpath("data", f"genome_{name}.npz")
        ) as path:
            return Path(path) if Path(path).is_file() else None
    except (FileNotFoundError, ModuleNotFoundError):
        return None


def _cache_bundle(name: str) -> Optional[Path]:
    """Bundle npz previously fetched into the user cache dir, if any."""

    path = cache_dir() / f"genome_{name}.npz"
    if path.is_file() and path.stat().st_size > 0:
        return path
    return None


def _download_file(url: str, dest: Path) -> None:
    """Stream ``url`` into ``dest`` with retries for flaky connections."""

    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Optional[Exception] = None
    for attempt in range(1, 8):
        try:
            with urllib.request.urlopen(url, timeout=180) as resp, open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh)
            break
        except Exception as exc:  # noqa: BLE001 - retry any network error
            last_error = exc
            print(
                f"pychipseeker: genome bundle download attempt {attempt} "
                f"failed ({exc}); retrying" if attempt < 7 else f" ({exc})",
                file=sys.stderr,
            )
            if attempt < 7:
                time.sleep(8 * attempt)
    else:
        raise RuntimeError(
            f"could not download genome bundle from {url}"
        ) from last_error
    tmp.replace(dest)


def bundle_url(name: str) -> str:
    """Download URL for a bundled genome table (npz)."""

    base = os.environ.get("PYCHIPSEEKER_GENOME_BASE") or DEFAULT_BUNDLE_BASE
    return f"{base}/genome_{name}.npz"


def _download_bundle(name: str) -> Path:
    """Fetch a genome bundle (npz + provenance JSON) into the cache dir.

    The provenance file is fetched first and the downloaded npz is verified
    against the recorded sha256; provenance travel together with the bundle
    so :func:`genome_provenance` also works offline after the first fetch.
    """

    base = os.environ.get("PYCHIPSEEKER_GENOME_BASE") or DEFAULT_BUNDLE_BASE
    try:
        _download_file(f"{base}/genome_{name}.json", cache_dir() / f"genome_{name}.json")
        with open(cache_dir() / f"genome_{name}.json", "r", encoding="utf-8") as fh:
            expected = json.load(fh).get("bundle_sha256")
    except Exception:  # noqa: BLE001 - provenance is a nicety, not a hard requirement
        expected = None

    dest = cache_dir() / f"genome_{name}.npz"
    _download_file(bundle_url(name), dest)
    if expected:
        digest = hashlib.sha256()
        with open(dest, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                digest.update(chunk)
        if digest.hexdigest() != expected:
            dest.unlink()
            raise RuntimeError(
                f"downloaded genome bundle {dest.name} failed its sha256 check"
            )
    return dest


def _bundle_path(name: str) -> Path:
    """Locate a bundle npz: installed data -> cache -> public download."""

    path = _installed_bundle(name) or _cache_bundle(name)
    if path is None:
        path = _download_bundle(name)
    return path


def _rebuild_frame(npz, table: str) -> pd.DataFrame:
    """Rebuild one frame from the stored column specs ``(name, dtype tag)``."""

    specs = json.loads(bytes(npz[f"{table}__cols"]).decode())
    columns: Dict[str, Any] = {}
    for name, tag in specs:
        key = f"{table}__{name}"
        if tag == "int64":
            integers = npz[key].astype(np.int64)
            if f"{key}__na" in npz.files:
                values = integers.astype(object)
                values[npz[f"{key}__na"].astype(bool)] = pd.NA
                columns[name] = values
            else:
                columns[name] = integers
            continue
        strings = npz[f"{key}__dict"][npz[f"{key}__idx"]]
        has_na = f"{key}__na" in npz.files
        mask = npz[f"{key}__na"].astype(bool) if has_na else None
        if tag.startswith("<M8"):
            columns[name] = strings
            continue
        if tag == "object":
            values = np.asarray(strings, dtype=object)
            if mask is not None:
                values[mask] = pd.NA
            columns[name] = values
            continue
        # extensible string dtypes ("string", "str", ...)
        arr = pd.array([str(v) for v in strings], dtype=tag)
        if mask is not None:
            arr[mask] = pd.NA
        columns[name] = arr
    return pd.DataFrame(columns)


def load_genome(
    assembly: Union[str, None],
    refresh: bool = False,
) -> GeneAnnotation:
    """Load a bundled assembly table and return a :class:`GeneAnnotation`.

    The loaded model is cached in-process; pass ``refresh=True`` to rebuild
    it.  Note that building the interval indexes for a full genome takes a
    few seconds on first use; later calls return instantly.

    Bundles currently ship for: ``hg38``/``GRCh38``, ``hg19``/``GRCh37``,
    ``mm10``/``GRCm38`` and ``mm9``/``GRCm37``.
    """

    name = canonical_assembly(assembly)
    if name in _loaded and not refresh:
        return _loaded[name]

    npz = np.load(str(_bundle_path(name)), allow_pickle=False)

    frames = [_rebuild_frame(npz, table) for table in TABLES]
    npz.close()

    model = GeneAnnotation(*frames)
    _loaded[name] = model
    return model


def genome_provenance(assembly: Union[str, None]) -> dict:
    """Return the provenance record of a bundled genome table.

    GENCODE release, source URL, checksums and counts; the JSON is looked up
    beside the bundle (installed data or download cache) and fetched on
    demand otherwise.
    """

    name = canonical_assembly(assembly)
    installed = _installed_bundle(name)
    candidates = []
    if installed is not None:
        candidates.append(installed.with_suffix(".json"))
    candidates.append(cache_dir() / f"genome_{name}.json")
    for candidate in candidates:
        if candidate.is_file():
            with open(candidate, "r", encoding="utf-8") as fh:
                return json.load(fh)
    _download_bundle(name)
    with open(cache_dir() / f"genome_{name}.json", "r", encoding="utf-8") as fh:
        return json.load(fh)


def supported_genomes() -> tuple[str, ...]:
    """Bundled assemblies that can be passed as ``TxDb``."""

    return SUPPORTED_GENOMES
