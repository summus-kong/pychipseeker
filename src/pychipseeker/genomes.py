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
from typing import Any, Dict, Optional, Union

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
        return {prefix: values.to_numpy(dtype=np.int64)}
    flat = values.to_numpy(dtype=object)
    missing = pd.isna(flat)
    strings = np.array(
        ["" if m else str(v) for v, m in zip(flat, missing)], dtype=object
    )
    uniq, inverse = np.unique(strings, return_inverse=True)
    payload: Dict[str, np.ndarray] = {
        f"{prefix}__dict": uniq.astype(str),
        f"{prefix}__idx": np.ascontiguousarray(inverse, dtype=np.int64),
    }
    if missing.any():
        payload[f"{prefix}__na"] = np.ascontiguousarray(missing, dtype=bool)
    return payload


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
    np.savez_compressed(dest, **payload)


def _bundle_path(assembly: str):
    from importlib.resources import as_file, files

    return as_file(
        files("pychipseeker").joinpath("data", f"genome_{assembly}.npz")
    )


def _rebuild_frame(npz, table: str) -> pd.DataFrame:
    """Rebuild one frame from the stored column specs ``(name, dtype tag)``."""

    specs = json.loads(bytes(npz[f"{table}__cols"]).decode())
    columns: Dict[str, Any] = {}
    for name, tag in specs:
        key = f"{table}__{name}"
        if tag == "int64" and f"{key}__na" not in npz.files:
            columns[name] = npz[key]
            continue
        strings = npz[f"{key}__dict"][npz[f"{key}__idx"]]
        has_na = f"{key}__na" in npz.files
        mask = npz[f"{key}__na"].astype(bool) if has_na else None
        if tag == "int64":
            # NA-ful integer column (e.g. nullable int dtype)
            values = np.empty(len(strings), dtype=object)
            values[:] = strings.astype(np.int64)
            if mask is not None:
                values[mask] = pd.NA
            columns[name] = values
            continue
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

    with _bundle_path(name) as path:
        npz = np.load(str(path), allow_pickle=False)

    frames = [_rebuild_frame(npz, table) for table in TABLES]
    npz.close()

    model = GeneAnnotation(*frames)
    _loaded[name] = model
    return model


def genome_provenance(assembly: Union[str, None]) -> dict:
    """Return the provenance record of a bundled genome table."""

    name = canonical_assembly(assembly)
    from importlib.resources import as_file, files

    with as_file(
        files("pychipseeker").joinpath("data", f"genome_{name}.json")
    ) as path:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)


def supported_genomes() -> tuple[str, ...]:
    """Bundled assemblies that can be passed as ``TxDb``."""

    return SUPPORTED_GENOMES
