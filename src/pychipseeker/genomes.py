"""Bundled genome annotation lookup — ChIPseeker's ``TxDb=`` analogue.

In R, ``ChIPseeker`` functions take a ``TxDb`` argument pointing at one of the
pre-built Bioconductor TxDb packages (``TxDb.Hsapiens.UCSC.hg19.knownGene``,
``TxDb.Mmusculus.UCSC.mm10.knownGene``, ...).  ``pychipseeker`` ships compact
bundled gene-annotation tables for the assemblies ChIPseeker users reach for
most often, so ``TxDb="hg38"`` (or an alias such as ``"GRCh38"``) works out of
the box, without any download or external dependency.

Supported assemblies and accepted aliases:

=============  ==================================
bundle          accepted names
=============  ==================================
``hg38``        hg38, GRCh38
``hg19``        hg19, GRCh37
``mm10``        mm10, GRCm38
``mm9``         mm9, GRCm37
``dm6``         dm6, BDGP6.46
``danRer11``    danRer11, GRCz11
``rn7``         rn7, Rnor_6.0
``mm39``        mm39, GRCm39
``sacCer3``     sacCer3, R64-1-1
=============  ==================================

The first four ship inside the package; the rest are fetched as Ensembl GTFs
on first use and cached as compact npz bundles.  Arbitrary annotations can be
loaded with ``ensembl:<species>[@release]``, ``gencode:<human|mouse>[@release]``
or ``url:<gtf>`` specs through :func:`load_annotation`.

Matching is case-insensitive; see :func:`load_genome` and
:func:`is_assembly_name`.  A bundle stores the six frames produced by
``GeneAnnotation.from_gtf`` (genes / transcripts / exons / introns / 5'UTR /
3'UTR) dict-encoded inside ``data/genome_<assembly>.npz``; provenance for
each bundle lives next to it in ``data/genome_<assembly>.json`` and can be
read with :func:`genome_provenance`.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd

from .annotation import GeneAnnotation
from .seqlevels import canonical_chrom, chromosome_style

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
    # downloadable on first use (Ensembl/UCSC assembly names)
    "dm6": "dm6",
    "bdgp6": "dm6",
    "bdgp6.46": "dm6",
    "drosophila_melanogaster": "dm6",
    "danrer11": "danRer11",
    "grcz11": "danRer11",
    "danio_rerio": "danRer11",
    "rn7": "rn7",
    "rnor6": "rn7",
    "rnor_6.0": "rn7",
    "rattus_norvegicus": "rn7",
    "mm39": "mm39",
    "grcm39": "mm39",
    "mus_musculus": "mm39",
    "saccer3": "sacCer3",
    "r64": "sacCer3",
    "r64-1-1": "sacCer3",
    "saccharomyces_cerevisiae": "sacCer3",
}

#: assemblies whose compact tables ship inside the wheel
BUNDLED_GENOMES: tuple[str, ...] = ("hg38", "hg19", "mm10", "mm9")

#: assemblies downloaded as a GTF and converted on first use
GENOME_SOURCES: Dict[str, Dict[str, str]] = {
    "dm6": {
        "provider": "ensembl",
        "release": "110",
        "url": (
            "https://ftp.ensembl.org/pub/release-110/gtf/"
            "drosophila_melanogaster/"
            "Drosophila_melanogaster.BDGP6.46.110.gtf.gz"
        ),
        "description": "Drosophila melanogaster BDGP6.46 (UCSC dm6), Ensembl 110",
    },
    "danRer11": {
        "provider": "ensembl",
        "release": "110",
        "url": (
            "https://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/"
            "Danio_rerio.GRCz11.110.gtf.gz"
        ),
        "description": "Danio rerio GRCz11 (UCSC danRer11), Ensembl 110",
    },
    "rn7": {
        "provider": "ensembl",
        "release": "110",
        "url": (
            "https://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/"
            "Rattus_norvegicus.Rnor_6.0.110.gtf.gz"
        ),
        "description": "Rattus norvegicus Rnor_6.0 (UCSC rn7), Ensembl 110",
    },
    "mm39": {
        "provider": "ensembl",
        "release": "110",
        "url": (
            "https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/"
            "Mus_musculus.GRCm39.110.gtf.gz"
        ),
        "description": "Mus musculus GRCm39 (UCSC mm39), Ensembl 110",
    },
    "sacCer3": {
        "provider": "ensembl",
        "release": "110",
        "url": (
            "https://ftp.ensembl.org/pub/release-110/gtf/"
            "saccharomyces_cerevisiae/"
            "Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz"
        ),
        "description": "Saccharomyces cerevisiae R64-1-1 (UCSC sacCer3), Ensembl 110",
    },
}

#: bundled assemblies first, then assemblies built from a downloaded GTF
SUPPORTED_GENOMES: tuple[str, ...] = BUNDLED_GENOMES + tuple(GENOME_SOURCES)

#: chromosome 1/2 sizes used by :func:`detect_assembly` (bp).  Sizes are for
#: the primary assembly and are exact for UCSC/Ensembl chromosome reports.
ASSEMBLY_CHROM_SIZES: Dict[str, Dict[str, int]] = {
    "hg38": {"1": 248956422, "2": 242193529},
    "hg19": {"1": 249250621, "2": 243199373},
    "mm10": {"1": 195471971, "2": 182113224},
    "mm9": {"1": 197195432, "2": 181976762},
    "mm39": {"1": 195154279, "2": 181755017},
    "dm6": {"1": 143726002, "2": 25286936},
    "danRer11": {"1": 59578282, "2": 45631630},
    "rn7": {"1": 282763074, "2": 267910886},
    "sacCer3": {"1": 230218, "2": 813184},
}

#: Ensembl species aliases used by ``TxDb="ensembl:<species>[@release]"``
ENSEMBL_SPECIES: Dict[str, tuple[str, str]] = {
    "human": ("homo_sapiens", "Homo_sapiens.GRCh38"),
    "mouse": ("mus_musculus", "Mus_musculus.GRCm39"),
    "rat": ("rattus_norvegicus", "Rattus_norvegicus.Rnor_6.0"),
    "zebrafish": ("danio_rerio", "Danio_rerio.GRCz11"),
    "fly": ("drosophila_melanogaster", "Drosophila_melanogaster.BDGP6.46"),
    "yeast": ("saccharomyces_cerevisiae", "Saccharomyces_cerevisiae.R64-1-1"),
}
DEFAULT_ENSEMBL_RELEASE = "110"

#: GENCODE aliases used by ``TxDb="gencode:<human|mouse>[@release]"``
GENCODE_SPECIES: Dict[str, tuple[str, str, str]] = {
    # name: (url directory, assembly, default release)
    "human": ("Gencode_human", "GRCh38", "44"),
    "mouse": ("Gencode_mouse", "GRCm39", "M33"),
}

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
    """True for a supported bundled or downloadable assembly name."""

    return isinstance(value, str) and value.strip().lower() in GENOME_ALIASES


def canonical_assembly(value: Any) -> str:
    """Map an assembly name or alias to the bundled assembly stem."""

    if not isinstance(value, str):
        raise ValueError(
            "TxDb must be a GeneAnnotation, a GTF path/DataFrame, a bundled "
            f"assembly name, or an ensembl:/gencode:/url: download spec; "
            f"supported assemblies: {', '.join(SUPPORTED_GENOMES)}"
        )
    name = value.strip().lower()
    if name not in GENOME_ALIASES:
        raise ValueError(
            f"unsupported TxDb assembly {value!r}; supported names: "
            f"{', '.join(SUPPORTED_GENOMES)}"
        )
    return GENOME_ALIASES[name]


def is_genome_spec(value: Any) -> bool:
    """True for ``ensembl:``, ``gencode:`` and ``url:`` annotation specs."""

    if not isinstance(value, str):
        return False
    return value.strip().lower().startswith(("ensembl:", "gencode:", "url:"))


def parse_genome_spec(spec: str) -> dict[str, str]:
    """Resolve a download spec to provider, release, URL and description.

    Accepted forms::

        ensembl:mouse
        ensembl:mus_musculus@110
        gencode:human@44
        url:https://example.org/annotation.gtf.gz
    """

    if not isinstance(spec, str):
        raise TypeError("genome spec must be a string")
    text = spec.strip()
    lower = text.lower()
    if lower.startswith("url:"):
        url = text[4:].strip()
        if not url.lower().startswith(("http://", "https://", "ftp://", "file://")):
            raise ValueError("url: specs must point at http(s), ftp or file URLs")
        return {
            "provider": "url",
            "release": "",
            "url": url,
            "description": f"GTF from {url}",
        }

    if lower.startswith("ensembl:"):
        body = text[len("ensembl:") :].strip()
        species, _, release = body.partition("@")
        key = species.strip().lower()
        if key not in ENSEMBL_SPECIES:
            matches = [
                name
                for name, (directory, _) in ENSEMBL_SPECIES.items()
                if directory == key
            ]
            if matches:
                key = matches[0]
            else:
                raise ValueError(
                    f"unknown Ensembl species {species!r}; choose from "
                    f"{', '.join(sorted(ENSEMBL_SPECIES))} or use url:<gtf>"
                )
        directory, stem = ENSEMBL_SPECIES[key]
        release = (release or DEFAULT_ENSEMBL_RELEASE).strip()
        url = (
            f"https://ftp.ensembl.org/pub/release-{release}/gtf/{directory}/"
            f"{stem}.{release}.gtf.gz"
        )
        return {
            "provider": "ensembl",
            "release": release,
            "url": url,
            "description": f"Ensembl {release} GTF for {directory}",
        }

    if lower.startswith("gencode:"):
        body = text[len("gencode:") :].strip()
        species, _, release = body.partition("@")
        key = species.strip().lower()
        if key not in GENCODE_SPECIES:
            matches = [
                name
                for name, (_, assembly, _) in GENCODE_SPECIES.items()
                if assembly.lower() == key
            ]
            if matches:
                key = matches[0]
            else:
                raise ValueError(
                    f"unknown GENCODE species {species!r}; choose from "
                    f"{', '.join(sorted(GENCODE_SPECIES))}"
                )
        directory, assembly, default_release = GENCODE_SPECIES[key]
        release = (release or default_release).strip()
        if key == "human":
            filename = f"gencode.v{release}.annotation.gtf.gz"
        else:
            filename = f"{assembly}.gencode.v{release}.annotation.gtf.gz"
        url = (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/"
            f"{directory}/release_{release}/{filename}"
        )
        return {
            "provider": "gencode",
            "release": release,
            "url": url,
            "description": f"GENCODE {release} GTF for {assembly}",
        }

    raise ValueError(
        "genome spec must start with ensembl:, gencode: or url: "
        "(for example ensembl:mouse@110 or url:https://.../genes.gtf.gz)"
    )


def _assembly_size_candidates(
    chrom_sizes: Dict[str, int],
    *,
    tolerance: float = 0.0005,
) -> list[str]:
    """Assemblies whose chromosome 1/2 sizes match ``chrom_sizes``."""

    lookup = {canonical_chrom(name): int(size) for name, size in chrom_sizes.items()}
    matches: list[str] = []
    for assembly, expected in ASSEMBLY_CHROM_SIZES.items():
        exact = [
            key
            for key, size in expected.items()
            if key in lookup
            and abs(int(lookup[key]) - int(size)) <= max(1.0, float(size) * tolerance)
        ]
        if len(exact) == len(expected):
            matches.append(assembly)
    return matches


def detect_assembly(
    peaks: Any = None,
    chrom_sizes: Optional[Dict[str, int]] = None,
    *,
    bigwig: Union[str, Path, None] = None,
    tolerance: float = 0.0005,
) -> str:
    """Infer a bundled assembly from chromosome names or chromosome sizes.

    ``chrom_sizes`` may come from a bigWig header (``bigwig=...`` also opens
    one directly) and is the most reliable signal.  With only a peak table the
    function is deliberately conservative: if the chromosome set is compatible
    with several assemblies (for example hg38 versus hg19) it raises a
    :class:`ValueError` listing the candidates instead of guessing.

    ``peaks`` accepts a DataFrame or any path understood by
    :func:`pychipseeker.io.read_peaks`.
    """

    names: list[str] = []
    max_end: Optional[int] = None

    if bigwig is not None:
        try:
            import pyBigWig
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "pyBigWig is required to inspect a bigWig header; install with "
                "`pip install pychipseeker[bigwig]`"
            ) from exc
        handle = pyBigWig.open(str(bigwig))
        if handle is None:
            raise ValueError(f"could not open bigWig: {bigwig}")
        try:
            chrom_sizes = dict(handle.chroms())
        finally:
            handle.close()

    if peaks is not None:
        if isinstance(peaks, pd.DataFrame):
            frame = peaks.copy()
        else:
            from .io import read_peaks

            frame = read_peaks(peaks)
        if "chrom" not in frame.columns:
            raise ValueError("peaks must contain a chrom column")
        names = [str(value) for value in pd.unique(frame["chrom"].astype(str))]
        if "end" in frame.columns:
            end_values = pd.to_numeric(frame["end"], errors="coerce")
            if end_values.notna().any():
                max_end = int(end_values.max())

    if chrom_sizes:
        candidates = _assembly_size_candidates(chrom_sizes, tolerance=tolerance)
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise ValueError(
                "chromosome sizes match several assemblies "
                f"({', '.join(candidates)}); pass an explicit TxDb/assembly"
            )
        names = names or list(chrom_sizes)
    if not names:
        raise ValueError(
            "cannot detect an assembly without peaks, chrom_sizes or bigwig"
        )

    keys = {canonical_chrom(name) for name in names}
    human = {str(i) for i in range(1, 23)} | {"x", "y", "m"}
    mouse = {str(i) for i in range(1, 20)} | {"x", "y", "m"}

    if keys & {"20", "21", "22"}:
        candidates = ["hg38", "hg19"]
        if max_end is not None and max_end > 249_000_000:
            candidates = ["hg19", "hg38"]
        raise ValueError(
            "chromosome names look human but cannot distinguish "
            f"{', '.join(candidates)} without chromosome sizes; pass "
            f"TxDb='{candidates[0]}' explicitly or use bigwig=... to detect it"
        )
    if keys <= mouse and (keys & {str(i) for i in range(1, 20)}):
        if keys <= human:
            raise ValueError(
                "chromosome names are compatible with both hg38/hg19 and "
                "mm10/mm39/mm9; pass an explicit TxDb or use bigwig=... to "
                "detect the assembly"
            )
        raise ValueError(
            "chromosome names look mouse-like, but mm10/mm39/mm9 cannot be "
            "distinguished without chromosome sizes; pass an explicit TxDb "
            "or use bigwig=... to detect it"
        )
    if keys <= human and any(key.isdigit() for key in keys):
        raise ValueError(
            "chromosome names are compatible with hg38/hg19; pass an explicit "
            "TxDb or use bigwig=... to detect the assembly"
        )
    if keys & {"i", "ii", "iii", "iv", "v"}:
        return "sacCer3"
    if keys & {"2l", "2r", "3l", "3r", "4", "x", "y", "m"} and "2l" in keys:
        return "dm6"
    raise ValueError(
        "could not detect an assembly from the supplied chromosome names "
        f"({chromosome_style(names)} style); pass TxDb=... explicitly"
    )


def _spec_digest(spec: str) -> str:
    import hashlib as _hashlib

    return _hashlib.sha1(spec.strip().encode("utf-8")).hexdigest()[:16]


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


def _write_npz(dest, payload: Dict[str, Any]) -> None:
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


def _encode_frame(frame: pd.DataFrame, table: str) -> Dict[str, Any]:
    specs = [[str(col), str(frame[col].dtype)] for col in frame.columns]
    out: Dict[str, Any] = {f"{table}__cols": json.dumps(specs).encode()}
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

    payload: Dict[str, Any] = {"__tables": json.dumps(list(TABLES)).encode()}
    for table in TABLES:
        frame = getattr(model, table)
        if drop:
            dropcols = [c for c in frame.columns if str(c) in drop]
            if dropcols:
                frame = frame.drop(columns=dropcols)
        payload.update(_encode_frame(frame, table))
    _write_npz(dest, payload)


#: Public download sources for genome bundles, tried in order.  Release
#: assets live on github.com (reliable worldwide); the raw.githubusercontent
#: URL is a fallback.  ``PYCHIPSEEKER_GENOME_BASE`` (if set) is tried first,
#: e.g. ``file:///D:/mirror`` or an internal mirror for air-gapped setups.
RELEASE_BUNDLE_BASE = (
    "https://github.com/summus-kong/pychipseeker/releases/download/"
    "genome-data-2026.10"
)
RAW_BUNDLE_BASE = (
    "https://raw.githubusercontent.com/summus-kong/pychipseeker/main/"
    "src/pychipseeker/data"
)


def _bundle_bases() -> tuple[str, ...]:
    bases = []
    env_base = os.environ.get("PYCHIPSEEKER_GENOME_BASE")
    if env_base:
        bases.append(env_base.rstrip("/"))
    bases += [RELEASE_BUNDLE_BASE, RAW_BUNDLE_BASE]
    return tuple(dict.fromkeys(bases))


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
            files("pychipseeker").joinpath("data").joinpath(f"genome_{name}.npz")
        ) as path:
            return Path(path) if Path(path).is_file() else None
    except (FileNotFoundError, ModuleNotFoundError):
        return None


def _installed_provenance(name: str) -> Optional[dict]:
    """Provenance JSON shipped inside the installed package, if any."""

    from importlib.resources import as_file, files

    try:
        with as_file(
            files("pychipseeker").joinpath("data").joinpath(f"genome_{name}.json")
        ) as path:
            if not Path(path).is_file():
                return None
            with open(path, encoding="utf-8") as fh:
                return json.load(fh)
    except (FileNotFoundError, ModuleNotFoundError):
        return None


def _cache_bundle(name: str) -> Optional[Path]:
    """Bundle npz previously fetched into the user cache dir, if any."""

    path = cache_dir() / f"genome_{name}.npz"
    if path.is_file() and path.stat().st_size > 0:
        return path
    return None


def _download_file(url: str, dest: Path, tries: int = 4) -> None:
    """Stream ``url`` into ``dest`` with retries for flaky connections."""

    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Optional[Exception] = None
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=180) as resp, open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh)
            break
        except Exception as exc:  # noqa: BLE001 - retry any network error
            last_error = exc
            if attempt < tries:
                print(
                    f"pychipseeker: genome bundle download attempt "
                    f"{attempt} failed ({exc}); retrying",
                    file=sys.stderr,
                )
                time.sleep(5 * attempt)
    else:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"could not download {url}") from last_error
    tmp.replace(dest)


def _download_bundle(name: str) -> Path:
    """Fetch a genome bundle (npz + provenance JSON) into the cache dir.

    Tries each source in :func:`_bundle_bases`.  The provenance JSON is
    fetched first (best effort) so the npz can be verified against its
    recorded sha256; keeping the JSON also makes
    :func:`genome_provenance` work offline afterwards.
    """

    dest = cache_dir() / f"genome_{name}.npz"
    json_dest = cache_dir() / f"genome_{name}.json"
    expected: Optional[str] = None
    last_error: Optional[Exception] = None
    for base in _bundle_bases():
        if expected is None:
            try:
                _download_file(f"{base}/genome_{name}.json", json_dest, tries=2)
                with open(json_dest, encoding="utf-8") as fh:
                    expected = json.load(fh).get("bundle_sha256")
            except Exception:  # noqa: BLE001 - provenance is a nicety
                expected = None
        try:
            _download_file(f"{base}/genome_{name}.npz", dest)
        except Exception as exc:  # noqa: BLE001 - fall through to next source
            last_error = exc
            continue
        if expected:
            digest = hashlib.sha256()
            with open(dest, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    digest.update(chunk)
            if digest.hexdigest() != expected:
                dest.unlink(missing_ok=True)
                last_error = RuntimeError(
                    f"genome bundle {dest.name} from {base} failed its sha256 check"
                )
                continue
        return dest
    raise RuntimeError(
        f"could not fetch genome bundle {name!r} from any source; download "
        f"genome_{name}.npz + genome_{name}.json manually into "
        f"{cache_dir()} (see NOTICE.md for the public URLs), or pass a local "
        "GTF instead of TxDb="
    ) from last_error


def _bundle_path(name: str) -> Path:
    """Locate a bundle npz: installed data -> cache -> public download."""

    path = _installed_bundle(name) or _cache_bundle(name)
    if path is None:
        path = _download_bundle(name)
    return path


# ---------------------------------------------------------------------------
# on-demand GTF downloads (new assemblies and ensembl:/gencode:/url: specs)
# ---------------------------------------------------------------------------

_loaded_specs: Dict[str, GeneAnnotation] = {}


def _write_json(path: Path, record: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
        fh.write("\n")
    tmp.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _gtf_suffix(url: str) -> str:
    path = url.lower().split("?", 1)[0]
    if path.endswith((".gtf.gz", ".gff.gz", ".gff3.gz")):
        return ".gtf.gz"
    return ".gtf"


def _gtf_cache_paths(spec: str) -> tuple[Path, Path, dict[str, str]]:
    info = parse_genome_spec(spec)
    digest = _spec_digest(spec)
    root = cache_dir() / "gtf"
    root.mkdir(parents=True, exist_ok=True)
    return (
        root / f"spec_{digest}{_gtf_suffix(info['url'])}",
        root / f"spec_{digest}.npz",
        info,
    )


def download_gtf(
    spec: str,
    dest: Optional[Union[str, Path]] = None,
) -> Path:
    """Download the GTF referenced by an ``ensembl:``/``gencode:``/``url:`` spec.

    The file is stored under ``<cache>/gtf`` and reused on later calls.
    ``dest`` overrides the cache location.
    """

    cache_path, _, info = _gtf_cache_paths(spec)
    path = cache_path if dest is None else Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.stat().st_size > 0:
        return path
    _download_file(info["url"], path)
    return path


def _npz_model(path: Path) -> GeneAnnotation:
    npz = np.load(str(path), allow_pickle=False)
    try:
        frames = [_rebuild_frame(npz, table) for table in TABLES]
    finally:
        npz.close()
    return GeneAnnotation(*frames, copy=False)


def load_annotation(spec: str, refresh: bool = False) -> GeneAnnotation:
    """Load a GTF from an ``ensembl:``/``gencode:``/``url:`` spec.

    The GTF is downloaded into the cache, parsed once and stored as a compact
    npz bundle, so later calls (including new processes) avoid re-parsing.
    """

    key = spec.strip()
    if key in _loaded_specs and not refresh:
        return _loaded_specs[key]

    _, npz_path, info = _gtf_cache_paths(key)
    if npz_path.is_file() and npz_path.stat().st_size > 0 and not refresh:
        model = _npz_model(npz_path)
        _loaded_specs[key] = model
        return model

    gtf_path = download_gtf(key)
    model = GeneAnnotation.from_gtf(gtf_path)
    tmp = npz_path.with_suffix(".npz.part")
    save_bundle(model, tmp)
    tmp.replace(npz_path)
    _write_json(
        npz_path.with_suffix(".json"),
        {
            "assembly": key,
            "provider": info["provider"],
            "release": info["release"],
            "url": info["url"],
            "gtf": str(gtf_path),
            "gtf_sha256": _sha256(gtf_path),
            "bundle_sha256": _sha256(npz_path),
            "downloaded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "counts": {table: len(getattr(model, table)) for table in TABLES},
        },
    )
    _loaded_specs[key] = model
    return model


def _build_manifest_bundle(name: str) -> Path:
    """Download a manifest GTF, parse it and cache ``genome_<name>.npz``."""

    source = GENOME_SOURCES[name]
    gtf_path = cache_dir() / f"genome_{name}{_gtf_suffix(source['url'])}"
    if not gtf_path.is_file() or gtf_path.stat().st_size == 0:
        _download_file(source["url"], gtf_path)
    model = GeneAnnotation.from_gtf(gtf_path)

    npz_path = cache_dir() / f"genome_{name}.npz"
    tmp = npz_path.with_suffix(".npz.part")
    save_bundle(model, tmp)
    tmp.replace(npz_path)
    _write_json(
        npz_path.with_suffix(".json"),
        {
            "assembly": name,
            "provider": source["provider"],
            "release": source["release"],
            "url": source["url"],
            "description": source["description"],
            "gtf": str(gtf_path),
            "gtf_sha256": _sha256(gtf_path),
            "bundle_sha256": _sha256(npz_path),
            "downloaded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "counts": {table: len(getattr(model, table)) for table in TABLES},
        },
    )
    return npz_path


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
        arr = pd.array(strings, dtype=tag)
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
    it.  ``hg38``/``hg19``/``mm10``/``mm9`` ship inside the package; ``dm6``,
    ``danRer11``, ``rn7``, ``mm39`` and ``sacCer3`` download their Ensembl
    GTF on first use and cache the parsed tables.
    """

    name = canonical_assembly(assembly)
    if name in _loaded and not refresh:
        return _loaded[name]

    if name in BUNDLED_GENOMES:
        path: Optional[Path] = _bundle_path(name)
    else:
        path = _cache_bundle(name)
        if path is None and name in GENOME_SOURCES:
            path = _build_manifest_bundle(name)
        elif path is None:
            # custom/test bundle registered through GENOME_ALIASES
            path = _bundle_path(name)

    if path is None:
        raise FileNotFoundError(f"no genome bundle found for {name!r}")
    model = _npz_model(path)
    _loaded[name] = model
    return model


def genome_provenance(assembly: Union[str, None]) -> dict:
    """Return the provenance record of a genome annotation.

    For bundled assemblies, the JSON beside the bundle records the GENCODE
    release, source URL and checksums.  For downloadable assemblies that have
    not been fetched yet, the manifest entry is returned with
    ``"downloaded": False``.
    """

    if is_genome_spec(assembly):
        spec = str(assembly).strip()
        _, npz_path, info = _gtf_cache_paths(spec)
        json_path = npz_path.with_suffix(".json")
        if json_path.is_file():
            with open(json_path, encoding="utf-8") as fh:
                return json.load(fh)
        return {
            "assembly": spec,
            "provider": info["provider"],
            "release": info["release"],
            "url": info["url"],
            "description": info["description"],
            "downloaded": False,
            "cache_dir": str(cache_dir()),
        }

    name = canonical_assembly(assembly)
    installed = _installed_bundle(name)
    candidates = []
    if installed is not None:
        candidates.append(installed.with_suffix(".json"))
    candidates.append(cache_dir() / f"genome_{name}.json")
    for candidate in candidates:
        if candidate.is_file():
            with open(candidate, encoding="utf-8") as fh:
                return json.load(fh)
    shipped = _installed_provenance(name)
    if shipped is not None:
        return shipped
    if name in GENOME_SOURCES:
        return {
            "assembly": name,
            **GENOME_SOURCES[name],
            "downloaded": False,
            "cache_dir": str(cache_dir()),
        }
    _download_bundle(name)
    with open(cache_dir() / f"genome_{name}.json", encoding="utf-8") as fh:
        return json.load(fh)


def supported_genomes() -> tuple[str, ...]:
    """Assemblies accepted by ``TxDb=``.

    The first four are bundled; the remaining entries are built from a
    download on first use.
    """

    return SUPPORTED_GENOMES


def download_genome(assembly: Union[str, None]) -> Path:
    """Fetch a genome annotation into the cache and return the npz path.

    An installed or already cached bundle is returned untouched.  For the
    downloadable assemblies the Ensembl GTF is fetched, parsed and cached as
    a compact npz next to the GTF.  Network failures raise
    :class:`RuntimeError` with the cache directory in the message.
    """

    name = canonical_assembly(assembly)
    path = _installed_bundle(name) or _cache_bundle(name)
    if path is not None:
        return path
    if name in BUNDLED_GENOMES:
        return _download_bundle(name)
    if name in GENOME_SOURCES:
        return _build_manifest_bundle(name)
    return _download_bundle(name)


def cached_genomes() -> list[dict[str, object]]:
    """Describe installed and downloaded bundles for ``pychipseeker cache``."""

    root = cache_dir()
    names = set(SUPPORTED_GENOMES)
    for bundle_path in root.glob("genome_*.npz"):
        stem = bundle_path.stem
        if stem.startswith("genome_"):
            names.add(stem[len("genome_") :])
    rows: list[dict[str, object]] = []

    # On-demand ensembl:/gencode:/url: specs are cached as gtf/spec_*.npz.
    gtf_root = root / "gtf"
    if gtf_root.is_dir():
        for npz_path in sorted(gtf_root.glob("spec_*.npz")):
            spec_record: dict[str, object] = {}
            json_path = npz_path.with_suffix(".json")
            if json_path.is_file():
                try:
                    with open(json_path, encoding="utf-8") as fh:
                        spec_record = json.load(fh)
                except (OSError, ValueError):
                    spec_record = {}
            rows.append(
                {
                    "assembly": spec_record.get("assembly", npz_path.stem),
                    "source": "spec-cache",
                    "provider": spec_record.get("provider", ""),
                    "release": spec_record.get("release", ""),
                    "description": spec_record.get("description", ""),
                    "size_mb": round(npz_path.stat().st_size / 1e6, 1),
                    "path": str(npz_path),
                    "cache_dir": str(root),
                }
            )
    for name in sorted(names):
        installed = _installed_bundle(name)
        cached = _cache_bundle(name)
        source_kind = (
            "installed" if installed is not None else ("cache" if cached else "")
        )
        path: Optional[Path] = installed or cached
        source = GENOME_SOURCES.get(name, {})
        record: dict[str, object] = {}
        if path is not None:
            provenance = path.with_suffix(".json")
            if provenance.is_file():
                try:
                    with open(provenance, encoding="utf-8") as fh:
                        record = json.load(fh)
                except (OSError, ValueError):
                    record = {}
        rows.append(
            {
                "assembly": name,
                "source": source_kind,
                "provider": record.get(
                    "provider",
                    source.get(
                        "provider",
                        "gencode" if name in BUNDLED_GENOMES else "",
                    ),
                ),
                "release": record.get("release", source.get("release", "")),
                "description": record.get(
                    "description", source.get("description", "")
                ),
                "size_mb": (
                    round(path.stat().st_size / 1e6, 1) if path is not None else None
                ),
                "path": str(path) if path is not None else "",
                "cache_dir": str(cache_dir()),
            }
        )
    return rows


def clear_genome_cache(assembly: Union[str, None] = None) -> list[Path]:
    """Delete downloaded bundles from the cache directory.

    Installed package data is never touched.  With ``assembly=None`` every
    cached bundle, provenance JSON, downloaded GTF and parsed spec is
    removed.  Returns the deleted paths.
    """

    root = cache_dir()
    gtf_root = root / "gtf"
    if assembly is None:
        targets = sorted(root.glob("genome_*.npz")) + sorted(root.glob("genome_*.json"))
        targets += sorted(root.glob("genome_*.gtf*"))
        targets += sorted(root.glob("parsed_*.pkl"))
        if gtf_root.is_dir():
            targets += sorted(gtf_root.glob("spec_*"))
        _loaded.clear()
        _loaded_specs.clear()
    elif is_genome_spec(assembly):
        spec = str(assembly).strip()
        _, npz_path, _ = _gtf_cache_paths(spec)
        targets = sorted((root / "gtf").glob(f"{npz_path.stem}.*"))
        _loaded_specs.pop(spec, None)
    else:
        name = canonical_assembly(assembly)
        targets = [
            root / f"genome_{name}.npz",
            root / f"genome_{name}.json",
            root / f"genome_{name}.gtf.gz",
            root / f"genome_{name}.gtf",
        ]
        targets += sorted(root.glob(f"parsed_{name}_*.pkl"))
        _loaded.pop(name, None)
    removed: list[Path] = []
    for path in targets:
        if path.exists():
            path.unlink()
            removed.append(path)
    if assembly is None and gtf_root.is_dir() and not any(gtf_root.iterdir()):
        gtf_root.rmdir()
    return removed


def prune_genome_cache() -> list[Path]:
    """Remove incomplete cache entries and orphaned provenance JSON files.

    Partial downloads (``*.part``), zero-byte bundles and JSON records whose
    npz bundle is missing are leftovers from interrupted downloads.  Returns
    the deleted paths.
    """

    root = cache_dir()
    gtf_root = root / "gtf"
    targets: list[Path] = []
    for pattern in ("*.part", "genome_*.npz", "genome_*.json", "genome_*.gtf", "genome_*.gtf.gz"):
        targets.extend(sorted(root.glob(pattern)))
    if gtf_root.is_dir():
        for pattern in ("*.part", "spec_*.npz", "spec_*.json", "spec_*.gtf", "spec_*.gtf.gz"):
            targets.extend(sorted(gtf_root.glob(pattern)))

    removed: list[Path] = []
    for path in dict.fromkeys(targets):
        if not path.exists():
            continue
        suffix = path.suffix.lower()
        incomplete = suffix == ".part" or path.stat().st_size == 0
        orphan = False
        if suffix == ".json":
            orphan = not path.with_suffix(".npz").is_file()
        if incomplete or orphan:
            path.unlink()
            removed.append(path)
    return removed
