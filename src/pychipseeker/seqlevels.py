"""Chromosome naming harmonisation across UCSC and Ensembl annotations.

Peak files usually use UCSC names (``chr1``, ``chrM``), while Ensembl GTF
annotations use bare names (``1``, ``MT``).  Matching them by string equality
silently drops otherwise valid peaks.  The helpers here canonicalise names for
comparison without rewriting unrelated chromosomes.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd

_MITOCHONDRIA = {"m", "mt", "mito", "mitochondria", "chrm", "chrmt"}


def canonical_chrom(value: object) -> str:
    """Return a comparison key for one chromosome name.

    The key removes a leading ``chr`` (case-insensitive) and maps common
    mitochondrial spellings to ``"m"``.  It is only used for matching; the
    original name is preserved in the data.
    """

    text = str(value).strip().lower()
    if text.startswith("chr"):
        text = text[3:]
    if text in _MITOCHONDRIA:
        return "m"
    return text


def chrom_aliases(name: object) -> tuple[str, ...]:
    """Return common spellings for one chromosome name, preferred first."""

    text = str(name).strip()
    key = canonical_chrom(text)
    bare = key.upper() if key in {"x", "y", "m"} else key
    prefixes = ("chr", "Chr", "CHR")
    aliases = [text]
    if text.lower().startswith("chr"):
        aliases.append(text[3:])
    else:
        aliases.extend(prefix + bare for prefix in prefixes)
    if key == "m":
        aliases.extend(["M", "MT", "chrM", "chrMT", "ChrM", "ChrMT"])
    else:
        aliases.extend([bare, bare.upper(), f"chr{bare}", f"Chr{bare}"])
    # Preserve order and remove duplicates.
    return tuple(dict.fromkeys(alias for alias in aliases if alias))


def alias_map(
    query: Iterable[object],
    target: Iterable[object],
) -> dict[str, str]:
    """Map query chromosome names onto target names where possible."""

    target_names = [str(name) for name in target]
    exact = set(target_names)
    by_key: dict[str, str] = {}
    for name in target_names:
        by_key.setdefault(canonical_chrom(name), name)

    mapping: dict[str, str] = {}
    for raw in query:
        name = str(raw)
        if name in exact:
            mapping[name] = name
            continue
        match = by_key.get(canonical_chrom(name))
        if match is not None:
            mapping[name] = match
    return mapping


def harmonize_frame(
    frame: pd.DataFrame,
    target: Iterable[object],
    column: str = "chrom",
    *,
    copy: bool = True,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Rename ``frame[column]`` values to matching target chromosome names.

    Returns the (possibly copied) frame and the applied mapping.  Values with
    no match are left untouched so callers can report them explicitly.
    """

    if column not in frame.columns:
        raise KeyError(f"frame has no {column!r} column")
    mapping = alias_map(frame[column].astype(str).unique(), target)
    if not mapping:
        return (frame.copy() if copy else frame), mapping
    out = frame.copy() if copy else frame
    original = out[column].astype(str)
    mapped = original.map(mapping)
    out[column] = mapped.where(mapped.notna(), original)
    return out, mapping


def chromosome_style(names: Sequence[object]) -> str:
    """Return ``"ucsc"``, ``"ensembl"`` or ``"mixed"`` for a name collection."""

    values = [str(name).strip() for name in names if str(name).strip()]
    if not values:
        return "mixed"
    prefixed = sum(value.lower().startswith("chr") for value in values)
    if prefixed == len(values):
        return "ucsc"
    if prefixed == 0:
        return "ensembl"
    return "mixed"
