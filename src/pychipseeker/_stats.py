"""Small statistical helpers shared by the compatibility and enrichment code.

SciPy is optional, so these routines are implemented directly.  Keeping a
single copy here avoids the compatibility and enrichment modules drifting
apart while they evolve independently.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np


def log_comb(n: int, k: int) -> float:
    """Return ``log(C(n, k))``, or ``-inf`` outside the valid range."""

    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """P(X >= k) for a hypergeometric distribution.

    ``M`` is the background size, ``n`` the query size and ``N`` the gene-set
    size.  The sum is evaluated in log space so that large backgrounds remain
    numerically stable.
    """

    if k <= 0:
        return 1.0
    upper = min(n, N)
    if k > upper:
        return 0.0
    terms = []
    for i in range(int(k), int(upper) + 1):
        terms.append(math.exp(log_comb(n, i) + log_comb(M - n, N - i) - log_comb(M, N)))
    return float(min(max(math.fsum(terms), 0.0), 1.0))


def p_adjust(pvalues: Sequence[float], method: str = "BH") -> np.ndarray:
    """Return Benjamini-Hochberg, Bonferroni or unadjusted p-values.

    ``method`` accepts the spellings used by ChIPseeker (``"BH"``,
    ``"bonferroni"``, ``"none"``) and by :mod:`pychipseeker.enrichment`
    (``"fdr"``, ``"bonf"``, ``"raw"``); matching is case-insensitive and
    whitespace is ignored.
    """

    p = np.asarray(pvalues, dtype=float)
    key = str(method).strip().lower()
    if key in {"none", "no", "nop", "raw"}:
        return p
    n = len(p)
    if n == 0:
        return p
    if key in {"bonferroni", "bonf"}:
        return np.minimum(p * n, 1.0)
    if key not in {"bh", "fdr", "benjamini-hochberg"}:
        raise ValueError(
            "p_adjust/pAdjustMethod must be one of BH/fdr, Bonferroni or none"
        )
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out
