"""Small process-pool helpers for the heavy analysis paths.

The library defaults to ``n_jobs=1`` so that importing and calling it never
spawns processes.  Passing ``n_jobs=-1`` uses every available core, and any
positive value selects an explicit worker count.

The helpers intentionally keep the work units coarse (whole chromosomes or
row chunks) because peak files can contain millions of rows and per-row tasks
would spend more time pickling than analysing.  All worker functions are
module-level so the pool also works with the ``spawn`` start method used on
Windows and macOS.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Callable, Optional, Sequence, TypeVar

T = TypeVar("T")
R = TypeVar("R")


def available_cores() -> int:
    """Return the number of usable CPUs, never less than one."""

    return max(1, int(os.cpu_count() or 1))


def process_start_method() -> str:
    """Return the multiprocessing start method used for worker pools."""

    return mp.get_start_method()


def effective_jobs(n_jobs: Optional[int], *, tasks: int = 1) -> int:
    """Resolve a user-facing ``n_jobs`` value into a worker count.

    ``None`` and ``1`` mean sequential execution; ``-1`` means all cores.
    The result is capped by the number of available tasks so a small input
    never starts more workers than it can use.
    """

    if n_jobs is None:
        return 1
    try:
        value = int(n_jobs)
    except (TypeError, ValueError) as exc:
        raise ValueError("n_jobs must be an integer") from exc
    if value == -1:
        workers = available_cores()
    elif value > 0:
        workers = value
    else:
        raise ValueError("n_jobs must be -1 or a positive integer")
    return max(1, min(workers, max(1, int(tasks))))


def chunk_bounds(n_items: int, n_chunks: int) -> list[tuple[int, int]]:
    """Split ``n_items`` rows into at most ``n_chunks`` contiguous ranges."""

    n_items = int(n_items)
    n_chunks = max(1, int(n_chunks))
    if n_items <= 0:
        return []
    n_chunks = min(n_chunks, n_items)
    size = (n_items + n_chunks - 1) // n_chunks
    return [
        (start, min(start + size, n_items))
        for start in range(0, n_items, size)
    ]


def parallel_map(
    func: Callable[[T], R],
    items: Sequence[T],
    *,
    n_jobs: Optional[int],
    initializer: Optional[Callable[..., None]] = None,
    initargs: Sequence[object] = (),
) -> list[R]:
    """Map ``func`` over ``items`` in order, using processes when requested.

    Results always follow the input order, which keeps annotation and matrix
    rows aligned with their peaks/features.  The sequential path is used for a
    single item, ``n_jobs<=1``, and inside daemonic processes (which cannot
    start children of their own).
    """

    workers = effective_jobs(n_jobs, tasks=len(items))
    if workers <= 1 or len(items) <= 1:
        if initializer is not None:
            initializer(*initargs)
        return [func(item) for item in items]

    if mp.current_process().daemon:
        # A daemonic worker cannot spawn children; stay correct rather than
        # failing deep inside multiprocessing.
        if initializer is not None:
            initializer(*initargs)
        return [func(item) for item in items]

    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initializer,
        initargs=tuple(initargs),
    ) as pool:
        return list(pool.map(func, items))
