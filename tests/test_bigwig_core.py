"""Unit tests for the shared bigWig integration helpers."""

from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from pychipseeker import _bigwig
from pychipseeker._bigwig import (
    TSS_MODE,
    WINDOW_MODE,
    BlockJob,
    IntervalIntegrator,
    plan_blocks,
    run_blocks,
    stats_matrix,
    tss_matrix,
    window_matrix,
)


def test_interval_integrator_rejects_bad_input() -> None:
    for value in ([], np.empty((0, 3)), np.asarray([0.0, 1.0, 2.0]), np.asarray([[0, 1]])):
        with pytest.raises(ValueError, match="non-empty"):
            IntervalIntegrator(value)


def test_interval_integrator_integrates_with_overlaps() -> None:
    integrator = IntervalIntegrator([(10, 20, 2.0), (0, 10, 1.0)])
    covered, integral = integrator.integrate(
        np.asarray([-5, 0, 5, 10, 15, 20, 25], dtype=np.int64)
    )
    np.testing.assert_allclose(covered, [0, 0, 5, 10, 15, 20, 20])
    np.testing.assert_allclose(integral, [0, 0, 5, 10, 20, 30, 30])


def test_tss_matrix_marks_uncovered_bins_and_applies_floor() -> None:
    integrator = IntervalIntegrator([(0, 100, 1.0)])
    arr = tss_matrix(integrator, np.asarray([50]), up=10, down=10, n_bins=4)
    np.testing.assert_allclose(arr, np.ones((1, 4)))

    floored = tss_matrix(
        integrator, np.asarray([50]), up=10, down=10, n_bins=4, min_value=2.0
    )
    np.testing.assert_allclose(floored, np.zeros((1, 4)))

    missing = tss_matrix(integrator, np.asarray([500]), up=10, down=10, n_bins=4)
    assert np.isnan(missing).all()
    zeroed = tss_matrix(
        integrator, np.asarray([500]), up=10, down=10, n_bins=4, min_value=0.0
    )
    assert np.isnan(zeroed).all()


def test_window_matrix_sum_mean_coverage_and_errors() -> None:
    integrator = IntervalIntegrator([(0, 10, 2.0)])
    starts = np.asarray([0])
    ends = np.asarray([20])

    np.testing.assert_allclose(
        window_matrix(integrator, starts, ends, 2, "sum"), [[20.0, 0.0]]
    )
    np.testing.assert_allclose(
        window_matrix(integrator, starts, ends, 2, "mean"), [[2.0, np.nan]]
    )
    np.testing.assert_allclose(
        window_matrix(integrator, starts, ends, 2, "coverage"), [[1.0, 0.0]]
    )

    zero_width = window_matrix(
        integrator, np.asarray([5]), np.asarray([5]), 2, "coverage"
    )
    assert np.isnan(zero_width).all()

    with pytest.raises(ValueError, match="unsupported value_type"):
        window_matrix(integrator, starts, ends, 2, "max")


def test_stats_matrix_handles_none_values() -> None:
    calls = []

    class Handle:
        def stats(self, chrom, start, end, nBins, type):
            calls.append((chrom, start, end, nBins, type))
            if start == 10:
                return None
            return [1.0, None, 3.0, 4.0]

    out = stats_matrix(
        Handle(), "chr1", np.asarray([0, 10]), np.asarray([10, 20]), 4, "max"
    )
    np.testing.assert_allclose(out[0], [1.0, 0.0, 3.0, 4.0])
    np.testing.assert_allclose(out[1], np.zeros(4))
    assert calls[1][4] == "max"


def test_plan_blocks_orders_dense_chromosomes_first() -> None:
    assert plan_blocks(
        np.asarray(["chr1"]),
        np.asarray([0]),
        np.asarray([0]),
        np.asarray([10]),
        np.asarray([], dtype=np.int64),
        n_blocks=2,
    ) == []

    chroms = np.asarray(["chr1", "chr1", "chr2", "chr2", "chr2"])
    order = np.asarray([4, 3, 2, 1, 0])
    starts = np.asarray([100, 0, 300, 200, 100])
    ends = np.asarray([110, 10, 310, 210, 110])
    tasks = plan_blocks(
        chroms,
        order,
        starts,
        ends,
        np.arange(5),
        n_blocks=2,
        min_block=1,
    )
    assert [(chrom, rows.tolist(), start, end) for chrom, rows, start, end in tasks] == [
        ("chr2", [4, 3, 2], 100, 310),
        ("chr1", [1, 0], 0, 110),
    ]


def _fake_handle(intervals, stats=None):
    class Handle:
        def intervals(self, chrom, start, end):
            return [(a, b, v) for a, b, v in intervals if b > start and a < end]

        def stats(self, chrom, start, end, nBins, type):
            return list(stats or [1.0] * int(nBins))

    return Handle()


def test_worker_run_tss_mode_orients_minus_strand(monkeypatch: pytest.MonkeyPatch) -> None:
    job = BlockJob(
        kind=TSS_MODE,
        bigwig="fake.bw",
        n_bins=4,
        up=10,
        down=10,
        dtype="float32",
        tss=np.asarray([5, 15]),
        strands=np.asarray(["+", "-"]),
    )
    monkeypatch.setattr(_bigwig, "_WORKER_JOB", job)
    monkeypatch.setattr(
        _bigwig, "_WORKER_HANDLE", _fake_handle([(0, 10, 1.0), (10, 20, 3.0)])
    )

    rows, matrix = _bigwig._worker_run(("chr1", np.asarray([0, 1]), 0, 20))
    assert rows.tolist() == [0, 1]
    assert matrix.dtype == np.float32
    # The minus-strand row is the forward profile flipped.
    np.testing.assert_allclose(matrix[0], [np.nan, 1.0, 1.0, 3.0])
    np.testing.assert_allclose(matrix[1], [np.nan, 3.0, 3.0, 1.0])

    bad_job = BlockJob(kind=TSS_MODE, bigwig="fake.bw", n_bins=4)
    monkeypatch.setattr(_bigwig, "_WORKER_JOB", bad_job)
    with pytest.raises(RuntimeError, match="tss/strands"):
        _bigwig._worker_run(("chr1", np.asarray([0]), 0, 20))


def test_worker_run_window_mode_and_empty_region(monkeypatch: pytest.MonkeyPatch) -> None:
    job = BlockJob(
        kind=WINDOW_MODE,
        bigwig="fake.bw",
        n_bins=2,
        value_type="sum",
        starts=np.asarray([0, 10]),
        ends=np.asarray([20, 20]),
        strands=np.asarray(["+", "+"]),
    )
    monkeypatch.setattr(_bigwig, "_WORKER_JOB", job)
    monkeypatch.setattr(_bigwig, "_WORKER_HANDLE", _fake_handle([(0, 20, 1.0)]))
    rows, matrix = _bigwig._worker_run(("chr1", np.asarray([0, 1]), 0, 20))
    np.testing.assert_allclose(matrix, [[10.0, 10.0], [5.0, 5.0]])

    stats_job = BlockJob(
        kind=WINDOW_MODE,
        bigwig="fake.bw",
        n_bins=2,
        value_type="max",
        starts=np.asarray([0]),
        ends=np.asarray([20]),
        strands=np.asarray(["+"]),
        dtype="float32",
    )
    monkeypatch.setattr(_bigwig, "_WORKER_JOB", stats_job)
    monkeypatch.setattr(
        _bigwig,
        "_WORKER_HANDLE",
        _fake_handle([(0, 20, 1.0)], stats=[1.0, 3.0]),
    )
    _, stats_matrix_out = _bigwig._worker_run(("chr1", np.asarray([0]), 0, 20))
    assert stats_matrix_out.dtype == np.float32
    np.testing.assert_allclose(stats_matrix_out, [[1.0, 3.0]])

    monkeypatch.setattr(
        _bigwig, "_WORKER_HANDLE", _fake_handle([])
    )
    rows, empty = _bigwig._worker_run(("chr2", np.asarray([3]), 0, 20))
    assert rows.tolist() == [3]
    assert empty is None


def test_worker_run_guards(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(_bigwig, "_WORKER_JOB", None)
    monkeypatch.setattr(_bigwig, "_WORKER_HANDLE", None)
    with pytest.raises(RuntimeError, match="not initialised"):
        _bigwig._worker_run(("chr1", np.asarray([0]), 0, 10))

    job = BlockJob(
        kind=WINDOW_MODE,
        bigwig="fake.bw",
        n_bins=2,
        value_type="sum",
    )
    monkeypatch.setattr(_bigwig, "_WORKER_JOB", job)
    monkeypatch.setattr(_bigwig, "_WORKER_HANDLE", _fake_handle([(0, 10, 1.0)]))
    with pytest.raises(RuntimeError, match="starts/ends/strands"):
        _bigwig._worker_run(("chr1", np.asarray([0]), 0, 10))


def test_worker_init_and_close_manage_handle(monkeypatch: pytest.MonkeyPatch) -> None:
    closed = []

    class Handle:
        def close(self):
            closed.append(True)

    fake = types.ModuleType("pyBigWig")
    fake.open = lambda path: Handle()
    monkeypatch.setitem(sys.modules, "pyBigWig", fake)

    job = BlockJob(kind=TSS_MODE, bigwig="fake.bw", n_bins=2)
    _bigwig._worker_init(job)
    assert _bigwig._WORKER_JOB is job
    assert _bigwig._WORKER_HANDLE is not None
    _bigwig._worker_close()
    assert closed == [True]
    assert _bigwig._WORKER_HANDLE is None
    _bigwig._worker_close()
    assert closed == [True]


def test_run_blocks_sequential_initialises_and_closes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events = []

    def fake_init(job):
        events.append(("init", job.bigwig))

    def fake_run(task):
        events.append(("run", task[0]))
        return task[1], None

    def fake_close():
        events.append(("close", None))

    monkeypatch.setattr(_bigwig, "_worker_init", fake_init)
    monkeypatch.setattr(_bigwig, "_worker_run", fake_run)
    monkeypatch.setattr(_bigwig, "_worker_close", fake_close)

    job = BlockJob(kind=TSS_MODE, bigwig="fake.bw", n_bins=2)
    tasks = [("chr1", np.asarray([0]), 0, 10), ("chr2", np.asarray([1]), 0, 10)]
    result = run_blocks(job, tasks, n_jobs=1)
    assert [rows.tolist() for rows, _ in result] == [[0], [1]]
    assert events == [
        ("init", "fake.bw"),
        ("run", "chr1"),
        ("run", "chr2"),
        ("close", None),
    ]


def test_run_blocks_uses_parallel_map(monkeypatch: pytest.MonkeyPatch) -> None:
    from pychipseeker import parallel

    calls = []

    monkeypatch.setattr(parallel, "effective_jobs", lambda n_jobs, tasks: 2)

    def fake_map(func, items, n_jobs, initializer, initargs):
        calls.append((n_jobs, initargs))
        initializer(*initargs)
        return [(item, None) for item in items]

    monkeypatch.setattr(parallel, "parallel_map", fake_map)
    monkeypatch.setattr(_bigwig, "_worker_init", lambda job: None)

    job = BlockJob(kind=TSS_MODE, bigwig="fake.bw", n_bins=2)
    tasks = [("chr1", np.asarray([0]), 0, 10), ("chr2", np.asarray([1]), 0, 10)]
    result = run_blocks(job, tasks, n_jobs=2)
    assert len(result) == 2
    assert calls == [(2, (job,))]
