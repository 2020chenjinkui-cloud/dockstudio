"""Unit tests for the serial job queue (mock engine, no real vina)."""

from __future__ import annotations

import asyncio
import threading
import time
from pathlib import Path

import pytest

import dockstudio.core.queue as queue_mod
from dockstudio.core.project import ProjectStore
from dockstudio.core.queue import EventBus, JobRunner
from dockstudio.engines.base import EngineManifest
from dockstudio.models import (
    Box,
    DockJobRequest,
    DockParams,
    DockResult,
    JobStatus,
    Pose,
)


class FakeEngine:
    manifest = EngineManifest(name="fake", scoring_functions=["vina"])

    def __init__(self, delay: float = 0.2):
        self.delay = delay
        self.current_concurrency = 0
        self.max_concurrency = 0
        self.cancel_calls = 0
        self._lock = threading.Lock()

    def is_available(self):
        return True, None

    def cancel_current(self):
        self.cancel_calls += 1

    def dock(self, receptor_pdbqt, ligand_pdbqt, box, params, work_dir, flex_pdbqt=None):
        with self._lock:
            self.current_concurrency += 1
            self.max_concurrency = max(self.max_concurrency, self.current_concurrency)
        try:
            work_dir.mkdir(parents=True, exist_ok=True)
            time.sleep(self.delay)
            (work_dir / "out.pdbqt").write_text(_FAKE_PDBQT, encoding="utf-8")
            return DockResult(
                engine="fake",
                engine_version="0.0",
                params=params,
                box=box,
                poses=[Pose(index=1, affinity=-7.5)],
                run_dir=str(work_dir),
            )
        finally:
            with self._lock:
                self.current_concurrency -= 1


_FAKE_PDBQT = """MODEL 1
REMARK VINA RESULT: -7.5 0.0 0.0
ROOT
ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O
ENDROOT
ENDMDL
"""


def _request() -> DockJobRequest:
    return DockJobRequest(
        engine="fake",
        receptor_id="r1",
        ligand_id="l1",
        box=Box(center=(0, 0, 0), size=(20, 20, 20)),
        params=DockParams(num_modes=1),
    )


class _Registry:
    class _Entry:
        pdbqt_path = str(Path(__file__))

    receptors = {"r1": _Entry()}  # noqa: RUF012
    ligands = {"l1": _Entry()}  # noqa: RUF012


_ORIGINAL_GET_ENGINE = queue_mod.get_engine


def _make_runner(settings, engine) -> JobRunner:
    store = ProjectStore(settings)
    runner = JobRunner(store, _Registry(), EventBus())
    queue_mod.get_engine = lambda name: engine if name == "fake" else _ORIGINAL_GET_ENGINE(name)
    return runner


@pytest.fixture(autouse=True)
def _restore_get_engine():
    yield
    queue_mod.get_engine = _ORIGINAL_GET_ENGINE


async def _wait_terminal(runner: JobRunner, job_id: str, timeout: float = 10.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = runner.get_job(job_id)
        if job is not None and job.is_terminal:
            return job
        await asyncio.sleep(0.05)
    raise TimeoutError(f"job {job_id} did not finish")


def test_serial_execution(settings) -> None:
    engine = FakeEngine(delay=0.3)

    async def main():
        runner = _make_runner(settings, engine)
        runner.start()
        store = runner.store
        project = store.create("p")
        j1 = store.new_job(project.id, _request())
        j2 = store.new_job(project.id, _request())
        runner.submit(j1)
        runner.submit(j2)
        done1, done2 = await asyncio.gather(
            _wait_terminal(runner, j1.id), _wait_terminal(runner, j2.id)
        )
        await runner.stop()
        return done1, done2

    done1, done2 = asyncio.run(main())
    assert engine.max_concurrency == 1  # strictly serial
    assert done1.status == JobStatus.completed
    assert done2.status == JobStatus.completed
    assert done1.result is not None and done1.result.poses[0].affinity == -7.5
    # job.json persisted with final state
    reloaded = ProjectStore(settings).find_job(done1.id)
    assert reloaded is not None and reloaded.status == JobStatus.completed


def test_cancel_queued_job(settings) -> None:
    engine = FakeEngine(delay=0.5)

    async def main():
        runner = _make_runner(settings, engine)
        runner.start()
        store = runner.store
        project = store.create("p")
        j1 = store.new_job(project.id, _request())
        j2 = store.new_job(project.id, _request())
        runner.submit(j1)
        runner.submit(j2)
        await asyncio.sleep(0.1)  # j1 starts running, j2 still queued
        canceled = runner.cancel(j2.id)
        assert canceled is not None and canceled.status == JobStatus.canceled
        done1 = await _wait_terminal(runner, j1.id)
        job2 = runner.get_job(j2.id)
        await runner.stop()
        return done1, job2

    done1, job2 = asyncio.run(main())
    assert done1.status in (JobStatus.completed, JobStatus.canceled)
    assert job2 is not None and job2.status == JobStatus.canceled


def test_cancel_running_job_calls_engine(settings) -> None:
    engine = FakeEngine(delay=2.0)

    async def main():
        runner = _make_runner(settings, engine)
        runner.start()
        store = runner.store
        project = store.create("p")
        j1 = store.new_job(project.id, _request())
        runner.submit(j1)
        await asyncio.sleep(0.1)
        job = runner.get_job(j1.id)
        assert job.status == JobStatus.running
        runner.cancel(j1.id)
        assert engine.cancel_calls >= 1
        final = await _wait_terminal(runner, j1.id)
        await runner.stop()
        return final

    final = asyncio.run(main())
    assert final.status in (JobStatus.canceled, JobStatus.completed, JobStatus.failed)


def test_worker_error_marks_failed(settings, monkeypatch) -> None:
    engine = FakeEngine(delay=0.05)
    monkeypatch.setattr(
        "dockstudio.core.queue.get_engine", lambda name: engine if name == "fake" else None
    )

    class BrokenEngine(FakeEngine):
        def dock(self, *args, **kwargs):
            raise RuntimeError("boom")

    broken = BrokenEngine()

    async def main():
        runner = _make_runner(settings, broken)
        runner.start()
        store = runner.store
        project = store.create("p")
        j1 = store.new_job(project.id, _request())
        runner.submit(j1)
        final = await _wait_terminal(runner, j1.id)
        await runner.stop()
        return final

    final = asyncio.run(main())
    assert final.status == JobStatus.failed
    assert "boom" in (final.error or "")
