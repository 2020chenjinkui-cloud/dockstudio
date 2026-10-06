"""Unit tests for ProjectStore."""

from __future__ import annotations

from dockstudio.core.project import ProjectStore
from dockstudio.models import Box, DockJobRequest, JobStatus


def _request() -> DockJobRequest:
    return DockJobRequest(
        receptor_id="rec1",
        ligand_id="lig1",
        box=Box(center=(0, 0, 0), size=(20, 20, 20)),
    )


def test_create_list_get(settings) -> None:
    store = ProjectStore(settings)
    p = store.create("测试项目")
    assert p.id and p.name == "测试项目"
    projects = store.list()
    assert [x.id for x in projects] == [p.id]
    got = store.get(p.id)
    assert got is not None and got.name == "测试项目"


def test_soft_delete(settings) -> None:
    from pathlib import Path

    store = ProjectStore(settings)
    p = store.create("gone")
    assert store.delete(p.id) is True
    assert store.get(p.id) is None
    assert store.list() == []
    assert (Path(p.path).parent / f".trash-{Path(p.path).name}").is_dir()
    assert store.delete(p.id) is False  # already deleted
    assert store.delete("nonexistent") is False


def test_detail_aggregates_jobs(settings) -> None:
    store = ProjectStore(settings)
    p = store.create("proj")
    job = store.new_job(p.id, _request())
    job.status = JobStatus.completed
    store.save_job(job)
    detail = store.detail(p.id)
    assert detail is not None
    assert len(detail["jobs"]) == 1
    summary = detail["jobs"][0]
    assert summary["id"] == job.id
    assert summary["status"] == "completed"
    assert summary["ligand_id"] == "lig1"
    assert summary["top_affinity"] is None


def test_mark_interrupted_on_restart(settings) -> None:
    store = ProjectStore(settings)
    p = store.create("proj")
    job = store.new_job(p.id, _request())
    job.status = JobStatus.running
    store.save_job(job)
    done = store.new_job(p.id, _request())
    done.status = JobStatus.completed
    store.save_job(done)

    store2 = ProjectStore(settings)
    n = store2.mark_interrupted()
    assert n == 1
    reloaded = store2.find_job(job.id)
    assert reloaded is not None and reloaded.status == JobStatus.interrupted
    assert store2.find_job(done.id).status == JobStatus.completed
