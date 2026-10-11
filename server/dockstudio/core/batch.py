"""Batch screening: submit a set of ligands as tagged jobs through the queue."""

from __future__ import annotations

import csv
import io
import json
import math
import uuid
from datetime import datetime
from pathlib import Path

from ..config import Settings
from ..logging import get_logger
from ..models import (
    BatchCreate,
    BatchRow,
    BatchSummary,
    DockJobRequest,
    Job,
    JobStatus,
)
from .project import ProjectStore

log = get_logger("core.batch")


class BatchStore:
    """Filesystem batch registry (batches.json) + summary/results aggregation."""

    def __init__(self, settings: Settings, store: ProjectStore):
        self.settings = settings
        self.store = store
        self.root = settings.data_dir / "projects"

    # ------------------------------------------------------------ persistence

    @property
    def registry_path(self) -> Path:
        return self.root / "batches.json"

    def _load(self) -> dict[str, dict]:
        if not self.registry_path.is_file():
            return {}
        try:
            return json.loads(self.registry_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            log.warning("could not read batches registry: %s", exc)
            return {}

    def _save(self, data: dict[str, dict]) -> None:
        tmp = self.registry_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
        tmp.replace(self.registry_path)

    def create(self, project_id: str, req: BatchCreate) -> BatchSummary:
        batch_id = uuid.uuid4().hex[:8]
        summary = BatchSummary(
            id=batch_id,
            project_id=project_id,
            name=req.name or f"batch-{batch_id}",
            total=len(req.ligand_ids),
            queued=len(req.ligand_ids),
            receptor_id=req.receptor_id,
            engine=req.engine,
            created_at=datetime.now(),
        )
        data = self._load()
        data[batch_id] = {
            "project_id": project_id,
            "name": summary.name,
            "receptor_id": req.receptor_id,
            "engine": req.engine,
            "created_at": summary.created_at.isoformat(),
        }
        self._save(data)
        return summary

    def get_meta(self, batch_id: str) -> dict | None:
        return self._load().get(batch_id)

    # ------------------------------------------------------------ aggregation

    def jobs_of(self, batch_id: str) -> list[Job]:
        meta = self.get_meta(batch_id)
        if meta is None:
            return []
        return [j for j in self.store.list_jobs(meta["project_id"]) if j.batch_id == batch_id]

    def summary(
        self, batch_id: str, runner_jobs: dict[str, Job] | None = None
    ) -> BatchSummary | None:
        meta = self.get_meta(batch_id)
        if meta is None:
            return None
        jobs = self.jobs_of(batch_id)
        if runner_jobs:
            jobs = [runner_jobs.get(j.id, j) for j in jobs]
        counts = {s: 0 for s in JobStatus}
        for j in jobs:
            counts[j.status] += 1
        terminal = all(j.is_terminal for j in jobs) if jobs else False
        return BatchSummary(
            id=batch_id,
            project_id=meta["project_id"],
            name=meta["name"],
            total=len(jobs),
            done=counts[JobStatus.completed],
            failed=counts[JobStatus.failed] + counts[JobStatus.interrupted],
            running=counts[JobStatus.running],
            queued=counts[JobStatus.queued],
            canceled=counts[JobStatus.canceled],
            receptor_id=meta["receptor_id"],
            engine=meta.get("engine", "vina"),
            created_at=datetime.fromisoformat(meta["created_at"]),
            is_terminal=terminal,
        )

    # ------------------------------------------------------------ results

    def results(
        self,
        batch_id: str,
        registry,
        sort: str = "affinity",
        order: str = "asc",
        limit: int = 50,
        offset: int = 0,
    ) -> dict | None:
        if self.get_meta(batch_id) is None:
            return None
        rows: list[BatchRow] = []
        for job in self.jobs_of(batch_id):
            ligand = registry.ligands.get(job.request.ligand_id)
            rows.append(
                BatchRow(
                    ligand_id=job.request.ligand_id,
                    ligand_name=ligand.name if ligand else None,
                    smiles=ligand.smiles if ligand else None,
                    job_id=job.id,
                    status=job.status,
                    top_affinity=(
                        job.result.poses[0].affinity
                        if job.result is not None and job.result.poses
                        else None
                    ),
                    error=job.error,
                )
            )
        reverse = order == "desc"
        if sort == "affinity":
            rows.sort(
                key=lambda r: (
                    r.top_affinity is None,
                    math.inf if r.top_affinity is None else r.top_affinity,
                ),
                reverse=reverse,
            )
        else:
            rows.sort(key=lambda r: r.ligand_name or r.ligand_id, reverse=reverse)
        total = len(rows)
        return {"total": total, "rows": rows[offset : offset + limit]}

    def results_csv(self, batch_id: str, registry) -> str | None:
        data = self.results(batch_id, registry, sort="affinity", order="asc", limit=10**9)
        if data is None:
            return None
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(
            ["rank", "ligand_name", "smiles", "top_affinity", "status", "job_id", "error"]
        )
        for rank, row in enumerate(data["rows"], start=1):
            writer.writerow(
                [
                    rank,
                    row.ligand_name or "",
                    row.smiles or "",
                    "" if row.top_affinity is None else f"{row.top_affinity:.3f}",
                    row.status.value,
                    row.job_id,
                    row.error or "",
                ]
            )
        return buf.getvalue()


def job_request_from_batch(req: BatchCreate, ligand_id: str) -> DockJobRequest:
    return DockJobRequest(
        engine=req.engine,
        receptor_id=req.receptor_id,
        ligand_id=ligand_id,
        box=req.box,
        params=req.params,
    )
