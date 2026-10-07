"""Project store: project directories under data_dir/projects/<id>/jobs/<jid>/."""

from __future__ import annotations

import json
import re
import uuid
from datetime import datetime
from pathlib import Path

from ..config import Settings
from ..logging import get_logger
from ..models import DockJobRequest, DockResult, Job, JobStatus, Project

log = get_logger("core.project")


class JobSummary(dict):
    """Lightweight job row for the project detail view."""


def _job_summary(job: Job) -> dict:
    return {
        "id": job.id,
        "status": job.status.value,
        "ligand_id": job.request.ligand_id,
        "engine": job.request.engine,
        "top_affinity": (
            job.result.poses[0].affinity if job.result is not None and job.result.poses else None
        ),
        "created_at": job.created_at,
        "finished_at": job.finished_at,
    }


class ProjectStore:
    """Filesystem-backed project registry (projects.json + per-job job.json)."""

    def __init__(self, settings: Settings):
        self.settings = settings
        self.root = settings.data_dir / "projects"
        self.root.mkdir(parents=True, exist_ok=True)
        self._trash_marker = ".trash-"

    @property
    def registry_path(self) -> Path:
        return self.root / "projects.json"

    def _load_registry(self) -> dict[str, dict]:
        if not self.registry_path.is_file():
            return {}
        try:
            return json.loads(self.registry_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            log.warning("could not read projects registry: %s", exc)
            return {}

    def _save_registry(self, data: dict[str, dict]) -> None:
        tmp = self.registry_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
        tmp.replace(self.registry_path)

    # ------------------------------------------------------------------ CRUD

    def create(self, name: str) -> Project:
        name = (name or "").strip() or "未命名项目"
        project_id = uuid.uuid4().hex[:8]
        project_dir = self.root / project_id
        (project_dir / "jobs").mkdir(parents=True, exist_ok=True)
        project = Project(
            id=project_id,
            name=name,
            path=str(project_dir),
            created_at=datetime.now(),  # noqa: DTZ005
        )
        data = self._load_registry()
        data[project_id] = project.model_dump(mode="json")
        self._save_registry(data)
        return project

    def list(self) -> list[Project]:
        projects: list[Project] = []
        for payload in self._load_registry().values():
            if payload.get("deleted"):
                continue
            projects.append(
                Project(
                    id=payload["id"],
                    name=payload["name"],
                    path=payload["path"],
                    created_at=payload["created_at"],
                )
            )
        projects.sort(key=lambda p: p.created_at, reverse=True)
        return projects

    def get(self, project_id: str) -> Project | None:
        payload = self._load_registry().get(project_id)
        if not payload or payload.get("deleted"):
            return None
        return Project(
            id=payload["id"],
            name=payload["name"],
            path=payload["path"],
            created_at=payload["created_at"],
        )

    def detail(self, project_id: str) -> dict | None:
        project = self.get(project_id)
        if project is None:
            return None
        jobs = [JobSummary(_job_summary(j)) for j in self.list_jobs(project_id)]
        jobs.sort(key=lambda j: j["created_at"], reverse=True)
        return {**project.model_dump(), "jobs": jobs}

    def delete(self, project_id: str) -> bool:
        data = self._load_registry()
        payload = data.get(project_id)
        if not payload or payload.get("deleted"):
            return False
        payload["deleted"] = True
        self._save_registry(data)
        project_dir = Path(payload["path"])
        if project_dir.is_dir():
            trash_dir = project_dir.parent / f"{self._trash_marker}{project_dir.name}"
            try:
                project_dir.replace(trash_dir)
            except OSError as exc:
                log.warning("could not move project dir to trash: %s", exc)
        return True

    # -------------------------------------------------------------------- jobs

    def job_dir(self, project_id: str, job_id: str) -> Path:
        return self.root / project_id / "jobs" / job_id

    def new_job(self, project_id: str, request: DockJobRequest, batch_id: str | None = None) -> Job:
        job_id = uuid.uuid4().hex[:8]
        job = Job(
            id=job_id,
            project_id=project_id,
            batch_id=batch_id,
            request=request,
            created_at=datetime.now(),  # noqa: DTZ005
        )
        job_dir = self.job_dir(project_id, job_id)
        job_dir.mkdir(parents=True, exist_ok=True)
        self.save_job(job)
        return job

    def save_job(self, job: Job) -> None:
        job_dir = self.job_dir(job.project_id, job.id)
        job_dir.mkdir(parents=True, exist_ok=True)
        path = job_dir / "job.json"
        tmp = path.with_suffix(".tmp")
        tmp.write_text(job.model_dump_json(indent=2), encoding="utf-8")
        tmp.replace(path)

    def load_job(self, path: Path) -> Job | None:
        try:
            return Job(**json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            log.warning("could not load job file %s: %s", path, exc)
            return None

    def find_job(self, job_id: str) -> Job | None:
        for path in self.root.glob("*/jobs/*/job.json"):
            if path.parent.name == job_id:
                job = self.load_job(path)
                if job is not None:
                    return job
        return None

    def list_jobs(self, project_id: str) -> list[Job]:
        jobs: list[Job] = []
        jobs_dir = self.root / project_id / "jobs"
        if jobs_dir.is_dir():
            for path in jobs_dir.glob("*/job.json"):
                job = self.load_job(path)
                if job is not None:
                    jobs.append(job)
        return jobs

    def set_job_result(self, job: Job, result: DockResult) -> None:
        job.result = result
        self.save_job(job)

    def mark_interrupted(self) -> int:
        """Any non-terminal job on disk belongs to a dead process -> interrupted."""
        count = 0
        for path in self.root.glob("*/jobs/*/job.json"):
            job = self.load_job(path)
            if job is not None and not job.is_terminal:
                job.status = JobStatus.interrupted
                job.progress.message = "服务重启打断"
                if job.finished_at is None:
                    job.finished_at = datetime.now()  # noqa: DTZ005
                self.save_job(job)
                count += 1
                log.info("job %s marked interrupted after restart", job.id)
        return count


_JOB_ID_RE = re.compile(r"^[0-9a-f]{8}$")


def valid_job_id(job_id: str) -> bool:
    return bool(_JOB_ID_RE.match(job_id))
