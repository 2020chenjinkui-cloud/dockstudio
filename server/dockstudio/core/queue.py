"""Serial asyncio job queue: run engine.dock in a thread, persist every state."""

from __future__ import annotations

import asyncio
import threading
from collections.abc import Callable
from datetime import datetime
from pathlib import Path

from ..analysis.export import export_poses
from ..engines import get_engine
from ..logging import get_logger
from ..models import (
    DockStudioError,
    EngineUnavailableError,
    Job,
    JobProgress,
    JobStatus,
    PreparationError,
)
from .orchestrator import PrepareRegistry
from .project import ProjectStore

log = get_logger("core.queue")

# Event types broadcast on the bus: "status" | "progress" | "finished"
JobEvent = tuple[str, object]
EventListener = Callable[[str, JobEvent], None]


class EventBus:
    """In-memory pub/sub for job events, consumed by the WebSocket endpoint."""

    def __init__(self) -> None:
        self._listeners: dict[str, list[EventListener]] = {}
        self._lock = threading.Lock()

    def subscribe(self, job_id: str, listener: EventListener) -> Callable[[], None]:
        with self._lock:
            self._listeners.setdefault(job_id, []).append(listener)

        def unsubscribe() -> None:
            with self._lock:
                listeners = self._listeners.get(job_id, [])
                if listener in listeners:
                    listeners.remove(listener)

        return unsubscribe

    def publish(self, job_id: str, event: JobEvent) -> None:
        with self._lock:
            listeners = list(self._listeners.get(job_id, []))
        for listener in listeners:
            try:
                listener(job_id, event)
            except Exception as exc:  # noqa: BLE001
                log.warning("event listener failed for job %s: %s", job_id, exc)


class JobRunner:
    """Runs jobs one at a time. State changes hit disk synchronously."""

    def __init__(
        self,
        store: ProjectStore,
        registry: PrepareRegistry,
        bus: EventBus | None = None,
        event_emitter: Callable[[Job, str, dict], None] | None = None,
    ):
        self.store = store
        self.registry = registry
        self.bus = bus or EventBus()
        self.event_emitter = event_emitter  # (job, action, payload), e.g. EventLog.emit
        self._queue: asyncio.Queue[str] = asyncio.Queue()
        self._jobs: dict[str, Job] = {}
        self._cancel_requested: set[str] = set()
        self._worker_task: asyncio.Task | None = None
        self._engine_lock = threading.Lock()
        self._active_engine: object | None = None
        self._active_job_id: str | None = None

    # ------------------------------------------------------------- enqueue

    def submit(self, job: Job) -> None:
        """Queue a job (called from a request handler on the event loop)."""
        self._jobs[job.id] = job
        self._queue.put_nowait(job.id)
        self._emit_status(job)

    # ------------------------------------------------------------- cancel

    def cancel(self, job_id: str) -> Job | None:
        job = self._jobs.get(job_id) or self.store.find_job(job_id)
        if job is None:
            return None
        if job.status == JobStatus.queued:
            self._cancel_requested.add(job_id)
            job.status = JobStatus.canceled
            job.progress.message = "已取消"
            job.finished_at = datetime.now()
            self.store.save_job(job)
            self._emit_finished(job)
            return job
        if job.status == JobStatus.running:
            self._cancel_requested.add(job_id)
            with self._engine_lock:
                engine = self._active_engine if self._active_job_id == job_id else None
            if engine is None:
                # The worker has not reached engine.dock() yet; the flag is
                # checked before launching and will cancel at that point.
                return job
            try:
                engine.cancel_current()  # type: ignore[attr-defined]
            except Exception as exc:  # noqa: BLE001
                log.warning("engine cancel_current failed: %s", exc)
            return job
        return job

    def _mark_canceled(self, job: Job) -> None:
        job.status = JobStatus.canceled
        job.progress.message = "已取消"
        job.finished_at = datetime.now()
        self.store.save_job(job)
        self._emit_finished(job)
        self._lifecycle(job, "job.canceled", {})

    # ------------------------------------------------------------- worker

    def start(self) -> None:
        if self._worker_task is None or self._worker_task.done():
            self._worker_task = asyncio.create_task(self._worker())

    async def stop(self) -> None:
        if self._worker_task is not None:
            self._worker_task.cancel()
            try:
                await self._worker_task
            except asyncio.CancelledError:
                pass
            self._worker_task = None

    async def _worker(self) -> None:
        while True:
            job_id = await self._queue.get()
            job = self._jobs.get(job_id) or self.store.find_job(job_id)
            try:
                if job is None:
                    log.warning("queued unknown job %s", job_id)
                    continue
                if job_id in self._cancel_requested or job.status == JobStatus.canceled:
                    self._emit_finished(job)
                    continue
                await self._run(job)
            except Exception as exc:  # noqa: BLE001
                log.error("job %s crashed the worker path: %s", job_id, exc)
                if job is not None and not job.is_terminal:
                    job.status = JobStatus.failed
                    job.error = f"内部错误: {exc}"
                    job.finished_at = datetime.now()
                    self.store.save_job(job)
                    self._emit_finished(job)
                    self._lifecycle(job, "job.failed", {"error": job.error})
            finally:
                self._cancel_requested.discard(job_id)
                self._queue.task_done()

    def _lifecycle(self, job: Job, action: str, payload: dict) -> None:
        if self.event_emitter is not None:
            try:
                self.event_emitter(job, action, payload)
            except Exception as exc:  # noqa: BLE001
                log.warning("event emitter failed for %s on job %s: %s", action, job.id, exc)

    async def _run(self, job: Job) -> None:
        loop = asyncio.get_running_loop()
        job.status = JobStatus.running
        job.started_at = datetime.now()
        self._set_progress(job, "docking", 5.0, "对接中")
        self._emit_status(job)
        self._lifecycle(job, "job.started", {"engine": job.request.engine})

        engine = get_engine(job.request.engine)
        if engine is None:
            raise PreparationError(f"unknown engine: {job.request.engine}")
        ok, reason = engine.is_available()
        if not ok:
            raise EngineUnavailableError(reason or "engine unavailable")
        receptor = self.registry.receptors.get(job.request.receptor_id)
        if receptor is None:
            raise PreparationError(f"unknown receptor_id: {job.request.receptor_id}")
        ligand = self.registry.ligands.get(job.request.ligand_id)
        if ligand is None:
            raise PreparationError(f"unknown ligand_id: {job.request.ligand_id}")
        if job.request.params.scoring not in engine.manifest.scoring_functions:
            raise PreparationError(
                f"scoring '{job.request.params.scoring}' not supported by {job.request.engine}"
            )

        run_dir = self.store.job_dir(job.project_id, job.id)
        # getattr: legacy registries/mocks may lack the M2 flex fields
        flex_pdbqt = getattr(receptor, "flex_pdbqt_path", None)
        flex_path: Path | None = None
        if flex_pdbqt:
            if not engine.manifest.supports_flex_residues:
                raise PreparationError(
                    f"engine {job.request.engine} 不支持柔性残基，"
                    "请使用支持 flex 的引擎（如 vina）或改用刚性受体。"
                )
            flex_path = Path(flex_pdbqt)

        with self._engine_lock:
            self._active_engine = engine
            self._active_job_id = job.id
        try:
            if job.id in self._cancel_requested:
                self._mark_canceled(job)
                return
            try:
                result = await loop.run_in_executor(
                    None,
                    lambda: engine.dock(
                        Path(receptor.pdbqt_path),
                        Path(ligand.pdbqt_path),
                        job.request.box,
                        job.request.params,
                        run_dir,
                        flex_path,
                    ),
                )
            except DockStudioError:
                if job.id in self._cancel_requested:
                    self._mark_canceled(job)
                    return
                raise
            except Exception as exc:
                if job.id in self._cancel_requested:
                    self._mark_canceled(job)
                    return
                raise PreparationError(f"dock crashed: {exc}") from exc
        finally:
            with self._engine_lock:
                if self._active_job_id == job.id:
                    self._active_engine = None
                    self._active_job_id = None

        if job.id in self._cancel_requested:
            self._mark_canceled(job)
            return

        result.run_dir = str(run_dir)
        self._set_progress(job, "exporting", 88.0, "导出姿态")
        self.store.save_job(job)
        pose_paths = await loop.run_in_executor(
            None, lambda: export_poses(run_dir / "out.pdbqt", run_dir)
        )
        for pose, pose_path in zip(result.poses, pose_paths, strict=False):
            pose.file = str(pose_path)

        if job.id in self._cancel_requested:
            self._mark_canceled(job)
            return

        self._set_progress(job, "analyzing", 95.0, "分析相互作用")
        self.store.save_job(job)

        job.result = result
        job.status = JobStatus.completed
        job.progress = JobProgress(stage="done", percent=100.0, message="完成")
        job.finished_at = datetime.now()
        self.store.save_job(job)
        self._emit_finished(job)
        self._lifecycle(
            job,
            "job.completed",
            {
                "top_affinity": (result.poses[0].affinity if result.poses else None),
                "n_poses": len(result.poses),
            },
        )
        log.info(
            "job %s completed: %d poses, top=%.2f",
            job.id,
            len(result.poses),
            result.poses[0].affinity if result.poses else float("nan"),
        )

    # ------------------------------------------------------------- helpers

    def _set_progress(self, job: Job, stage: str, percent: float | None, message: str) -> None:
        job.progress = JobProgress(stage=stage, percent=percent, message=message)
        self.store.save_job(job)
        self._emit("progress", job, job.progress)

    def _emit_status(self, job: Job) -> None:
        self._emit("status", job, job)

    def _emit_finished(self, job: Job) -> None:
        self._emit("finished", job, job)

    def _emit(self, kind: str, job: Job, payload: object) -> None:
        self._jobs[job.id] = job
        self.bus.publish(job.id, (kind, payload))

    def get_job(self, job_id: str) -> Job | None:
        return self._jobs.get(job_id) or self.store.find_job(job_id)
