"""FastAPI application for DockStudio (M0 sync chain + M1 projects/jobs/WS)."""

from __future__ import annotations

import asyncio
import json
import os
import sys
import traceback
import uuid
from contextlib import asynccontextmanager
from datetime import datetime
from pathlib import Path
from typing import Annotated

import uvicorn
from fastapi import (
    FastAPI,
    File,
    Form,
    Query,
    Request,
    Response,
    UploadFile,
    WebSocket,
    WebSocketDisconnect,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from pydantic import BaseModel

from .. import __version__
from ..analysis.export import export_poses
from ..analysis.interactions import analyze_job_poses
from ..analysis.report import generate_report
from ..config import Settings, get_settings
from ..core.batch import BatchStore
from ..core.events import EventLog, env_fingerprint, sha256_file
from ..core.orchestrator import PrepareRegistry, new_run_dir
from ..core.project import ProjectStore
from ..core.queue import EventBus, JobRunner
from ..engines import available_engines, get_engine
from ..logging import get_logger, setup_logging
from ..models import (
    AiRescoreRequest,
    AiRescoreTask,
    BatchCreate,
    DockJobRequest,
    DockResult,
    DockStudioError,
    EngineUnavailableError,
    Job,
    JobStatus,
    PreparationError,
    PreparedLigand,
    PreparedReceptor,
    Project,
)
from ..prep.ligand import prepare_ligand_from_file, prepare_ligand_from_smiles
from ..prep.receptor import prepare_receptor

log = get_logger("api.app")

# 放行来源：vite 开发服务器 + Tauri 生产 webview
# （Windows/Linux 为 http://tauri.localhost，macOS 为 tauri://localhost）
_ALLOWED_ORIGINS = [
    "http://localhost:5173",
    "http://127.0.0.1:5173",
    "http://tauri.localhost",
    "tauri://localhost",
]


class CreateProjectRequest(BaseModel):
    name: str


class PocketsRequest(BaseModel):
    receptor_id: str
    top_n: int = 5


class FlexRequest(BaseModel):
    flex_residues: list[str]


MAX_UPLOAD_BYTES = 50 * 1024 * 1024  # 50 MB


async def _read_upload(file: UploadFile, limit: int = MAX_UPLOAD_BYTES) -> bytes:
    """Read an uploaded file enforcing the size cap (413 with readable detail)."""
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = await file.read(1 << 20)
        if not chunk:
            break
        total += len(chunk)
        if total > limit:
            raise PayloadTooLargeError(
                f"文件超过大小上限 {limit // (1024 * 1024)} MB，请压缩或裁切后重试。"
            )
        chunks.append(chunk)
    return b"".join(chunks)


class PayloadTooLargeError(Exception):
    """Uploaded payload exceeded the size cap (maps to 413)."""


def _error(status: int, detail: str) -> JSONResponse:
    return JSONResponse(status_code=status, content={"detail": detail})


def _not_found(what: str) -> JSONResponse:
    return _error(404, f"{what} not found")


def _job_protein_pdb(registry: PrepareRegistry, job: Job) -> Path:
    receptor = registry.receptors.get(job.request.receptor_id)
    if receptor is None:
        raise PreparationError(f"unknown receptor_id: {job.request.receptor_id}")
    rec_dir = Path(receptor.pdbqt_path).parent
    protein = rec_dir / "receptor_protein.pdb"
    return protein if protein.is_file() else rec_dir / "receptor_clean.pdb"


def _job_pose_paths(job: Job) -> list[Path]:
    if job.result is None:
        return []
    return [Path(p.file) for p in job.result.poses if p.file]


def create_app(settings: Settings | None = None) -> FastAPI:
    settings = settings or get_settings()
    setup_logging(settings.data_dir / "logs")
    registry = PrepareRegistry(settings)
    store = ProjectStore(settings)
    interrupted = store.mark_interrupted()
    if interrupted:
        log.info("marked %d stale jobs interrupted", interrupted)
    bus = EventBus()

    def _runner_emitter(job: Job, action: str, payload: dict) -> None:
        EventLog(store.root / job.project_id).emit(action, "job", job.id, payload, actor="system")

    runner = JobRunner(store, registry, bus, event_emitter=_runner_emitter)
    batches = BatchStore(settings, store)
    from ..analysis.ai_rescore import RescoreRunner

    rescore = RescoreRunner(store, settings, registry)

    @asynccontextmanager
    async def lifespan(_app: FastAPI):
        runner.start()
        rescore.start()
        if os.environ.get("DOCKSTUDIO_ENABLE_WATCHDOG") == "1":
            from ..core.lifecycle import install_parent_watchdog

            # Start after the heavy app imports; creating this thread during
            # uvicorn's factory import can deadlock the frozen onefile child.
            install_parent_watchdog()
        try:
            yield
        finally:
            rescore.stop()
            await runner.stop()

    app = FastAPI(title="DockStudio", version=__version__, lifespan=lifespan)
    app.state.settings = settings
    app.state.registry = registry
    app.state.store = store
    app.state.runner = runner
    app.state.rescore = rescore

    def events_for(project_id: str) -> EventLog:
        return EventLog(store.root / project_id)

    def _input_hashes(job: Job) -> dict:
        receptor = registry.receptors.get(job.request.receptor_id)
        ligand = registry.ligands.get(job.request.ligand_id)
        out: dict[str, str] = {}
        for name, entry in (("receptor", receptor), ("ligand", ligand)):
            if entry is None:
                continue
            path = Path(getattr(entry, "pdbqt_path", ""))
            if not path.is_file():
                path = Path(getattr(entry, "original_path", "") or "")
            if path.is_file():
                out[name] = sha256_file(path)
        return out

    def _save_input_hashes(job: Job) -> dict:
        hashes = _input_hashes(job)
        job_dir = store.job_dir(job.project_id, job.id)
        job_dir.mkdir(parents=True, exist_ok=True)
        (job_dir / "inputs.sha256.json").write_text(json.dumps(hashes, indent=2), encoding="utf-8")
        return hashes

    def _load_input_hashes(job: Job) -> dict:
        sidecar = store.job_dir(job.project_id, job.id) / "inputs.sha256.json"
        if sidecar.is_file():
            try:
                return json.loads(sidecar.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                return {}
        return {}

    app.add_middleware(
        CORSMiddleware,
        allow_origins=_ALLOWED_ORIGINS,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.exception_handler(DockStudioError)
    async def dockstudio_error_handler(_request: Request, exc: DockStudioError) -> JSONResponse:
        status = 503 if isinstance(exc, EngineUnavailableError) else 400
        return _error(status, str(exc))

    @app.exception_handler(PayloadTooLargeError)
    async def payload_too_large_handler(
        _request: Request, exc: PayloadTooLargeError
    ) -> JSONResponse:
        return _error(413, str(exc))

    @app.exception_handler(Exception)
    async def unknown_error_handler(_request: Request, exc: Exception) -> JSONResponse:
        log.error("unhandled error: %s\n%s", exc, traceback.format_exc())
        return _error(500, f"internal error: {type(exc).__name__}: {exc}")

    # ------------------------------------------------------------- M0 (unchanged)

    @app.get("/api/health")
    async def health() -> dict[str, str]:
        return {"status": "ok", "version": __version__}

    @app.get("/api/engines")
    async def engines() -> list[dict]:
        out = []
        for engine in available_engines():
            ok, reason = engine.is_available()
            out.append(
                {
                    **engine.manifest.model_dump(),
                    "available": ok,
                    "unavailable_reason": reason,
                }
            )
        return out

    @app.post("/api/prepare/receptor")
    async def prepare_receptor_endpoint(
        file: Annotated[UploadFile, File()],
    ) -> PreparedReceptor:
        content = await _read_upload(file)
        if not content.strip():
            raise PreparationError("上传文件为空，请选择有效的 PDB 文件。")
        receptor = await asyncio.get_running_loop().run_in_executor(
            None, prepare_receptor, content, file.filename or "receptor.pdb", registry.receptors_dir
        )
        registry.add_receptor(receptor)
        events_for("_prepared").emit(
            "receptor.prepared",
            "receptor",
            receptor.id,
            {
                "filename": file.filename,
                "original_sha256": sha256_file(Path(receptor.original_path)),
                "n_chains": len(receptor.chains),
                "chains": receptor.chains,
                "n_cocrystals": len(receptor.cocrystal_ligands),
                "cocrystals": [c.model_dump(mode="json") for c in receptor.cocrystal_ligands],
                "warnings": receptor.warnings,
            },
        )
        return receptor

    @app.post("/api/prepare/ligand")
    async def prepare_ligand_endpoint(
        file: Annotated[UploadFile | None, File()] = None,
        smiles: Annotated[str | None, Form()] = None,
        name: Annotated[str | None, Form()] = None,
    ) -> PreparedLigand:
        has_file = file is not None and bool(file.filename)
        has_smiles = smiles is not None and smiles.strip() != ""
        if has_file == has_smiles:
            raise PreparationError("provide exactly one of 'file' or 'smiles'")
        loop = asyncio.get_running_loop()
        if has_file:
            assert file is not None and file.filename is not None
            content = await _read_upload(file)
            if not content.strip():
                raise PreparationError("上传文件为空，请选择有效的配体文件。")
            safe_name = Path(file.filename).name or "ligand.sdf"
            tmp_path = registry.ligands_dir / f"upload_{safe_name}"
            tmp_path.write_bytes(content)
            try:
                ligand = await loop.run_in_executor(
                    None, prepare_ligand_from_file, tmp_path, registry.ligands_dir
                )
            finally:
                tmp_path.unlink(missing_ok=True)
        else:
            ligand = await loop.run_in_executor(
                None, prepare_ligand_from_smiles, smiles or "", name, registry.ligands_dir
            )
        registry.add_ligand(ligand)
        events_for("_prepared").emit(
            "ligand.prepared",
            "ligand",
            ligand.id,
            {
                "source": ligand.source,
                "smiles": ligand.smiles,
                "n_rotatable": ligand.n_rotatable_bonds,
                "name": ligand.name,
            },
        )
        return ligand

    @app.post("/api/dock")
    def dock_endpoint(req: DockJobRequest) -> DockResult:
        engine = get_engine(req.engine)
        if engine is None:
            raise PreparationError(f"unknown engine: {req.engine}")
        ok, reason = engine.is_available()
        if not ok:
            raise EngineUnavailableError(reason or f"engine {req.engine} is unavailable")
        if req.params.scoring not in engine.manifest.scoring_functions:
            raise PreparationError(
                f"scoring '{req.params.scoring}' not supported by {req.engine} "
                f"(supported: {engine.manifest.scoring_functions})"
            )
        receptor = registry.receptors.get(req.receptor_id)
        if receptor is None:
            raise PreparationError(f"unknown receptor_id: {req.receptor_id}")
        ligand = registry.ligands.get(req.ligand_id)
        if ligand is None:
            raise PreparationError(f"unknown ligand_id: {req.ligand_id}")
        flex_path = Path(receptor.flex_pdbqt_path) if receptor.flex_pdbqt_path else None
        if flex_path is not None and not engine.manifest.supports_flex_residues:
            raise PreparationError(
                f"engine {req.engine} 不支持柔性残基受体变体，请选择支持 flex 的引擎。"
            )

        run_dir = new_run_dir(settings)
        result = engine.dock(
            Path(receptor.pdbqt_path),
            Path(ligand.pdbqt_path),
            req.box,
            req.params,
            run_dir,
            flex_path,
        )
        result.run_dir = str(run_dir)
        pose_paths = export_poses(run_dir / "out.pdbqt", run_dir)
        for pose, pose_path in zip(result.poses, pose_paths, strict=False):
            pose.file = str(pose_path)
        log.info(
            "dock run %s: %d poses, top=%.2f",
            run_dir.name,
            len(result.poses),
            result.poses[0].affinity,
        )
        return result

    @app.get("/api/file", response_class=PlainTextResponse)
    async def read_file(path: str) -> str:
        data_root = settings.data_dir.resolve()
        try:
            target = Path(path).resolve()
        except OSError:
            raise PreparationError(f"invalid path: {path!r}") from None
        if not target.is_relative_to(data_root):
            return _error(403, "path is outside the data directory")
        if not target.is_file():
            raise PreparationError(f"file not found: {target}")
        return PlainTextResponse(target.read_text(encoding="utf-8", errors="replace"))

    # ------------------------------------------------------------- M1: projects

    @app.post("/api/projects", status_code=201)
    async def create_project(body: CreateProjectRequest) -> dict:
        project = store.create(body.name)
        events_for(project.id).emit(
            "project.create", "project", project.id, {"name": project.name}, actor="user"
        )
        return project.model_dump(mode="json")

    @app.get("/api/projects")
    async def list_projects() -> list[dict]:
        return [p.model_dump(mode="json") for p in store.list()]

    @app.get("/api/projects/{pid}")
    async def project_detail(pid: str):
        detail = store.detail(pid)
        if detail is None:
            return _not_found("project")
        return detail

    @app.delete("/api/projects/{pid}", status_code=204)
    async def delete_project(pid: str) -> Response:
        project = store.get(pid)
        if not store.delete(pid):
            return _not_found("project")
        # project dir moved to .trash-*; record the event in the trash copy too
        trash = store.root / f".trash-{pid}"
        EventLog(trash if trash.is_dir() else store.root / pid).emit(
            "project.delete",
            "project",
            pid,
            {"name": project.name if project else None},
        )
        return Response(status_code=204)

    @app.get("/api/prepare/receptors")
    async def list_receptors() -> list[PreparedReceptor]:
        return list(registry.receptors.values())

    @app.get("/api/prepare/ligands")
    async def list_ligands() -> list[PreparedLigand]:
        return list(registry.ligands.values())

    # ------------------------------------------------------------- M1: jobs

    @app.post("/api/projects/{pid}/jobs", status_code=201)
    async def submit_job(pid: str, body: DockJobRequest):
        if store.get(pid) is None:
            return _not_found("project")
        if body.receptor_id not in registry.receptors:
            raise PreparationError(f"unknown receptor_id: {body.receptor_id}")
        if body.ligand_id not in registry.ligands:
            raise PreparationError(f"unknown ligand_id: {body.ligand_id}")
        engine = get_engine(body.engine)
        if engine is not None:
            if body.params.scoring not in engine.manifest.scoring_functions:
                raise PreparationError(
                    f"scoring '{body.params.scoring}' not supported by {body.engine} "
                    f"(supported: {engine.manifest.scoring_functions})"
                )
            ok, reason = engine.is_available()
            if not ok:
                raise EngineUnavailableError(reason or f"engine {body.engine} is unavailable")
        job = store.new_job(pid, body)
        job.env = env_fingerprint()
        hashes = _save_input_hashes(job)
        store.save_job(job)
        runner.submit(job)
        events_for(pid).emit(
            "job.created",
            "job",
            job.id,
            {
                "engine": body.engine,
                "receptor_id": body.receptor_id,
                "ligand_id": body.ligand_id,
                "box": body.box.model_dump(),
                "params": body.params.model_dump(),
                "batch_id": job.batch_id,
                "env": job.env.model_dump(),
                "input_sha256": hashes,
            },
        )
        return {"job_id": job.id}

    @app.get("/api/jobs/{jid}")
    async def get_job(jid: str) -> JSONResponse:
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        return JSONResponse(job.model_dump(mode="json"), headers={"X-Poll-Fallback": "true"})

    @app.post("/api/jobs/{jid}/cancel")
    async def cancel_job(jid: str):
        job = runner.cancel(jid)
        if job is None:
            return _not_found("job")
        return job

    # ------------------------------------------------------------- M1: analysis

    @app.get("/api/jobs/{jid}/analysis")
    def job_analysis(jid: str, pose_index: int = Query(default=1, ge=1)):
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        if job.status != JobStatus.completed or job.result is None:
            return _error(400, f"job is not completed (status={job.status.value})")
        pose_paths = _job_pose_paths(job)
        if not (1 <= pose_index <= len(pose_paths)):
            return _error(400, f"pose_index out of range (1..{len(pose_paths)})")
        protein_pdb = _job_protein_pdb(registry, job)
        job_dir = store.job_dir(job.project_id, job.id)
        analyses = analyze_job_poses(
            protein_pdb,
            [pose_paths[pose_index - 1]],
            job_dir,
            pose_indices=[pose_index],
            box=job.result.box,
        )
        return analyses[0]

    def _diagram_for_job(job: Job, pose_index: int) -> str:
        """SVG of the 2D interaction diagram, cached under job dir."""
        pose_paths = _job_pose_paths(job)
        protein_pdb = _job_protein_pdb(registry, job)
        job_dir = store.job_dir(job.project_id, job.id)
        cache_path = job_dir / "analysis" / f"diagram_pose_{pose_index}.svg"
        if cache_path.is_file():
            return cache_path.read_text(encoding="utf-8")
        analysis = analyze_job_poses(
            protein_pdb,
            [pose_paths[pose_index - 1]],
            job_dir,
            pose_indices=[pose_index],
            box=job.result.box if job.result else None,
        )[0]
        from ..analysis.diagram import generate_diagram_svg

        ligand = registry.ligands.get(job.request.ligand_id)
        svg = generate_diagram_svg(
            pose_paths[pose_index - 1],
            analysis.interactions,
            title=(ligand.name or ligand.id) if ligand else job.request.ligand_id,
        )
        cache_path.write_text(svg, encoding="utf-8")
        return svg

    # ------------------------------------------------------------- M3: diagram

    @app.get("/api/jobs/{jid}/diagram")
    def job_diagram(jid: str, pose_index: int = Query(default=1, ge=1)):
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        if job.status != JobStatus.completed or job.result is None:
            return _error(400, f"job is not completed (status={job.status.value})")
        if not (1 <= pose_index <= len(_job_pose_paths(job))):
            return _error(400, f"pose_index out of range (1..{len(_job_pose_paths(job))})")
        return Response(_diagram_for_job(job, pose_index), media_type="image/svg+xml")

    # ------------------------------------------------------------- M1: report

    @app.get("/api/jobs/{jid}/report", response_class=HTMLResponse)
    def job_report(
        jid: str,
        sections: str | None = Query(default=None),
        title: str | None = Query(default=None, max_length=80),
    ) -> Response:
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        if job.status != JobStatus.completed or job.result is None:
            return _error(400, f"job is not completed (status={job.status.value})")
        from ..analysis.report import ALL_SECTIONS

        selected: list[str] | None = None
        if sections:
            selected = [s.strip() for s in sections.split(",") if s.strip()]
            invalid = [s for s in selected if s not in ALL_SECTIONS]
            if invalid:
                return _error(400, f"非法 section: {invalid}（可选: {', '.join(ALL_SECTIONS)}）")
        project = store.get(job.project_id)
        pose_paths = _job_pose_paths(job)
        job_dir = store.job_dir(job.project_id, job.id)
        protein_pdb = _job_protein_pdb(registry, job)
        box = job.result.box
        analyses = analyze_job_poses(
            protein_pdb, pose_paths[:1], job_dir, pose_indices=[1], box=box
        )
        diagram_svg = ""
        if selected is None or "diagram" in selected:
            try:
                diagram_svg = _diagram_for_job(job, 1)
            except Exception as exc:  # noqa: BLE001
                log.warning("diagram generation failed (non-fatal): %s", exc)
        ligand = registry.ligands.get(job.request.ligand_id)
        if project is None:  # project was soft-deleted; build a transient shell
            project = Project(
                id=job.project_id,
                name="（已删除项目）",
                path=str(store.job_dir(job.project_id, job.id).parent.parent),
                created_at=datetime.now(),
            )
        report_html = generate_report(
            project=project,
            request=job.request,
            result=job.result,
            pose_analyses=analyses,
            receptor_pdb_path=protein_pdb,
            pose_sdf_paths=pose_paths,
            receptor_name=job.request.receptor_id,
            ligand_name=(ligand.name or ligand.id) if ligand else job.request.ligand_id,
            sections=selected,
            title=title,
            diagram_svg=diagram_svg,
        )
        return HTMLResponse(
            report_html,
            headers={"Content-Disposition": f'attachment; filename="dockstudio-report-{jid}.html"'},
        )

    # ------------------------------------------------------------- M3: compare & stats

    @app.get("/api/compare")
    async def compare_jobs(job_ids: str = Query(default="")):
        from ..models import ComparisonResult, ComparisonRow

        rows: list[ComparisonRow] = []
        for jid in [s.strip() for s in job_ids.split(",") if s.strip()]:
            job = runner.get_job(jid)
            if job is None:
                continue
            ligand = registry.ligands.get(job.request.ligand_id)
            poses = job.result.poses if job.result else []
            rows.append(
                ComparisonRow(
                    job_id=job.id,
                    ligand_id=job.request.ligand_id,
                    ligand_name=ligand.name if ligand else None,
                    smiles=ligand.smiles if ligand else None,
                    engine=job.request.engine,
                    scoring=job.request.params.scoring,
                    exhaustiveness=job.request.params.exhaustiveness,
                    top_affinity=poses[0].affinity if poses else None,
                    n_poses=len(poses),
                    status=job.status,
                )
            )
        return ComparisonResult(rows=rows)

    @app.get("/api/stats")
    async def stats() -> dict:
        projects = store.list()
        all_jobs: list[Job] = []
        for p in projects:
            all_jobs.extend(store.list_jobs(p.id))
        batches_path = settings.data_dir / "projects" / "batches.json"
        n_batches = 0
        if batches_path.is_file():
            try:
                n_batches = len(json.loads(batches_path.read_text(encoding="utf-8")))
            except (OSError, json.JSONDecodeError):
                n_batches = 0
        completed = [j for j in all_jobs if j.status == JobStatus.completed]
        durations = [
            (j.finished_at - j.started_at).total_seconds()
            for j in completed
            if j.started_at and j.finished_at
        ]
        by_engine: dict[str, int] = {}
        for j in all_jobs:
            by_engine[j.request.engine] = by_engine.get(j.request.engine, 0) + 1
        return {
            "projects": len(projects),
            "jobs": len(all_jobs),
            "batches": n_batches,
            "completed_jobs": len(completed),
            "avg_dock_seconds": (round(sum(durations) / len(durations), 2) if durations else None),
            "by_engine": by_engine,
        }

    # ------------------------------------------------------------- M4: AI rescore (Boltz-2)

    @app.get("/api/capabilities")
    async def capabilities() -> dict:
        from ..analysis.ai_rescore import get_availability

        return {"ai_rescore": get_availability().model_dump()}

    @app.post("/api/jobs/{jid}/rescore/ai", status_code=202)
    async def submit_rescore(jid: str, body: AiRescoreRequest | None = None):
        from ..analysis.ai_rescore import get_availability

        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        if job.status != JobStatus.completed or job.result is None:
            return _error(400, f"job is not completed (status={job.status.value})")
        available = get_availability()
        if not available.available:
            return _error(503, available.reason or "Boltz-2 不可用")
        body = body or AiRescoreRequest()
        n_poses = len(job.result.poses)
        pose_indices = body.pose_indices or list(range(1, n_poses + 1))
        if any(i < 1 or i > n_poses for i in pose_indices):
            return _error(400, f"pose_indices 超出范围（1..{n_poses}）")
        if not pose_indices:
            return _error(400, "pose_indices 为空")
        task = AiRescoreTask(
            id=uuid.uuid4().hex[:8],
            job_id=job.id,
            pose_indices=pose_indices,
            created_at=datetime.now(),
            msa_server_used=body.use_msa_server,
        )
        rescore.submit(task, job)
        events_for(job.project_id).emit(
            "rescore.created",
            "rescore",
            task.id,
            {
                "job_id": job.id,
                "pose_indices": pose_indices,
                "msa_server_used": task.msa_server_used,
            },
        )
        return {"task_id": task.id}

    @app.get("/api/jobs/{jid}/rescore")
    async def list_rescore(jid: str):
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        return [t.model_dump(mode="json") for t in rescore.tasks_of_job(job)]

    @app.get("/api/rescore/{task_id}")
    async def get_rescore(task_id: str):
        task = rescore.get_task(task_id)
        if task is None:
            return _not_found("rescore task")
        return task

    @app.post("/api/rescore/{task_id}/cancel")
    async def cancel_rescore(task_id: str):
        task = rescore.get_task(task_id)
        if task is None:
            return _not_found("rescore task")
        job = runner.get_job(task.job_id)
        if job is None:
            return _not_found("rescore task job")
        result = rescore.cancel(task_id, job)
        return result if result is not None else task

    # ------------------------------------------------------------- M4: audit & provenance (2.0)

    @app.get("/api/projects/{pid}/events.jsonl", response_class=PlainTextResponse)
    async def project_events(pid: str, since: str | None = Query(default=None)):
        # `_prepared` is the pseudo-project that stores receptor/ligand prep
        # events.  It is not in the project registry, but audit needs to read it.
        if pid != "_prepared" and store.get(pid) is None:
            return _not_found("project")
        events = list(events_for(pid).read(since=since))
        body = "\n".join(json.dumps(e, ensure_ascii=False) for e in events)
        return PlainTextResponse(body + ("\n" if body else ""), media_type="text/plain")

    @app.get("/api/jobs/{jid}/audit.zip")
    def job_audit_zip(jid: str):
        from ..analysis.audit import build_audit_zip

        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        project = store.get(job.project_id)
        related_ids = {
            job.id,
            job.project_id,
            job.request.receptor_id,
            job.request.ligand_id,
        }
        if job.batch_id:
            related_ids.add(job.batch_id)

        def relevant(event: dict) -> bool:
            subject = event.get("subject") or {}
            return subject.get("id") in related_ids or jid in json.dumps(event.get("payload", {}))

        collected = [e for e in events_for(job.project_id).read() if relevant(e)] + [
            e for e in events_for("_prepared").read() if relevant(e)
        ]
        unique: dict[tuple, dict] = {}
        for event in collected:
            key = (
                event.get("ts"),
                event.get("action"),
                (event.get("subject") or {}).get("kind"),
                (event.get("subject") or {}).get("id"),
            )
            unique[key] = event
        events = sorted(unique.values(), key=lambda e: str(e.get("ts", "")))
        input_hashes = _load_input_hashes(job) or _input_hashes(job)
        receptor = registry.receptors.get(job.request.receptor_id)
        ligand = registry.ligands.get(job.request.ligand_id)
        prepared = {
            "receptor": receptor.model_dump(mode="json") if receptor is not None else None,
            "ligand": ligand.model_dump(mode="json") if ligand is not None else None,
        }
        zip_path = store.job_dir(job.project_id, job.id) / "audit.zip"
        build_audit_zip(
            job,
            store.root / job.project_id,
            events,
            zip_path,
            project_name=(project.name if project else job.project_id),
            input_hashes=input_hashes,
            prepared=prepared,
        )
        return Response(
            zip_path.read_bytes(),
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="audit-{jid}.zip"'},
        )

    # ------------------------------------------------------------- M4: demo + co-crystal (2.0)

    _DEMO_FLAG = "_demo_loaded.json"

    @app.get("/api/demo/example")
    def demo_example():
        flag_path = registry.receptors_dir / _DEMO_FLAG
        demo_dir = Path(__file__).resolve().parent.parent / "demo_data"
        note = (
            "示例数据：PDB 1HSG（HIV-1 蛋白酶 + 共晶配体 MK1），"
            "配体为 indinavir（SMILES 来自 PDB 共晶配体）；盒子取 MK1 质心 22.5Å 立方。"
        )
        flag = None
        if flag_path.is_file():
            try:
                flag = json.loads(flag_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                flag = None
        if (
            flag
            and flag.get("receptor_id") in registry.receptors
            and flag.get("ligand_id") in registry.ligands
        ):
            receptor = registry.receptors[flag["receptor_id"]]
            ligand = registry.ligands[flag["ligand_id"]]
        else:
            receptor = prepare_receptor(
                (demo_dir / "1HSG.pdb").read_bytes(), "1HSG.pdb", registry.receptors_dir
            )
            registry.add_receptor(receptor)
            smiles_text = (demo_dir / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()
            ligand = prepare_ligand_from_smiles(smiles_text, "indinavir", registry.ligands_dir)
            registry.add_ligand(ligand)
            flag = {"receptor_id": receptor.id, "ligand_id": ligand.id}
            flag_path.write_text(json.dumps(flag), encoding="utf-8")
            prep_events = events_for("_prepared")
            prep_events.emit(
                "receptor.prepared",
                "receptor",
                receptor.id,
                {
                    "filename": "1HSG.pdb",
                    "original_sha256": sha256_file(Path(receptor.original_path)),
                    "n_chains": len(receptor.chains),
                    "chains": receptor.chains,
                    "n_cocrystals": len(receptor.cocrystal_ligands),
                    "cocrystals": [c.model_dump(mode="json") for c in receptor.cocrystal_ligands],
                    "warnings": receptor.warnings,
                },
            )
            prep_events.emit(
                "ligand.prepared",
                "ligand",
                ligand.id,
                {
                    "source": ligand.source,
                    "smiles": ligand.smiles,
                    "n_rotatable": ligand.n_rotatable_bonds,
                    "name": ligand.name,
                },
            )
        mk1 = next((c for c in receptor.cocrystal_ligands if c.resname == "MK1"), None)
        if mk1 is None:
            raise PreparationError("示例受体未识别到共晶配体 MK1")
        from ..models import Box, DemoExample

        return DemoExample(
            receptor=receptor,
            ligand=ligand,
            box=Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
            note=note,
        )

    @app.get("/api/demo/library")
    async def demo_library() -> dict:
        demo_path = Path(__file__).resolve().parent.parent / "demo_data" / "demo_library.smi"
        items = []
        for line in demo_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            items.append(
                {"name": parts[1].strip() if len(parts) > 1 else parts[0], "smiles": parts[0]}
            )
        return {"name": "demo-library", "items": items}

    def _job_cocrystal_meta(job: Job) -> list[dict]:
        receptor = registry.receptors.get(job.request.receptor_id)
        if receptor is None:
            return []
        meta_path = Path(receptor.pdbqt_path).parent / "cocrystals.json"
        if not meta_path.is_file():
            return []
        try:
            return json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return []

    @app.get("/api/jobs/{jid}/crystal_ligand", response_class=PlainTextResponse)
    async def crystal_ligand(jid: str):
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        meta = _job_cocrystal_meta(job)
        if not meta:
            return _error(404, "该受体没有可用的共晶配体")
        sdf_path = Path(meta[0]["sdf_path"])
        if not sdf_path.is_file():
            return _error(404, "共晶配体 SDF 缺失")
        return PlainTextResponse(sdf_path.read_text(encoding="utf-8", errors="replace"))

    @app.get("/api/jobs/{jid}/crystal_rmsd")
    def crystal_rmsd(jid: str, pose_index: int = Query(default=1, ge=1)):
        from ..models import CrystalRmsd

        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        pose_paths = _job_pose_paths(job)
        if not (1 <= pose_index <= len(pose_paths)):
            return _error(404, f"pose {pose_index} 不存在")
        meta = _job_cocrystal_meta(job)
        if not meta:
            return CrystalRmsd(pose_index=pose_index, rmsd=None, reference=None)
        reference = f"{meta[0]['resname']}:{meta[0]['chain']}"
        sdf_path = Path(meta[0]["sdf_path"])
        if not sdf_path.is_file():
            return CrystalRmsd(pose_index=pose_index, rmsd=None, reference=reference)
        from ..analysis.crystal import heavy_atom_rmsd

        rmsd = heavy_atom_rmsd(sdf_path, pose_paths[pose_index - 1])
        return CrystalRmsd(pose_index=pose_index, rmsd=rmsd, reference=reference)

    # ------------------------------------------------------------- M1: WebSocket

    @app.websocket("/ws/jobs/{jid}")
    async def job_ws(websocket: WebSocket, jid: str) -> None:
        await websocket.accept()
        queue: asyncio.Queue[tuple[str, object]] = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def listener(job_id: str, event: tuple[str, object]) -> None:
            kind, payload = event

            def enqueue() -> None:
                queue.put_nowait((kind, payload))

            loop.call_soon_threadsafe(enqueue)

        unsubscribe = bus.subscribe(jid, listener)
        try:
            job = runner.get_job(jid)
            if job is not None:
                await websocket.send_text(
                    json.dumps(
                        {"type": "status", "data": job.model_dump(mode="json")},
                        ensure_ascii=False,
                    )
                )
            while True:
                kind, payload = await queue.get()
                data = (
                    payload.model_dump(mode="json") if isinstance(payload, BaseModel) else payload
                )
                await websocket.send_text(
                    json.dumps({"type": kind, "data": data}, ensure_ascii=False)
                )
        except WebSocketDisconnect:
            pass
        finally:
            unsubscribe()

    # ------------------------------------------------------------- M2: pockets

    @app.post("/api/pockets")
    def pockets_endpoint(body: PocketsRequest):
        from ..analysis.pockets import pocket_response

        receptor = registry.receptors.get(body.receptor_id)
        if receptor is None:
            return _error(400, f"unknown receptor_id: {body.receptor_id}")
        rec_dir = Path(receptor.pdbqt_path).parent
        protein = rec_dir / "receptor_protein.pdb"
        if not protein.is_file():
            protein = rec_dir / "receptor_clean.pdb"
        return pocket_response(body.receptor_id, protein, rec_dir, top_n=body.top_n)

    # ------------------------------------------------------------- M2: flex

    @app.post("/api/prepare/receptor/{rid}/flex", status_code=201)
    def flex_variant_endpoint(rid: str, body: FlexRequest):
        from ..prep.receptor import prepare_flex_variant

        base = registry.receptors.get(rid)
        if base is None:
            return _not_found("receptor")
        if base.parent_id:  # flex variants nest off the rigid base only
            base = registry.receptors.get(base.parent_id, base)
        variant = prepare_flex_variant(base, body.flex_residues, registry.receptors_dir)
        registry.add_receptor(variant)
        events_for("_prepared").emit(
            "receptor.flex_variant",
            "receptor",
            variant.id,
            {"parent_id": rid, "flex_residues": body.flex_residues},
        )
        return variant

    # ------------------------------------------------------------- M2: library

    @app.post("/api/prepare/library", status_code=201)
    async def library_endpoint(
        file: Annotated[UploadFile, File()],
        name: Annotated[str | None, Form()] = None,
    ):
        from ..prep.library import prepare_library

        content = await _read_upload(file)
        if not content.strip():
            raise PreparationError("库文件为空。")
        lib = await asyncio.get_running_loop().run_in_executor(
            None, prepare_library, content, file.filename or "library.smi", registry, name
        )
        events_for("_prepared").emit(
            "library.prepared",
            "library",
            lib.id,
            {
                "name": lib.name,
                "n_ligands": len(lib.ligand_ids),
                "n_failed": len(lib.failed),
                "filename": file.filename,
            },
        )
        return lib

    # ------------------------------------------------------------- M2: batches

    @app.post("/api/projects/{pid}/batches", status_code=201)
    async def submit_batch(pid: str, body: BatchCreate):
        from ..core.batch import job_request_from_batch

        if store.get(pid) is None:
            return _not_found("project")
        for ligand_id in body.ligand_ids:
            if ligand_id not in registry.ligands:
                return _error(400, f"unknown ligand_id: {ligand_id}")
        if body.receptor_id not in registry.receptors:
            return _error(400, f"unknown receptor_id: {body.receptor_id}")
        engine = get_engine(body.engine)
        if engine is not None:
            if body.params.scoring not in engine.manifest.scoring_functions:
                return _error(
                    400,
                    f"scoring '{body.params.scoring}' not supported by {body.engine} "
                    f"(supported: {engine.manifest.scoring_functions})",
                )
            ok, reason = engine.is_available()
            if not ok:
                return _error(503, reason or f"engine {body.engine} is unavailable")
        summary = batches.create(pid, body)
        for ligand_id in body.ligand_ids:
            job = store.new_job(pid, job_request_from_batch(body, ligand_id), batch_id=summary.id)
            job.env = env_fingerprint()
            _save_input_hashes(job)
            store.save_job(job)
            runner.submit(job)
        events_for(pid).emit(
            "batch.created",
            "batch",
            summary.id,
            {
                "name": summary.name,
                "total": summary.total,
                "receptor_id": body.receptor_id,
                "engine": body.engine,
                "ligand_ids": body.ligand_ids,
            },
        )
        return {"batch_id": summary.id}

    def _batch_summary_or_404(bid: str):
        jobs_mem = {jid: runner._jobs[jid] for jid in runner._jobs}
        summary = batches.summary(bid, runner_jobs=jobs_mem)
        if summary is None:
            return None
        return summary

    @app.get("/api/batches/{bid}")
    async def get_batch(bid: str):
        summary = _batch_summary_or_404(bid)
        if summary is None:
            return _not_found("batch")
        return summary

    @app.post("/api/batches/{bid}/cancel")
    async def cancel_batch(bid: str):
        if batches.get_meta(bid) is None:
            return _not_found("batch")
        for job in batches.jobs_of(bid):
            if not job.is_terminal:
                runner.cancel(job.id)
        meta = batches.get_meta(bid) or {}
        events_for(meta.get("project_id", "_prepared")).emit("batch.canceled", "batch", bid, {})
        return _batch_summary_or_404(bid)

    @app.get("/api/batches/{bid}/results")
    async def batch_results(
        bid: str,
        sort: str = Query(default="affinity", pattern="^(affinity|name)$"),
        order: str = Query(default="asc", pattern="^(asc|desc)$"),
        limit: int = Query(default=50, ge=1, le=10000),
        offset: int = Query(default=0, ge=0),
    ):
        result = batches.results(bid, registry, sort=sort, order=order, limit=limit, offset=offset)
        if result is None:
            return _not_found("batch")
        return {
            "total": result["total"],
            "rows": [r.model_dump(mode="json") for r in result["rows"]],
        }

    @app.get("/api/batches/{bid}/results.csv", response_class=PlainTextResponse)
    async def batch_results_csv(bid: str) -> Response:
        csv_text = batches.results_csv(bid, registry)
        if csv_text is None:
            return _not_found("batch")
        return PlainTextResponse(
            csv_text,
            media_type="text/csv",
            headers={"Content-Disposition": f'attachment; filename="batch-{bid}-results.csv"'},
        )

    @app.websocket("/ws/batches/{bid}")
    async def batch_ws(websocket: WebSocket, bid: str) -> None:
        await websocket.accept()
        queue: asyncio.Queue[tuple[str, object]] = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def listener(_job_id: str, _event: tuple[str, object]) -> None:
            def enqueue() -> None:
                summary = batches.summary(bid)
                if summary is not None:
                    kind = "finished" if summary.is_terminal else "batch_progress"
                    queue.put_nowait((kind, summary))

            loop.call_soon_threadsafe(enqueue)

        # subscribe to every job of the batch (discovery is lazy: jobs may be
        # added between page load and ws connect, so also send one snapshot now)
        unsubscribers = []
        try:
            for job in batches.jobs_of(bid):
                unsubscribers.append(bus.subscribe(job.id, listener))
            summary = batches.summary(bid)
            if summary is not None:
                await websocket.send_text(
                    json.dumps(
                        {"type": "batch_progress", "data": summary.model_dump(mode="json")},
                        ensure_ascii=False,
                    )
                )
            while True:
                kind, payload = await queue.get()
                data = (
                    payload.model_dump(mode="json") if isinstance(payload, BaseModel) else payload
                )
                await websocket.send_text(
                    json.dumps({"type": kind, "data": data}, ensure_ascii=False)
                )
        except WebSocketDisconnect:
            pass
        finally:
            for unsub in unsubscribers:
                unsub()

    return app


def main() -> None:
    from ..core.lifecycle import acquire_single_instance_lock

    settings = get_settings()
    if not acquire_single_instance_lock(settings.data_dir):
        print(
            "DockStudio backend is already running for this data directory; "
            "reusing the existing instance.",
            file=sys.stderr,
        )
        return
    # Started from the FastAPI lifespan, after all application imports complete.
    os.environ["DOCKSTUDIO_ENABLE_WATCHDOG"] = "1"
    uvicorn.run("dockstudio.api.app:create_app", factory=True, host="127.0.0.1", port=8000)
