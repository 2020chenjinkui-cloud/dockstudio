"""FastAPI application for DockStudio (M0 sync chain + M1 projects/jobs/WS)."""

from __future__ import annotations

import asyncio
import json
import traceback
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
from ..core.orchestrator import PrepareRegistry, new_run_dir
from ..core.project import ProjectStore
from ..core.queue import EventBus, JobRunner
from ..engines import available_engines, get_engine
from ..logging import get_logger, setup_logging
from ..models import (
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

_ALLOWED_ORIGINS = ["http://localhost:5173", "http://127.0.0.1:5173"]


class CreateProjectRequest(BaseModel):
    name: str


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
    runner = JobRunner(store, registry, bus)

    app = FastAPI(title="DockStudio", version=__version__)
    app.state.settings = settings
    app.state.registry = registry
    app.state.store = store
    app.state.runner = runner

    app.add_middleware(
        CORSMiddleware,
        allow_origins=_ALLOWED_ORIGINS,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.on_event("startup")
    async def _startup() -> None:
        runner.start()

    @app.on_event("shutdown")
    async def _shutdown() -> None:
        await runner.stop()

    @app.exception_handler(DockStudioError)
    async def dockstudio_error_handler(_request: Request, exc: DockStudioError) -> JSONResponse:
        status = 503 if isinstance(exc, EngineUnavailableError) else 400
        return _error(status, str(exc))

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
        content = await file.read()
        receptor = prepare_receptor(
            content, file.filename or "receptor.pdb", registry.receptors_dir
        )
        registry.add_receptor(receptor)
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
        if has_file:
            assert file is not None and file.filename is not None
            content = await file.read()
            tmp_path = registry.ligands_dir / f"upload_{file.filename}"
            tmp_path.write_bytes(content)
            try:
                ligand = prepare_ligand_from_file(tmp_path, registry.ligands_dir)
            finally:
                tmp_path.unlink(missing_ok=True)
        else:
            ligand = prepare_ligand_from_smiles(smiles or "", name, registry.ligands_dir)
        registry.add_ligand(ligand)
        return ligand

    @app.post("/api/dock")
    async def dock_endpoint(req: DockJobRequest) -> DockResult:
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

        run_dir = new_run_dir(settings)
        result = engine.dock(
            Path(receptor.pdbqt_path),
            Path(ligand.pdbqt_path),
            req.box,
            req.params,
            run_dir,
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
        return store.create(body.name).model_dump(mode="json")

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
        if not store.delete(pid):
            return _not_found("project")
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
        job = store.new_job(pid, body)
        runner.submit(job)
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
    async def job_analysis(jid: str, pose_index: int = Query(default=1, ge=1)):
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
        )
        return analyses[0]

    # ------------------------------------------------------------- M1: report

    @app.get("/api/jobs/{jid}/report", response_class=HTMLResponse)
    async def job_report(jid: str) -> Response:
        job = runner.get_job(jid)
        if job is None:
            return _not_found("job")
        if job.status != JobStatus.completed or job.result is None:
            return _error(400, f"job is not completed (status={job.status.value})")
        project = store.get(job.project_id)
        pose_paths = _job_pose_paths(job)
        job_dir = store.job_dir(job.project_id, job.id)
        protein_pdb = _job_protein_pdb(registry, job)
        analyses = analyze_job_poses(protein_pdb, pose_paths[:1], job_dir, pose_indices=[1])
        ligand = registry.ligands.get(job.request.ligand_id)
        if project is None:  # project was soft-deleted; build a transient shell
            project = Project(
                id=job.project_id,
                name="（已删除项目）",
                path=str(store.job_dir(job.project_id, job.id).parent.parent),
                created_at=datetime.now(),  # noqa: DTZ005
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
        )
        return HTMLResponse(
            report_html,
            headers={"Content-Disposition": f'attachment; filename="dockstudio-report-{jid}.html"'},
        )

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

    return app


def main() -> None:
    uvicorn.run("dockstudio.api.app:create_app", factory=True, host="127.0.0.1", port=8000)


app = create_app()
