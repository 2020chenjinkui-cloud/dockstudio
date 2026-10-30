"""Boltz-2 AI affinity rescore: discovery, YAML generation, task execution.

Boltz-2 is NOT bundled (multi-GB torch + weights). Users install it themselves;
DockStudio discovers it via DOCKSTUDIO_BOLTZ (venv python) or
DOCKSTUDIO_BOLTZ_CMD (arbitrary command prefix, e.g. `wsl boltz` / docker).

The rescore semantics: input is (protein sequence + ligand SMILES), pose
geometry is NOT fed to the model — Boltz-2 predicts affinity independently,
giving a cross-check signal against the docking score.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import threading
from pathlib import Path

from ..config import Settings
from ..core.project import ProjectStore
from ..logging import get_logger
from ..models import (
    AiAvailability,
    AiRescorePoseResult,
    AiRescoreTask,
    Job,
    JobStatus,
    PreparationError,
)

log = get_logger("analysis.ai_rescore")

DOCKSTUDIO_BOLTZ_ENV = "DOCKSTUDIO_BOLTZ"

_INSTALL_HINT = (
    "未检测到 Boltz-2。请任选其一：\n"
    "1) 独立 venv：`python -m venv boltz-venv && pip install boltz`，"
    "然后设置环境变量 DOCKSTUDIO_BOLTZ=<该 venv 的 python 路径>；\n"
    "2) 已装好的命令前缀：设置 DOCKSTUDIO_BOLTZ_CMD（如 `wsl boltz` 或"
    " docker run --rm --gpus all ghcr.io/jwohlwend/boltz）。\n"
    "注意 Boltz-2 本体 + 权重数十 GB，强烈建议 NVIDIA GPU；CPU 可跑但很慢。"
    "详见 docs/ai.md。"
)

_AA3 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


# ---------------------------------------------------------------------------
# availability
# ---------------------------------------------------------------------------


def get_availability(env: dict[str, str] | None = None) -> AiAvailability:
    """cmd prefix > venv python > unavailable (Chinese install hints)."""
    env = env if env is not None else os.environ
    cmd = (env.get("DOCKSTUDIO_BOLTZ_CMD") or "").strip()
    if cmd:
        return AiAvailability(available=True, mode="cmd")
    python = (env.get(DOCKSTUDIO_BOLTZ_ENV) or "").strip()
    if python:
        ok, detail = _probe_boltz_python(python, env)
        if ok:
            return AiAvailability(available=True, mode="venv")
        return AiAvailability(
            available=False,
            mode="venv",
            reason=(
                f"DOCKSTUDIO_BOLTZ={python} 指向的 Python 中无法导入 boltz"
                f"（{detail}）。请确认在该环境 `pip install boltz` 成功，或改用 DOCKSTUDIO_BOLTZ_CMD。"
            ),
        )
    return AiAvailability(available=False, reason=_INSTALL_HINT)


def _probe_boltz_python(python: str, env: dict[str, str] | None = None) -> tuple[bool, str]:
    full_env = dict(os.environ)
    if env:
        full_env.update({k: v for k, v in env.items() if k.startswith(("PYTHON", "HF_", "CONDA"))})
    try:
        proc = subprocess.run(
            [python, "-c", "import boltz; print(boltz.__version__)"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
            check=False,
            env=full_env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, str(exc)
    if proc.returncode == 0 and (proc.stdout or "").strip():
        return True, ""
    tail = ((proc.stderr or "") + (proc.stdout or "")).strip().splitlines()
    return False, (tail[-1] if tail else f"exit {proc.returncode}")[:200]


# ---------------------------------------------------------------------------
# YAML generation (pure-text sequence extraction, no Biopython)
# ---------------------------------------------------------------------------


def extract_chain_sequence(pdb_path: Path, prefer_chain: str | None = None) -> tuple[str, str]:
    """Pull the one-letter sequence of one chain from a raw PDB.

    Chain choice: prefer_chain if present, else the chain with the most residues.
    Returns (chain_id, sequence)."""
    chain_res: dict[str, dict[str, str]] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27:
            continue
        resname = line[17:20].strip()
        chain = line[21:22].strip() or "_"
        key = line[22:27].strip()
        if resname in _AA3:
            chain_res.setdefault(chain, {})[key] = resname
    if not chain_res:
        raise PreparationError(f"{pdb_path.name} 中没有可用的蛋白残基序列。")

    def natural(key: str) -> tuple[int, str]:
        m = re.match(r"(\d+)([A-Za-z]?)", key)
        return (int(m.group(1)) if m else 0, key)

    chain = prefer_chain if prefer_chain in chain_res else None
    if chain is None:
        chain = max(chain_res, key=lambda c: len(chain_res[c]))
    residues = chain_res[chain]
    sequence = "".join(_AA3[residues[k]] for k in sorted(residues, key=natural))
    return chain, sequence


def build_boltz_yaml(
    sequence: str, smiles: str, protein_id: str = "A", ligand_id: str = "L"
) -> str:
    """Boltz YAML for sequence + smiles affinity prediction."""
    smiles_q = smiles.replace("'", "''")
    return (
        "version: 1\n"
        "sequences:\n"
        f"  - protein:\n"
        f"      id: {protein_id}\n"
        f"      sequence: {sequence}\n"
        f"  - ligand:\n"
        f"      id: {ligand_id}\n"
        f"      smiles: '{smiles_q}'\n"
        "properties:\n"
        "  - affinity:\n"
        f"      binder: {ligand_id}\n"
    )


def parse_affinity_json(ai_dir: Path) -> tuple[float | None, float | None]:
    """Find predictions/*/affinity_*.json under the boltz output dir."""
    files = sorted(ai_dir.rglob("affinity_*.json"))
    if not files:
        raise PreparationError(f"Boltz-2 输出目录中没有 affinity JSON：{ai_dir}")
    try:
        data = json.loads(files[0].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PreparationError(f"无法解析 Boltz-2 输出 {files[0].name}: {exc}") from exc
    log_ic50 = data.get("affinity_pred_value")
    prob = data.get("affinity_probability_binary")
    return (
        float(log_ic50) if log_ic50 is not None else None,
        float(prob) if prob is not None else None,
    )


# ---------------------------------------------------------------------------
# task execution
# ---------------------------------------------------------------------------


def _error_hint(stderr: str) -> str:
    tail = (
        " | ".join([line.strip() for line in stderr.splitlines() if line.strip()][-5:])
        or "<no output>"
    )
    hint = ""
    if re.search(
        r"huggingface|hf.co|HF_ENDPOINT|Connection|timeout|SSLError", stderr, re.IGNORECASE
    ):
        hint = "（疑似 Hugging Face 权重下载失败：请设置 HF_ENDPOINT=https://hf-mirror.com 后重试）"
    return f"{tail}{hint}"


class RescoreRunner:
    """Serial AI-rescore executor: one task at a time, cancellable, crash-safe
    (task.json persisted on every state change; stale tasks -> interrupted)."""

    def __init__(self, store: ProjectStore, settings: Settings, registry=None):
        self.store = store
        self.settings = settings
        self.registry = registry
        self._tasks: dict[str, AiRescoreTask] = {}
        self._current_proc: subprocess.Popen | None = None
        self._lock = threading.Lock()
        self._wake = threading.Event()
        self._cancel_requested: set[str] = set()
        self._thread: threading.Thread | None = None

    # ------------------------------------------------------------ discovery

    def _ai_dir(self, job: Job, task_id: str) -> Path:
        return self.store.job_dir(job.project_id, job.id) / "ai" / task_id

    def _task_path(self, job: Job, task_id: str) -> Path:
        return self._ai_dir(job, task_id) / "task.json"

    def _load_task(self, job: Job, task_id: str) -> AiRescoreTask | None:
        path = self._task_path(job, task_id)
        if not path.is_file():
            return None
        try:
            return AiRescoreTask(**json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError, ValueError):
            return None

    # ------------------------------------------------------------ lifecycle

    def start(self) -> None:
        self.mark_interrupted()
        if self._thread is None or not self._thread.is_alive():
            self._wake.clear()
            self._thread = threading.Thread(target=self._worker, daemon=True)
            self._thread.start()

    def stop(self) -> None:
        self._wake.set()

    def mark_interrupted(self) -> int:
        """Non-terminal task.json files left by a crashed process -> interrupted."""
        from datetime import datetime

        count = 0
        for path in self.store.root.glob("*/jobs/*/ai/*/task.json"):
            try:
                task = AiRescoreTask(**json.loads(path.read_text(encoding="utf-8")))
            except (OSError, json.JSONDecodeError, ValueError):
                continue
            if task.status in (JobStatus.queued, JobStatus.running):
                task.status = JobStatus.interrupted
                task.finished_at = task.finished_at or datetime.now()
                tmp = path.with_suffix(".tmp")
                tmp.write_text(task.model_dump_json(indent=2), encoding="utf-8")
                tmp.replace(path)
                self._tasks[task.id] = task
                count += 1
        return count

    def get_task(self, task_id: str) -> AiRescoreTask | None:
        if task_id in self._tasks:
            return self._tasks[task_id]
        for path in self.store.root.glob("*/jobs/*/ai/*/task.json"):
            if path.parent.name == task_id:
                try:
                    task = AiRescoreTask(**json.loads(path.read_text(encoding="utf-8")))
                    self._tasks[task_id] = task
                    return task
                except (OSError, json.JSONDecodeError, ValueError):
                    return None
        return None

    def tasks_of_job(self, job: Job) -> list[AiRescoreTask]:
        out = []
        ai_root = self.store.job_dir(job.project_id, job.id) / "ai"
        if ai_root.is_dir():
            for path in ai_root.glob("*/task.json"):
                try:
                    out.append(AiRescoreTask(**json.loads(path.read_text(encoding="utf-8"))))
                except (OSError, json.JSONDecodeError, ValueError):
                    continue
        for task in self._tasks.values():
            if task.job_id == job.id and all(t.id != task.id for t in out):
                out.append(task)
        out.sort(key=lambda t: t.created_at, reverse=True)
        return out

    def submit(self, task: AiRescoreTask, job: Job) -> None:
        self._tasks[task.id] = task
        task_dir = self._ai_dir(job, task.id)
        task_dir.mkdir(parents=True, exist_ok=True)
        self._save(task, job)
        self._wake.set()

    def cancel(self, task_id: str, job: Job) -> AiRescoreTask | None:
        task = self._tasks.get(task_id) or self._load_task(job, task_id)
        if task is None:
            return None
        if task.status == JobStatus.queued:
            self._cancel_requested.add(task_id)
        elif task.status == JobStatus.running:
            with self._lock:
                proc = self._current_proc
            if proc is not None and proc.poll() is None:
                try:
                    proc.kill()
                except OSError:
                    pass
        return task

    # ------------------------------------------------------------ worker

    def _save(self, task: AiRescoreTask, job: Job) -> None:
        self._tasks[task.id] = task
        path = self._task_path(job, task.id)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(task.model_dump_json(indent=2), encoding="utf-8")
        tmp.replace(path)

    def _worker(self) -> None:
        while True:
            self._wake.wait(timeout=2.0)
            self._wake.clear()
            pending = [t for t in self._tasks.values() if t.status == JobStatus.queued]
            if not pending:
                continue
            task = min(pending, key=lambda t: t.created_at)
            if task.id in self._cancel_requested:
                self._cancel_requested.discard(task.id)
                job_obj = self.store.find_job(task.job_id)
                if job_obj is not None:
                    task.status = JobStatus.canceled
                    self._save(task, job_obj)
                continue
            job = self.store.find_job(task.job_id)
            if job is None:
                continue
            try:
                self._run(task, job)
            except Exception as exc:  # noqa: BLE001
                if task.id in self._cancel_requested:
                    log.info("rescore task %s canceled", task.id)
                    task.status = JobStatus.canceled
                    task.progress.message = "已取消"
                else:
                    log.error("rescore task %s crashed: %s", task.id, exc)
                    task.status = JobStatus.failed
                    task.error = str(exc)
                from datetime import datetime

                task.finished_at = task.finished_at or datetime.now()
                self._save(task, job)

    def _run(self, task: AiRescoreTask, job: Job) -> None:
        from datetime import datetime

        available = get_availability()
        if not available.available:
            raise PreparationError(available.reason or "Boltz-2 不可用")

        task.status = JobStatus.running
        task.progress.stage = "preparing"
        task.progress.message = "构建 Boltz-2 输入"
        task.progress.percent = 5.0
        self._save(task, job)

        task_dir = self._ai_dir(job, task.id)
        yaml_path, _seq, _smiles = self._write_inputs(job, task_dir)
        task.progress.stage = "scoring"
        task.progress.message = "Boltz-2 推理（配体级预测，与姿态几何无关）"
        task.progress.percent = 20.0
        self._save(task, job)
        if task.id in self._cancel_requested:
            task.status = JobStatus.canceled
            task.finished_at = datetime.now()
            self._save(task, job)
            return

        # Boltz-2's affinity head consumes protein sequence + ligand SMILES; pose
        # coordinates are not part of the input.  Running it once per pose would
        # burn identical GPU work N times and produce identical numbers, so run it
        # exactly once and expose the result as a ligand-level prediction.
        out_dir = task_dir / "out_ligand"
        out_dir.mkdir(parents=True, exist_ok=True)
        self._invoke_boltz(available, yaml_path, out_dir, task.msa_server_used)
        log_ic50, prob = parse_affinity_json(out_dir)
        task.results = [
            AiRescorePoseResult(
                pose_index=pose_index,
                log_ic50=log_ic50,
                binder_probability=prob,
                rank_by_score=None,
            )
            for pose_index in task.pose_indices
        ]
        self._save(task, job)

        if task.id in self._cancel_requested:
            task.status = JobStatus.canceled
        else:
            task.status = JobStatus.completed
            task.progress.stage = "done"
            task.progress.percent = 100.0
            task.progress.message = "完成（配体级预测）"
        task.finished_at = datetime.now()
        self._save(task, job)
        self._lifecycle(
            job,
            task,
            "rescore.completed" if task.status == JobStatus.completed else "rescore.canceled",
            {
                "n_results": len(task.results),
                "results": [r.model_dump(mode="json") for r in task.results]
                if task.status == JobStatus.completed
                else [],
            },
        )

    def _lifecycle(self, job: Job, task: AiRescoreTask, action: str, payload: dict) -> None:
        try:
            from ..core.events import EventLog

            EventLog(self.store.root / job.project_id).emit(
                action, "rescore", task.id, payload, actor="system"
            )
        except Exception as exc:  # noqa: BLE001
            log.warning("rescore event emit failed (%s): %s", action, exc)

    def _write_inputs(self, job: Job, task_dir: Path) -> tuple[Path, str, str]:
        from ..core.orchestrator import get_registry

        registry = self.registry or get_registry(self.settings)
        receptor = registry.receptors.get(job.request.receptor_id)
        if receptor is None:
            raise PreparationError(f"unknown receptor_id: {job.request.receptor_id}")
        ligand = registry.ligands.get(job.request.ligand_id)
        if ligand is None or not ligand.smiles:
            raise PreparationError("配体缺少 SMILES，无法做 AI 复核。")
        original = Path(receptor.original_path)
        preferred = receptor.cocrystal_ligands[0].chain if receptor.cocrystal_ligands else None
        chain, sequence = extract_chain_sequence(original, prefer_chain=preferred)
        yaml_text = build_boltz_yaml(sequence, ligand.smiles)
        yaml_path = task_dir / "input.yaml"
        yaml_path.write_text(yaml_text, encoding="utf-8")
        (task_dir / "meta.json").write_text(
            json.dumps({"chain": chain, "n_residues": len(sequence), "smiles": ligand.smiles}),
            encoding="utf-8",
        )
        return yaml_path, sequence, ligand.smiles

    def _invoke_boltz(
        self, available: AiAvailability, yaml_path: Path, out_dir: Path, use_msa: bool
    ) -> None:
        if available.mode == "cmd":
            prefix = shlex.split(os.environ["DOCKSTUDIO_BOLTZ_CMD"], posix=True)
        else:
            python = os.environ.get("DOCKSTUDIO_BOLTZ", "")
            prefix = [python, "-m", "boltz"]
        cmd = prefix + ["predict", str(yaml_path), "--out_dir", str(out_dir)]
        if use_msa:
            cmd.append("--use_msa_server")
        log.info("running boltz: %s", " ".join(cmd))
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                cwd=out_dir,
            )
        except OSError as exc:
            raise PreparationError(f"启动 Boltz-2 失败：{exc}") from exc
        with self._lock:
            self._current_proc = proc
        try:
            try:
                stdout, stderr = proc.communicate(timeout=self.settings.job_timeout_s)
            except subprocess.TimeoutExpired as exc:
                proc.kill()
                raise PreparationError(f"Boltz-2 超时（>{self.settings.job_timeout_s}s）") from exc
        finally:
            with self._lock:
                if self._current_proc is proc:
                    self._current_proc = None
        out_log = out_dir / "boltz.log"
        out_log.write_text((stdout or "") + "\n--- stderr ---\n" + (stderr or ""), encoding="utf-8")
        if proc.returncode != 0:
            raise PreparationError(f"Boltz-2 退出码 {proc.returncode}：{_error_hint(stderr or '')}")
