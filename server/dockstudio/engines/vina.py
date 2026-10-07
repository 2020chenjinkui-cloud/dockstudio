"""AutoDock Vina engine adapter: subprocess wrapper around the official binary."""

from __future__ import annotations

import os
import platform
import re
import shutil
import signal
import subprocess
import threading
from pathlib import Path

from ..config import Settings, get_settings, resolve_vina_binary
from ..logging import get_logger
from ..models import Box, DockingError, DockParams, DockResult, EngineUnavailableError, Pose
from .base import EngineManifest

log = get_logger("engines.vina")

_RESULT_RE = re.compile(
    r"^\s*REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)",
    re.IGNORECASE,
)

# Error knowledge base: common vina/meeko stderr patterns -> Chinese advice.
_ERROR_KB: list[tuple[re.Pattern[str], str]] = [
    (
        re.compile(r"PDBQT parsing error|parse error.*pdbqt|could not parse", re.IGNORECASE),
        "PDBQT 文件解析失败：请重新准备受体/配体（meeko 输出可能不完整），并检查文件未被手工编辑。",
    ),
    (
        re.compile(r"missing residues|not enough residues|cannot find residue", re.IGNORECASE),
        "受体中存在缺失残基或非标准残基：请在准备阶段修复 PDB（补全缺失残基或移除非标准残基）后重试。",
    ),
    (
        re.compile(r"charge|Gasteiger|electrostatic", re.IGNORECASE),
        "原子电荷计算失败：可能含未知原子类型，请在配体准备时检查元素组成，必要时删除金属/稀有原子。",
    ),
    (
        re.compile(
            r"grid|outside the box|search space is too|too large|box.*exceed", re.IGNORECASE
        ),
        "搜索盒设置超限：请缩小 box 尺寸（建议各边 ≤ 30 Å）或将盒子中心移到结合位点。",
    ),
    (
        re.compile(r"ligand is outside|outside the search space", re.IGNORECASE),
        "配体初始坐标不在搜索盒内：请核对 box 中心是否覆盖结合口袋（可用共晶配体质心自动定盒）。",
    ),
    (
        re.compile(r"no torsion|0 torsions|rigid ligand", re.IGNORECASE),
        "配体无可旋转键：对接退化为刚体放置，打分可能不具区分度，属正常现象。",
    ),
    (
        re.compile(r"exhaustiveness.*out of range|invalid.*exhaustiveness", re.IGNORECASE),
        "exhaustiveness 参数超出合法范围：请使用 1-64 之间的整数。",
    ),
    (
        re.compile(r"atom type.*unknown|unknown atom type|unsupported element", re.IGNORECASE),
        "存在 Vina 不支持的原子类型：请检查配体元素（Vina 仅支持常见有机元素），必要时删除重原子。",
    ),
]


def annotate_error(text: str) -> str:
    """Append a Chinese hint from the error knowledge base, if one matches."""
    for pattern, hint in _ERROR_KB:
        if pattern.search(text):
            return f"{hint}（原始错误片段：{_key_lines(text, 3)}）"
    return ""


def parse_vina_out(out_path: Path) -> list[Pose]:
    """Parse a vina multi-MODEL PDBQT output into ordered poses."""
    poses: list[Pose] = []
    for line in out_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _RESULT_RE.match(line)
        if m:
            poses.append(
                Pose(
                    index=len(poses) + 1,
                    affinity=float(m.group(1)),
                    rmsd_lb=float(m.group(2)),
                    rmsd_ub=float(m.group(3)),
                )
            )
    if not poses:
        raise DockingError(f"vina output contained no binding poses: {out_path}")
    return poses


def _key_lines(text: str, limit: int = 5) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return " | ".join(lines[-limit:]) if lines else "<no output>"


class VinaEngine:
    """Runs vina.exe via subprocess (no PyPI wheel on Windows)."""

    manifest = EngineManifest(
        name="vina",
        scoring_functions=["vina", "vinardo"],
        supports_flex_residues=True,
        homepage="https://vina.scripps.edu",
    )

    def __init__(self, settings: Settings | None = None):
        self.settings = settings or get_settings()
        self._current: subprocess.Popen | None = None
        self._lock = threading.Lock()

    def _binary(self) -> Path:
        if platform.system().lower() not in self.manifest.platforms:
            raise EngineUnavailableError(f"vina is not supported on {platform.system()}")
        binary = resolve_vina_binary(self.settings)
        if binary is None:
            binary = self._bundled_fallback()
        if binary is None:
            raise EngineUnavailableError(
                "vina binary not found; set DOCKSTUDIO_VINA or bundle vina under resources/bin"
            )
        return binary

    @staticmethod
    def _bundled_fallback() -> Path | None:
        """Repo-root fallback: this file lives at repo/server/dockstudio/engines/vina.py,
        so parents[3] is the repo root (config.bundled_bin_dir assumes a flatter layout)."""
        exe = "vina.exe" if platform.system() == "Windows" else "vina"
        root = Path(__file__).resolve().parents[3]
        for name in ("windows", "linux-x86_64", "linux-aarch64", "mac-x86_64", "mac-aarch64"):
            candidate = root / "resources" / "bin" / name / exe
            if candidate.is_file():
                return candidate
        return None

    def _version(self, binary: Path) -> str:
        try:
            proc = subprocess.run(
                [str(binary), "--version"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise EngineUnavailableError(f"failed to run vina --version: {exc}") from exc
        output = (proc.stdout or "") + (proc.stderr or "")
        if proc.returncode != 0:
            raise EngineUnavailableError(
                f"vina --version exited {proc.returncode}: {_key_lines(output)}"
            )
        m = re.search(r"(\d+\.\d+(?:\.\d+)?)", output)
        return m.group(1) if m else "unknown"

    def is_available(self) -> tuple[bool, str | None]:
        try:
            binary = self._binary()
            self.manifest.version = self._version(binary)
        except EngineUnavailableError as exc:
            return False, str(exc)
        return True, None

    def _subprocess_env(self) -> dict[str, str] | None:
        """Extra environment for the engine subprocess (smina needs BABEL_LIBDIR)."""
        return None

    def _parse_poses(self, out_path: Path, run_log: str) -> list[Pose]:
        return parse_vina_out(out_path)

    def cancel_current(self) -> None:
        """Best-effort kill of the currently running vina process tree."""
        with self._lock:
            proc = self._current
        if proc is None or proc.poll() is not None:
            return
        log.info("cancel requested, killing vina pid=%s", proc.pid)
        try:
            if platform.system() == "Windows":
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                    capture_output=True,
                    check=False,
                )
            else:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    proc.kill()
        except OSError as exc:
            log.warning("taskkill failed: %s", exc)

    def _wait(self, proc: subprocess.Popen) -> tuple[str, str]:
        """Wait for the subprocess with a timeout; raise DockingError on timeout."""
        try:
            stdout, stderr = proc.communicate(timeout=self.settings.job_timeout_s)
        except subprocess.TimeoutExpired as exc:
            self.cancel_current()
            raise DockingError(
                f"vina timed out after {self.settings.job_timeout_s}s: "
                f"{_key_lines((exc.stderr or b'') if isinstance(exc.stderr, bytes) else (exc.stderr or ''))}"
            ) from exc
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf-8", errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace")
        return stdout or "", stderr or ""

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        binary = self._binary()
        version = self._version(binary)
        self.manifest.version = version
        work_dir.mkdir(parents=True, exist_ok=True)
        out_path = work_dir / "out.pdbqt"

        # vina.exe (MinGW build) cannot reliably open non-ASCII paths; stage
        # inputs in work_dir and pass bare filenames with cwd=work_dir.
        local_receptor = work_dir / "receptor.pdbqt"
        local_ligand = work_dir / "ligand.pdbqt"
        if receptor_pdbqt.resolve() != local_receptor.resolve():
            shutil.copy2(receptor_pdbqt, local_receptor)
        if ligand_pdbqt.resolve() != local_ligand.resolve():
            shutil.copy2(ligand_pdbqt, local_ligand)
        flex_name: str | None = None
        if flex_pdbqt is not None:
            local_flex = work_dir / "flex.pdbqt"
            if flex_pdbqt.resolve() != local_flex.resolve():
                shutil.copy2(flex_pdbqt, local_flex)
            flex_name = local_flex.name

        cmd = [
            str(binary),
            "--receptor",
            local_receptor.name,
            "--ligand",
            local_ligand.name,
            "--center_x",
            str(box.center[0]),
            "--center_y",
            str(box.center[1]),
            "--center_z",
            str(box.center[2]),
            "--size_x",
            str(box.size[0]),
            "--size_y",
            str(box.size[1]),
            "--size_z",
            str(box.size[2]),
            "--exhaustiveness",
            str(params.exhaustiveness),
            "--num_modes",
            str(params.num_modes),
            "--energy_range",
            str(params.energy_range),
            "--cpu",
            str(params.cpu),
            "--scoring",
            params.scoring,
            "--out",
            out_path.name,
        ]
        if params.seed is not None:
            cmd += ["--seed", str(params.seed)]
        if flex_name is not None:
            cmd += ["--flex", flex_name]

        log.info("running vina: %s", " ".join(cmd))
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                cwd=work_dir,
                env=self._subprocess_env(),
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
                if platform.system() == "Windows"
                else 0,
                start_new_session=platform.system() != "Windows",
            )
        except OSError as exc:
            raise EngineUnavailableError(f"failed to launch vina: {exc}") from exc

        with self._lock:
            self._current = proc
        try:
            stdout, stderr = self._wait(proc)
        finally:
            with self._lock:
                if self._current is proc:
                    self._current = None

        run_log = (stdout or "") + "\n--- stderr ---\n" + (stderr or "")
        if proc.returncode != 0:
            hint = annotate_error(run_log)
            msg = f"vina exited with code {proc.returncode}: {_key_lines(run_log)}"
            raise DockingError(f"{msg}\n{hint}" if hint else msg)
        if not out_path.is_file():
            hint = annotate_error(run_log)
            msg = f"vina produced no output file: {_key_lines(run_log)}"
            raise DockingError(f"{msg}\n{hint}" if hint else msg)

        poses = self._parse_poses(out_path, run_log)
        return DockResult(
            engine=self.manifest.name,
            engine_version=version,
            params=params,
            box=box,
            poses=poses,
            log=run_log,
            run_dir=str(work_dir),
        )
