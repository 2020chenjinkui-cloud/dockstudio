"""gnina engine adapter (Linux + NVIDIA GPU only; unavailable elsewhere with a
Chinese explanation). Reuses the vina subprocess adapter; parses gnina PDBQT
output (`REMARK minimizedAffinity` per MODEL plus CNNscore/CNNaffinity rows in
stdout)."""

from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

from ..models import DockingError, EngineUnavailableError, Pose
from .base import EngineManifest
from .vina import VinaEngine, _key_lines

_AFF_RE = re.compile(r"^\s*REMARK\s+minimizedAffinity\s+(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_CNN_RE = re.compile(
    r"CNNscore[:=]\s*(-?\d+(?:\.\d+)?).*?CNNaffinity[:=]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE
)


def parse_gnina_out(out_path: Path, run_log: str) -> list[Pose]:
    """gnina out.pdbqt: `REMARK minimizedAffinity <x>` per model."""
    poses: list[Pose] = []
    if out_path.is_file():
        for line in out_path.read_text(encoding="utf-8", errors="replace").splitlines():
            m = _AFF_RE.match(line)
            if m:
                poses.append(Pose(index=len(poses) + 1, affinity=float(m.group(1))))
    # fallback: stdout mode table (vina-style) if present
    if not poses:
        for line in run_log.splitlines():
            m = re.match(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)", line)
            if m:
                poses.append(Pose(index=int(m.group(1)), affinity=float(m.group(2))))
    if not poses:
        raise DockingError(f"gnina output contained no binding poses: {_key_lines(run_log)}")
    return poses


class GninaEngine(VinaEngine):
    """gnina: deep-learning rescoring fork of Vina."""

    manifest = EngineManifest(
        name="gnina",
        scoring_functions=["vina", "vinardo", "ad4_scoring", "dkoes_fast"],
        supports_flex_residues=True,
        requires_gpu=True,
        platforms=["linux"],
        homepage="https://github.com/gnina/gnina",
    )

    def _binary(self) -> Path:
        if platform.system() != "Linux":
            raise EngineUnavailableError(
                "gnina 需要 Linux + NVIDIA GPU，Windows 用户请用 WSL2"
                "（指引见 docs/deploy.md#wsl2）。"
            )
        exe = "gnina"
        env = os.environ.get("DOCKSTUDIO_GNINA")
        if env and Path(env).is_file():
            return Path(env)
        found = shutil.which(exe)
        if found:
            return Path(found)
        candidate = Path(__file__).resolve().parents[3] / "resources" / "bin" / "linux-x86_64" / exe
        if candidate.is_file():
            return candidate
        raise EngineUnavailableError(
            "gnina 二进制未找到：设置 DOCKSTUDIO_GNINA 或放入 resources/bin/linux-x86_64。"
        )

    def _gpu_ok(self) -> tuple[bool, str]:
        if shutil.which("nvidia-smi"):
            try:
                proc = subprocess.run(
                    ["nvidia-smi", "-L"],
                    capture_output=True,
                    timeout=15,
                    check=False,
                )
                if proc.returncode == 0 and b"GPU" in (proc.stdout or b""):
                    return True, ""
            except (OSError, subprocess.TimeoutExpired):
                pass
        for candidate in ("libcuda.so", "libcuda.so.1"):
            for libdir in ("/usr/lib/x86_64-linux-gnu", "/usr/lib64", "/usr/local/cuda/lib64"):
                if (Path(libdir) / candidate).exists():
                    return True, ""
        return False, "检测到 gnina 但无 NVIDIA GPU，可用 --no_gpu 强制 CPU（极慢）"

    def is_available(self) -> tuple[bool, str | None]:
        try:
            binary = self._binary()
            self.manifest.version = self._version(binary)
        except EngineUnavailableError as exc:
            return False, str(exc)
        ok, reason = self._gpu_ok()
        if not ok:
            return False, reason
        return True, None

    def _parse_poses(self, out_path: Path, run_log: str) -> list[Pose]:
        return parse_gnina_out(out_path, run_log)
