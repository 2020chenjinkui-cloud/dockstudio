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
from typing import Any

from ..config import bundled_engine_binary
from ..models import DockingError, EngineUnavailableError, Pose
from .base import EngineManifest
from .vina import VinaEngine, _key_lines

_AFF_RE = re.compile(r"^\s*REMARK\s+minimizedAffinity\s+(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_CNN_RE = re.compile(
    r"CNNscore[:=\s]\s*(-?\d+(?:\.\d+)?).*?CNNaffinity[:=\s]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE
)
_CNN_SCORE_RE = re.compile(r"CNNscore[:=\s]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_CNN_AFF_RE = re.compile(r"CNNaffinity[:=\s]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_MODEL_RE = re.compile(r"^\s*MODEL\s+(\d+)", re.IGNORECASE)
_GNINA_STDOUT_TABLE_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)"
)


def parse_gnina_out(out_path: Path, run_log: str) -> list[Pose]:
    """gnina out.pdbqt: `REMARK minimizedAffinity <x>`, `REMARK CNNscore <y>`,
    `REMARK CNNaffinity <z>` per model, plus stdout mode table / log support."""
    poses: list[Pose] = []
    if out_path.is_file():
        current_pose: dict[str, Any] | None = None
        for line in out_path.read_text(encoding="utf-8", errors="replace").splitlines():
            m_mod = _MODEL_RE.match(line)
            if m_mod:
                if current_pose is not None and current_pose.get("affinity") is not None:
                    poses.append(Pose(**current_pose))
                current_pose = {
                    "index": int(m_mod.group(1)),
                    "affinity": None,
                    "cnn_score": None,
                    "cnn_affinity": None,
                }
                continue

            m_aff = _AFF_RE.match(line)
            if m_aff:
                if current_pose is None:
                    current_pose = {
                        "index": len(poses) + 1,
                        "affinity": float(m_aff.group(1)),
                        "cnn_score": None,
                        "cnn_affinity": None,
                    }
                else:
                    current_pose["affinity"] = float(m_aff.group(1))

            m_comb = _CNN_RE.search(line)
            if m_comb and current_pose is not None:
                current_pose["cnn_score"] = float(m_comb.group(1))
                current_pose["cnn_affinity"] = float(m_comb.group(2))
            else:
                m_cs = _CNN_SCORE_RE.search(line)
                if m_cs and current_pose is not None:
                    current_pose["cnn_score"] = float(m_cs.group(1))
                m_ca = _CNN_AFF_RE.search(line)
                if m_ca and current_pose is not None:
                    current_pose["cnn_affinity"] = float(m_ca.group(1))

        if current_pose is not None and current_pose.get("affinity") is not None:
            poses.append(Pose(**current_pose))

    # fallback: stdout mode table (gnina 4-column or vina-style) if present
    if not poses:
        has_gnina_table_header = bool(
            re.search(r"CNNscore.*?CNNaffinity", run_log, re.IGNORECASE)
        )
        for line in run_log.splitlines():
            if has_gnina_table_header:
                m_gnina = _GNINA_STDOUT_TABLE_RE.match(line)
                if m_gnina:
                    poses.append(
                        Pose(
                            index=int(m_gnina.group(1)),
                            affinity=float(m_gnina.group(2)),
                            cnn_score=float(m_gnina.group(3)),
                            cnn_affinity=float(m_gnina.group(4)),
                        )
                    )
                    continue
            m_vina = re.match(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)", line)
            if m_vina:
                poses.append(Pose(index=int(m_vina.group(1)), affinity=float(m_vina.group(2))))

    # Supplement/backfill CNN score & CNN affinity from run_log if missing on poses
    cnn_pairs: list[tuple[float, float]] = []
    for m in _CNN_RE.finditer(run_log):
        cnn_pairs.append((float(m.group(1)), float(m.group(2))))
    if not cnn_pairs:
        sc = _CNN_SCORE_RE.findall(run_log)
        af = _CNN_AFF_RE.findall(run_log)
        if sc and af and len(sc) == len(af):
            cnn_pairs = [(float(s), float(a)) for s, a in zip(sc, af)]

    for i, (cs, ca) in enumerate(cnn_pairs):
        if i < len(poses):
            if poses[i].cnn_score is None:
                poses[i].cnn_score = cs
            if poses[i].cnn_affinity is None:
                poses[i].cnn_affinity = ca

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
        bundled = bundled_engine_binary(exe)
        if bundled is not None:
            return bundled
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
        return False, (
            "检测到 gnina 但无 NVIDIA GPU；设置 DOCKSTUDIO_GNINA_CPU=1 可强制 CPU（--no_gpu，极慢）"
        )

    def is_available(self) -> tuple[bool, str | None]:
        try:
            binary = self._binary()
            self.manifest.version = self._version(binary)
        except EngineUnavailableError as exc:
            return False, str(exc)
        ok, reason = self._gpu_ok()
        if not ok and os.environ.get("DOCKSTUDIO_GNINA_CPU") != "1":
            return False, reason
        return True, None

    def _build_cmd(
        self,
        binary: str,
        box,
        params,
        receptor_name: str,
        ligand_name: str,
        out_name: str,
    ) -> list[str]:
        """gnina CLI 不含 --energy_range（vina 专有），从继承的命令里剔除。"""
        cmd = super()._build_cmd(binary, box, params, receptor_name, ligand_name, out_name)
        if "--energy_range" in cmd:
            i = cmd.index("--energy_range")
            del cmd[i : i + 2]
        if os.environ.get("DOCKSTUDIO_GNINA_CPU") == "1":
            cmd.append("--no_gpu")
        return cmd

    def _parse_poses(self, out_path: Path, run_log: str) -> list[Pose]:
        return parse_gnina_out(out_path, run_log)
