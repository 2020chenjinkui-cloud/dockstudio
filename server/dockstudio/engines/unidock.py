"""Uni-Dock engine adapter (Linux + NVIDIA GPU; batch capable). CLI mirrors
vina with `--ligands` and `--search_mode fast|balance|detail`; output PDBQT
carries the same `REMARK VINA RESULT` rows as vina."""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
from pathlib import Path

from ..models import Box, DockParams, EngineUnavailableError
from .base import EngineManifest
from .vina import VinaEngine


class UniDockEngine(VinaEngine):
    """Uni-Dock: GPU-batched molecular docking."""

    manifest = EngineManifest(
        name="unidock",
        scoring_functions=["vina", "vinardo", "ad4"],
        supports_flex_residues=False,
        supports_batch_ligands=True,
        requires_gpu=True,
        platforms=["linux"],
        homepage="https://github.com/dptech-corp/Uni-Dock",
    )

    def _binary(self) -> Path:
        if platform.system() != "Linux":
            raise EngineUnavailableError(
                "Uni-Dock 需要 Linux + NVIDIA GPU，Windows 用户请用 WSL2"
                "（指引见 docs/deploy.md#wsl2）。"
            )
        env = os.environ.get("DOCKSTUDIO_UNIDOCK")
        if env and Path(env).is_file():
            return Path(env)
        found = shutil.which("unidock")
        if found:
            return Path(found)
        candidate = (
            Path(__file__).resolve().parents[3] / "resources" / "bin" / "linux-x86_64" / "unidock"
        )
        if candidate.is_file():
            return candidate
        raise EngineUnavailableError(
            "Uni-Dock 二进制未找到：设置 DOCKSTUDIO_UNIDOCK 或放入 resources/bin/linux-x86_64。"
        )

    def _gpu_ok(self) -> tuple[bool, str]:
        if shutil.which("nvidia-smi"):
            try:
                proc = subprocess.run(
                    ["nvidia-smi", "-L"], capture_output=True, timeout=15, check=False
                )
                if proc.returncode == 0 and b"GPU" in (proc.stdout or b""):
                    return True, ""
            except (OSError, subprocess.TimeoutExpired):
                pass
        for candidate in ("libcuda.so", "libcuda.so.1"):
            for libdir in ("/usr/lib/x86_64-linux-gnu", "/usr/lib64", "/usr/local/cuda/lib64"):
                if (Path(libdir) / candidate).exists():
                    return True, ""
        return False, "检测到 Uni-Dock 但无 NVIDIA GPU（找不到 nvidia-smi 或 libcuda）"

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

    def _build_cmd(
        self,
        binary: str,
        box: Box,
        params: DockParams,
        receptor_name: str,
        ligand_name: str,
        out_name: str,
    ) -> list[str]:
        # Uni-Dock: --ligands (batch spelling), --search_mode in place of
        # --exhaustiveness, no --cpu/--energy_range.
        mode = (
            "fast"
            if params.exhaustiveness <= 4
            else ("balance" if params.exhaustiveness <= 16 else "detail")
        )
        return [
            binary,
            "--receptor",
            receptor_name,
            "--ligands",
            ligand_name,
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
            "--search_mode",
            mode,
            "--num_modes",
            str(params.num_modes),
            "--scoring",
            params.scoring,
            "--out",
            out_name,
        ]
