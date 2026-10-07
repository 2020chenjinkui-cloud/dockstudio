"""smina engine adapter: reuses the vina subprocess wrapper with smina-specific
binary resolution (DOCKSTUDIO_SMINA -> PATH -> bundled) and result parsing
(smina writes `REMARK minimizedAffinity`, not `REMARK VINA RESULT`)."""

from __future__ import annotations

import os
import platform
import re
import shutil
from pathlib import Path

from ..models import DockingError, EngineUnavailableError, Pose
from .base import EngineManifest
from .vina import VinaEngine, _key_lines

_AFF_RE = re.compile(r"^\s*REMARK\s+minimizedAffinity\s+(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_TABLE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$")


def parse_smina_out(out_path: Path, stdout: str) -> list[Pose]:
    """Parse smina output: affinities from out.pdbqt REMARK lines, RMSD from the
    stdout score table when available."""
    rmsd_by_idx: dict[int, tuple[float, float]] = {}
    for line in stdout.splitlines():
        m = _TABLE_RE.match(line)
        if m:
            rmsd_by_idx[int(m.group(1))] = (float(m.group(3)), float(m.group(4)))
    poses: list[Pose] = []
    if out_path.is_file():
        for line in out_path.read_text(encoding="utf-8", errors="replace").splitlines():
            m = _AFF_RE.match(line)
            if m:
                idx = len(poses) + 1
                rmsd = rmsd_by_idx.get(idx, (0.0, 0.0))
                poses.append(
                    Pose(
                        index=idx,
                        affinity=float(m.group(1)),
                        rmsd_lb=rmsd[0],
                        rmsd_ub=rmsd[1],
                    )
                )
    if not poses:
        # last resort: score the stdout table alone (out.pdbqt may be empty if
        # the OpenBabel backend failed to write)
        for line in stdout.splitlines():
            m = _TABLE_RE.match(line)
            if m and abs(float(m.group(2))) > 1e-6:
                poses.append(
                    Pose(
                        index=int(m.group(1)),
                        affinity=float(m.group(2)),
                        rmsd_lb=float(m.group(3)),
                        rmsd_ub=float(m.group(4)),
                    )
                )
    if not poses:
        raise DockingError(f"smina output contained no binding poses: {_key_lines(stdout)}")
    return poses


class SminaEngine(VinaEngine):
    """smina (Vina 1.1.2 fork) via bundled binary; OpenBabel plugins live next to it."""

    manifest = EngineManifest(
        name="smina",
        scoring_functions=["vina", "vinardo", "ad4_scoring"],
        supports_flex_residues=True,
        homepage="https://github.com/mwojcikowski/smina",
    )

    def _binary(self) -> Path:
        exe = "smina.exe" if platform.system() == "Windows" else "smina"
        env = os.environ.get("DOCKSTUDIO_SMINA")
        if env:
            p = Path(env)
            if p.is_file():
                return p
        found = shutil.which("smina")
        if found:
            return Path(found)
        root = Path(__file__).resolve().parents[3]
        for name in ("windows", "linux-x86_64", "linux-aarch64", "mac-x86_64", "mac-aarch64"):
            candidate = root / "resources" / "bin" / name / exe
            if candidate.is_file():
                return candidate
        raise EngineUnavailableError(
            "smina 二进制未找到：请设置 DOCKSTUDIO_SMINA 环境变量或将其放入 resources/bin。<br>"
            "smina is not available: set DOCKSTUDIO_SMINA or place the binary under resources/bin"
        )

    @staticmethod
    def _bundled_fallback() -> Path | None:  # unused; resolution done in _binary
        return None

    def _subprocess_env(self) -> dict[str, str]:
        env = dict(os.environ)
        binary_dir = str(self._binary().parent)
        env.setdefault("BABEL_LIBDIR", binary_dir)
        return env

    def _parse_poses(self, out_path: Path, run_log: str) -> list[Pose]:
        return parse_smina_out(out_path, run_log)
