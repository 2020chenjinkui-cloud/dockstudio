"""Application configuration: paths, concurrency, engine resolution."""

from __future__ import annotations

import os
import platform
import shutil
from functools import lru_cache
from pathlib import Path

from pydantic import BaseModel


def default_data_dir() -> Path:
    """Per-user data directory for runs and projects."""
    system = platform.system()
    if system == "Windows":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
        return base / "DockStudio"
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "DockStudio"
    return Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share")) / "dockstudio"


class Settings(BaseModel):
    data_dir: Path
    max_workers: int = max(1, (os.cpu_count() or 2) - 1)
    job_timeout_s: int = 3600
    vina_binary: str | None = None  # explicit path override

    @property
    def runs_dir(self) -> Path:
        return self.data_dir / "runs"

    @property
    def prepared_dir(self) -> Path:
        return self.data_dir / "prepared"


def _platform_bin_dirnames() -> list[str]:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "windows":
        return ["windows"]
    if system == "darwin":
        return ["mac-aarch64" if machine in ("arm64", "aarch64") else "mac-x86_64"]
    return ["linux-aarch64" if machine in ("arm64", "aarch64") else "linux-x86_64"]


def bundled_bin_dir() -> Path:
    """repo resources/bin/<platform>/ next to the source tree (works in dev
    layout and in the PyInstaller-frozen layout under _MEIPASS)."""
    import sys

    roots: list[Path] = [Path(__file__).resolve().parents[2]]  # repo root
    if getattr(sys, "frozen", False):
        meipass = Path(getattr(sys, "_MEIPASS", sys.executable))
        # onedir 布局：datas 落在 _MEIPASS(resources/bin) 与其父目录两种都试
        roots.insert(0, meipass)
        roots.insert(0, meipass.parent)
    for root in roots:
        for name in _platform_bin_dirnames():
            candidate = root / "resources" / "bin" / name
            if candidate.is_dir():
                return candidate
    return Path(__file__).resolve().parents[2] / "resources" / "bin"


def resolve_vina_binary(settings: "Settings") -> Path | None:
    """Locate a vina executable: explicit setting -> PATH -> bundled binary."""
    if settings.vina_binary:
        p = Path(settings.vina_binary)
        return p if p.is_file() else None
    found = shutil.which("vina")
    if found:
        return Path(found)
    exe = "vina.exe" if platform.system() == "Windows" else "vina"
    candidate = bundled_bin_dir() / exe
    return candidate if candidate.is_file() else None


@lru_cache
def get_settings() -> Settings:
    env_dir = os.environ.get("DOCKSTUDIO_DATA_DIR")
    env_vina = os.environ.get("DOCKSTUDIO_VINA")
    settings = Settings(
        data_dir=Path(env_dir) if env_dir else default_data_dir(),
        vina_binary=env_vina or None,
    )
    settings.data_dir.mkdir(parents=True, exist_ok=True)
    settings.runs_dir.mkdir(parents=True, exist_ok=True)
    settings.prepared_dir.mkdir(parents=True, exist_ok=True)
    return settings
