"""DockEngine protocol: the plugin contract every docking engine implements.

Adding a new engine = one module in this package exposing `get_engine() -> DockEngine`.
No frontend or API changes needed; capabilities are read from the manifest.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from pydantic import BaseModel, Field

from ..models import Box, DockParams, DockResult


class EngineManifest(BaseModel):
    """Static capability declaration; the API serves it so the UI can render
    engine-specific parameter forms without hardcoding anything."""

    name: str
    version: str = "unknown"
    scoring_functions: list[str] = Field(default_factory=lambda: ["vina"])
    supports_flex_residues: bool = False
    supports_batch_ligands: bool = False
    requires_gpu: bool = False
    platforms: list[str] = Field(
        default_factory=lambda: ["windows", "linux", "darwin"]
    )  # platform.system().lower() values
    homepage: str = ""


@runtime_checkable
class DockEngine(Protocol):
    """Engine adapters are stateless; each dock() call is self-contained."""

    manifest: EngineManifest

    def is_available(self) -> tuple[bool, str | None]:
        """(ok, reason_if_not). Reason is shown to the user when greyed out."""
        ...

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        """Run one docking job. Must raise DockingError on failure, never leak
        raw subprocess errors. work_dir is created if missing; all intermediate
        files stay inside it."""
        ...

    def cancel_current(self) -> None:
        """Best-effort termination of a currently running dock() call.
        M1 queue is serial, so at most one run per engine instance is active.
        Default implementations may no-op; cancellation then just waits."""
        ...
