"""DockEngine protocol: the plugin contract every docking engine implements.

Adding a new engine = one module in this package exposing `get_engine() -> DockEngine`.
No frontend or API changes needed; capabilities are read from the manifest.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from pydantic import BaseModel, Field

from ..logging import get_logger
from ..models import Box, DockParams, DockResult

log = get_logger("engines.base")


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

    def dock_batch(
        self,
        receptor_pdbqt: Path,
        batch: list[tuple[str, Path]] | list[Path],
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
        **kwargs: Any,
    ) -> dict[str, DockResult]:
        """Run batch docking for multiple ligands against a single receptor.
        Returns a mapping of ligand_id -> DockResult for successfully docked ligands."""
        ...

    def cancel_current(self) -> None:
        """Best-effort termination of a currently running dock() call.
        M1 queue is serial, so at most one run per engine instance is active.
        Default implementations may no-op; cancellation then just waits."""
        ...


class BaseEngine:
    """Base class for docking engines providing common helpers and default fallback
    batch implementation for engines without native batch capabilities."""

    manifest: EngineManifest

    def is_available(self) -> tuple[bool, str | None]:
        return False, "Not implemented"

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        raise NotImplementedError("dock() must be implemented by subclasses")

    def dock_batch(
        self,
        receptor_pdbqt: Path,
        batch: list[tuple[str, Path]] | list[Path],
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
        **kwargs: Any,
    ) -> dict[str, DockResult]:
        """Default sequential fallback batch implementation.

        Iterates over the batch of ligands, allocating an isolated sub-directory
        inside `work_dir` for each ligand, and invokes `self.dock(...)`.
        Handles failures on individual ligands gracefully without failing the entire batch.

        Also supports argument-order polymorphism where batch is passed as the first parameter:
        `dock_batch(batch, receptor, box, params, work_dir)`.
        """
        # Argument-order polymorphism support
        if isinstance(receptor_pdbqt, (list, tuple)) and not isinstance(batch, (list, tuple)):
            raw_batch = receptor_pdbqt
            receptor_path = Path(batch)
        else:
            raw_batch = batch
            receptor_path = Path(receptor_pdbqt)

        # Normalize batch items to list of (ligand_id, Path)
        normalized: list[tuple[str, Path]] = []
        for item in raw_batch:
            if isinstance(item, tuple) and len(item) == 2:
                lig_id, p = item
                normalized.append((str(lig_id), Path(p)))
            elif isinstance(item, (str, Path)):
                p = Path(item)
                normalized.append((p.stem, p))
            elif hasattr(item, "id") and hasattr(item, "pdbqt_path"):
                normalized.append((str(item.id), Path(item.pdbqt_path)))
            else:
                log.warning("Unknown item format in dock_batch: %r", item)

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        results: dict[str, DockResult] = {}
        for lig_id, lig_path in normalized:
            if not lig_path.is_file():
                log.warning("Ligand file not found for %s: %s", lig_id, lig_path)
                continue
            sub_work_dir = work_dir / f"dock_{lig_id}"
            sub_work_dir.mkdir(parents=True, exist_ok=True)
            try:
                res = self.dock(
                    receptor_path,
                    lig_path,
                    box,
                    params,
                    sub_work_dir,
                    flex_pdbqt=flex_pdbqt,
                )
                results[lig_id] = res
            except Exception as exc:
                log.warning("dock_batch: failed to dock ligand %s: %s", lig_id, exc)
                continue

        return results

    def cancel_current(self) -> None:
        """Default cancellation hook; subclasses override with process termination."""
        pass
