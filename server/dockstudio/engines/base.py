"""DockEngine protocol: the plugin contract every docking engine implements.

Adding a new engine = one module in this package exposing `get_engine() -> DockEngine`.
No frontend or API changes needed; capabilities are read from the manifest.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from pydantic import BaseModel, Field

from ..logging import get_logger
from ..models import Box, DockingError, DockParams, DockResult

log = get_logger("engines.base")

# ---------------------------------------------------------------------------
# Structured degradation markers
# ---------------------------------------------------------------------------
# The 2.0.5 audit found scientific signals being dropped *silently*: an
# ``except Exception: pass`` around the GNINA CNN rescore, ``cnn_score`` turning
# into ``None`` without a warning, and an empty receptor producing an
# ``affinity == 0.0`` "pose" that read as success (docs/qa/gnina-real-machine-
# 2026-11-18.md, GN-03/GN-04/GN-08).  Every such site now emits one of the
# markers below.  They are plain strings on purpose: they travel through
# existing channels (``ScreeningCandidateRow.alerts`` -> ``screening_results.csv``
# "Alerts" column -> CLI table -> ``checkpoint.jsonl`` -> ``DockResult.log``)
# without changing any persisted schema, so old checkpoints stay readable.

#: In-place GNINA rescore could not be run (engine error, missing pose, ...).
GNINA_RESCORE_UNAVAILABLE = "GNINA_RESCORE_UNAVAILABLE"
#: GNINA ran, but the output carried no CNN score / CNN affinity.
GNINA_CNN_MISSING = "GNINA_CNN_MISSING"

#: Prefix used for the machine-greppable lines appended to ``DockResult.log``.
DEGRADED_PREFIX = "DEGRADED"

#: CSV/log cells are single-line; keep reasons short so a run stays readable.
_MAX_DEGRADATION_REASON_CHARS = 240


def marker(name: str, reason: str = "") -> str:
    """Render ``NAME:one-line-reason`` (safe for a CSV cell and for grep)."""
    text = " ".join(str(reason).split())
    if len(text) > _MAX_DEGRADATION_REASON_CHARS:
        text = text[: _MAX_DEGRADATION_REASON_CHARS - 1] + "\u2026"
    return f"{name}:{text}" if text else f"{name}:unknown"


def log_line(name: str, reason: str = "") -> str:
    """Render a ``DEGRADED <NAME>:<reason>`` line for ``DockResult.log``."""
    return f"{DEGRADED_PREFIX} {marker(name, reason)}"


#: Ligand file suffixes every adapter may stage and hand to an engine verbatim.
#: Engines dispatch their parser on the *extension*, so the staged copy must keep
#: it: renaming an exported ``pose_1.sdf`` to ``ligand.pdbqt`` makes gnina's PDBQT
#: reader choke on the SDF header (GN-08).
SUPPORTED_LIGAND_SUFFIXES = (".pdbqt", ".sdf", ".mol2", ".pdb", ".mol")

#: Records that make a receptor file usable. A receptor without any of them can
#: only produce meaningless scores (GN-04).
_RECEPTOR_ATOM_RECORDS = ("ATOM", "HETATM")


@lru_cache(maxsize=64)
def _count_receptor_atoms(path_str: str, mtime_ns: int, size: int) -> int:
    """Count ATOM/HETATM records in a receptor file (memoised on stat identity)."""
    del mtime_ns, size  # part of the cache key only
    try:
        text = Path(path_str).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 0
    return sum(1 for line in text.splitlines() if line.startswith(_RECEPTOR_ATOM_RECORDS))


def validate_receptor_file(receptor_path: Path) -> int:
    """Reject receptor inputs that cannot yield a meaningful score; return atom count.

    A 0-byte (or atom-less) PDBQT makes gnina and vina *succeed* with
    ``affinity == 0.0`` plus random CNN values, which downstream code cannot tell
    apart from a real hit (GN-04, confirmed on a P4).  The failure therefore
    surfaces as a readable Chinese ``DockingError`` instead of a pseudo pose.

    Call sites (accurate as of 2.0.5 - this is *not* a base-class hook, so an
    engine added later must call it explicitly):

    * ``VinaEngine.dock`` - inherited by ``GninaEngine``/``UniDockEngine`` for the
      single-ligand path (``UniDockEngine`` overrides only the batch path);
    * ``UniDockEngine.dock_batch`` - the native GPU batch entry point.

    ``BaseEngine.dock``/``dock_batch`` are Protocols/fallbacks and do **not**
    call it, so a custom engine can still bypass the guard (see
    docs/fixes/gnina-rescore-2026-11-18.md, finding T9-F2).
    """
    try:
        stat = receptor_path.stat()
    except OSError as exc:
        raise DockingError(
            f"受体文件不存在或无法访问：{receptor_path}（{exc.strerror or exc}）"
        ) from exc
    if not receptor_path.is_file():
        raise DockingError(f"受体路径不是文件：{receptor_path}")
    atoms = _count_receptor_atoms(str(receptor_path.resolve()), stat.st_mtime_ns, stat.st_size)
    if atoms == 0:
        raise DockingError(
            f"受体文件不含任何原子记录（ATOM/HETATM）：{receptor_path}（{stat.st_size} 字节）。"
            "空受体无法用于对接，重新准备受体后再试。"
        )
    return atoms


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
