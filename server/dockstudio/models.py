"""Core domain models shared across engines, API and frontend."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class Box(BaseModel):
    """Docking search box in angstroms."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]


class DockParams(BaseModel):
    """Engine-agnostic docking parameters (subset honored per engine manifest)."""

    exhaustiveness: int = Field(default=8, ge=1, le=64)
    num_modes: int = Field(default=9, ge=1, le=100)
    energy_range: float = Field(default=3.0, gt=0)
    cpu: int = Field(default=1, ge=1)
    seed: int | None = None
    scoring: str = "vina"  # must be in engine manifest scoring_functions


class Pose(BaseModel):
    """One docked binding mode."""

    index: int
    affinity: float  # kcal/mol, lower is better
    rmsd_lb: float | None = None
    rmsd_ub: float | None = None
    file: str | None = None  # path to this pose's SDF/PDBQT (filled after export)


class DockResult(BaseModel):
    """Uniform result object every engine adapter must produce."""

    engine: str
    engine_version: str
    params: DockParams
    box: Box
    poses: list[Pose]
    log: str = ""
    run_dir: str


class PreparedLigand(BaseModel):
    id: str
    source: Literal["file", "smiles"]
    pdbqt_path: str
    n_rotatable_bonds: int | None = None
    smiles: str | None = None
    name: str | None = None


class CocrystalLigand(BaseModel):
    """HETATM group detected in the receptor, usable for autobox."""

    resname: str
    chain: str
    n_atoms: int
    center: tuple[float, float, float]


class PreparedReceptor(BaseModel):
    id: str
    pdbqt_path: str
    original_path: str
    chains: list[str] = Field(default_factory=list)
    cocrystal_ligands: list[CocrystalLigand] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# M1: projects, async jobs, interaction analysis

from datetime import datetime
from enum import Enum


class Project(BaseModel):
    """A project = one directory under data_dir/projects/<id>/ with jobs/ inside."""

    id: str
    name: str
    path: str
    created_at: datetime


class JobStatus(str, Enum):
    queued = "queued"
    running = "running"
    completed = "completed"
    failed = "failed"
    canceled = "canceled"
    interrupted = "interrupted"  # found unfinished after server restart


class JobProgress(BaseModel):
    stage: str = "queued"  # queued|docking|exporting|analyzing
    percent: float | None = None
    message: str = ""


class DockJobRequest(BaseModel):
    """Persisted request snapshot (params + input ids)."""

    engine: str = "vina"
    receptor_id: str
    ligand_id: str
    box: Box
    params: DockParams = Field(default_factory=DockParams)


class Job(BaseModel):
    id: str
    project_id: str
    status: JobStatus = JobStatus.queued
    progress: JobProgress = Field(default_factory=JobProgress)
    request: DockJobRequest
    result: DockResult | None = None
    error: str | None = None
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None

    @property
    def is_terminal(self) -> bool:
        return self.status in {
            JobStatus.completed,
            JobStatus.failed,
            JobStatus.canceled,
            JobStatus.interrupted,
        }


class Interaction(BaseModel):
    """One protein-ligand contact for the result page."""

    type: str  # hbond|hydrophobic|pi_stacking|salt_bridge|halogen|water_bridge
    residue: str = ""  # e.g. "ASP25:A"
    distance: float | None = None
    ligand_atom_ids: list[int] = Field(default_factory=list)
    detail: str = ""


class PoseAnalysis(BaseModel):
    pose_index: int
    source: str = "builtin"  # plip | builtin
    interactions: list[Interaction] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Errors


class DockStudioError(Exception):
    """Base error with a user-readable message."""


class EngineUnavailableError(DockStudioError):
    """Engine binary missing or unsupported on this platform."""


class PreparationError(DockStudioError):
    """Receptor/ligand preparation failed."""


class DockingError(DockStudioError):
    """Engine run failed (nonzero exit, unparseable output, timeout)."""
