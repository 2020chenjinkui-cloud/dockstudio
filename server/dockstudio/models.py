"""Core domain models shared across engines, API and frontend."""

from __future__ import annotations

import os
from typing import Literal

from pydantic import BaseModel, Field


def _default_cpu() -> int:
    # cpu 只影响速度不影响结果：默认吃到 8 线程封顶，避免新手向导动辄 7 分钟的等待
    return max(2, min(os.cpu_count() or 4, 8))


class Box(BaseModel):
    """Docking search box in angstroms."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]


class DockParams(BaseModel):
    """Engine-agnostic docking parameters (subset honored per engine manifest)."""

    exhaustiveness: int = Field(default=8, ge=1, le=64)
    num_modes: int = Field(default=9, ge=1, le=100)
    energy_range: float = Field(default=3.0, gt=0)
    cpu: int = Field(default_factory=_default_cpu, ge=1)
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
    # M2: 柔性变体支持（在受体准备时定义，供 dock 时走 --flex）
    flex_residues: list[str] = Field(default_factory=list)  # 形如 ["A:42"]
    flex_pdbqt_path: str | None = None
    receptor_json_path: str | None = None
    parent_id: str | None = None  # 柔性变体指向其基础受体 id


class PreparedLibrary(BaseModel):
    """配体库（批量筛选用）：一个已逐条准备的配体清单。"""

    id: str
    name: str
    path: str  # 库目录
    ligand_ids: list[str] = Field(default_factory=list)
    failed: list[dict] = Field(default_factory=list)  # [{index, smiles|source, reason}]


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
    batch_id: str | None = None  # M2：批量筛选归属
    env: EnvFingerprint | None = None  # 2.0：任务创建时的软件环境指纹
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
    pose_check: PoseCheck | None = None  # M3，质检结果随分析返回（可空）


# ---------------------------------------------------------------------------
# M2: 口袋预测与批量筛选


class Pocket(BaseModel):
    id: int
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    score: float  # 相对体积评分，越大越可能是口袋
    n_points: int  # 网格点数（体积代理）


class PocketResponse(BaseModel):
    receptor_id: str
    method: str = "builtin_grid"
    pockets: list[Pocket] = Field(default_factory=list)


class BatchCreate(BaseModel):
    name: str | None = None
    engine: str = "vina"
    receptor_id: str
    ligand_ids: list[str] = Field(min_length=1)
    box: Box
    params: DockParams = Field(default_factory=DockParams)


class BatchSummary(BaseModel):
    id: str
    project_id: str
    name: str
    total: int
    done: int = 0  # completed
    failed: int = 0
    running: int = 0
    queued: int = 0
    canceled: int = 0
    receptor_id: str
    engine: str = "vina"
    created_at: datetime
    is_terminal: bool = False


class BatchRow(BaseModel):
    ligand_id: str
    ligand_name: str | None = None
    smiles: str | None = None
    job_id: str
    status: JobStatus
    top_affinity: float | None = None
    error: str | None = None


# ---------------------------------------------------------------------------
# M3: 姿态质检与结果对比


class PoseCheckIssue(BaseModel):
    rule: str  # e.g. "clash" | "bond_length" | "posebusters:..."
    severity: str = "warning"  # info|warning|error
    message: str = ""


class PoseCheck(BaseModel):
    pose_index: int
    ok: bool
    source: str = "posebusters"  # posebusters | builtin
    issues: list[PoseCheckIssue] = Field(default_factory=list)


class ComparisonRow(BaseModel):
    job_id: str
    ligand_id: str
    ligand_name: str | None = None
    smiles: str | None = None
    engine: str
    scoring: str
    exhaustiveness: int
    top_affinity: float | None = None
    n_poses: int = 0
    status: JobStatus


class ComparisonResult(BaseModel):
    rows: list[ComparisonRow] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# M4: Boltz-2 AI 亲和性复核


class AiRescorePoseResult(BaseModel):
    pose_index: int
    log_ic50: float | None = None  # Boltz-2 affinity_pred_value（log10 IC50，μM 基准）
    binder_probability: float | None = None  # affinity_probability_binary 0-1
    rank_by_score: int | None = None


class AiRescoreTask(BaseModel):
    id: str
    job_id: str
    status: JobStatus = JobStatus.queued
    progress: JobProgress = Field(default_factory=JobProgress)
    pose_indices: list[int] = Field(default_factory=list)
    results: list[AiRescorePoseResult] = Field(default_factory=list)
    error: str | None = None
    created_at: datetime
    finished_at: datetime | None = None
    msa_server_used: bool = True  # True=序列经公共 ColabFold MSA 服务器（有提示义务）


class AiRescoreRequest(BaseModel):
    pose_indices: list[int] | None = None  # None=全部姿态
    use_msa_server: bool = True


class AiAvailability(BaseModel):
    available: bool
    reason: str | None = None  # 不可用时的中文说明/安装指引
    mode: str | None = None  # "venv" | "cmd" | None


# ---------------------------------------------------------------------------
# 2.0: 溯源与演示


class CrystalRmsd(BaseModel):
    pose_index: int
    rmsd: float | None = None  # 与共晶配体重原子 RMSD（Å）；无共晶时 null
    reference: str | None = None  # 参考构象来源，如 "MK1:A"


class DemoExample(BaseModel):
    """内置示例一键载入的返回：受体/配体已就绪 + 推荐盒子。"""

    receptor: PreparedReceptor
    ligand: PreparedLigand
    box: Box
    note: str = ""


class EnvFingerprint(BaseModel):
    """记录任务时的软件环境（复现用）。"""

    python: str = ""
    rdkit: str = ""
    meeko: str = ""
    openbabel: str = ""
    platform: str = ""


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
