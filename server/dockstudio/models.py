"""Core domain models shared across engines, API and frontend."""

from __future__ import annotations

import os
from datetime import datetime
from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, Field


def _default_cpu() -> int:
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

    # --- Orthogonal rescoring / in-place refinement (GNINA --minimize / --score_only) ---
    # When `minimize` is True the engine must locally minimize the *supplied* pose
    # in place instead of running a fresh global search; `score_only` evaluates the
    # supplied pose without any coordinate change.
    minimize: bool = False
    score_only: bool = False


class Pose(BaseModel):
    """One docked binding mode."""

    index: int
    affinity: float  # kcal/mol, lower is better
    rmsd_lb: float | None = None
    rmsd_ub: float | None = None
    file: str | None = None  # path to this pose's SDF/PDBQT (filled after export)
    cnn_score: float | None = None  # GNINA CNN pose score (0~1)
    cnn_affinity: float | None = None  # GNINA CNN predicted affinity (pKd)


class DockResult(BaseModel):
    """Uniform result object every engine adapter must produce."""

    engine: str
    engine_version: str
    params: DockParams
    box: Box
    poses: list[Pose]
    log: str = ""
    run_dir: str


class LigandQCProfile(BaseModel):
    """配体理化性质、天然产物特征与假阳性警示画像。"""

    # Core QC Flags & Alert Names
    has_pains: bool = False
    pains_names: list[str] = Field(default_factory=list)
    has_brenk: bool = False
    brenk_names: list[str] = Field(default_factory=list)
    tcm_alerts: list[str] = Field(default_factory=list)  # e.g. ["Catechol_Polyphenol", "Quinone_Redox"]
    pass_all: bool = True

    # Physicochemical & descriptor properties (defaults for full forward-compatibility)
    heavy_atoms: int = 0
    mw: float = 0.0
    clogp: float = 0.0
    tpsa: float = 0.0
    hbd: int = 0
    hba: int = 0
    rotatable_bonds: int = 0
    fsp3: float = 0.0
    qed: float = 0.0
    inchikey: str = ""
    murcko_scaffold: str = ""
    is_glycoside: bool = False
    parent_glycoside_id: str | None = None
    derivative_type: Literal["original", "aglycone", "stereoisomer"] = "original"
    min_energy_kcal: float | None = None
    pains_alerts: list[str] = Field(default_factory=list)
    tcm_warnings: list[str] = Field(default_factory=list)

    def model_post_init(self, __context: Any) -> None:
        """Harmonize dual naming conventions and evaluate pass_all."""
        if self.pains_names and not self.pains_alerts:
            self.pains_alerts = list(self.pains_names)
        elif self.pains_alerts and not self.pains_names:
            self.pains_names = list(self.pains_alerts)
        if self.pains_names:
            self.has_pains = True

        if self.brenk_names:
            self.has_brenk = True

        if self.tcm_alerts and not self.tcm_warnings:
            self.tcm_warnings = list(self.tcm_alerts)
        elif self.tcm_warnings and not self.tcm_alerts:
            self.tcm_alerts = list(self.tcm_warnings)

        if self.has_pains or self.has_brenk or self.tcm_alerts or self.tcm_warnings:
            self.pass_all = False


class PreparedLigand(BaseModel):
    id: str
    source: Literal["file", "smiles"]
    pdbqt_path: str
    n_rotatable_bonds: int | None = None
    smiles: str | None = None
    name: str | None = None
    qc_profile: LigandQCProfile | None = None
    min_energy_kcal: float | None = None
    derivative_type: Literal["original", "aglycone", "stereoisomer"] = "original"
    parent_id: str | None = None


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
    flex_residues: list[str] = Field(default_factory=list)  # ["A:42"]
    flex_pdbqt_path: str | None = None
    receptor_json_path: str | None = None
    parent_id: str | None = None


class PreparedLibrary(BaseModel):
    """配体库（批量筛选用）：一个已逐条准备的配体清单。"""

    id: str
    name: str
    path: str  # 库目录
    ligand_ids: list[str] = Field(default_factory=list)
    failed: list[dict] = Field(default_factory=list)  # [{index, smiles|source, reason}]
    n_raw: int = 0
    n_desalt_fail: int = 0
    n_dedup: int = 0
    n_chiral_fail: int = 0
    n_prep_fail: int = 0
    n_ready: int = 0
    source_herbs: dict[str, list[str]] = Field(default_factory=dict)
    manifest: dict[str, Any] = Field(default_factory=dict)

    def to_funnel_dict(self) -> dict[str, int]:
        """Convert stage counts into standard funnel dict matching /verify T5."""
        return {
            "raw": self.n_raw,
            "desalt": self.n_desalt_fail,
            "dedup": self.n_dedup,
            "chiral_failed": self.n_chiral_fail,
            "prep_failed": self.n_prep_fail,
            "dock_failed": 0,
            "docked_success": self.n_ready,
        }


# ---------------------------------------------------------------------------
# M1: projects, async jobs, interaction analysis
# ---------------------------------------------------------------------------


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
    interrupted = "interrupted"


class JobProgress(BaseModel):
    stage: str = "queued"
    percent: float | None = None
    message: str = ""


class DockJobRequest(BaseModel):
    engine: str = "vina"
    receptor_id: str
    ligand_id: str
    box: Box
    params: DockParams = Field(default_factory=DockParams)


class Job(BaseModel):
    id: str
    project_id: str
    batch_id: str | None = None
    env: EnvFingerprint | None = None
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
    type: str
    residue: str = ""
    distance: float | None = None
    ligand_atom_ids: list[int] = Field(default_factory=list)
    detail: str = ""


class PoseAnalysis(BaseModel):
    pose_index: int
    source: str = "builtin"
    interactions: list[Interaction] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    pose_check: PoseCheck | None = None


# ---------------------------------------------------------------------------
# M2: 口袋预测与批量筛选
# ---------------------------------------------------------------------------


class Pocket(BaseModel):
    id: int
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    score: float
    n_points: int


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
    done: int = 0
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

    # HTVS & Scoring Rigor extensions
    le: float | None = None
    sile: float | None = None
    lle: float | None = None
    strain_energy: float | None = None
    cnn_score: float | None = None
    cnn_affinity: float | None = None
    derivative_type: str = "original"
    parent_id: str | None = None
    cluster_id: int | None = None
    is_cluster_rep: bool = False
    ecr_score: float | None = None
    alerts: list[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# M3: 姿态质检与结果对比
# ---------------------------------------------------------------------------


class PoseCheckIssue(BaseModel):
    rule: str
    severity: str = "warning"
    message: str = ""


class PoseCheck(BaseModel):
    pose_index: int
    ok: bool
    source: str = "posebusters"
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
# ---------------------------------------------------------------------------


class AiRescorePoseResult(BaseModel):
    pose_index: int
    log_ic50: float | None = None
    binder_probability: float | None = None
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
    msa_server_used: bool = True


class AiRescoreRequest(BaseModel):
    pose_indices: list[int] | None = None
    use_msa_server: bool = True


class AiAvailability(BaseModel):
    available: bool
    reason: str | None = None
    mode: str | None = None


# ---------------------------------------------------------------------------
# 2.0: 溯源与演示
# ---------------------------------------------------------------------------


class CrystalRmsd(BaseModel):
    pose_index: int
    rmsd: float | None = None
    reference: str | None = None


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
# HTVS: Hierarchical Virtual Screening & Rigor Models
# ---------------------------------------------------------------------------


class ScreeningConfig(BaseModel):
    """三级漏斗高通量虚拟筛选配置。"""

    receptor_id: str = ""
    library_id: str = ""
    box: Box | None = None
    tcm_mode: bool = True
    tcm_deglycosylation: bool = True  # 自动衍生中药苷元同屏比对
    enable_tiered_funnel: bool = True  # 启用 Tier1 -> Tier2 -> Tier3 级联漏斗
    tier1_keep_ratio: float = 0.15  # Tier 1 晋级比例 (按 SILE 排序)
    tier1_max_candidates: int = 1000  # Tier 1 最多保留数
    tier2_seeds: list[int] = Field(default_factory=lambda: [42, 1024, 2026])
    tier2_keep_ratio: float = 0.02  # Tier 2 晋级比例
    tier2_max_candidates: int = 100  # Tier 2 最多保留数
    max_strain_kcal: float = 8.0  # 配体构象张力能上限
    max_rmsd_seed: float = 2.0  # 跨种子构象自洽 RMSD 上限 (Å)
    max_sigma_dG: float = 0.35  # 跨种子亲和力标准差上限 (kcal/mol)
    anchor_residues: list[str] = Field(default_factory=list)  # 关键药效团锚定残基，如 ["ASP25:A"]
    enable_gnina_rescore: bool = True  # 启用 GNINA CNN 正交重打分 (GPU 可用时)
    engine: str = "vina"  # 首选引擎
    cpu: int = Field(default_factory=_default_cpu)


class ScreeningCandidateRow(BaseModel):
    """高通量筛选最终多维评估行。"""

    rank: int = 0
    ligand_id: str
    ligand_name: str = ""
    derivative_type: str = "original"  # original | aglycone | stereoisomer
    parent_id: str | None = None
    smiles: str = ""
    affinity_mean: float = 0.0  # kcal/mol (多 seed 均值)
    affinity_std: float = 0.0  # kcal/mol (多 seed 标准差)
    affinity_seed0: float | None = None  # kcal/mol (主种子 seed_0 单次对接值, 供 /verify T1 重算对齐)
    pose_rmsd_seed: float | None = None  # 跨种子构象自洽 RMSD (Å)
    le: float = 0.0  # Ligand Efficiency (-dG / N_heavy)
    sile: float = 0.0  # Size-Independent LE (-dG / N_heavy^0.3)
    lle: float = 0.0  # Lipophilic LE (pKd - cLogP)
    strain_energy: float = 0.0  # 配体内部构象张力能 (kcal/mol)
    affinity_corrected: float | None = None  # ΔG_corr = ΔG_raw + 0.5*max(0, ΔE_strain-3.0)
    cnn_score: float | None = None  # GNINA CNN Pose Score (0~1)
    cnn_affinity: float | None = None  # GNINA CNN pKd
    anchor_match_ratio: float = 1.0  # 关键残基相互作用满足率 (0~1)
    posebusters_ok: bool = True  # 物理几何质检是否通过
    ecr_score: float = 0.0  # 指数共识综合得分 (越高越好)
    ecr_rank: int = 0  # ECR 排名
    cluster_id: int = 0  # Murcko/Butina 骨架簇编号
    is_cluster_rep: bool = False  # 是否为该骨架簇代表分子
    best_pose_path: str | None = None  # 真实最优对接构象文件绝对路径 (交付 poses/ 时使用)
    alerts: list[str] = Field(default_factory=list)  # PAINS / 中药假阳性警示标签


class ScreeningFunnelReport(BaseModel):
    """三级漏斗全链路损耗统计与交付报告。

    Mass-conservation invariants (both must hold, see :attr:`is_conserved`):

    1. Preparation stage:
       ``n_raw - (n_desalt_fail + n_dedup + n_chiral_fail + n_prep_fail) == n_ready``
    2. Docking + funnel stage:
       ``n_ready + n_aglycone_generated - (n_dock_failed + n_funnel_filtered) == n_docked_success``

    ``n_aglycone_generated`` records the extra TCM dual-track aglycone derivatives
    added to the docking pool (they are derivatives, not new library entries, so
    they must be disclosed explicitly instead of silently inflating the pool).
    ``n_funnel_filtered`` records molecules eliminated by the hierarchical funnel
    thresholds (Tier-1/Tier-2 ratio cut-offs, convergence/strain QC rejection) and
    is deliberately kept separate from ``n_dock_failed`` (engine crashes), so that
    "正常淘汰" and "对接程序失败" never share the same accounting bucket.
    """

    n_raw: int = 0
    n_aglycone_generated: int = 0
    n_desalt_fail: int = 0
    n_dedup: int = 0
    n_chiral_fail: int = 0
    n_prep_fail: int = 0
    n_ready: int = 0
    n_tier1_in: int = 0
    n_tier1_out: int = 0
    n_tier2_in: int = 0
    n_tier2_out: int = 0
    n_tier3_in: int = 0
    n_tier3_out: int = 0
    n_tier1_filtered: int = 0
    n_tier2_filtered: int = 0
    n_funnel_filtered: int = 0
    n_dock_failed: int = 0
    n_docked_success: int = 0
    funnel_csv_path: str | None = None
    checkpoint_path: str | None = None
    results_csv_path: str | None = None
    candidates: list[ScreeningCandidateRow] = Field(default_factory=list)

    @property
    def prep_losses(self) -> int:
        return self.n_desalt_fail + self.n_dedup + self.n_chiral_fail + self.n_prep_fail

    @property
    def n_docking_pool(self) -> int:
        """Total ligands actually submitted to docking (incl. derived aglycones)."""
        return self.n_ready + self.n_aglycone_generated

    @property
    def is_conserved(self) -> bool:
        prep_conserved = self.n_raw - self.prep_losses == self.n_ready
        dock_conserved = (
            self.n_docking_pool - (self.n_dock_failed + self.n_funnel_filtered)
            == self.n_docked_success
        )
        return prep_conserved and dock_conserved

    def to_funnel_rows(self) -> list[tuple[str, int]]:
        return [
            ("raw", self.n_raw),
            ("aglycone_generated", self.n_aglycone_generated),
            ("desalt", self.n_desalt_fail),
            ("dedup", self.n_dedup),
            ("chiral_failed", self.n_chiral_fail),
            ("prep_failed", self.n_prep_fail),
            ("tier1_filtered", self.n_tier1_filtered),
            ("tier2_filtered", self.n_tier2_filtered),
            ("funnel_filtered", self.n_funnel_filtered),
            ("dock_failed", self.n_dock_failed),
            ("docked_success", self.n_docked_success),
        ]


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class DockStudioError(Exception):
    """Base error with a user-readable message."""


class EngineUnavailableError(DockStudioError):
    """Engine binary missing or unsupported on this platform."""


class PreparationError(DockStudioError):
    """Receptor/ligand preparation failed."""


class DockingError(DockStudioError):
    """Engine run failed (nonzero exit, unparseable output, timeout)."""
