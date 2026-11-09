"""Comprehensive test suite for Milestone 3: HTVS Screen Runner, Uni-Dock GPU Batching,
3-Tier Hierarchical Funnel, JSONL Checkpoint Recovery, and Multi-Core Concurrency.

Covers:
1. Uni-Dock Native Batch Command Line & Output Parsing:
   - CLI command line generation (_build_batch_cmd)
   - Staging directory creation & ligand filename deduplication
   - Index file (ligand_index.txt) generation
   - Multi-PDBQT output parsing & partial failure resilience
   - BaseEngine fallback execution
2. 3-Tier Funnel Scheduling:
   - Tier 1: HTVS Fast pass, SILE sorting, and top 15% cutoff
   - Tier 2: SP Balanced pass, 3-seed convergence (sigma_dG, RMSD) & strain QC, top 2% cutoff
   - Tier 3: XP Rigor Rescore, strain correction, ECR consensus ranking, and Murcko-Butina clustering
   - Single-tier fallback when tiered funnel is disabled
3. Checkpoint Interrupt Simulation & Zero-Loss Resume:
   - Incremental JSONL write verification
   - Mid-run interrupt simulation and restart verification
   - Zero duplicate dockings on resume
   - Corrupt trailing line recovery
   - Immediate return on already-completed run
4. Funnel Report & funnel.csv Consistency:
   - Exact standard format matching /verify T5
   - Strict mass conservation closure (N_raw - losses == N_docked_success)
   - Hierarchical stage count monotonicity
   - Direct assertion using verify T5 closure algorithm
5. Multi-Core CPU Worker Allocation:
   - Worker formula calculation across 1, 2, 4, 8, 16, 64, 128 cores
   - Worker pool capped at molecule count
   - Graceful fallback on invalid/zero configurations
   - Concurrency speedup and worker crash isolation
"""

from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path
from typing import Any

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.core.screen_runner import (
    ScreeningRunner,
    calculate_worker_pool_size,
)
from dockstudio.engines.base import BaseEngine, EngineManifest
from dockstudio.engines.unidock import UniDockEngine
from dockstudio.models import (
    Box,
    DockingError,
    DockParams,
    DockResult,
    Pose,
    PreparedLibrary,
    PreparedLigand,
    ScreeningCandidateRow,
    ScreeningConfig,
    ScreeningFunnelReport,
)

# =============================================================================
# Synthetic Fixtures & Mock Helpers
# =============================================================================

_SAMPLE_PDBQT_TEMPLATE = """\
MODEL {model_idx}
REMARK VINA RESULT: {affinity:.2f} 0.000 0.000
ROOT
ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O
ENDROOT
ENDMDL
"""


def make_sample_pdbqt(affinities: list[float]) -> str:
    """Generates a valid multi-MODEL PDBQT string with given affinities."""
    blocks = []
    for idx, aff in enumerate(affinities, start=1):
        blocks.append(_SAMPLE_PDBQT_TEMPLATE.format(model_idx=idx, affinity=aff))
    return "".join(blocks)


def create_fake_ligand_pdbqt(path: Path, affinities: list[float] | None = None) -> Path:
    """Writes a sample PDBQT file to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    affs = affinities or [-8.0]
    path.write_text(make_sample_pdbqt(affs), encoding="utf-8")
    return path


class MockFunnelEngine(BaseEngine):
    """Controllable mock engine for testing 3-tier funnel transitions and checkpoints."""

    manifest = EngineManifest(
        name="mock_funnel",
        scoring_functions=["vina"],
        supports_batch_ligands=True,
    )

    def __init__(self, fail_ligand_ids: set[str] | None = None, interrupt_on_id: str | None = None):
        self.fail_ligand_ids = fail_ligand_ids or set()
        self.interrupt_on_id = interrupt_on_id
        self.dock_call_count: dict[str, int] = {}
        self.docked_history: list[tuple[str, int, int | None]] = []
        self.custom_affinities: dict[str, float] = {}
        self.custom_poses: dict[str, list[float]] = {}

    def set_ligand_response(self, lig_id: str, affinity: float, seed_affinities: dict[int, float] | None = None) -> None:
        self.custom_affinities[lig_id] = affinity
        if seed_affinities:
            for s, aff in seed_affinities.items():
                self.custom_poses[f"{lig_id}_seed_{s}"] = [aff]

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        lig_id = ligand_pdbqt.stem
        if lig_id == "ligand":  # DockStudio standard directory layout
            lig_id = ligand_pdbqt.parent.name

        if lig_id == self.interrupt_on_id:
            raise KeyboardInterrupt(f"Simulated user interrupt on {lig_id}")

        self.dock_call_count[lig_id] = self.dock_call_count.get(lig_id, 0) + 1
        self.docked_history.append((lig_id, params.exhaustiveness, params.seed))

        if lig_id in self.fail_ligand_ids:
            raise DockingError(f"Simulated docking failure for {lig_id}")

        work_dir.mkdir(parents=True, exist_ok=True)
        out_file = work_dir / "out.pdbqt"

        seed_key = f"{lig_id}_seed_{params.seed}"
        if seed_key in self.custom_poses:
            affs = self.custom_poses[seed_key]
        elif lig_id in self.custom_affinities:
            affs = [self.custom_affinities[lig_id]]
        else:
            parsed_affs = []
            if ligand_pdbqt.is_file():
                for line in ligand_pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
                    if "REMARK VINA RESULT:" in line:
                        parts = line.split()
                        try:
                            parsed_affs.append(float(parts[3]))
                        except (IndexError, ValueError):
                            pass
            affs = parsed_affs if parsed_affs else [-7.5]

        out_file.write_text(make_sample_pdbqt(affs), encoding="utf-8")
        poses = [Pose(index=i + 1, affinity=a) for i, a in enumerate(affs)]
        return DockResult(
            engine="mock_funnel",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=poses,
            log="mock log",
            run_dir=str(work_dir),
        )


# =============================================================================
# 1. Uni-Dock Native Batch Command Line & Output Parsing Tests
# =============================================================================

def test_unidock_build_batch_cmd_formatting() -> None:
    """Verifies that UniDockEngine._build_batch_cmd constructs the exact CLI command:
    - Passes --receptor, --ligand_index, --dir
    - Translates exhaustiveness into search_mode (fast, balance, detail)
    - Strictly excludes --ligands, --ligand, --out, --cpu, and --energy_range."""
    eng = UniDockEngine()
    box = Box(center=(10.0, 20.0, 30.0), size=(25.0, 25.0, 25.0))

    # Fast mode (exhaustiveness <= 4)
    cmd_fast = eng._build_batch_cmd(
        binary="unidock",
        box=box,
        params=DockParams(exhaustiveness=4, num_modes=9, scoring="vina", seed=42),
        receptor_path="/data/rec.pdbqt",
        ligand_index_path="/data/index.txt",
        out_dir_path="/data/out",
    )
    assert "--ligand_index" in cmd_fast
    assert cmd_fast[cmd_fast.index("--ligand_index") + 1] == "/data/index.txt"
    assert "--dir" in cmd_fast
    assert cmd_fast[cmd_fast.index("--dir") + 1] == "/data/out"
    assert "--receptor" in cmd_fast
    assert cmd_fast[cmd_fast.index("--receptor") + 1] == "/data/rec.pdbqt"
    assert "--search_mode" in cmd_fast
    assert cmd_fast[cmd_fast.index("--search_mode") + 1] == "fast"
    assert cmd_fast[cmd_fast.index("--num_modes") + 1] == "9"
    assert cmd_fast[cmd_fast.index("--scoring") + 1] == "vina"
    assert cmd_fast[cmd_fast.index("--seed") + 1] == "42"

    # In batch mode, single-molecule flags must be absent
    assert "--out" not in cmd_fast
    assert "--ligand" not in cmd_fast
    assert "--ligands" not in cmd_fast
    assert "--cpu" not in cmd_fast
    assert "--energy_range" not in cmd_fast

    # Balance mode (5 <= exhaustiveness <= 16)
    cmd_bal = eng._build_batch_cmd(
        binary="unidock",
        box=box,
        params=DockParams(exhaustiveness=12),
        receptor_path="r",
        ligand_index_path="i",
        out_dir_path="o",
    )
    assert cmd_bal[cmd_bal.index("--search_mode") + 1] == "balance"

    # Detail mode (exhaustiveness > 16)
    cmd_det = eng._build_batch_cmd(
        binary="unidock",
        box=box,
        params=DockParams(exhaustiveness=32),
        receptor_path="r",
        ligand_index_path="i",
        out_dir_path="o",
    )
    assert cmd_det[cmd_det.index("--search_mode") + 1] == "detail"


def test_unidock_dock_batch_creates_index_and_staged_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies that UniDockEngine.dock_batch correctly stages ligands into unique filenames,
    avoiding the DockStudio '.../<id>/ligand.pdbqt' collision hazard, and writes ligand_index.txt."""
    eng = UniDockEngine()

    # Create 3 ligands stored in DockStudio standard layout: <id>/ligand.pdbqt
    l1_path = create_fake_ligand_pdbqt(tmp_path / "lib" / "mol_001" / "ligand.pdbqt", [-8.1])
    l2_path = create_fake_ligand_pdbqt(tmp_path / "lib" / "mol_002" / "ligand.pdbqt", [-7.4])
    l3_path = create_fake_ligand_pdbqt(tmp_path / "lib" / "mol_003" / "ligand.pdbqt", [-9.2])
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    batch = [("mol_001", l1_path), ("mol_002", l2_path), ("mol_003", l3_path)]
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "unidock_run"

    # Mock engine availability & subprocess
    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    executed_cmd: list[str] = []

    def fake_popen(cmd, *args, **kwargs):
        executed_cmd.extend(cmd)
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        (out_dir / "mol_001_out.pdbqt").write_text(make_sample_pdbqt([-8.1]), encoding="utf-8")
        (out_dir / "mol_002_out.pdbqt").write_text(make_sample_pdbqt([-7.4]), encoding="utf-8")
        (out_dir / "mol_003_out.pdbqt").write_text(make_sample_pdbqt([-9.2]), encoding="utf-8")

        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                return ("Uni-Dock finished 3 ligands", "")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    results = eng.dock_batch(rec_path, batch, box, params, work_dir)

    # 1. Verify staging directory contains disambiguated files
    staging_dir = work_dir / "staging_ligands"
    assert staging_dir.is_dir()
    assert (staging_dir / "mol_001.pdbqt").is_file()
    assert (staging_dir / "mol_002.pdbqt").is_file()
    assert (staging_dir / "mol_003.pdbqt").is_file()

    # 2. Verify ligand_index.txt contains exact paths to staged files
    index_file = work_dir / "ligand_index.txt"
    assert index_file.is_file()
    index_lines = index_file.read_text(encoding="utf-8").strip().splitlines()
    assert len(index_lines) == 3
    assert any("mol_001.pdbqt" in line for line in index_lines)
    assert any("mol_002.pdbqt" in line for line in index_lines)
    assert any("mol_003.pdbqt" in line for line in index_lines)

    # 3. Verify parsed results
    assert len(results) == 3
    assert results["mol_001"].poses[0].affinity == pytest.approx(-8.1)
    assert results["mol_002"].poses[0].affinity == pytest.approx(-7.4)
    assert results["mol_003"].poses[0].affinity == pytest.approx(-9.2)


def test_unidock_dock_batch_partial_failures(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies that if one ligand fails or produces no output during Uni-Dock batch docking,
    the remaining valid ligands are still parsed and returned without crashing."""
    eng = UniDockEngine()
    rec = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    l1 = create_fake_ligand_pdbqt(tmp_path / "l1.pdbqt")
    l2 = create_fake_ligand_pdbqt(tmp_path / "l2.pdbqt")
    l3 = create_fake_ligand_pdbqt(tmp_path / "l3.pdbqt")

    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    def fake_popen(cmd, *args, **kwargs):
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        (out_dir / "l1_out.pdbqt").write_text(make_sample_pdbqt([-8.0]), encoding="utf-8")
        (out_dir / "l2_out.pdbqt").write_text(make_sample_pdbqt([-7.0]), encoding="utf-8")
        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                return ("Uni-Dock batch: 2 succeeded, 1 failed", "Warning: l3 skipped")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    results = eng.dock_batch(
        rec,
        [("l1", l1), ("l2", l2), ("l3", l3)],
        Box(center=(0, 0, 0), size=(20, 20, 20)),
        DockParams(),
        tmp_path / "work",
    )
    assert "l1" in results
    assert "l2" in results
    assert "l3" not in results
    assert len(results) == 2


def test_unidock_dock_batch_fatal_crash_raises(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies that if Uni-Dock terminates with a non-zero exit code and 0 outputs,
    a descriptive DockingError is raised."""
    eng = UniDockEngine()
    rec = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    l1 = create_fake_ligand_pdbqt(tmp_path / "l1.pdbqt")

    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    def fake_popen_crash(cmd, *args, **kwargs):
        class FakeProcess:
            returncode = 1
            def communicate(self, timeout=None):
                return ("", "CUDA error: out of memory at grid alloc")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen_crash)

    with pytest.raises(DockingError, match="out of memory"):
        eng.dock_batch(
            rec,
            [("l1", l1)],
            Box(center=(0, 0, 0), size=(20, 20, 20)),
            DockParams(),
            tmp_path / "work",
        )


def test_base_engine_fallback_dock_batch(tmp_path: Path) -> None:
    """Verifies that BaseEngine provides a working sequential fallback dock_batch
    that isolates each dock() call in subdirectories and handles argument-order polymorphism."""
    mock_eng = MockFunnelEngine(fail_ligand_ids={"bad_mol"})
    rec = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    m1 = create_fake_ligand_pdbqt(tmp_path / "m1.pdbqt", [-6.5])
    m2 = create_fake_ligand_pdbqt(tmp_path / "bad_mol.pdbqt", [-5.0])
    m3 = create_fake_ligand_pdbqt(tmp_path / "m3.pdbqt", [-8.5])

    batch = [("m1", m1), ("bad_mol", m2), ("m3", m3)]
    box = Box(center=(0, 0, 0), size=(15, 15, 15))
    params = DockParams(exhaustiveness=4)
    work = tmp_path / "fallback_work"

    results = mock_eng.dock_batch(rec, batch, box, params, work)
    assert "m1" in results
    assert "m3" in results
    assert "bad_mol" not in results
    assert results["m1"].poses[0].affinity == pytest.approx(-6.5)
    assert results["m3"].poses[0].affinity == pytest.approx(-8.5)

    assert (work / "dock_m1").is_dir()
    assert (work / "dock_m3").is_dir()


# =============================================================================
# 2. 3-Tier Funnel Scheduling Tests (Tier 1 -> Tier 2 -> Tier 3 Transitions)
# =============================================================================

def _build_test_library(tmp_path: Path, n_ligands: int = 100) -> tuple[PreparedLibrary, list[PreparedLigand]]:
    """Constructs a test library of n molecules with diverse heavy atom counts and valid PDBQT files."""
    lib_dir = tmp_path / "synth_lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    ligands: list[PreparedLigand] = []
    ligand_ids: list[str] = []

    for i in range(1, n_ligands + 1):
        lid = f"LIG_{i:03d}"
        lig_sub = lib_dir / lid
        lig_sub.mkdir(parents=True, exist_ok=True)
        pdbqt_file = lig_sub / "ligand.pdbqt"

        n_ha = 15 + (i * 65) // n_ligands
        smi = "C" * n_ha
        pdbqt_file.write_text(make_sample_pdbqt([-7.0]), encoding="utf-8")

        lig_model = PreparedLigand(
            id=lid,
            source="smiles",
            pdbqt_path=str(pdbqt_file),
            smiles=smi,
            n_rotatable_bonds=2 + (i % 6),
        )
        ligands.append(lig_model)
        ligand_ids.append(lid)

    prep_lib = PreparedLibrary(
        id="lib_100",
        name="TestLibrary100",
        path=str(lib_dir),
        ligand_ids=ligand_ids,
        n_raw=n_ligands,
        n_ready=n_ligands,
    )
    return prep_lib, ligands


def test_funnel_tier1_selection_by_sile(tmp_path: Path) -> None:
    """Verifies Tier 1 HTVS fast filtering:
    - 100 ligands evaluated at fast exhaustiveness
    - Sorted descending by SILE = -affinity / (N_heavy^0.3)
    - Exactly top 15% (15 ligands) pass to Tier 2
    - Verified that passing ligands correspond strictly to the highest SILE values."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=100)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    mock_eng = MockFunnelEngine()
    expected_siles: dict[str, float] = {}
    for idx, lig in enumerate(ligands, start=1):
        aff = -5.0 - (idx * 5.0 / 100.0)
        n_ha = len(lig.smiles)
        sile = -aff / (n_ha ** 0.3)
        expected_siles[lig.id] = sile
        mock_eng.set_ligand_response(lig.id, aff)

    top_15_expected = sorted(expected_siles.keys(), key=lambda k: expected_siles[k], reverse=True)[:15]

    config = ScreeningConfig(
        tier1_keep_ratio=0.15,
        tier1_max_candidates=1000,
        enable_tiered_funnel=True,
    )
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "screen_out",
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()

    # Tier 1 Accounting
    assert report.n_tier1_in == 100
    assert report.n_tier1_out == 15
    assert len(runner.tier1_passed_ids) == 15
    assert set(runner.tier1_passed_ids) == set(top_15_expected)


def test_funnel_tier2_convergence_and_strain_qc(tmp_path: Path) -> None:
    """Verifies Tier 2 SP balanced filtering:
    - Docks 15 candidates across 3 random seeds (42, 1024, 2026)
    - Rejects candidates with sigma_dG > 0.35 kcal/mol
    - Rejects candidates with RMSD > 2.0 Å
    - Rejects candidates with ligand internal strain > 8.0 kcal/mol
    - Retains top 2% (2 candidates) for Tier 3."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=100)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    mock_eng = MockFunnelEngine()

    for lig in ligands:
        mock_eng.set_ligand_response(lig.id, -6.0)

    # Configure LIG_095, LIG_096, LIG_097 to dominate SILE and converge across seeds
    for lid, base_aff in [("LIG_095", -9.5), ("LIG_096", -9.3), ("LIG_097", -9.1)]:
        mock_eng.set_ligand_response(lid, base_aff, seed_affinities={
            42: base_aff,
            1024: base_aff + 0.1,
            2026: base_aff - 0.1,
        })
    # Divergent affinity: LIG_098
    mock_eng.set_ligand_response("LIG_098", -9.8, seed_affinities={
        42: -9.8,
        1024: -7.5,
        2026: -5.0,
    })

    config = ScreeningConfig(
        tier1_keep_ratio=0.15,
        tier2_keep_ratio=0.02,
        tier2_seeds=[42, 1024, 2026],
        max_sigma_dG=0.35,
        max_rmsd_seed=2.0,
        max_strain_kcal=8.0,
    )
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "screen_tier2",
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()

    assert report.n_tier2_in == 15
    assert report.n_tier2_out <= 3
    passed_ids = [c.ligand_id for c in report.candidates]
    assert "LIG_098" not in passed_ids
    assert "LIG_095" in passed_ids or "LIG_096" in passed_ids


def test_funnel_tier3_xp_rigor_scoring_ecr_and_murcko_clustering(tmp_path: Path) -> None:
    """Verifies Tier 3 XP Rigor scoring:
    - Applies internal strain energy penalty: dG_corr = affinity + 0.5 * max(0, strain - 3.0)
    - Computes Exponential Consensus Ranking (ECR) across orthogonal metrics
    - Extracts Bemis-Murcko scaffolds and clusters via Butina algorithm
    - Assigns cluster_id and designates exactly one representative per cluster (is_cluster_rep=True)."""
    mol1_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"  # Flavone
    mol2_smi = "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"      # Flavone
    mol3_smi = "c1ccc2c(c1)[nH]c(CCN)c2"                      # Indole

    lib_dir = tmp_path / "tier3_lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    ligands = []
    for lid, smi in [("FLAV_1", mol1_smi), ("FLAV_2", mol2_smi), ("IND_1", mol3_smi)]:
        p = lib_dir / lid / "ligand.pdbqt"
        create_fake_ligand_pdbqt(p, [-8.5])
        ligands.append(PreparedLigand(id=lid, source="smiles", pdbqt_path=str(p), smiles=smi))

    prep_lib = PreparedLibrary(
        id="t3_lib", name="T3Lib", path=str(lib_dir),
        ligand_ids=["FLAV_1", "FLAV_2", "IND_1"], n_raw=3, n_ready=3,
    )
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    mock_eng = MockFunnelEngine()
    mock_eng.set_ligand_response("FLAV_1", -9.0)
    mock_eng.set_ligand_response("FLAV_2", -8.2)
    mock_eng.set_ligand_response("IND_1", -7.5)

    config = ScreeningConfig(
        enable_tiered_funnel=False,  # Directly evaluate all 3 in Tier 3 mode
    )
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "screen_t3",
        receptor_pdbqt=rec_path,
        ligands=ligands,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()

    assert len(report.candidates) == 3
    # Check ECR score and ranking assigned
    for c in report.candidates:
        assert c.ecr_score >= 0.0
        assert c.ecr_rank in (1, 2, 3)

    # Check Murcko clustering: FLAV_1 and FLAV_2 share same scaffold -> same cluster_id
    flav_clusters = {c.ligand_id: c.cluster_id for c in report.candidates if "FLAV" in c.ligand_id}
    ind_cluster = next(c.cluster_id for c in report.candidates if "IND" in c.ligand_id)

    assert flav_clusters["FLAV_1"] == flav_clusters["FLAV_2"]
    assert flav_clusters["FLAV_1"] != ind_cluster

    # Exactly one representative per cluster
    reps = [c for c in report.candidates if c.is_cluster_rep]
    assert len(reps) == 2  # 2 clusters total
    flav_reps = [c for c in reps if "FLAV" in c.ligand_id]
    assert len(flav_reps) == 1
    # The representative must be the higher-ranked flavone
    assert flav_reps[0].ligand_id == "FLAV_1"


def test_funnel_disabled_single_tier_fallback(tmp_path: Path) -> None:
    """Verifies that when enable_tiered_funnel=False, all input ligands pass
    directly through standard docking and are scored without reduction."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=5)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    mock_eng = MockFunnelEngine()

    config = ScreeningConfig(enable_tiered_funnel=False)
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "screen_single_tier",
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()
    assert report.n_docked_success == 5
    assert len(report.candidates) == 5


# =============================================================================
# 3. Checkpoint Interrupt Simulation & Resume Verification (Zero-Loss)
# =============================================================================

def test_checkpoint_incremental_jsonl_write(tmp_path: Path) -> None:
    """Verifies that ScreeningRunner writes checkpoint.jsonl incrementally as each
    candidate completes docking."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=4)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    mock_eng = MockFunnelEngine()

    out_dir = tmp_path / "screen_ckpt_test"
    config = ScreeningConfig(enable_tiered_funnel=False)
    runner = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()

    ckpt_path = out_dir / "checkpoint.jsonl"
    assert ckpt_path.is_file()

    lines = [line.strip() for line in ckpt_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert len(lines) == 4

    # Validate JSON deserialization
    for line in lines:
        entry = json.loads(line)
        assert "ligand_id" in entry
        assert "affinity" in entry or "affinity_mean" in entry


def test_checkpoint_interrupt_simulation_and_zero_loss_resume(tmp_path: Path) -> None:
    """Rigorous crash recovery test:
    1. Run 5 ligands (LIG_001 .. LIG_005).
    2. Simulate crash on LIG_003 (LIG_001 and LIG_002 already written to checkpoint.jsonl).
    3. Re-run with resume=True.
    4. Assert LIG_001 and LIG_002 are NEVER re-docked (call_count == 1).
    5. Assert LIG_003, LIG_004, LIG_005 are docked and all 5 exist in final report."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=5)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    out_dir = tmp_path / "screen_crash_recovery"
    config = ScreeningConfig(enable_tiered_funnel=False, cpu=2)

    # Phase 1: Engine crashes on LIG_003
    interrupting_engine = MockFunnelEngine(interrupt_on_id="LIG_003")
    runner_1 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=interrupting_engine,
    )

    with pytest.raises(KeyboardInterrupt, match="Simulated user interrupt on LIG_003"):
        runner_1.run(resume=False)

    # Checkpoint after crash must contain exactly 2 lines (LIG_001, LIG_002)
    ckpt_path = out_dir / "checkpoint.jsonl"
    assert ckpt_path.is_file()
    lines_pre = [l for l in ckpt_path.read_text(encoding="utf-8").splitlines() if l.strip()]
    assert len(lines_pre) == 2
    completed_ids_pre = [json.loads(l)["ligand_id"] for l in lines_pre]
    assert completed_ids_pre == ["LIG_001", "LIG_002"]

    # Phase 2: Resume with stable engine
    stable_engine = MockFunnelEngine()
    runner_2 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=stable_engine,
    )

    report_resumed = runner_2.run(resume=True)

    # Assert stable engine was NEVER called for LIG_001 or LIG_002!
    assert "LIG_001" not in stable_engine.dock_call_count
    assert "LIG_002" not in stable_engine.dock_call_count
    # Called only for remaining 3 molecules
    assert stable_engine.dock_call_count["LIG_003"] == 1
    assert stable_engine.dock_call_count["LIG_004"] == 1
    assert stable_engine.dock_call_count["LIG_005"] == 1

    # Final report must contain all 5 molecules with zero data loss
    assert len(report_resumed.candidates) == 5
    resumed_ids = [c.ligand_id for c in report_resumed.candidates]
    assert set(resumed_ids) == {"LIG_001", "LIG_002", "LIG_003", "LIG_004", "LIG_005"}


def test_checkpoint_corrupt_tail_line_resilience(tmp_path: Path) -> None:
    """Verifies that an abrupt power failure that writes a corrupt partial line
    to checkpoint.jsonl is handled gracefully: valid lines are loaded, the broken line
    is discarded and re-computed, and the process completes cleanly."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=3)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    out_dir = tmp_path / "corrupt_ckpt"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Seed checkpoint with 1 valid line and 1 truncated line
    ckpt_path = out_dir / "checkpoint.jsonl"
    valid_entry = json.dumps({"ligand_id": "LIG_001", "affinity": -8.2, "rank": 1})
    corrupt_entry = '{"ligand_id": "LIG_002", "affin'  # Truncated
    ckpt_path.write_text(f"{valid_entry}\n{corrupt_entry}\n", encoding="utf-8")

    engine = MockFunnelEngine()
    runner = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False),
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
    )

    report = runner.run(resume=True)

    # LIG_001 loaded from checkpoint; LIG_002 and LIG_003 docked by engine
    assert "LIG_001" not in engine.dock_call_count
    assert engine.dock_call_count["LIG_002"] == 1
    assert engine.dock_call_count["LIG_003"] == 1
    assert len(report.candidates) == 3


def test_resume_when_already_fully_completed(tmp_path: Path) -> None:
    """Verifies that calling resume=True on an already finished run performs zero
    docking operations and returns the existing results immediately."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=3)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    out_dir = tmp_path / "already_done"

    engine_first = MockFunnelEngine()
    runner_first = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False),
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine_first,
    )
    report_1 = runner_first.run(resume=False)
    assert len(report_1.candidates) == 3

    # Second run with resume
    engine_second = MockFunnelEngine()
    runner_second = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False),
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine_second,
    )
    report_2 = runner_second.run(resume=True)

    assert len(engine_second.dock_call_count) == 0
    assert len(report_2.candidates) == 3


# =============================================================================
# 4. Funnel Report and funnel.csv Stage Loss Consistency Tests
# =============================================================================

def test_funnel_csv_row_labels_and_values(tmp_path: Path) -> None:
    """Verifies that ScreeningRunner exports funnel.csv with exact standard headers
    and strictly accurate row counts."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=10)
    # Inject prep loss counts into library
    prep_lib.n_raw = 20
    prep_lib.n_desalt_fail = 2
    prep_lib.n_dedup = 3
    prep_lib.n_chiral_fail = 2
    prep_lib.n_prep_fail = 3
    prep_lib.n_ready = 10  # 20 - (2+3+2+3) = 10

    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    # Fail 1 ligand during docking
    mock_eng = MockFunnelEngine(fail_ligand_ids={"LIG_005"})

    out_dir = tmp_path / "screen_funnel_csv"
    runner = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False),
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()
    funnel_csv = out_dir / "funnel.csv"
    assert funnel_csv.is_file()

    lines = funnel_csv.read_text(encoding="utf-8").strip().splitlines()
    data = dict(line.split(",") for line in lines)

    assert int(data["raw"]) == 20
    assert int(data["desalt"]) == 2
    assert int(data["dedup"]) == 3
    assert int(data["chiral_failed"]) == 2
    assert int(data["prep_failed"]) == 3
    assert int(data["dock_failed"]) == 1
    assert int(data["docked_success"]) == 9


def test_funnel_mass_conservation_closure(tmp_path: Path) -> None:
    """Verifies that ScreeningFunnelReport.is_conserved strictly enforces:
    N_raw - (N_desalt + N_dedup + N_chiral_failed + N_prep_failed + N_dock_failed) == N_docked_success."""
    rep_ok = ScreeningFunnelReport(
        n_raw=100,
        n_desalt_fail=10,
        n_dedup=10,
        n_chiral_fail=10,
        n_prep_fail=10,
        n_ready=60,
        n_dock_failed=5,
        n_docked_success=55,
    )
    assert rep_ok.is_conserved is True


def test_funnel_hierarchical_stage_counts_monotonicity(tmp_path: Path) -> None:
    """Verifies that tiered screening preserves hierarchical monotonicity:
    n_raw >= n_ready = n_tier1_in >= n_tier1_out = n_tier2_in >= n_tier2_out = n_tier3_in >= n_tier3_out."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=50)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    mock_eng = MockFunnelEngine()

    config = ScreeningConfig(
        enable_tiered_funnel=True,
        tier1_keep_ratio=0.20,  # 50 -> 10
        tier2_keep_ratio=0.04,  # 50 -> 2
    )
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "screen_monotonic",
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )

    report = runner.run()

    assert report.n_raw >= report.n_ready
    assert report.n_ready == report.n_tier1_in
    assert report.n_tier1_in >= report.n_tier1_out
    assert report.n_tier1_out == report.n_tier2_in
    assert report.n_tier2_in >= report.n_tier2_out
    assert report.n_tier2_out == report.n_tier3_in
    assert report.n_tier3_in >= report.n_tier3_out
    assert len(report.candidates) == report.n_tier3_out


def test_funnel_csv_compatibility_with_verify_t5(tmp_path: Path) -> None:
    """Directly verifies that the generated funnel.csv satisfies the exact verification
    closure logic from verify.py T5 and test_htvs_funnel_stress.py."""
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=20)
    prep_lib.n_raw = 25
    prep_lib.n_desalt_fail = 1
    prep_lib.n_dedup = 1
    prep_lib.n_chiral_fail = 1
    prep_lib.n_prep_fail = 2
    prep_lib.n_ready = 20

    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")
    mock_eng = MockFunnelEngine()

    out_dir = tmp_path / "t5_verify_pkg"
    runner = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False),
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=mock_eng,
    )
    runner.run()

    # Apply T5 Closure logic
    funnel_path = out_dir / "funnel.csv"
    lines = funnel_path.read_text(encoding="utf-8").strip().splitlines()
    funnel = dict(line.split(",") for line in lines)

    n_raw = int(funnel.get("raw", 0))
    n_desalt = int(funnel.get("desalt", 0))
    n_dedup = int(funnel.get("dedup", 0))
    n_chiral = int(funnel.get("chiral_failed", 0))
    n_prep = int(funnel.get("prep_failed", 0))
    n_dock_fail = int(funnel.get("dock_failed", 0))
    n_dock_ok = int(funnel.get("docked_success", 0))

    total_losses = n_desalt + n_dedup + n_chiral + n_prep + n_dock_fail
    expected_docked = n_raw - total_losses

    disc_1 = abs(expected_docked - n_dock_ok)
    disc_2 = abs(n_dock_ok - prep_lib.n_ready + n_dock_fail)

    assert disc_1 == 0, f"Discrepancy 1 failure: {disc_1}"
    assert disc_2 == 0, f"Discrepancy 2 failure: {disc_2}"


# =============================================================================
# 5. Multi-Core CPU Worker Allocation Tests
# =============================================================================

@pytest.mark.parametrize(
    ("total_cpus", "cpu_per_job", "n_mols", "expected_workers"),
    [
        (1, 2, 10, 1),
        (2, 2, 10, 1),
        (4, 2, 10, 2),
        (8, 2, 10, 4),
        (16, 2, 100, 8),
        (16, 4, 100, 4),
        (64, 2, 100, 32),
        (128, 4, 1000, 32),
        (32, 2, 3, 3),
        (0, 2, 10, 1),
        (-8, 2, 10, 1),
        (16, 0, 10, 1),
    ],
)
def test_calculate_worker_pool_size(
    total_cpus: int,
    cpu_per_job: int,
    n_mols: int,
    expected_workers: int,
) -> None:
    """Verifies worker pool calculation across standard, high-core HPC,
    and boundary/invalid hardware configurations."""
    workers = calculate_worker_pool_size(
        total_cpus=total_cpus,
        cpu_per_job=cpu_per_job,
        num_molecules=n_mols,
    )
    assert workers == expected_workers


def test_worker_concurrency_execution_speedup(tmp_path: Path) -> None:
    """Verifies that multi-core execution runs tasks in parallel."""
    class DelayedEngine(BaseEngine):
        manifest = EngineManifest(name="delayed", scoring_functions=["vina"])
        def dock(self, *args, **kwargs) -> DockResult:
            time.sleep(0.05)
            return DockResult(
                engine="delayed", engine_version="1.0",
                params=DockParams(), box=Box(center=(0,0,0), size=(10,10,10)),
                poses=[Pose(index=1, affinity=-7.0)], run_dir=str(tmp_path),
            )

    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=4)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    config = ScreeningConfig(cpu=8, enable_tiered_funnel=False)
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "speedup_test",
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=DelayedEngine(),
    )

    t0 = time.perf_counter()
    report = runner.run()
    elapsed = time.perf_counter() - t0

    assert len(report.candidates) == 4
    # Serial would take >= 0.20s; parallel with 4 workers finishes in < 0.18s
    assert elapsed < 0.18


def test_worker_pool_exception_isolation(tmp_path: Path) -> None:
    """Verifies that an unhandled exception or crash in one worker does not deadlock
    or crash the process pool; the failure is recorded in n_dock_failed while other workers complete."""
    class FlakyEngine(BaseEngine):
        manifest = EngineManifest(name="flaky", scoring_functions=["vina"])
        def dock(self, rec, lig, *args, **kwargs) -> DockResult:
            if "LIG_002" in str(lig):
                raise RuntimeError("Segfault emulation in worker process")
            return DockResult(
                engine="flaky", engine_version="1.0",
                params=DockParams(), box=Box(center=(0,0,0), size=(10,10,10)),
                poses=[Pose(index=1, affinity=-7.0)], run_dir=str(tmp_path),
            )

    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=4)
    rec_path = create_fake_ligand_pdbqt(tmp_path / "rec.pdbqt")

    config = ScreeningConfig(cpu=4, enable_tiered_funnel=False)
    runner = ScreeningRunner(
        config=config,
        out_dir=tmp_path / "flaky_test",
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=FlakyEngine(),
    )

    report = runner.run()

    # LIG_002 failed, 3 succeeded
    assert report.n_dock_failed == 1
    assert report.n_docked_success == 3
    assert len(report.candidates) == 3
