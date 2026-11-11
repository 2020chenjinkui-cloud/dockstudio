"""Adversarial stress-test suite for Milestone 3 HTVS ScreenRunner.

Engineered by Challenger 2 to empirically stress-test:
1. Mid-run crash simulation at arbitrary compounds (compounds 7 and 19) on 30 ligands,
   restarting with resume=True, verifying zero duplicate dockings and zero data loss.
2. Checkpoint corruption recovery: trailing partial lines, non-JSON bytes, empty lines,
   verifying graceful recovery of all valid records without data loss.
3. Funnel conservation mass accounting across diverse input configurations (100% pass,
   50% pass, 0% pass) against verify.py T5 validation logic (total_disc == 0).
4. CSV formula injection protection: inputs starting with '=', '+', '-', '@', '\t', '\r'
   are safely escaped in exported CSV files.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import tempfile
from pathlib import Path
from typing import Any

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.cli.commands.verify import (
    DeliveryPackage,
    _verify_t5,
    load_energy_table,
    load_funnel_data,
)
from dockstudio.core.screen_runner import (
    CheckpointManager,
    ScreeningRunner,
    _csv_safe,
    export_funnel_csv,
    export_screening_results_csv,
)
from dockstudio.engines.base import BaseEngine, EngineManifest
from dockstudio.models import (
    Box,
    DockParams,
    DockResult,
    Pose,
    PreparedLibrary,
    PreparedLigand,
    ScreeningCandidateRow,
    ScreeningConfig,
    ScreeningFunnelReport,
)


class MockAdversarialCrashEngine(BaseEngine):
    """Mock engine that tracks call counts and can inject crashes at specific call indices or ligand IDs."""

    def __init__(self, crash_at_ids: set[str] | None = None, fail_dock_ids: set[str] | None = None):
        super().__init__()
        self.crash_at_ids = set(crash_at_ids or [])
        self.fail_dock_ids = set(fail_dock_ids or [])
        self.dock_call_count: dict[str, int] = {}
        self.call_history: list[str] = []

    @property
    def manifest(self) -> EngineManifest:
        return EngineManifest(
            name="mock_crash_engine",
            version="1.0.0",
            description="Mock engine with deterministic crash injection",
            supported_scoring=["vina"],
            supports_batch_ligands=False,
            supports_gpu=False,
        )

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
    ) -> DockResult:
        lig_id = ligand_pdbqt.stem
        if lig_id == "ligand":
            lig_id = ligand_pdbqt.parent.name

        self.dock_call_count[lig_id] = self.dock_call_count.get(lig_id, 0) + 1
        self.call_history.append(lig_id)

        # Crash injection simulation
        if lig_id in self.crash_at_ids:
            # Remove so it won't crash again on resume
            self.crash_at_ids.remove(lig_id)
            raise KeyboardInterrupt(f"Simulated abrupt process crash at ligand {lig_id}")

        if lig_id in self.fail_dock_ids:
            return DockResult(
                engine="mock_crash_engine",
                engine_version="1.0.0",
                params=params,
                box=box,
                poses=[],
                log=f"Docking failed for {lig_id}",
                run_dir=str(work_dir),
            )

        work_dir.mkdir(parents=True, exist_ok=True)
        out_file = work_dir / "out.pdbqt"
        sample_pdbqt = (
            "REMARK VINA RESULT: -7.50 0.00 0.00\n"
            "ATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00          C\n"
        )
        out_file.write_text(sample_pdbqt, encoding="utf-8")

        pose = Pose(
            index=1,
            affinity=-7.5,
            rmsd_lb=0.0,
            rmsd_ub=0.0,
        )
        return DockResult(
            engine="mock_crash_engine",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=[pose],
            log="OK",
            run_dir=str(work_dir),
        )


def _make_dummy_pdbqt(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00          N\n", encoding="utf-8")
    return path


def _build_test_library(tmp_path: Path, n_ligands: int = 30) -> tuple[PreparedLibrary, list[PreparedLigand]]:
    ligands_dir = tmp_path / "test_ligands"
    ligands_dir.mkdir(parents=True, exist_ok=True)
    ligands: list[PreparedLigand] = []

    for i in range(1, n_ligands + 1):
        lig_id = f"LIG_{i:03d}"
        lig_sub = ligands_dir / lig_id
        lig_sub.mkdir(parents=True, exist_ok=True)
        p = lig_sub / "ligand.pdbqt"
        p.write_text(
            f"REMARK VINA RESULT: -7.50 0.00 0.00\nATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00          C\n",
            encoding="utf-8"
        )
        lig = PreparedLigand(
            id=lig_id,
            source="smiles",
            name=f"Compound_{i:03d}",
            pdbqt_path=str(p),
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            n_rotatable_bonds=2,
        )
        ligands.append(lig)

    prep_lib = PreparedLibrary(
        id=f"lib_{n_ligands}",
        name=f"TestLib_{n_ligands}",
        path=str(ligands_dir),
        ligand_ids=[lig.id for lig in ligands],
        n_raw=n_ligands,
        n_ready=n_ligands,
        summary_csv_path=str(tmp_path / "summary.csv"),
    )
    return prep_lib, ligands


# =============================================================================
# 1. Mid-Run Crash Simulation & Zero-Loss Resumption Tests
# =============================================================================

def test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume(tmp_path: Path):
    """Crash simulation on 30 ligands:
    - Crash at compound 7 (LIG_007)
    - Resume with resume=True -> Crash at compound 19 (LIG_019)
    - Resume with resume=True -> Finishes remaining compounds
    Asserts:
    - Zero duplicate completed dockings across runs.
    - Exactly 30 processed entries in checkpoint.jsonl.
    - Zero data loss in final candidates report (all 30 present).
    """
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=30)
    rec_path = _make_dummy_pdbqt(tmp_path / "receptor.pdbqt")
    out_dir = tmp_path / "crash_resume_run"

    # Configure engine to crash at LIG_007 and LIG_019
    engine = MockAdversarialCrashEngine(crash_at_ids={"LIG_007", "LIG_019"})

    # Run 1: Should fail at LIG_007
    config = ScreeningConfig(enable_tiered_funnel=False, cpu=2)
    runner1 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,  # 1 worker -> deterministic sequential execution
    )

    with pytest.raises(KeyboardInterrupt, match="Simulated abrupt process crash at ligand LIG_007"):
        runner1.run(resume=False)

    ckpt1 = out_dir / "checkpoint.jsonl"
    assert ckpt1.is_file()
    lines1 = [json.loads(line) for line in ckpt1.read_text(encoding="utf-8").strip().splitlines() if line.strip()]
    # Up to 6 compounds (LIG_001..LIG_006) should be checkpointed
    assert len(lines1) == 6
    assert {d["ligand_id"] for d in lines1} == {f"LIG_{i:03d}" for i in range(1, 7)}

    # Run 2: Resume -> Should crash at LIG_019
    runner2 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )

    with pytest.raises(KeyboardInterrupt, match="Simulated abrupt process crash at ligand LIG_019"):
        runner2.run(resume=True)

    lines2 = [json.loads(line) for line in ckpt1.read_text(encoding="utf-8").strip().splitlines() if line.strip()]
    # 18 compounds (LIG_001..LIG_018) should now be checkpointed
    assert len(lines2) == 18
    assert {d["ligand_id"] for d in lines2} == {f"LIG_{i:03d}" for i in range(1, 19)}

    # Run 3: Resume -> Should finish all remaining compounds to 30
    runner3 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )
    report = runner3.run(resume=True)

    # Assertions
    # 1. Total lines in checkpoint.jsonl must be exactly 30
    lines3 = [json.loads(line) for line in ckpt1.read_text(encoding="utf-8").strip().splitlines() if line.strip()]
    assert len(lines3) == 30
    ids_in_ckpt = [d["ligand_id"] for d in lines3]
    assert len(set(ids_in_ckpt)) == 30, "Duplicate entries found in checkpoint.jsonl!"

    # 2. Docking call counts: exactly 1 per completed ligand (retried ligands were retried once after crash)
    for i in range(1, 31):
        lig_id = f"LIG_{i:03d}"
        if lig_id in ("LIG_007", "LIG_019"):
            assert engine.dock_call_count[lig_id] == 2, f"Crash compound {lig_id} should be attempted then retried"
        else:
            assert engine.dock_call_count[lig_id] == 1, f"Compound {lig_id} docked more than once: {engine.dock_call_count[lig_id]}"

    # 3. Final report candidates must be exactly 30 with zero data loss
    assert len(report.candidates) == 30
    candidate_ids = [c.ligand_id for c in report.candidates]
    assert len(set(candidate_ids)) == 30
    assert set(candidate_ids) == {f"LIG_{i:03d}" for i in range(1, 31)}


def test_concurrent_crash_and_resumption(tmp_path: Path):
    """Tests resume behavior when running under higher concurrency (workers > 1).
    Ensures that when a crash occurs in a worker thread, pending jobs are cancelled,
    valid completed jobs are flushed, and resuming reaches full completion with zero duplicates.
    """
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=20)
    rec_path = _make_dummy_pdbqt(tmp_path / "receptor.pdbqt")
    out_dir = tmp_path / "concurrent_crash_run"

    engine = MockAdversarialCrashEngine(crash_at_ids={"LIG_010"})

    # Run with 4 cpus (2 workers)
    config = ScreeningConfig(enable_tiered_funnel=False, cpu=4)
    runner1 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )

    with pytest.raises(KeyboardInterrupt):
        runner1.run(resume=False)

    # Resume to completion
    runner2 = ScreeningRunner(
        config=config,
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )
    report = runner2.run(resume=True)

    assert len(report.candidates) == 20
    candidate_ids = {c.ligand_id for c in report.candidates}
    assert candidate_ids == {f"LIG_{i:03d}" for i in range(1, 21)}

    ckpt = out_dir / "checkpoint.jsonl"
    lines = [json.loads(line) for line in ckpt.read_text(encoding="utf-8").strip().splitlines() if line.strip()]
    assert len(lines) == 20
    assert len({d["ligand_id"] for d in lines}) == 20


# =============================================================================
# 2. Checkpoint JSONL Corruption Recovery Tests
# =============================================================================

def test_checkpoint_corruption_recovery_trailing_partial_and_garbage(tmp_path: Path):
    """Adversarial stress test on CheckpointManager:
    - Multiple valid lines
    - Middle empty lines
    - Truncated half-written trailing line
    - Binary non-JSON bytes
    - Trailing whitespace
    Verifies that all valid lines are loaded without error, corrupt lines safely ignored,
    and subsequent writes/resumption succeed seamlessly.
    """
    ckpt_path = tmp_path / "corrupted_checkpoint.jsonl"
    ckpt = CheckpointManager(ckpt_path)

    # Write 5 valid entries
    valid_entries = [
        {"ligand_id": f"LIG_{i:03d}", "affinity": -8.0 - i * 0.1, "tier": 1}
        for i in range(1, 6)
    ]
    raw_content = ""
    for e in valid_entries[:3]:
        raw_content += json.dumps(e) + "\n"

    # Inject middle corruption and empty lines
    raw_content += "\n\n"
    raw_content += "{'invalid': single_quoted_json}\n"
    raw_content += "NOT_JSON_AT_ALL_12345\n"
    raw_content += "\n"

    # Add 4th and 5th valid entries
    for e in valid_entries[3:]:
        raw_content += json.dumps(e) + "\n"

    # Inject trailing truncated line and binary garbage bytes
    raw_content += '{"ligand_id": "LIG_006", "affinity": -9.'  # Half-written
    ckpt_path.write_text(raw_content, encoding="utf-8", errors="replace")

    # Append raw non-JSON binary bytes
    with open(ckpt_path, "ab") as f:
        f.write(b"\x00\xff\xfe\xca\xfe\xba\xbe\n")

    # Load completed
    state = ckpt.load_completed()
    by_id = state.by_id

    # Must recover exactly the 5 valid entries
    assert len(by_id) == 5
    for i in range(1, 6):
        lig_id = f"LIG_{i:03d}"
        assert lig_id in by_id
        assert by_id[lig_id]["affinity"] == pytest.approx(-8.0 - i * 0.1)

    assert "LIG_006" not in by_id

    # Verify that recording new entries appends safely and can be read back
    ckpt.record_single("LIG_006", -9.5)
    new_state = ckpt.load_completed()
    assert len(new_state.by_id) == 6
    assert "LIG_006" in new_state.by_id
    assert new_state.by_id["LIG_006"]["affinity"] == -9.5


# =============================================================================
# 3. Funnel Conservation Mass Accounting & /verify T5 Tests
# =============================================================================

@pytest.mark.parametrize(
    ("n_raw", "n_desalt", "n_dedup", "n_chiral", "n_prep", "n_dock_fail", "scenario_name"),
    [
        (50, 0, 0, 0, 0, 0, "100% Pass (0 losses)"),
        (100, 10, 10, 10, 10, 10, "50% Pass (50 losses, 50 success)"),
        (20, 5, 5, 5, 5, 0, "0% Pass via Preparation (0 ready)"),
        (20, 0, 0, 0, 0, 20, "0% Pass via Docking (20 dock failures)"),
        (1, 0, 0, 0, 0, 0, "Single compound 100% Pass"),
        (1, 1, 0, 0, 0, 0, "Single compound 0% Pass"),
        (1000, 150, 100, 50, 100, 200, "Large scale 400 success"),
    ],
)
def test_funnel_mass_accounting_across_diverse_scenarios_verify_t5(
    tmp_path: Path,
    n_raw: int,
    n_desalt: int,
    n_dedup: int,
    n_chiral: int,
    n_prep: int,
    n_dock_fail: int,
    scenario_name: str,
):
    """Directly verifies that funnel mass conservation holds and passes /verify T5 logic:
    total_disc == disc_1 + disc_2 == 0.
    """
    n_prep_losses = n_desalt + n_dedup + n_chiral + n_prep
    n_ready = max(0, n_raw - n_prep_losses)
    n_dock_success = max(0, n_ready - n_dock_fail)

    report = ScreeningFunnelReport(
        n_raw=n_raw,
        n_desalt_fail=n_desalt,
        n_dedup=n_dedup,
        n_chiral_fail=n_chiral,
        n_prep_fail=n_prep,
        n_ready=n_ready,
        n_dock_failed=n_dock_fail,
        n_docked_success=n_dock_success,
    )

    funnel_csv = tmp_path / f"funnel_{scenario_name.replace(' ', '_').replace('%', 'pct')}.csv"
    export_funnel_csv(report, funnel_csv)

    # 1. Parse via load_funnel_data
    funnel = load_funnel_data(funnel_csv)
    raw_val = funnel.get("raw", 0)
    total_losses = (
        funnel.get("desalt", 0)
        + funnel.get("dedup", 0)
        + funnel.get("chiral_failed", 0)
        + funnel.get("prep_failed", 0)
        + funnel.get("dock_failed", 0)
    )
    docked_success = funnel.get("docked_success", 0)

    expected_docked = raw_val - total_losses
    disc_1 = abs(expected_docked - docked_success)
    assert disc_1 == 0, f"Scenario '{scenario_name}' failed disc_1: {disc_1} != 0"

    # 2. Check directly against _verify_t5
    # Build simulated delivery package
    pkg_dir = tmp_path / f"pkg_{scenario_name.replace(' ', '_').replace('%', 'pct')}"
    pkg_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(funnel_csv, pkg_dir / "funnel.csv")

    pkg = DeliveryPackage(pkg_dir)

    # Mock energy rows matching docked_success
    mock_energy_rows = [{"id": f"MOL_{i}", "affinity": -8.0} for i in range(n_dock_success)]

    res = _verify_t5(pkg, mock_energy_rows, None)
    if n_raw > 0:
        assert res.passed is True, f"Scenario '{scenario_name}' failed _verify_t5: {res.message}, {res.observed}"
        assert res.status == "PASS"
        assert "差值=0" in res.observed
    else:
        assert res.passed is False


def test_verify_t5_with_actual_screening_results_csv(tmp_path: Path):
    """Tests what happens when _verify_t5 is run on a package containing BOTH
    funnel.csv and screening_results.csv generated by ScreeningRunner.
    Empirically verifies whether load_energy_table correctly parses screening_results.csv
    or if there is a column header mismatch (e.g. Affinity_Mean_kcal_mol vs affinity).
    """
    prep_lib, ligands = _build_test_library(tmp_path, n_ligands=5)
    rec_path = _make_dummy_pdbqt(tmp_path / "receptor.pdbqt")
    out_dir = tmp_path / "screen_delivery_pkg"

    engine = MockAdversarialCrashEngine()
    runner = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False),
        out_dir=out_dir,
        receptor_pdbqt=rec_path,
        prepared_library=prep_lib,
        engine=engine,
    )
    report = runner.run()

    # Discover package
    pkg = DeliveryPackage(out_dir)
    assert pkg.funnel_path is not None
    assert pkg.energy_path is not None

    # T5 verification with candidate energy rows
    candidate_energy_rows = [
        {"id": c.ligand_id, "affinity": c.affinity_mean, "smiles": c.smiles}
        for c in report.candidates
    ]
    res_t5 = _verify_t5(pkg, candidate_energy_rows, None)
    assert res_t5.passed is True
    assert res_t5.status == "PASS"
    assert "差值=0" in res_t5.observed

    # Cross-Milestone Finding: load_energy_table on screening_results.csv
    # Documents that load_energy_table currently expects 'affinity' or 'vina_score',
    # while screening_results.csv headers use 'Affinity_Mean_kcal_mol' and 'Ligand_ID'.
    raw_energy_rows = load_energy_table(pkg.energy_path)
    # Highlight finding for M4 acceptance alignment
    print(f"NOTICE [Cross-Milestone M3/M4]: load_energy_table parsed {len(raw_energy_rows)} rows from "
          f"{pkg.energy_path.name} due to Affinity_Mean_kcal_mol vs affinity header convention.")


# =============================================================================
# 4. CSV Formula Injection Protection Tests
# =============================================================================

def test_csv_formula_injection_defense():
    """Unit test for _csv_safe sanitizing malicious formula triggers:
    '=', '+', '-', '@', '\t', '\r'
    """
    dangerous_payloads = [
        "=SUM(A1:A10)",
        "=cmd|' /C calc'!A0",
        "+12345",
        "-5+5",
        "@SUM(B1:B2)",
        "\tmalicious_tab",
        "\rmalicious_cr",
        '=HYPERLINK("http://evil.com","Click")',
    ]
    for payload in dangerous_payloads:
        safe = _csv_safe(payload)
        assert safe.startswith("'"), f"Payload {payload!r} was not escaped with leading single quote: {safe!r}"
        assert safe[1:] == payload

    safe_payloads = [
        "Normal_Name",
        "Compound_123",
        "CCO",
        "benzene",
        "",
        None,
    ]
    for payload in safe_payloads:
        safe = _csv_safe(payload)
        assert not safe.startswith("'") or payload == "", f"Safe string was improperly escaped: {payload!r} -> {safe!r}"


def test_screening_results_csv_formula_injection_adversarial_e2e(tmp_path: Path):
    """End-to-end adversarial test: creates candidate rows containing formula injection
    strings in ligand_id, ligand_name, derivative_type, smiles, alerts.
    Verifies that export_screening_results_csv outputs all strings safely escaped with single quote.
    """
    csv_file = tmp_path / "screening_results.csv"

    malicious_candidates = [
        ScreeningCandidateRow(
            rank=1,
            ligand_id="=1+1",
            ligand_name="+MALICIOUS_NAME",
            derivative_type="-TAUTOMER",
            parent_id="@PARENT_EVIL",
            affinity_mean=-9.5,
            affinity_std=0.1,
            sile=1.8,
            le=0.45,
            lle=2.1,
            strain_energy=1.2,
            alerts=["=ALERT1", "+ALERT2"],
            smiles="=SMILES_INJECTION",
        ),
        ScreeningCandidateRow(
            rank=2,
            ligand_id="@EVIL_ID",
            ligand_name="\tTAB_NAME",
            derivative_type="original",
            affinity_mean=-8.0,
            affinity_std=0.2,
            sile=1.5,
            le=0.4,
            lle=1.9,
            strain_energy=0.8,
            alerts=[],
            smiles="CCO",
        ),
    ]

    export_screening_results_csv(malicious_candidates, csv_file)
    assert csv_file.is_file()

    with open(csv_file, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) == 2

    row1 = rows[0]
    assert row1["Ligand_ID"] == "'=1+1"
    assert row1["Ligand_Name"] == "'+MALICIOUS_NAME"
    assert row1["Derivative_Type"] == "'-TAUTOMER"
    assert row1["Parent_ID"] == "'@PARENT_EVIL"
    assert row1["Alerts"] == "'=ALERT1;+ALERT2"
    assert row1["SMILES"] == "'=SMILES_INJECTION"

    row2 = rows[1]
    assert row2["Ligand_ID"] == "'@EVIL_ID"
    assert row2["Ligand_Name"] == "'\tTAB_NAME"
    assert row2["Derivative_Type"] == "original"
    assert row2["SMILES"] == "CCO"
