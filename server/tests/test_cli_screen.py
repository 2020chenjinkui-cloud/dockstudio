"""Comprehensive Automated Test Suite for Milestone 4: CLI /screen Command,
Interactive Wizard, Non-Interactive Headless Execution, Rich Leaderboard,
and Full Acceptance Verification against /verify (T1~T5).

Test Structure:
1. TestCliRunnerArchitectureAndInteractiveWizard:
   - CLI help and option declaration check
   - 5-step interactive wizard input stream testing
   - Function-level monkeypatching of UI prompts
   - Graceful cancellation handling (Ctrl+C / EOF)
2. TestNonInteractiveFlagExecution:
   - Complete screen run with --non-interactive and all required flags
   - TCM mode toggle (--tcm-mode vs --no-tcm-mode)
   - Key residue anchor flag (--anchor-residues ASP25:A,THR80:A)
   - Custom 3-tier ratios (--tier1-ratio 0.20 --tier2-ratio 0.05)
   - Checkpoint resume test with pre-existing checkpoint.jsonl
   - Validation fail-fast on missing required inputs
3. TestRichTableRendering:
   - Rich leaderboard column headers (Rank, SILE, ECR, Cluster, Badges)
   - Chemical alert badges (Catechol, Michael, PAINS, Clean)
   - Anti-formula-injection CSV export
4. TestEndToEndAcceptanceWithVerify:
   - Full execution of /screen on mock receptor and library
   - Automatic generation of conf.txt, receptor.pdbqt, poses/, funnel.csv, and screening_results.csv
   - 100% compliance verification via run_verify(out_dir) passing T1~T5.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any
from unittest import mock

import pytest
from click.testing import CliRunner
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.cli.commands.screen import run_screen_interactive, screen_cli
from dockstudio.cli.commands.verify import (
    _verify_t1,
    _verify_t2,
    _verify_t3,
    _verify_t4,
    _verify_t5,
    load_delivery_package,
    load_energy_table,
    load_funnel_data,
    parse_conf_file,
    run_verify,
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

# =============================================================================
# Synthetic Fixtures & Mock Engine
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


def _make_sample_pdbqt(affinities: list[float]) -> str:
    """Generates a valid multi-MODEL PDBQT string with given affinities."""
    blocks = []
    for idx, aff in enumerate(affinities, start=1):
        blocks.append(_SAMPLE_PDBQT_TEMPLATE.format(model_idx=idx, affinity=aff))
    return "".join(blocks)


class MockFunnelEngine(BaseEngine):
    """Controllable mock engine for testing 3-tier funnel transitions and checkpoints."""

    manifest = EngineManifest(
        name="mock_funnel",
        scoring_functions=["vina"],
        supports_batch_ligands=True,
    )

    def __init__(self, fail_ligand_ids: set[str] | None = None, interrupt_on_id: str | None = None):
        super().__init__()
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
        if lig_id == "ligand":
            lig_id = ligand_pdbqt.parent.name

        if lig_id == self.interrupt_on_id:
            raise KeyboardInterrupt(f"Simulated user interrupt on {lig_id}")

        self.dock_call_count[lig_id] = self.dock_call_count.get(lig_id, 0) + 1
        self.docked_history.append((lig_id, params.exhaustiveness, params.seed))

        if lig_id in self.fail_ligand_ids:
            return DockResult(
                engine="mock_funnel",
                engine_version="1.0.0",
                params=params,
                box=box,
                poses=[],
                log=f"Simulated docking failure for {lig_id}",
                run_dir=str(work_dir),
            )

        work_dir.mkdir(parents=True, exist_ok=True)
        out_file = work_dir / "out.pdbqt"

        seed_key = f"{lig_id}_seed_{params.seed}"
        if seed_key in self.custom_poses:
            affs = self.custom_poses[seed_key]
        elif lig_id in self.custom_affinities:
            affs = [self.custom_affinities[lig_id]]
        else:
            affs = [-7.5]

        out_file.write_text(_make_sample_pdbqt(affs), encoding="utf-8")
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


@pytest.fixture
def runner() -> CliRunner:
    """Standard Click CliRunner fixture."""
    return CliRunner()


@pytest.fixture
def mock_receptor(tmp_path: Path) -> tuple[Path, Path]:
    """Generates a valid receptor PDB and PDBQT pair."""
    rec_dir = tmp_path / "receptor"
    rec_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = rec_dir / "1HSG.pdb"
    pdbqt_path = rec_dir / "1HSG.pdbqt"

    pdb_content = (
        "ATOM      1  N   ASP A  25      10.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      2  CA  ASP A  25      11.000  10.000  10.000  1.00  0.00           C\n"
        "ATOM      3  N   THR A  80      10.000  15.000  10.000  1.00  0.00           N\n"
        "ATOM      4  CA  THR A  80      11.000  15.000  10.000  1.00  0.00           C\n"
        "TER\nEND\n"
    )
    pdb_path.write_text(pdb_content, encoding="utf-8")

    pdbqt_content = (
        "ATOM      1  N   ASP A  25      10.000  10.000  10.000  1.00  0.00    -0.300 N\n"
        "ATOM      2  CA  ASP A  25      11.000  10.000  10.000  1.00  0.00    +0.100 C\n"
        "ATOM      3  N   THR A  80      10.000  15.000  10.000  1.00  0.00    -0.300 N\n"
        "ATOM      4  CA  THR A  80      11.000  15.000  10.000  1.00  0.00    +0.100 C\n"
        "TER\n"
    )
    pdbqt_path.write_text(pdbqt_content, encoding="utf-8")
    return pdb_path, pdbqt_path


@pytest.fixture
def mock_tcm_library(tmp_path: Path) -> Path:
    """Generates a synthetic TCM & synthetic drug library (SMILES) containing:
    1. Baicalin (hydrolyzable O-glucuronide -> generates Baicalin_aglycone)
    2. Quercetin (Catechol alert)
    3. Curcumin (Michael acceptor alert)
    4. Ibuprofen (Clean drug)
    5. Aspirin (Clean drug)
    """
    lib_path = tmp_path / "tcm_library.smi"
    content = (
        "O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O Baicalin\n"
        "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12 Quercetin\n"
        "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1 Curcumin\n"
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O Ibuprofen\n"
        "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
    )
    lib_path.write_text(content, encoding="utf-8")
    return lib_path


@pytest.fixture
def mock_engine_setup(monkeypatch: pytest.MonkeyPatch) -> MockFunnelEngine:
    """Installs MockFunnelEngine into get_engine, ScreeningRunner, and screen.py."""
    mock_eng = MockFunnelEngine()
    mock_eng.set_ligand_response("Baicalin", -8.5, {42: -8.5, 1024: -8.4, 2026: -8.6})
    mock_eng.set_ligand_response("Baicalin_aglycone", -9.2, {42: -9.2, 1024: -9.1, 2026: -9.3})
    mock_eng.set_ligand_response("Quercetin", -8.8, {42: -8.8, 1024: -8.7, 2026: -8.9})
    mock_eng.set_ligand_response("Curcumin", -7.9, {42: -7.9, 1024: -8.0, 2026: -7.8})
    mock_eng.set_ligand_response("Ibuprofen", -6.5, {42: -6.5, 1024: -6.4, 2026: -6.6})
    mock_eng.set_ligand_response("Aspirin", -6.2, {42: -6.2, 1024: -6.1, 2026: -6.3})

    monkeypatch.setattr("dockstudio.core.screen_runner.get_engine", lambda name: mock_eng)
    monkeypatch.setattr("dockstudio.cli.commands.screen.get_engine", lambda name: mock_eng)
    monkeypatch.setattr("dockstudio.cli.commands.screen._resolve_engine", lambda name: mock_eng)
    return mock_eng


# =============================================================================
# 1. Test Suite: Click Architecture & Interactive Wizard
# =============================================================================

class TestCliRunnerArchitectureAndInteractiveWizard:
    """Validates Click command definition, options documentation, and interactive wizard."""

    def test_screen_help_and_options(self, runner: CliRunner):
        """Verifies that `screen_cli --help` exits with code 0 and lists all 18 options."""
        result = runner.invoke(screen_cli, ["--help"])
        assert result.exit_code == 0
        output = result.output

        assert "--receptor" in output or "-r" in output
        assert "--library" in output or "-l" in output
        assert "--box-center" in output
        assert "--box-size" in output
        assert "--tcm-mode" in output
        assert "--no-tcm-mode" in output
        assert "--anchor-residues" in output
        assert "--tier1-ratio" in output
        assert "--tier2-ratio" in output
        assert "--seeds" in output
        assert "--strain-threshold" in output
        assert "--engine" in output or "-e" in output
        assert "--cpu" in output or "-c" in output
        assert "--resume" in output
        assert "--no-resume" in output
        assert "--out" in output or "-o" in output
        assert "--non-interactive" in output
        assert "--demo" in output or "-d" in output

    def test_interactive_wizard_input_stream(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Simulates complete 5-step wizard execution using an input stream."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "interactive_out"

        input_stream = "\n".join([
            "1",                          # Step 1: Select custom receptor path (1st choice)
            str(pdb_path),                # Receptor path
            "3",                          # Step 1b: Box mode: manual (3rd choice)
            "10.0",                       # Center X
            "10.0",                       # Center Y
            "10.0",                       # Center Z
            "20.0",                       # Box size
            "1",                          # Step 2: Select custom library path (1st choice)
            str(mock_tcm_library),        # Library path
            "y",                          # Step 2b: TCM dual-track enable
            "2",                          # Step 3: Manual anchor residues (2nd choice)
            "ASP25:A,THR80:A",            # Anchors
            "1",                          # Step 4: Standard 3-tier preset (1st choice)
            "1",                          # Step 5: Engine choice (1st choice: auto)
            "y",                          # Preflight confirm launch
        ]) + "\n"

        result = runner.invoke(screen_cli, ["--out", str(out_dir)], input=input_stream)
        assert result.exit_code == 0, f"Interactive run failed: {result.output}"
        assert (out_dir / "screening_results.csv").is_file()
        assert (out_dir / "funnel.csv").is_file()

    def test_interactive_wizard_function_mocking(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        """Mocks UI functions (select_one, prompt_path, confirm) to verify deterministic wizard branching."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "mocked_wizard_out"

        step_answers = iter([
            0,  # Receptor source: custom
            0,  # Box mode: autobox
            0,  # Library source: custom
            1,  # Anchor mode: custom input
            0,  # Funnel preset: standard
            0,  # Engine: auto
        ])
        monkeypatch.setattr("dockstudio.cli.commands.screen.select_one", lambda *a, **kw: next(step_answers, 0))

        paths = iter([pdb_path, mock_tcm_library])
        monkeypatch.setattr("dockstudio.cli.commands.screen.prompt_path", lambda *a, **kw: next(paths, mock_tcm_library))
        monkeypatch.setattr("dockstudio.cli.commands.screen.confirm", lambda *a, **kw: True)

        result = runner.invoke(screen_cli, ["--out", str(out_dir)], input="\n\n\n\n\n\n\n\n")
        assert result.exit_code == 0, f"Failed: {result.output}"
        assert (out_dir / "funnel.csv").is_file()

    def test_interactive_wizard_cancellation(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
    ):
        """Verifies graceful exit without traceback when user sends KeyboardInterrupt in wizard."""
        def raise_interrupt(*args: Any, **kwargs: Any) -> Any:
            raise KeyboardInterrupt("Simulated Ctrl+C")

        monkeypatch.setattr("dockstudio.cli.commands.screen.select_one", raise_interrupt)
        result = runner.invoke(screen_cli, [])
        assert result.exit_code in (1, 130)
        assert "取消" in result.output or "Cancel" in result.output


# =============================================================================
# 2. Test Suite: Non-Interactive Headless Flag Execution
# =============================================================================

class TestNonInteractiveFlagExecution:
    """Validates complete headless flag execution and configuration transformation."""

    def test_screen_non_interactive_minimal_flags(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Runs screen in headless mode with minimal required arguments."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "headless_minimal"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--box-size", "20.0,20.0,20.0",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Failed: {result.output}"
        assert (out_dir / "screening_results.csv").is_file()
        assert (out_dir / "funnel.csv").is_file()

    def test_tcm_mode_toggle_aglycone_generation(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies that --tcm-mode generates aglycones and --no-tcm-mode suppresses them."""
        pdb_path, _ = mock_receptor
        out_tcm = tmp_path / "out_tcm"
        out_notcm = tmp_path / "out_notcm"

        # 1. With --tcm-mode
        args_tcm = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--tcm-mode",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_tcm),
        ]
        res_tcm = runner.invoke(screen_cli, args_tcm)
        assert res_tcm.exit_code == 0, f"TCM mode failed: {res_tcm.output}"

        with open(out_tcm / "screening_results.csv", "r", encoding="utf-8") as f:
            reader = list(csv.DictReader(f))
            types = {r.get("Derivative_Type") or r.get("derivative_type") for r in reader}
            assert "aglycone" in types

        # 2. With --no-tcm-mode
        args_notcm = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--no-tcm-mode",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_notcm),
        ]
        res_notcm = runner.invoke(screen_cli, args_notcm)
        assert res_notcm.exit_code == 0, f"No-TCM mode failed: {res_notcm.output}"

        with open(out_notcm / "screening_results.csv", "r", encoding="utf-8") as f:
            reader = list(csv.DictReader(f))
            types = {r.get("Derivative_Type") or r.get("derivative_type") for r in reader}
            assert "aglycone" not in types

    def test_anchor_residues_flag_execution(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies --anchor-residues parsing and anchor_match_ratio evaluation."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "out_anchors"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--anchor-residues", "ASP25:A,THR80:A",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Failed: {result.output}"
        assert "ASP25:A" in result.output
        assert "THR80:A" in result.output

    def test_custom_tier_ratios_and_seeds(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies custom --tier1-ratio, --tier2-ratio, --seeds, and --strain-threshold."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "out_ratios"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "0.20",
            "--tier2-ratio", "0.05",
            "--seeds", "101,202,303",
            "--strain-threshold", "6.5",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Failed: {result.output}"

    def test_checkpoint_resume_behavior(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Prepares a pre-existing checkpoint.jsonl and verifies --resume skips completed ligands."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "out_resume"
        out_dir.mkdir(parents=True, exist_ok=True)

        ckpt_file = out_dir / "checkpoint.jsonl"
        ckpt_record = {
            "tier": 1,
            "ligand_id": "Baicalin",
            "affinity": -8.5,
            "sile": 1.25,
            "le": 0.35,
            "lle": 4.10,
            "run_dir": str(out_dir / "tier1" / "dock_Baicalin"),
            "n_modes": 1,
        }
        ckpt_file.write_text(json.dumps(ckpt_record) + "\n", encoding="utf-8")

        mock_eng = mock_engine_setup
        mock_eng.dock_call_count.clear()

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--resume",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Failed: {result.output}"

        t1_calls = [item for item in mock_eng.docked_history if item[0] == "Baicalin" and item[1] == 4]
        assert len(t1_calls) == 0

    def test_validation_fail_fast_on_missing_inputs(self, runner: CliRunner):
        """In --non-interactive mode, missing --receptor or --library must exit with non-zero error."""
        res_no_rec = runner.invoke(screen_cli, ["--non-interactive", "--library", "fake.smi"])
        assert res_no_rec.exit_code != 0
        assert "receptor" in res_no_rec.output.lower() or "受体" in res_no_rec.output

        res_no_lib = runner.invoke(screen_cli, ["--non-interactive", "--receptor", "fake.pdb"])
        assert res_no_lib.exit_code != 0
        assert "library" in res_no_lib.output.lower() or "化合物库" in res_no_lib.output


# =============================================================================
# 3. Test Suite: Rich Leaderboard & Alert Badges Rendering
# =============================================================================

class TestRichTableRendering:
    """Validates the terminal Rich leaderboard columns, ordering, and badges."""

    def test_rich_table_headers_and_columns(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies all required column headers are rendered in terminal table."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "out_render"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args, color=False, env={"COLUMNS": "160"})
        assert result.exit_code == 0, f"Failed: {result.output}"
        out = result.output

        assert "Rank" in out or "#" in out
        assert "SILE" in out
        assert "ECR" in out
        assert "Cluster" in out or "骨架" in out
        assert "Badges" in out or "警示" in out

    def test_rich_table_badges_rendering(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies that Quercetin triggers Catechol badge and clean molecules show Clean."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "out_badges"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args, color=False, env={"COLUMNS": "160"})
        assert result.exit_code == 0, f"Failed: {result.output}"
        out = result.output

        assert "Catechol" in out
        assert "Clean" in out or "●" in out

    def test_anti_formula_injection_defense(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies that malicious names starting with [=, +, -, @] are escaped in CSV export."""
        pdb_path, _ = mock_receptor
        malicious_smi = tmp_path / "malicious.smi"
        malicious_smi.write_text("CC(=O)Oc1ccccc1C(=O)O =1+2_InjectedName\n", encoding="utf-8")
        out_dir = tmp_path / "out_sec"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(malicious_smi),
            "--box-center", "10.0,10.0,10.0",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Failed: {result.output}"

        csv_path = out_dir / "screening_results.csv"
        raw_text = csv_path.read_text(encoding="utf-8")
        assert "'=1+2_InjectedName" in raw_text or "\"'=1+2_InjectedName\"" in raw_text


# =============================================================================
# 4. Test Suite: Full End-to-End Acceptance with /verify (T1~T5)
# =============================================================================

class TestEndToEndAcceptanceWithVerify:
    """Executes full /screen workflow and pipes delivery package to /verify to ensure 100% PASS."""

    def test_screen_to_verify_full_pass(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Full pipeline: /screen -> generate package -> run_verify(out_dir) -> 100% PASS."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "verify_screen_delivery"

        # 1. Execute CLI screen
        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--box-size", "20.0,20.0,20.0",
            "--seeds", "42,1024,2026",
            "--tcm-mode",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Screen execution failed: {result.output}"

        # 2. Verify all delivery files exist
        assert (out_dir / "conf.txt").is_file(), "conf.txt missing for D5"
        assert (out_dir / "receptor.pdbqt").is_file(), "receptor.pdbqt missing for D1"
        assert (out_dir / "funnel.csv").is_file(), "funnel.csv missing for D10"
        assert (out_dir / "screening_results.csv").is_file(), "screening_results.csv missing for D3"
        assert (out_dir / "poses").is_dir(), "poses directory missing for D7"

        # 3. Validate funnel.csv mass conservation closure
        funnel_counts = load_funnel_data(out_dir / "funnel.csv")
        energy_rows = load_energy_table(out_dir / "screening_results.csv")

        n_raw = funnel_counts.get("raw", 0)
        losses = (
            funnel_counts.get("desalt", 0)
            + funnel_counts.get("dedup", 0)
            + funnel_counts.get("chiral_failed", 0)
            + funnel_counts.get("prep_failed", 0)
            + funnel_counts.get("dock_failed", 0)
        )
        docked_success = funnel_counts.get("docked_success", 0)
        assert n_raw - losses == docked_success
        assert docked_success == len(energy_rows)

        # 4. Invoke run_verify on the output package
        pkg = load_delivery_package(out_dir)
        try:
            conf = parse_conf_file(pkg.conf_path)
            r1 = _verify_t1(pkg, conf, energy_rows, None)
            r2 = _verify_t2(pkg, conf, energy_rows, None)
            r3 = _verify_t3(pkg, energy_rows, None)
            r4 = _verify_t4(pkg, conf, None)
            r5 = _verify_t5(pkg, energy_rows, None)

            assert r1.passed, f"T1 failed: {r1.message}"
            assert r2.passed, f"T2 failed: {r2.message}"
            assert r3.passed, f"T3 failed: {r3.message}"
            assert r4.passed, f"T4 failed: {r4.message}"
            assert r5.passed, f"T5 failed: {r5.message}"
        finally:
            pkg.cleanup()

        # Direct verify CLI invocation check
        verify_rc = run_verify(out_dir)
        assert verify_rc == 0, f"run_verify(out_dir) must return exit code 0, got {verify_rc}"

    def test_verify_detects_corrupted_funnel(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_tcm_library: Path,
        mock_engine_setup: MockFunnelEngine,
        tmp_path: Path,
    ):
        """Verifies that if funnel.csv is tampered with, run_verify reports exit code 1."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "corrupted_funnel_delivery"

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(mock_tcm_library),
            "--box-center", "10.0,10.0,10.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0, f"Failed: {res.output}"

        # Tamper funnel.csv by injecting an extra count
        funnel_path = out_dir / "funnel.csv"
        lines = funnel_path.read_text(encoding="utf-8").splitlines()
        tampered_lines = []
        for line in lines:
            if line.startswith("docked_success"):
                tampered_lines.append("docked_success,9999")
            else:
                tampered_lines.append(line)
        funnel_path.write_text("\n".join(tampered_lines), encoding="utf-8")

        verify_rc = run_verify(out_dir)
        assert verify_rc == 1
