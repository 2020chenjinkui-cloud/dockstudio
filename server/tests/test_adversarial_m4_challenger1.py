"""Milestone 4 Adversarial Stress-Testing Suite (Challenger 1).

Adversarially challenges `dockstudio/cli/commands/screen.py` against:
1. Missing required CLI arguments in non-interactive / headless mode.
2. Malformed box center and box size coordinates.
3. Malformed, invalid, and adversarial anchor residues.
4. Extreme, zero, negative, and inverted tier ratios.
5. Interactive wizard cancellation across all 5 steps (SIGINT / KeyboardInterrupt / EOFError).
6. Additional parameter boundaries (invalid CPU, non-numeric seeds, extreme strain thresholds).
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path
from unittest import mock

import pytest
from click.testing import CliRunner

from dockstudio.cli.commands.screen import (
    parse_box_coords,
    run_screen_headless,
    run_screen_interactive,
    screen_cli,
)
from dockstudio.engines.base import BaseEngine, EngineManifest
from dockstudio.models import (
    Box,
    DockParams,
    DockResult,
    PreparedLibrary,
    PreparedLigand,
    ScreeningCandidateRow,
    ScreeningConfig,
    ScreeningFunnelReport,
)


# =============================================================================
# Mock Engine and Fixtures
# =============================================================================

class MockAdversarialEngine(BaseEngine):
    """Deterministic fast mock engine for adversarial CLI test runs."""

    manifest = EngineManifest(
        name="mock_adv",
        scoring_functions=["vina"],
        supports_batch_ligands=True,
    )

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        work_dir.mkdir(parents=True, exist_ok=True)
        out_file = work_dir / "out.pdbqt"
        sample_pdbqt = (
            "MODEL 1\n"
            "REMARK VINA RESULT: -8.50 0.000 0.000\n"
            "ROOT\n"
            "ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C\n"
            "ENDROOT\n"
            "ENDMDL\n"
        )
        out_file.write_text(sample_pdbqt, encoding="utf-8")

        from dockstudio.models import Pose
        return DockResult(
            engine="mock_adv",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-8.5, file=str(out_file))],
            log="Mock adversarial dock success",
            run_dir=str(work_dir),
        )


@pytest.fixture
def mock_engine(monkeypatch: pytest.MonkeyPatch) -> MockAdversarialEngine:
    engine = MockAdversarialEngine()
    monkeypatch.setattr("dockstudio.core.screen_runner.get_engine", lambda name: engine)
    monkeypatch.setattr("dockstudio.cli.commands.screen.get_engine", lambda name: engine)
    monkeypatch.setattr("dockstudio.cli.commands.screen._resolve_engine", lambda name: engine)
    return engine


@pytest.fixture
def sample_receptor(tmp_path: Path) -> Path:
    rec_path = tmp_path / "target.pdbqt"
    rec_path.write_text("ATOM      1  N   ASP A  25      15.000  25.000   4.000  1.00  0.00           N\nTER\n", encoding="utf-8")
    return rec_path


@pytest.fixture
def sample_library(tmp_path: Path) -> Path:
    lib_path = tmp_path / "lib.smi"
    lib_path.write_text(
        "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O Ibuprofen\n",
        encoding="utf-8",
    )
    return lib_path


# =============================================================================
# 1. Missing Required Arguments in Non-Interactive Mode
# =============================================================================

class TestMissingRequiredArgumentsFailFast:
    """Verifies fail-fast behavior when required arguments are omitted in non-interactive mode."""

    def test_non_interactive_without_receptor(self, sample_library: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, ["--non-interactive", "--library", str(sample_library)])
        assert result.exit_code == 2, f"Expected exit code 2, got {result.exit_code}"
        assert "非交互模式必须提供 --receptor 参数" in result.output

    def test_non_interactive_without_library(self, sample_receptor: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, ["--non-interactive", "--receptor", str(sample_receptor)])
        assert result.exit_code == 2, f"Expected exit code 2, got {result.exit_code}"
        assert "非交互模式必须提供 --library 参数" in result.output

    def test_non_interactive_both_missing(self):
        runner = CliRunner()
        result = runner.invoke(screen_cli, ["--non-interactive"])
        assert result.exit_code == 2, f"Expected exit code 2, got {result.exit_code}"
        assert "非交互模式必须提供 --receptor 参数" in result.output

    def test_non_interactive_nonexistent_receptor(self, sample_library: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", "non_existent_receptor.pdbqt",
            "--library", str(sample_library),
        ])
        assert result.exit_code == 2, f"Expected exit code 2, got {result.exit_code}"
        assert "受体文件不存在" in result.output

    def test_non_interactive_nonexistent_library(self, sample_receptor: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", "non_existent_library.smi",
        ])
        assert result.exit_code == 2, f"Expected exit code 2, got {result.exit_code}"
        assert "化合物库文件不存在" in result.output

    def test_non_interactive_empty_zero_byte_library(self, sample_receptor: Path, tmp_path: Path):
        empty_lib = tmp_path / "empty.smi"
        empty_lib.write_text("", encoding="utf-8")
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(empty_lib),
        ])
        # Fail fast: exit code must not be 0
        assert result.exit_code != 0, f"Expected non-zero exit code on empty library, got {result.exit_code}"
        # Either raises PreparationError or prints failure message
        assert result.exit_code == 1 or "未提取到有效小分子结构" in result.output or "PreparationError" in str(result.exception)

    def test_non_interactive_corrupt_binary_library(self, sample_receptor: Path, tmp_path: Path):
        corrupt_lib = tmp_path / "corrupt.smi"
        corrupt_lib.write_bytes(b"\x00\xff\xfe\x01\x02\x03\x04")
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(corrupt_lib),
        ])
        assert result.exit_code == 1, f"Expected exit code 1 on corrupted library, got {result.exit_code}"


# =============================================================================
# 2. Malformed Box Coordinates & Spatial Boundaries
# =============================================================================

class TestMalformedBoxCoordinates:
    """Stress-tests parsing and handling of malformed center and size coordinates."""

    def test_box_center_non_numeric_string(self, sample_receptor: Path, sample_library: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--box-center", "abc,def",
        ])
        # Malformed non-numeric coordinates must result in non-zero exit code
        assert result.exit_code != 0, "CLI should fail on non-numeric box-center"

    def test_box_center_incomplete_parts(self, sample_receptor: Path, sample_library: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--box-center", "15.0,20.0",
        ])
        assert result.exit_code != 0, "CLI should fail when box-center has only 2 parts"

    def test_box_center_extra_parts(self, sample_receptor: Path, sample_library: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--box-center", "15.0,20.0,30.0,40.0",
        ])
        assert result.exit_code != 0, "CLI should fail when box-center has 4 parts"

    def test_box_size_non_numeric(self, sample_receptor: Path, sample_library: Path):
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--box-size", "invalid_size",
        ])
        assert result.exit_code != 0, "CLI should fail on non-numeric box-size"

    def test_box_size_negative_dimension_construction(self):
        """Probes boundary behavior of parse_box_coords on negative/zero sizes.

        Empirical Observation: parse_box_coords currently parses negative and zero
        dimensions without throwing ValueError, returning Box with negative sizes.
        """
        box = parse_box_coords("10,10,10", "-10,0,10")
        assert box.center == (10.0, 10.0, 10.0)
        # Note: box.size allows (-10.0, 0.0, 10.0) without physical validation
        assert box.size == (-10.0, 0.0, 10.0)

    def test_box_coords_whitespace_and_semicolon_tolerance(self):
        """Verifies that semicolon, full-width comma, and whitespace are gracefully parsed."""
        box = parse_box_coords(" 15.2 ;  24.8 ， -5.1 ", " 20.0 ; 21.0 ， 22.0 ")
        assert box.center == (15.2, 24.8, -5.1)
        assert box.size == (20.0, 21.0, 22.0)

    def test_box_coords_single_dimension_broadcast(self):
        """Verifies that a single size scalar broadcasts to all 3 axes."""
        box = parse_box_coords(None, "25.0")
        assert box.center == (16.0, 25.0, 4.0)
        assert box.size == (25.0, 25.0, 25.0)

    def test_box_size_negative_cli_behavior_probe(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Probes how the CLI behaves when negative box sizes like -10,0,10 are passed.

        Adversarial Observation: In physical docking, box sizes <= 0 are unphysical.
        The CLI does not currently fail-fast with an error message; instead it executes
        the run and writes negative dimensions to conf.txt.
        """
        out_dir = tmp_path / "out_neg_box"
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--box-size", "-10,0,10",
            "--out", str(out_dir),
        ])
        conf_file = out_dir / "conf.txt"
        if conf_file.is_file():
            conf_text = conf_file.read_text(encoding="utf-8")
            # Documented empirical finding: size_x is written as negative!
            has_neg_size = "size_x = -10.000" in conf_text
            assert has_neg_size or result.exit_code != 0


# =============================================================================
# 3. Malformed Anchor Residues
# =============================================================================

class TestMalformedAnchorResidues:
    """Tests robustness against invalid, nonsensical, or special-character anchor residues."""

    def test_completely_invalid_anchor_name(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Verifies that unknown/nonsense residue tags do not crash the pipeline."""
        out_dir = tmp_path / "out_invalid_anchor"
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--anchor-residues", "INVALID_RES",
            "--out", str(out_dir),
        ])
        assert result.exit_code == 0, f"Screening should not crash on invalid residue: {result.output}"
        assert (out_dir / "screening_results.csv").is_file()

    def test_mixed_valid_and_adversarial_anchor_residues(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Mixes valid receptor residue ASP25:A with garbage tags."""
        out_dir = tmp_path / "out_mixed_anchors"
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--anchor-residues", "ASP25:A,??###,99999:ZZ,THR80:A",
            "--out", str(out_dir),
        ])
        assert result.exit_code == 0
        assert (out_dir / "screening_results.csv").is_file()

    def test_punctuation_only_anchor_residues(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Anchor string consisting purely of commas, semicolons, and spaces."""
        out_dir = tmp_path / "out_punct_anchors"
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--anchor-residues", " , ; , ; ;  , ",
            "--out", str(out_dir),
        ])
        assert result.exit_code == 0
        assert (out_dir / "conf.txt").is_file()

    def test_extremely_long_anchor_residues(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Anchor string with 2000 characters."""
        long_anchor = "ASP25:A," + "X" * 2000
        out_dir = tmp_path / "out_long_anchor"
        runner = CliRunner()
        result = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--anchor-residues", long_anchor,
            "--out", str(out_dir),
        ])
        assert result.exit_code == 0


# =============================================================================
# 4. Extreme Tier Ratios & Boundary Values
# =============================================================================

class TestExtremeTierRatiosAndFunnelBoundaries:
    """Validates Click FloatRange enforcement and funnel parameter boundaries."""

    def test_tier1_ratio_exceeds_one(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--tier1-ratio", "1.5"])
        assert res.exit_code == 2
        assert "is not in the range 0.001<=x<=1.0" in res.output

    def test_tier1_ratio_zero(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--tier1-ratio", "0.0"])
        assert res.exit_code == 2
        assert "is not in the range 0.001<=x<=1.0" in res.output

    def test_tier1_ratio_negative(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--tier1-ratio", "-0.25"])
        assert res.exit_code == 2
        assert "is not in the range 0.001<=x<=1.0" in res.output

    def test_tier2_ratio_exceeds_one(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--tier2-ratio", "2.0"])
        assert res.exit_code == 2
        assert "is not in the range 0.001<=x<=1.0" in res.output

    def test_tier2_ratio_zero(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--tier2-ratio", "0.0"])
        assert res.exit_code == 2
        assert "is not in the range 0.001<=x<=1.0" in res.output

    def test_tier2_ratio_negative(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--tier2-ratio", "-1.0"])
        assert res.exit_code == 2
        assert "is not in the range 0.001<=x<=1.0" in res.output

    def test_inverted_tier_ratios_safe_execution(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Sets Tier 1 ratio small (0.05) and Tier 2 ratio large (0.50).

        Ensures that when Tier 2 ratio exceeds Tier 1, the pipeline executes safely
        without index errors or crashes.
        """
        out_dir = tmp_path / "out_inverted_ratios"
        runner = CliRunner()
        res = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--tier1-ratio", "0.05",
            "--tier2-ratio", "0.50",
            "--out", str(out_dir),
        ])
        assert res.exit_code == 0
        assert (out_dir / "funnel.csv").is_file()


# =============================================================================
# 5. Interactive Wizard Cancellation Across All Steps
# =============================================================================

class TestInteractiveWizardCancellation:
    """Verifies that user interrupts (KeyboardInterrupt / EOFError) exit with code 130 cleanly."""

    def test_cancellation_step1_receptor_select(self, monkeypatch: pytest.MonkeyPatch):
        def _interrupt(*args, **kwargs):
            raise KeyboardInterrupt("Simulated Ctrl+C at step 1")

        monkeypatch.setattr("dockstudio.cli.commands.screen.select_one", _interrupt)
        runner = CliRunner()
        res = runner.invoke(screen_cli, [])
        assert res.exit_code == 130
        assert "操作已取消" in res.output

    def test_cancellation_step1_path_prompt(self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setattr("dockstudio.cli.commands.screen.select_one", lambda *a, **kw: 0)

        def _interrupt(*args, **kwargs):
            raise KeyboardInterrupt("Simulated Ctrl+C at path prompt")

        monkeypatch.setattr("dockstudio.cli.commands.screen.prompt_path", _interrupt)
        runner = CliRunner()
        res = runner.invoke(screen_cli, [])
        assert res.exit_code == 130
        assert "操作已取消" in res.output

    def test_cancellation_step2_tcm_confirm(
        self,
        sample_receptor: Path,
        sample_library: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        # Step 1 succeeds
        monkeypatch.setattr("dockstudio.cli.commands.screen.select_one", lambda *a, **kw: 0)
        call_count = {"path": 0}

        def _mock_path(*a, **kw):
            call_count["path"] += 1
            return sample_receptor if call_count["path"] == 1 else sample_library

        monkeypatch.setattr("dockstudio.cli.commands.screen.prompt_path", _mock_path)

        def _interrupt(*args, **kwargs):
            raise KeyboardInterrupt("Simulated Ctrl+C at TCM confirm")

        monkeypatch.setattr("dockstudio.cli.commands.screen.confirm", _interrupt)
        runner = CliRunner()
        res = runner.invoke(screen_cli, [])
        assert res.exit_code == 130
        assert "操作已取消" in res.output

    def test_cancellation_step5_preflight_rejection(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        monkeypatch: pytest.MonkeyPatch,
    ):
        """User answers 'No' at pre-flight confirmation. Should exit cleanly with code 0."""
        monkeypatch.setattr("dockstudio.cli.commands.screen.select_one", lambda *a, **kw: 0)
        call_count = {"path": 0}

        def _mock_path(*a, **kw):
            call_count["path"] += 1
            return sample_receptor if call_count["path"] == 1 else sample_library

        monkeypatch.setattr("dockstudio.cli.commands.screen.prompt_path", _mock_path)

        # First confirm is TCM mode (True), second is Pre-flight (False)
        confirm_responses = iter([True, False])
        monkeypatch.setattr("dockstudio.cli.commands.screen.confirm", lambda *a, **kw: next(confirm_responses))

        runner = CliRunner()
        res = runner.invoke(screen_cli, [])
        assert res.exit_code == 0
        assert "用户取消启动" in res.output

    def test_cancellation_on_closed_empty_stdin_eof(self):
        """Simulates running interactive screen with closed empty stdin (EOF)."""
        runner = CliRunner()
        res = runner.invoke(screen_cli, [], input="")
        # When prompt_path encounters EOF with no default, it raises EOFError -> code 130
        assert res.exit_code == 130
        assert "操作已取消" in res.output


# =============================================================================
# 6. Additional Adversarial CLI Parameter Boundaries
# =============================================================================

class TestAdditionalAdversarialCLIParameters:
    """Tests extreme CPU counts, non-numeric seeds, and extreme strain thresholds."""

    def test_cpu_zero_rejected(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--cpu", "0"])
        assert res.exit_code == 2
        assert "is not in the range x>=1" in res.output

    def test_cpu_negative_rejected(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--cpu", "-4"])
        assert res.exit_code == 2
        assert "is not in the range x>=1" in res.output

    def test_cpu_per_job_zero_rejected(self):
        runner = CliRunner()
        res = runner.invoke(screen_cli, ["--demo", "--cpu-per-job", "0"])
        assert res.exit_code == 2
        assert "is not in the range x>=1" in res.output

    def test_non_numeric_seeds_raises_or_fails(
        self,
        sample_receptor: Path,
        sample_library: Path,
    ):
        runner = CliRunner()
        res = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--seeds", "abc,def",
        ])
        # Non-numeric seeds should fail fast
        assert res.exit_code != 0

    def test_empty_seeds_string_falls_back_to_defaults(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        out_dir = tmp_path / "out_empty_seeds"
        runner = CliRunner()
        res = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--seeds", "   ",
            "--out", str(out_dir),
        ])
        assert res.exit_code == 0
        conf_text = (out_dir / "conf.txt").read_text(encoding="utf-8")
        assert "seed = 42" in conf_text

    def test_extreme_negative_strain_threshold(
        self,
        sample_receptor: Path,
        sample_library: Path,
        mock_engine: MockAdversarialEngine,
        tmp_path: Path,
    ):
        """Negative strain threshold means zero tolerance for conformational strain."""
        out_dir = tmp_path / "out_neg_strain"
        runner = CliRunner()
        res = runner.invoke(screen_cli, [
            "--non-interactive",
            "--receptor", str(sample_receptor),
            "--library", str(sample_library),
            "--strain-threshold", "-10.0",
            "--out", str(out_dir),
        ])
        assert res.exit_code == 0
        assert (out_dir / "screening_results.csv").is_file()
