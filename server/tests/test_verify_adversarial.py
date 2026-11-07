"""Adversarial stress tests for dockstudio.cli.commands.verify (T1~T5 acceptance verifier).

Written by empirical_challenger to stress-test boundary conditions, corrupted inputs,
numerical invariants, and failure modes.
"""

from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.cli.commands.verify import run_verify


def _build_valid_package(root: Path, omit: list[str] | None = None) -> None:
    """Helper to build a structurally valid delivery package."""
    omit = omit or []

    # 1. conf.txt
    if "conf.txt" not in omit:
        (root / "conf.txt").write_text(
            "center_x = 10.0\n"
            "center_y = 10.0\n"
            "center_z = 10.0\n"
            "size_x = 20.0\n"
            "size_y = 20.0\n"
            "size_z = 20.0\n"
            "seed = 42\n"
            "exhaustiveness = 8\n"
            "num_modes = 9\n"
            "energy_range = 3.0\n",
            encoding="utf-8",
        )

    # 2. receptor.pdbqt
    if "receptor.pdbqt" not in omit:
        (root / "receptor.pdbqt").write_text(
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00          N\n",
            encoding="utf-8",
        )

    # 3. crystal_ligand.sdf (centered exactly at 10.0, 10.0, 10.0)
    if "crystal_ligand.sdf" not in omit:
        m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(m, randomSeed=42)
        conf = m.GetConformer()
        pts = conf.GetPositions()
        c = pts.mean(axis=0)
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, pts[i] - c + [10.0, 10.0, 10.0])
        (root / "crystal_ligand.sdf").write_text(Chem.MolToMolBlock(m), encoding="utf-8")

    # 4. poses/ and energy.csv
    poses_dir = root / "poses"
    if "poses" not in omit:
        poses_dir.mkdir(exist_ok=True)

    smiles_list = [
        "CC(=O)O", "CCO", "CCC", "CCCC", "c1ccccc1",
        "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS",
    ]
    energy_lines = ["id,affinity,smiles,inchikey"]
    for idx, smi in enumerate(smiles_list, 1):
        cid = f"MOL_{idx}"
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        ik = Chem.MolToInchiKey(mol)
        score = -10.0 + idx * 0.2
        energy_lines.append(f"{cid},{score:.2f},{smi},{ik}")
        if "poses" not in omit:
            (poses_dir / f"pose_{cid}.sdf").write_text(Chem.MolToMolBlock(mol), encoding="utf-8")

    if "energy.csv" not in omit:
        (root / "energy.csv").write_text("\n".join(energy_lines), encoding="utf-8")

    # 5. funnel.csv (strictly closed)
    if "funnel.csv" not in omit:
        (root / "funnel.csv").write_text(
            "raw,100\n"
            "desalt,10\n"
            "dedup,10\n"
            "chiral_failed,10\n"
            "prep_failed,10\n"
            "dock_failed,50\n"
            "docked_success,10\n",
            encoding="utf-8",
        )


def _zip_directory(source_dir: Path, zip_file: Path) -> None:
    """Helper to zip a directory."""
    with zipfile.ZipFile(zip_file, "w", zipfile.ZIP_DEFLATED) as z:
        for file in source_dir.rglob("*"):
            if file.is_file():
                z.write(file, file.relative_to(source_dir))


# ===========================================================================
# 1. Input Robustness & Return Code 2 Tests
# ===========================================================================

def test_nonexistent_package_path(tmp_path: Path):
    """Target: Nonexistent package path -> assert return code 2."""
    bogus_path = tmp_path / "nonexistent_delivery_package_12345"
    assert not bogus_path.exists()

    code = run_verify(bogus_path)
    assert code == 2


def test_missing_conf_file(tmp_path: Path):
    """Target: Package missing conf.txt -> assert return code 2."""
    pkg_dir = tmp_path / "pkg_no_conf"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["conf.txt"])

    code = run_verify(pkg_dir)
    assert code == 2


def test_unsupported_file_format(tmp_path: Path):
    """Target: Unsupported package format (.tar.gz) -> assert return code 2."""
    unsupported_file = tmp_path / "package.tar.gz"
    unsupported_file.write_bytes(b"dummy tar content")

    code = run_verify(unsupported_file)
    assert code == 2


def test_corrupted_zip_package(tmp_path: Path):
    """Target: Corrupted/truncated ZIP package -> assert return code 2.

    CRITICAL CHALLENGE:
    Currently verify.py does not catch zipfile.BadZipFile in load_delivery_package
    or run_verify, resulting in an unhandled crash instead of returning 2.
    """
    corrupted_zip = tmp_path / "corrupted_package.zip"
    corrupted_zip.write_bytes(b"PK\x03\x04CORRUPTED_GARBAGE_DATA_XYZ_1234567890")

    # The spec requires exit code 2 when a package cannot be parsed.
    # An unhandled zipfile.BadZipFile exception is a bug.
    code = run_verify(corrupted_zip)
    assert code == 2


def test_malformed_conf_missing_center(tmp_path: Path):
    """Target: conf.txt missing center coordinates -> assert return code 2.

    CRITICAL CHALLENGE:
    parse_conf_file raises ValueError('conf.txt 缺少必要的盒子坐标定义: ...'),
    which is not caught in run_verify, resulting in an unhandled crash.
    """
    pkg_dir = tmp_path / "pkg_bad_conf"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["conf.txt"])
    (pkg_dir / "conf.txt").write_text("size_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\nseed = 42\n")

    code = run_verify(pkg_dir)
    assert code == 2


def test_malformed_conf_non_numeric(tmp_path: Path):
    """Target: conf.txt containing non-numeric coordinates -> assert return code 2.

    CRITICAL CHALLENGE:
    float() conversion failure raises ValueError, uncaught in run_verify.
    """
    pkg_dir = tmp_path / "pkg_str_conf"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["conf.txt"])
    (pkg_dir / "conf.txt").write_text("center_x = NOT_A_FLOAT\ncenter_y = 10.0\ncenter_z = 10.0\nsize_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\n")

    code = run_verify(pkg_dir)
    assert code == 2


def test_energy_csv_sparse_row_handling(tmp_path: Path):
    """Target: energy.csv containing row with missing columns -> assert return code 0 or 1, not crash.

    CRITICAL CHALLENGE:
    When a row has fewer fields than the header, csv.DictReader sets v = None.
    verify.py calls v.strip() without null check, raising AttributeError: 'NoneType' object has no attribute 'strip'.
    """
    pkg_dir = tmp_path / "pkg_sparse_csv"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["energy.csv"])
    # Row 1 has fewer columns than the header
    (pkg_dir / "energy.csv").write_text(
        "id,affinity,smiles,inchikey\n"
        "MOL_1,-9.5\n",
        encoding="utf-8",
    )

    # Must complete verification and return 1 (or 2), NOT crash with uncaught AttributeError
    code = run_verify(pkg_dir)
    assert code in (1, 2)


# ===========================================================================
# 2. Acceptance Criteria Violation Tests (Return Code 1)
# ===========================================================================

def test_funnel_discrepancy_t5_failure(tmp_path: Path):
    """Target: Funnel arithmetic discrepancy -> assert return code 1 (T5 failure)."""
    pkg_dir = tmp_path / "pkg_funnel_fail"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["funnel.csv"])

    # Discrepancy: raw 100 - losses 90 = 10 expected, but docked_success = 8 (loss mismatch)
    (pkg_dir / "funnel.csv").write_text(
        "raw,100\n"
        "desalt,10\n"
        "dedup,10\n"
        "chiral_failed,10\n"
        "prep_failed,10\n"
        "dock_failed,50\n"
        "docked_success,8\n",  # discrepancy!
        encoding="utf-8",
    )

    code = run_verify(pkg_dir)
    assert code == 1


def test_box_centroid_deviation_t4_failure(tmp_path: Path):
    """Target: Box centroid deviation > 2.0 Å -> assert return code 1 (T4 failure)."""
    pkg_dir = tmp_path / "pkg_box_fail"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["conf.txt"])

    # Shift box center to (15.0, 10.0, 10.0), distance to crystal ligand (10, 10, 10) is 5.0 Å > 2.0 Å
    (pkg_dir / "conf.txt").write_text(
        "center_x = 15.0\n"
        "center_y = 10.0\n"
        "center_z = 10.0\n"
        "size_x = 20.0\n"
        "size_y = 20.0\n"
        "size_z = 20.0\n"
        "seed = 42\n"
        "exhaustiveness = 8\n"
        "num_modes = 9\n",
        encoding="utf-8",
    )

    code = run_verify(pkg_dir)
    assert code == 1


def test_inchikey_mismatch_t3_failure(tmp_path: Path):
    """Target: InChIKey mismatch between energy table and pose -> assert return code 1 (T3 failure)."""
    pkg_dir = tmp_path / "pkg_inchikey_fail"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["energy.csv"])

    # Corrupt InChIKey of MOL_1
    lines = ["id,affinity,smiles,inchikey"]
    lines.append("MOL_1,-9.80,CC(=O)O,WRONGINCHIKEY1234567890123")
    smiles_list = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS"]
    for idx, smi in enumerate(smiles_list, 2):
        cid = f"MOL_{idx}"
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        ik = Chem.MolToInchiKey(mol)
        lines.append(f"{cid},{-10.0 + idx * 0.2:.2f},{smi},{ik}")

    (pkg_dir / "energy.csv").write_text("\n".join(lines), encoding="utf-8")

    code = run_verify(pkg_dir)
    assert code == 1


def test_missing_seed_t1_failure(tmp_path: Path):
    """Target: Missing seed in conf.txt -> assert return code 1 (T1 failure)."""
    pkg_dir = tmp_path / "pkg_no_seed"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["conf.txt"])

    # Valid coordinates but omit seed
    (pkg_dir / "conf.txt").write_text(
        "center_x = 10.0\n"
        "center_y = 10.0\n"
        "center_z = 10.0\n"
        "size_x = 20.0\n"
        "size_y = 20.0\n"
        "size_z = 20.0\n",
        encoding="utf-8",
    )

    code = run_verify(pkg_dir)
    assert code == 1


# ===========================================================================
# 3. Acceptance Criteria Passed Tests (Return Code 0)
# ===========================================================================

def test_normal_valid_package_directory(tmp_path: Path):
    """Target: Normal valid package directory -> assert return code 0."""
    pkg_dir = tmp_path / "pkg_valid_dir"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir)

    code = run_verify(pkg_dir)
    assert code == 0


def test_normal_valid_package_zip(tmp_path: Path):
    """Target: Normal valid package ZIP archive -> assert return code 0."""
    pkg_dir = tmp_path / "pkg_valid_for_zip"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir)

    zip_file = tmp_path / "valid_package.zip"
    _zip_directory(pkg_dir, zip_file)

    code = run_verify(zip_file)
    assert code == 0


def test_valid_package_without_crystal_ligand_t4_skipped(tmp_path: Path):
    """Target: Valid package without co-crystal ligand (apo receptor) -> T4 SKIP, overall code 0."""
    pkg_dir = tmp_path / "pkg_apo_valid"
    pkg_dir.mkdir()
    _build_valid_package(pkg_dir, omit=["crystal_ligand.sdf"])

    code = run_verify(pkg_dir)
    assert code == 0
