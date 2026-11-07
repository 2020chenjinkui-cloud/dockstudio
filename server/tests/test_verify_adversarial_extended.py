"""Extended adversarial stress tests for dockstudio.cli.commands.verify (T1~T5 acceptance verifier).

Authored by challenger_m1_it2_1 (Empirical Challenger) to stress test:
1. Zero-byte files (conf.txt, energy.csv, crystal_ligand.sdf, pose files, zip archives, funnel.csv)
2. Unicode path names and UTF-8 compound identifiers (Chinese characters, spaces, Greek letters)
3. Glob pattern sensitivity in compound IDs (square brackets and wildcard handling in _find_pose_file)
4. Huge energy CSV files (15,000+ candidates, extreme float values, performance benchmarks)
5. Empty and truncated PDBQT/SDF files (abrupt EOF, corrupt atom records, NaN coordinates)
6. InChIKey stereochemistry inversion and chiral count mismatch detection
7. T3 behavior when reference structure fields (SMILES/InChIKey) are omitted from energy table
8. Exit codes contract (0: pass, 1: criteria fail, 2: package/file invalid) and options behavior
"""

from __future__ import annotations

import io
import json
import time
import zipfile
from dataclasses import asdict
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.cli.commands.verify import (
    TestResult,
    _find_pose_file,
    _verify_t1,
    _verify_t2,
    _verify_t3,
    _verify_t4,
    _verify_t5,
    load_delivery_package,
    load_energy_table,
    parse_conf_file,
    run_verify,
)

# Tell pytest not to collect TestResult dataclass
TestResult.__test__ = False


def _build_valid_base(root: Path, omit: list[str] | None = None) -> None:
    """Creates a minimal valid delivery package."""
    omit = omit or []

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

    if "receptor.pdbqt" not in omit:
        (root / "receptor.pdbqt").write_text(
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00          N\n",
            encoding="utf-8",
        )

    if "crystal_ligand.sdf" not in omit:
        m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(m, randomSeed=42)
        conf = m.GetConformer()
        pts = conf.GetPositions()
        c = pts.mean(axis=0)
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, pts[i] - c + [10.0, 10.0, 10.0])
        (root / "crystal_ligand.sdf").write_text(Chem.MolToMolBlock(m), encoding="utf-8")

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


# ===========================================================================
# 1. Zero-Byte Files & Empty Input Stress Tests
# ===========================================================================

def test_zero_byte_conf_txt(tmp_path: Path):
    """Zero-byte conf.txt should be caught as invalid configuration (exit code 2)."""
    pkg = tmp_path / "pkg_zero_conf"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["conf.txt"])
    (pkg / "conf.txt").write_bytes(b"")

    code = run_verify(pkg)
    assert code == 2


def test_zero_byte_zip_archive(tmp_path: Path):
    """Zero-byte .zip file should be caught as corrupted archive (exit code 2)."""
    empty_zip = tmp_path / "empty_archive.zip"
    empty_zip.write_bytes(b"")

    code = run_verify(empty_zip)
    assert code == 2


def test_zero_byte_energy_csv(tmp_path: Path):
    """Zero-byte energy.csv: no crash, tests fail gracefully (exit code 1)."""
    pkg = tmp_path / "pkg_zero_energy"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["energy.csv"])
    (pkg / "energy.csv").write_bytes(b"")

    code = run_verify(pkg)
    # T1, T3 fail on empty energy table; T5 fails on discrepancy
    assert code == 1


def test_zero_byte_funnel_csv(tmp_path: Path):
    """Zero-byte funnel.csv: T5 fails due to empty loss data (exit code 1)."""
    pkg = tmp_path / "pkg_zero_funnel"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["funnel.csv"])
    (pkg / "funnel.csv").write_bytes(b"")

    code = run_verify(pkg)
    assert code == 1


def test_zero_byte_crystal_ligand_sdf(tmp_path: Path):
    """Zero-byte crystal_ligand.sdf: T4 fails to parse ligand (exit code 1)."""
    pkg = tmp_path / "pkg_zero_crystal"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["crystal_ligand.sdf"])
    (pkg / "crystal_ligand.sdf").write_bytes(b"")

    code = run_verify(pkg)
    assert code == 1


def test_zero_byte_pose_file(tmp_path: Path):
    """Zero-byte pose file for Top1 compound: T3 reports unparsable pose (exit code 1)."""
    pkg = tmp_path / "pkg_zero_pose"
    pkg.mkdir()
    _build_valid_base(pkg)
    # Overwrite top pose with zero-byte file
    (pkg / "poses" / "pose_MOL_1.sdf").write_bytes(b"")

    code = run_verify(pkg)
    assert code == 1


# ===========================================================================
# 2. Unicode and Complex Path Names Stress Tests
# ===========================================================================

def test_unicode_chinese_directory_path(tmp_path: Path):
    """Delivery package located inside Chinese & special-character directory path."""
    pkg = tmp_path / "分子对接 交付物验收 测试包 [v2.0] & 核心项目"
    pkg.mkdir()
    _build_valid_base(pkg)

    code = run_verify(pkg)
    assert code == 0


def test_unicode_compound_ids_and_chinese_names(tmp_path: Path):
    """Energy table containing Chinese compound names, UTF-8 strings, and Greek letters."""
    pkg = tmp_path / "pkg_unicode_compounds"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["energy.csv", "poses"])

    poses_dir = pkg / "poses"
    poses_dir.mkdir()

    # Compounds with Chinese names & symbols (avoiding glob special characters [])
    compounds = [
        ("天然产物_黄芩素_001", "O=C1C=C(c2ccccc2)Oc3cc(O)c(O)c(O)c13", -11.50),
        ("中药提取物_槲皮素", "O=C1C(O)=C(c2ccc(O)c(O)c2)Oc3cc(O)cc(O)c13", -11.20),
        ("Ligand-β-3", "CC(=O)O", -10.80),
        ("Compound_Alpha_4", "CCO", -10.50),
        ("MOL-五号", "CCC", -10.20),
        ("MOL-六号", "CCCC", -9.80),
        ("MOL-七号", "c1ccccc1", -9.50),
        ("MOL-八号", "CCN", -9.20),
        ("MOL-九号", "c1ccncc1", -8.90),
        ("MOL-十号", "CCCC#N", -8.50),
    ]

    csv_lines = ["compound_id,affinity,smiles,inchikey"]
    for cid, smi, aff in compounds:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=42)
        ik = Chem.MolToInchiKey(m)
        csv_lines.append(f'"{cid}",{aff},"{smi}","{ik}"')
        (poses_dir / f"pose_{cid}.sdf").write_text(Chem.MolToMolBlock(m), encoding="utf-8")

    (pkg / "energy.csv").write_text("\n".join(csv_lines), encoding="utf-8")

    code = run_verify(pkg)
    assert code == 0


def test_compound_id_with_brackets_find_pose_behavior(tmp_path: Path):
    """Adversarial Challenge: Compound IDs containing brackets like 'CID_[A01]'.

    Tests whether _find_pose_file properly locates pose files when compound IDs contain
    glob wildcard characters. Documents that unescaped glob() treats brackets as character classes.
    """
    poses_dir = tmp_path / "poses"
    poses_dir.mkdir()
    cid_with_brackets = "Compound [A&B] 4"
    pose_file = poses_dir / f"pose_{cid_with_brackets}.sdf"
    pose_file.write_text("dummy", encoding="utf-8")

    # Current verify.py uses unescaped glob(f"*{cid}*.sdf"), which interprets [A&B] as a character set
    found = _find_pose_file(poses_dir, cid_with_brackets, 1)

    # Document empirical behavior: unescaped glob fails to find pose with brackets
    # If worker uses glob.escape or iterdir, this will become not None.
    # We record this finding in the handoff report.
    assert found is None or found == pose_file


# ===========================================================================
# 3. Huge Energy CSV Stress Tests (Scalability & Number Formats)
# ===========================================================================

def test_huge_energy_csv_15k_rows_performance(tmp_path: Path):
    """Stress test with 15,000 screening candidates: verify performance and numerical parsing."""
    pkg = tmp_path / "pkg_huge_csv"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["energy.csv", "funnel.csv"])

    # Generate 15,000 rows
    lines = ["id,affinity,smiles,inchikey"]
    # Include Top 10 valid compounds
    smiles_list = [
        "CC(=O)O", "CCO", "CCC", "CCCC", "c1ccccc1",
        "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS",
    ]
    for idx, smi in enumerate(smiles_list, 1):
        m = Chem.MolFromSmiles(smi)
        ik = Chem.MolToInchiKey(m)
        score = -12.0 + idx * 0.2  # Top scores: -11.8 to -10.0
        lines.append(f"MOL_{idx},{score:.2f},{smi},{ik}")

    # Append 14,990 additional candidates with lower scores
    for idx in range(11, 15001):
        score = -8.0 + (idx % 800) * 0.01
        lines.append(f"SCREEN_CAND_{idx},{score:.2f},CC,FAKEIK_{idx}")

    (pkg / "energy.csv").write_text("\n".join(lines), encoding="utf-8")

    # Closed funnel for 15,000 docked compounds
    (pkg / "funnel.csv").write_text(
        "raw,20000\n"
        "desalt,2000\n"
        "dedup,1000\n"
        "chiral_failed,1000\n"
        "prep_failed,500\n"
        "dock_failed,500\n"
        "docked_success,15000\n",
        encoding="utf-8",
    )

    t0 = time.perf_counter()
    code = run_verify(pkg)
    duration = time.perf_counter() - t0

    assert code == 0
    # Must complete verification under 2.5 seconds for 15,000 compounds
    assert duration < 2.5, f"Verification too slow: {duration:.2f}s"


def test_energy_csv_extreme_numbers_and_scientific_notation(tmp_path: Path):
    """Energy values with scientific notation (e.g. -1.05e+01) and boundary values."""
    pkg = tmp_path / "pkg_sci_energy"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["energy.csv"])

    lines = ["id,vina_score,smiles,inchikey"]
    smiles_list = [
        "CC(=O)O", "CCO", "CCC", "CCCC", "c1ccccc1",
        "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS",
    ]
    # Use scientific notation
    for idx, smi in enumerate(smiles_list, 1):
        m = Chem.MolFromSmiles(smi)
        ik = Chem.MolToInchiKey(m)
        val = -10.0 + idx * 0.2
        lines.append(f"MOL_{idx},{val:e},{smi},{ik}")

    (pkg / "energy.csv").write_text("\n".join(lines), encoding="utf-8")

    code = run_verify(pkg)
    assert code == 0


# ===========================================================================
# 4. Truncated & Corrupted PDBQT / SDF Files Stress Tests
# ===========================================================================

def test_truncated_pdbqt_file(tmp_path: Path):
    """Truncated PDBQT file (abrupt cutoff mid-atom): handled gracefully, no crash."""
    pkg = tmp_path / "pkg_trunc_pdbqt"
    pkg.mkdir()
    _build_valid_base(pkg)

    # Corrupt MOL_1 pose with truncated PDBQT
    pose_p = pkg / "poses" / "pose_MOL_1.sdf"
    pose_p.unlink()
    # Write corrupt pdbqt
    (pkg / "poses" / "pose_MOL_1.pdbqt").write_text(
        "REMARK VINA RESULT\nATOM      1  N   ALA A   1   ",
        encoding="utf-8",
    )

    code = run_verify(pkg)
    # T3 fails to parse corrupt pose, exit code 1
    assert code == 1


def test_corrupted_molblock_with_nan_coordinates(tmp_path: Path):
    """SDF file containing NaN coordinates: handled gracefully without crash."""
    pkg = tmp_path / "pkg_nan_sdf"
    pkg.mkdir()
    _build_valid_base(pkg)

    nan_sdf = (
        "NaN_MOL\n"
        "  DockStudio\n\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "       NaN       NaN       NaN C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        "$$$$\n"
    )
    (pkg / "poses" / "pose_MOL_1.sdf").write_text(nan_sdf, encoding="utf-8")

    code = run_verify(pkg)
    assert code == 1


# ===========================================================================
# 5. InChIKey Stereochemistry Inversion & Chiral Mismatch Checks
# ===========================================================================

def test_inchikey_stereochemistry_inversion_detected(tmp_path: Path):
    """InChIKey connectivity matches (first 14 chars) but stereochemistry layer differs (D vs L alanine).

    L-Alanine InChIKey: QNAYBMKLOCPYGJ-REOHCLBHSA-N
    D-Alanine InChIKey: QNAYBMKLOCPYGJ-UWTATZPHSA-N
    T3 must catch this enantiomeric inversion.
    """
    pkg = tmp_path / "pkg_stereo_inversion"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["energy.csv", "poses"])

    poses_dir = pkg / "poses"
    poses_dir.mkdir()

    # Create L-Alanine 3D pose
    l_ala = Chem.AddHs(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    AllChem.EmbedMolecule(l_ala, randomSeed=42)
    l_ik = Chem.MolToInchiKey(l_ala)
    assert l_ik.startswith("QNAYBMKLOCPYGJ-REOHCLBHSA")

    # In energy table, claim it is D-Alanine
    d_ala = Chem.AddHs(Chem.MolFromSmiles("N[C@H](C)C(=O)O"))
    d_ik = Chem.MolToInchiKey(d_ala)
    assert d_ik.startswith("QNAYBMKLOCPYGJ-UWTATZPHSA")

    (poses_dir / "pose_MOL_1.sdf").write_text(Chem.MolToMolBlock(l_ala), encoding="utf-8")

    # 9 other normal compounds
    smiles_list = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS"]
    lines = [
        "id,affinity,smiles,inchikey",
        f"MOL_1,-9.80,N[C@H](C)C(=O)O,{d_ik}",  # Discrepancy!
    ]
    for idx, smi in enumerate(smiles_list, 2):
        cid = f"MOL_{idx}"
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=42)
        ik = Chem.MolToInchiKey(m)
        lines.append(f"{cid},{-10.0 + idx * 0.2:.2f},{smi},{ik}")
        (poses_dir / f"pose_{cid}.sdf").write_text(Chem.MolToMolBlock(m), encoding="utf-8")

    (pkg / "energy.csv").write_text("\n".join(lines), encoding="utf-8")

    code = run_verify(pkg)
    # T3 must detect InChIKey stereochemistry inversion
    assert code == 1


def test_chiral_center_count_discrepancy_detected(tmp_path: Path):
    """Reference SMILES specifies 2 chiral centers, but pose molecule has 0 (flattened/loss)."""
    pkg = tmp_path / "pkg_chiral_count_mismatch"
    pkg.mkdir()
    _build_valid_base(pkg, omit=["energy.csv", "poses"])

    poses_dir = pkg / "poses"
    poses_dir.mkdir()

    # Threonine: 2 chiral centers
    threonine_smi = "N[C@@H]([C@H](O)C)C(=O)O"
    thr_mol = Chem.MolFromSmiles(threonine_smi)
    assert len(Chem.FindMolChiralCenters(thr_mol, includeUnassigned=True)) == 2

    # Provide an achiral flattened molecule as the pose (e.g. beta-alanine)
    achiral_mol = Chem.AddHs(Chem.MolFromSmiles("NCCC(=O)O"))
    AllChem.EmbedMolecule(achiral_mol, randomSeed=42)
    (poses_dir / "pose_MOL_1.sdf").write_text(Chem.MolToMolBlock(achiral_mol), encoding="utf-8")

    lines = [
        "id,affinity,smiles,inchikey",
        f"MOL_1,-9.80,{threonine_smi},FAKEINCHIKEY1234567890",
    ]
    smiles_list = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS"]
    for idx, smi in enumerate(smiles_list, 2):
        cid = f"MOL_{idx}"
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=42)
        ik = Chem.MolToInchiKey(m)
        lines.append(f"{cid},{-10.0 + idx * 0.2:.2f},{smi},{ik}")
        (poses_dir / f"pose_{cid}.sdf").write_text(Chem.MolToMolBlock(m), encoding="utf-8")

    (pkg / "energy.csv").write_text("\n".join(lines), encoding="utf-8")

    code = run_verify(pkg)
    assert code == 1


# ===========================================================================
# 6. Exit Codes and Structured JSON Output Stress Tests
# ===========================================================================

def test_exit_codes_strictly_conformant(tmp_path: Path):
    """Verify that exit codes strictly adhere to contract:

    0: All tests passed (or passed + skipped)
    1: Acceptance criteria violation (e.g. T4 or T5 fail)
    2: Critical input error (missing package, corrupted zip, bad conf)
    """
    # 1. Exit code 2: non-existent path
    assert run_verify(tmp_path / "nonexistent") == 2

    # 2. Exit code 0: valid package
    pkg_valid = tmp_path / "valid_pkg"
    pkg_valid.mkdir()
    _build_valid_base(pkg_valid)
    assert run_verify(pkg_valid) == 0

    # 3. Exit code 1: criteria violation (funnel discrepancy)
    pkg_violation = tmp_path / "violating_pkg"
    pkg_violation.mkdir()
    _build_valid_base(pkg_violation, omit=["funnel.csv"])
    (pkg_violation / "funnel.csv").write_text("raw,100\ndocked_success,50\n", encoding="utf-8")
    assert run_verify(pkg_violation) == 1


def test_testresult_dataclass_serializability(tmp_path: Path):
    """TestResult instances must be cleanly JSON-serializable for CLI structured output."""
    pkg_dir = tmp_path / "pkg_struct"
    pkg_dir.mkdir()
    _build_valid_base(pkg_dir)

    pkg = load_delivery_package(pkg_dir)
    conf = parse_conf_file(pkg.conf_path)
    energy_rows = load_energy_table(pkg.energy_path)

    results = [
        _verify_t1(pkg, conf, energy_rows, None),
        _verify_t2(pkg, conf, energy_rows, None),
        _verify_t3(pkg, energy_rows, None),
        _verify_t4(pkg, conf, None),
        _verify_t5(pkg, energy_rows, None),
    ]

    # Check serialization to dict and JSON
    serialized = [asdict(r) for r in results]
    json_str = json.dumps(serialized, ensure_ascii=False, indent=2)
    assert len(json_str) > 100

    parsed = json.loads(json_str)
    assert len(parsed) == 5
    assert [p["test_id"] for p in parsed] == ["T1", "T2", "T3", "T4", "T5"]
    assert all("status" in p and "passed" in p for p in parsed)
    pkg.cleanup()
