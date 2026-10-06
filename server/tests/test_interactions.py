"""Unit tests for pose interaction analysis (builtin geometry + PLIP)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from dockstudio.analysis.export import export_poses
from dockstudio.analysis.interactions import _run_builtin, analyze_pose
from dockstudio.engines.vina import VinaEngine
from dockstudio.models import Box, DockParams
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.prep.receptor import prepare_receptor


def _engine_available() -> bool:
    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


@pytest.fixture(scope="module")
def pose_files(tmp_path_factory, indinavir_smiles: str, hsg_pdb_bytes: bytes):
    """One quick vina run reused by all tests in this module."""
    if not _engine_available():
        pytest.skip("vina binary not available")
    base = tmp_path_factory.mktemp("inter")
    lig = prepare_ligand_from_smiles(indinavir_smiles, "indinavir", base / "lig")
    rec = prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", base / "rec")
    eng = VinaEngine()
    mk1 = rec.cocrystal_ligands[0]
    run = base / "run"
    eng.dock(
        Path(rec.pdbqt_path),
        Path(lig.pdbqt_path),
        Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=2, num_modes=3, cpu=max(1, (os.cpu_count() or 2) - 1), seed=42),
        run,
    )
    poses = export_poses(run / "out.pdbqt", run)
    protein = Path(rec.pdbqt_path).parent / "receptor_protein.pdb"
    return protein, poses


def test_builtin_detects_hbonds(pose_files) -> None:
    protein, poses = pose_files
    interactions = _run_builtin(protein, poses[0])
    assert interactions, "no interactions detected at all"
    hbonds = [i for i in interactions if i.type == "hbond"]
    assert hbonds, f"no hbonds among {len(interactions)} interactions"
    # active-site catalytic aspartates (ASP25/ASP25') contact the ligand
    residues = {i.residue for i in interactions}
    assert any(r.startswith("ASP25") or r.startswith("ASP 25") for r in residues), residues  # noqa: PIE810
    for i in interactions:
        assert ":" in i.residue
        if i.distance is not None:
            assert 0 < i.distance < 6.0


def test_analyze_pose_runs_and_caches(pose_files, tmp_path) -> None:
    protein, poses = pose_files
    analysis = analyze_pose(protein, poses[0], tmp_path, use_plip=False)
    assert analysis.source == "builtin"
    assert analysis.pose_index == 1
    assert analysis.interactions


def test_analyze_pose_prefers_plip_when_available(pose_files, tmp_path) -> None:
    pytest.importorskip("plip")
    protein, poses = pose_files
    analysis = analyze_pose(protein, poses[0], tmp_path, use_plip=True)
    assert analysis.source == "plip"
    assert analysis.interactions
