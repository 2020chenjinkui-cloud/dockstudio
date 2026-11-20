"""PoseCheck: clashing SDF must report an error; a real 1HSG pose must be ok."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from dockstudio.analysis.posecheck import check_pose
from dockstudio.models import Box

_CLASHING_SDF = """collision
  DockStudio

  4  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1000    0.3000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.0000    9.0000    9.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""

_SMALL_SDF = (
    "e\n  DockStudio\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0\n"
    "    1.4500    0.0000    0.0000 C   0  0\n"
    "    1.9000    1.2700    0.0000 O   0  0\n"
    "  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n$$$$\n"
)


def _engine_available() -> bool:
    from dockstudio.engines.vina import VinaEngine

    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


@pytest.fixture(scope="module")
def pose_files(tmp_path_factory, indinavir_smiles: str, hsg_pdb_bytes: bytes):
    """One quick vina run reused by pose-dependent tests in this file."""
    if not _engine_available():
        pytest.skip("vina binary not available")
    from dockstudio.analysis.export import export_poses
    from dockstudio.engines.vina import VinaEngine
    from dockstudio.models import DockParams
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    base = tmp_path_factory.mktemp("posecheck")
    lig = prepare_ligand_from_smiles(indinavir_smiles, "indinavir", base / "lig")
    rec = prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", base / "rec")
    eng = VinaEngine()
    mk1 = rec.cocrystal_ligands[0]
    box = Box(center=mk1.center, size=(22.5, 22.5, 22.5))
    run = base / "run"
    # exhaustiveness=1 was too low once the fixture started using the *real*
    # stereoisomer: measured on the P4, exh=1 returns +8.408 kcal/mol with an
    # internal_steric_clash (a physically impossible pose), while exh>=4 returns
    # -11.27 with no errors.  The old flat SMILES happened to land in the pocket by
    # luck, which is why the low setting looked fine before.
    eng.dock(
        Path(rec.pdbqt_path),
        Path(lig.pdbqt_path),
        box,
        DockParams(exhaustiveness=4, num_modes=3, cpu=max(1, (os.cpu_count() or 2) - 1), seed=42),
        run,
    )
    return rec, export_poses(run / "out.pdbqt", run), box


def test_clashing_pose_is_error(tmp_path: Path) -> None:
    sdf = tmp_path / "pose_1.sdf"
    sdf.write_text(_CLASHING_SDF, encoding="utf-8")
    check = check_pose(sdf, use_posebusters=False)
    assert check.source == "builtin"
    assert not check.ok
    assert any(i.severity == "error" for i in check.issues)
    assert any("碰撞" in i.message or i.rule == "clash" for i in check.issues)


def test_valid_sdf_ok(tmp_path: Path) -> None:
    sdf = tmp_path / "pose_1.sdf"
    sdf.write_text(_SMALL_SDF, encoding="utf-8")
    check = check_pose(
        sdf, box=Box(center=(0.6, 0.5, 0.0), size=(20, 20, 20)), use_posebusters=False
    )
    assert check.ok
    assert check.issues == []


def test_out_of_box_warning(tmp_path: Path) -> None:
    sdf = tmp_path / "pose_2.sdf"
    sdf.write_text(
        "far\n  DockStudio\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "   50.0000   50.0000   50.0000 C   0  0\nM  END\n$$$$\n",
        encoding="utf-8",
    )
    check = check_pose(sdf, box=Box(center=(0, 0, 0), size=(20, 20, 20)), use_posebusters=False)
    assert any(i.rule == "out_of_box" and i.severity == "warning" for i in check.issues)
    assert check.ok  # warning only, not an error
    assert check.pose_index == 2


def test_1hsg_real_pose_ok(pose_files) -> None:
    _, poses, box = pose_files
    check = check_pose(poses[0], box=box)
    # real indinavir pose has no geometry errors (warnings tolerable)
    assert not any(i.severity == "error" for i in check.issues), check.issues


def test_docking_score_is_physically_plausible(pose_files) -> None:
    """A positive affinity means the pose is worse than not binding at all.

    This guard exists because exhaustiveness=1 used to hand back +8.408 kcal/mol with
    an internal_steric_clash and nothing complained: the pose check only looked at the
    geometry of the SDF, never at the score.  A dock reported as successful must be
    favourable, otherwise downstream ranking silently promotes a non-binder.
    """
    _, _, _ = pose_files  # ensures the shared vina run happened and stayed valid
    from dockstudio.analysis.export import export_poses  # noqa: F401  (fixture parity)

    import tempfile
    from pathlib import Path as _P

    from dockstudio.engines.vina import VinaEngine
    from dockstudio.models import Box, DockParams
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    fixtures = _P(__file__).parent / "fixtures"
    td = _P(tempfile.mkdtemp())
    rec = prepare_receptor((fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", td / "rec")
    lig = prepare_ligand_from_smiles(
        (fixtures / "indinavir_smiles.txt").read_text().strip(), "indinavir", td / "lig"
    )
    mk1 = rec.cocrystal_ligands[0]
    res = VinaEngine().dock(
        _P(rec.pdbqt_path),
        _P(lig.pdbqt_path),
        Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=4, num_modes=3, cpu=2, seed=42),
        td / "run",
    )
    top = res.top_affinity
    assert top is not None, "no pose produced"
    assert top < 0.0, (
        f"top affinity {top:.3f} kcal/mol is not favourable -- the pose is worse than "
        "not binding, so it must not be reported as a successful dock"
    )
