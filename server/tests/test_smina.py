"""smina engine: real e2e when binary is available, plus parser unit tests."""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

import pytest

from dockstudio.engines.smina import SminaEngine, parse_smina_out


def _engine_available() -> bool:
    try:
        ok, _ = SminaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


SMINA_STDOUT = textwrap.dedent(
    """\
    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -8.9       0.000      0.000
    2       -7.6       10.280     15.977
    3       -7.4       2.112      4.447
    """
)

SMINA_PDBQT = textwrap.dedent(
    """\
    MODEL 1
    REMARK minimizedAffinity -8.88390541
    REMARK minimizedRMSD -1
    ROOT
    ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
    ENDROOT
    ENDMDL
    MODEL 2
    REMARK minimizedAffinity -7.60000000
    REMARK minimizedRMSD -1
    ROOT
    ATOM      1  C   UNL     1       2.000   2.000   3.000  1.00  0.00    +0.000 C
    ENDROOT
    ENDMDL
    """
)


def test_parse_remark_affinities_with_stdout_table(tmp_path: Path) -> None:
    out = tmp_path / "out.pdbqt"
    out.write_text(SMINA_PDBQT, encoding="utf-8")
    poses = parse_smina_out(out, SMINA_STDOUT)
    assert len(poses) == 2
    assert poses[0].affinity == pytest.approx(-8.88390541)
    # rmsd paired from stdout table rows
    assert poses[1].rmsd_lb == pytest.approx(10.28)


def test_parse_stdout_only_fallback(tmp_path: Path) -> None:
    # out.pdbqt missing/empty (e.g. OpenBabel write failed) -> table fallback
    out = tmp_path / "out.pdbqt"
    out.write_text("", encoding="utf-8")
    poses = parse_smina_out(out, SMINA_STDOUT)
    assert len(poses) == 3
    assert poses[2].affinity == pytest.approx(-7.4)
    assert poses[2].rmsd_ub == pytest.approx(4.447)


def test_parse_garbage_raises(tmp_path: Path) -> None:
    from dockstudio.models import DockingError

    out = tmp_path / "out.pdbqt"
    out.write_text("nothing\n", encoding="utf-8")
    with pytest.raises(DockingError):
        parse_smina_out(out, "junk\n")


def test_manifest() -> None:
    eng = SminaEngine()
    assert {"vina", "vinardo", "ad4_scoring"} == set(eng.manifest.scoring_functions)


@pytest.mark.skipif(not _engine_available(), reason="smina binary not available")
def test_smina_real_dock(tmp_path: Path, indinavir_smiles: str, hsg_pdb_bytes: bytes) -> None:
    from dockstudio.models import Box, DockParams
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    lig = prepare_ligand_from_smiles(indinavir_smiles, "indinavir", tmp_path / "lig")
    rec = prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", tmp_path / "rec")
    eng = SminaEngine()
    mk1 = rec.cocrystal_ligands[0]
    result = eng.dock(
        Path(rec.pdbqt_path),
        Path(lig.pdbqt_path),
        Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=1, num_modes=3, cpu=max(1, (os.cpu_count() or 2) - 1), seed=42),
        tmp_path / "run",
    )
    assert len(result.poses) >= 1
    # same order of magnitude as vina's ~-10 kcal/mol on 1HSG/indinavir
    assert -14.0 <= result.poses[0].affinity <= -5.0
    assert result.engine_version.startswith("2020")
