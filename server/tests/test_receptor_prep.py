"""Unit tests for receptor preparation on 1HSG."""

from __future__ import annotations

from pathlib import Path

import pytest

from dockstudio.models import PreparationError
from dockstudio.prep.receptor import prepare_receptor


def test_prepare_1hsg(hsg_pdb_bytes: bytes, tmp_path: Path) -> None:
    rec = prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", tmp_path)
    assert rec.chains == ["A", "B"]
    assert Path(rec.pdbqt_path).is_file()
    assert Path(rec.pdbqt_path).stat().st_size > 0
    assert Path(rec.original_path).is_file()

    mk1 = [c for c in rec.cocrystal_ligands if c.resname == "MK1"]
    assert len(mk1) == 1
    c = mk1[0]
    assert c.n_atoms >= 5
    # MK1 centroid near the known 1HSG active site
    assert c.center[0] == pytest.approx(13.07, abs=0.5)
    assert c.center[1] == pytest.approx(22.47, abs=0.5)
    assert c.center[2] == pytest.approx(5.56, abs=0.5)

    # waters must be stripped from the cleaned PDB
    clean = Path(rec.original_path).parent / "receptor_clean.pdb"
    assert clean.is_file()
    assert "HOH" not in clean.read_text(encoding="utf-8")


def test_mmcif_rejected(tmp_path: Path) -> None:
    with pytest.raises(PreparationError, match="mmCIF"):
        prepare_receptor(b"data_x", "x.cif", tmp_path)


def test_no_atom_records(tmp_path: Path) -> None:
    with pytest.raises(PreparationError):
        prepare_receptor(b"HEADER    nothing here\nEND\n", "x.pdb", tmp_path)
