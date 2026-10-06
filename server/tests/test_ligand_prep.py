"""Unit tests for ligand preparation."""

from __future__ import annotations

from pathlib import Path

import pytest

from dockstudio.models import PreparationError
from dockstudio.prep.ligand import prepare_ligand_from_smiles


def test_prepare_indinavir(indinavir_smiles: str, tmp_path: Path) -> None:
    lig = prepare_ligand_from_smiles(indinavir_smiles, "indinavir", tmp_path)
    assert lig.id and len(lig.id) == 8
    assert lig.n_rotatable_bonds is not None and lig.n_rotatable_bonds > 0
    assert lig.smiles
    pdbqt = Path(lig.pdbqt_path)
    assert pdbqt.is_file()
    text = pdbqt.read_text(encoding="utf-8")
    assert "TORSDOF" in text  # meeko torsion footer proves valid ligand PDBQT


def test_empty_smiles(tmp_path: Path) -> None:
    with pytest.raises(PreparationError):
        prepare_ligand_from_smiles("   ", None, tmp_path)


def test_invalid_smiles(tmp_path: Path) -> None:
    with pytest.raises(PreparationError):
        prepare_ligand_from_smiles("not_a_smiles(((", None, tmp_path)
