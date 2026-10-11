"""Co-crystal reference comparison: heavy-atom RMSD between crystal and docked pose."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from ..logging import get_logger

log = get_logger("analysis.crystal")


def _read_sdf(path: Path) -> Chem.Mol:
    mol = Chem.MolFromMolBlock(path.read_text(encoding="utf-8"), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not read SDF {path}")
    return mol


def _strip_bond_types(mol: Chem.Mol) -> Chem.Mol:
    """Copy the mol with all bonds as SINGLE — RMSD is about geometry/topology,
    not bond order; mixed typing would break substructure matching."""
    mol = Chem.Mol(mol)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    return mol


def heavy_atom_rmsd(crystal_sdf: Path, pose_sdf: Path) -> float | None:
    """Heavy-atom RMSD after best-fit alignment (GetBestRMS substructure search
    after normalizing bond types, which is what makes a PDB-derived crystal SDF
    and a vina pose comparable)."""
    crystal = _read_sdf(crystal_sdf)
    pose = _read_sdf(pose_sdf)
    ref = _strip_bond_types(Chem.RemoveHs(crystal))
    probe = _strip_bond_types(Chem.RemoveHs(pose))
    try:
        rmsd = rdMolAlign.GetBestRMS(probe, ref)
        return round(float(rmsd), 3)
    except Exception as exc:  # noqa: BLE001
        log.warning("RMSD alignment failed: %s", exc)
        return None
