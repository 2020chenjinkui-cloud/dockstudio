"""Ligand preparation: RDKit 3D embedding + meeko PDBQT writing."""

from __future__ import annotations

import uuid
from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from ..logging import get_logger
from ..models import PreparationError, PreparedLigand

log = get_logger("prep.ligand")


def _embed_3d(mol: Chem.Mol) -> Chem.Mol:
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise PreparationError("failed to generate 3D conformer for ligand")
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)
    except Exception as exc:  # noqa: BLE001  # optimization failure is not fatal
        log.warning("force-field optimization failed: %s", exc)
    return mol


def _mol_to_pdbqt(mol: Chem.Mol, name: str) -> str:
    try:
        setups = MoleculePreparation().prepare(mol)
    except Exception as exc:
        raise PreparationError(f"meeko failed to prepare ligand: {exc}") from exc
    if not setups:
        raise PreparationError("meeko produced no molecule setup for ligand")
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok or not pdbqt:
        raise PreparationError(f"meeko PDBQT writing failed: {err}")
    return pdbqt


def _write_ligand(mol: Chem.Mol, source: str, name: str | None, out_dir: Path) -> PreparedLigand:
    ligand_id = uuid.uuid4().hex[:8]
    lig_dir = out_dir / ligand_id
    lig_dir.mkdir(parents=True, exist_ok=True)
    pdbqt = _mol_to_pdbqt(mol, name or ligand_id)
    pdbqt_path = lig_dir / "ligand.pdbqt"
    pdbqt_path.write_text(pdbqt, encoding="utf-8")
    return PreparedLigand(
        id=ligand_id,
        source=source,  # type: ignore[arg-type]
        pdbqt_path=str(pdbqt_path),
        n_rotatable_bonds=Chem.rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(mol)),
        smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),
        name=name,
    )


def prepare_ligand_from_smiles(smiles: str, name: str | None, out_dir: Path) -> PreparedLigand:
    """Build a ligand PDBQT from a SMILES string."""
    smiles = (smiles or "").strip()
    if not smiles:
        raise PreparationError("SMILES 字符串为空，请粘贴有效的 SMILES。")
    if "." in smiles:
        raise PreparationError("SMILES 含有多组分（'.' 分隔），请先拆分，每次只提交单个组分。")
    if any(ch.isspace() for ch in smiles):
        smiles = smiles.split()[0]
        log.info("SMILES contained whitespace; truncated to %r", smiles)
    bad = [c for c in smiles if ord(c) > 127]
    if bad:
        raise PreparationError(f"SMILES 含非法字符 {bad[0]!r}（只接受 ASCII），请检查输入法。")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PreparationError(f"RDKit 无法解析 SMILES：{smiles[:60]!r}")
    mol = _embed_3d(mol)
    return _write_ligand(mol, "smiles", name, out_dir)


def prepare_ligand_from_file(path: Path, out_dir: Path) -> PreparedLigand:
    """Build a ligand PDBQT from an .sdf/.mol/.mol2 file.

    Parsing goes through in-memory mol blocks so unicode/non-ASCII paths work
    (RDKit's file-based suppliers fail on such paths on Windows)."""
    suffix = path.suffix.lower()
    text = path.read_text(encoding="utf-8", errors="replace")
    mol: Chem.Mol | None = None
    if suffix == ".sdf":
        block = text.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
        if mol is not None:
            try:
                mol.UpdatePropertyCache()
                Chem.FastFindRings(mol)
                Chem.SanitizeMol(mol)
            except Exception as exc:
                raise PreparationError(f"RDKit could not sanitize {path.name}: {exc}") from exc
    elif suffix == ".mol":
        mol = Chem.MolFromMolBlock(text, sanitize=True, removeHs=False)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2Block(text, sanitize=True, removeHs=False)
    else:
        raise PreparationError(
            f"unsupported ligand file type: {suffix or '<none>'} (use .sdf/.mol/.mol2)"
        )
    if mol is None:
        raise PreparationError(f"RDKit could not parse ligand file: {path.name}")
    has_3d = mol.GetNumConformers() > 0 and any(
        abs(mol.GetConformer().GetAtomPosition(i).z) > 1e-4 for i in range(mol.GetNumAtoms())
    )
    if not has_3d:
        mol = _embed_3d(mol)
    elif mol.GetNumAtoms() != Chem.AddHs(mol).GetNumAtoms():
        mol = _embed_3d(mol)  # missing hydrogens
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else path.stem
    return _write_ligand(mol, "file", name or None, out_dir)
