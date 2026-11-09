"""Ligand preparation: RDKit 3D embedding + meeko PDBQT writing with TCM optimization."""

from __future__ import annotations

import uuid
from pathlib import Path
from typing import Literal

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from ..logging import get_logger
from ..models import LigandQCProfile, PreparationError, PreparedLigand
from .tcm import (
    desalt_and_filter,
    generate_conformers_and_min_energy,
    scan_pains_and_tcm_alerts,
)

log = get_logger("prep.ligand")


def _embed_3d(mol: Chem.Mol) -> tuple[Chem.Mol, float]:
    """Generates ETKDGv3 conformers, MMFF94s minimization, and returns (best_mol, min_energy)."""
    return generate_conformers_and_min_energy(mol)


def _mol_to_pdbqt(mol: Chem.Mol, name: str) -> str:
    """Prepares and writes ligand PDBQT using Meeko."""
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


def _write_ligand(
    mol: Chem.Mol,
    source: Literal["file", "smiles"],
    name: str | None,
    out_dir: Path,
    qc_profile: LigandQCProfile | None = None,
    min_energy_kcal: float | None = None,
    derivative_type: Literal["original", "aglycone", "stereoisomer"] = "original",
    parent_id: str | None = None,
    ligand_id: str | None = None,
) -> PreparedLigand:
    lid = ligand_id or uuid.uuid4().hex[:8]
    lig_dir = out_dir / lid
    lig_dir.mkdir(parents=True, exist_ok=True)
    pdbqt = _mol_to_pdbqt(mol, name or lid)
    pdbqt_path = lig_dir / "ligand.pdbqt"
    pdbqt_path.write_text(pdbqt, encoding="utf-8")
    mol_no_h = Chem.RemoveHs(mol)
    return PreparedLigand(
        id=lid,
        source=source,
        pdbqt_path=str(pdbqt_path),
        n_rotatable_bonds=Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_no_h),
        smiles=Chem.MolToSmiles(mol_no_h),
        name=name,
        qc_profile=qc_profile,
        min_energy_kcal=min_energy_kcal,
        derivative_type=derivative_type,
        parent_id=parent_id,
    )


def prepare_ligand_from_smiles(
    smiles: str,
    name: str | None,
    out_dir: Path,
    derivative_type: Literal["original", "aglycone", "stereoisomer"] = "original",
    parent_id: str | None = None,
    ligand_id: str | None = None,
    min_heavy: int = 3,
    max_heavy: int = 120,
) -> PreparedLigand:
    """Build a ligand PDBQT from a SMILES string with intelligent desalting,
    elemental/size filtering, MMFF94s conformer ensemble energy baseline,
    and PAINS/TCM risk profile."""
    smiles = (smiles or "").strip()
    if not smiles:
        raise PreparationError("SMILES 字符串为空，请粘贴有效的 SMILES。")
    if any(ch.isspace() for ch in smiles):
        smiles = smiles.split()[0]
        log.info("SMILES contained whitespace; truncated to %r", smiles)
    bad = [c for c in smiles if ord(c) > 127]
    if bad:
        raise PreparationError(f"SMILES 含非法字符 {bad[0]!r}（只接受 ASCII），请检查输入法。")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PreparationError(f"RDKit 无法解析 SMILES：{smiles[:60]!r}")

    # Detect genuine multi-component organic mixtures (excluding salts and solvates)
    if "." in smiles:
        try:
            from rdkit.Chem.SaltRemover import SaltRemover
            remover = SaltRemover()
            stripped = remover.StripMol(mol, dontRemoveEverything=True)
            frags = Chem.GetMolFrags(stripped, asMols=True)
            substantial_frags = [f for f in frags if f.GetNumHeavyAtoms() >= min_heavy]
            unique_species = {Chem.MolToSmiles(f) for f in substantial_frags}
            if len(unique_species) > 1:
                raise PreparationError("SMILES 含有多个独立组分（'.' 分隔），请先拆分，每次只提交单个组分。")
        except PreparationError:
            raise
        except Exception as exc:
            log.debug("Multi-component check bypassed due to error: %s", exc)

    # Intelligent desalting & organic element / size whitelist
    clean_mol, fail_reason = desalt_and_filter(mol, min_heavy=min_heavy, max_heavy=max_heavy)
    if clean_mol is None:
        raise PreparationError(fail_reason or "配体脱盐与元素/尺寸过滤失败")

    # Generate PAINS & TCM False-Positive Warning Profile
    qc_profile = scan_pains_and_tcm_alerts(clean_mol)
    qc_profile.derivative_type = derivative_type
    qc_profile.parent_glycoside_id = parent_id

    # Conformer Ensemble Minimization & Global Minimum Energy (E_0)
    mol_3d, min_energy = generate_conformers_and_min_energy(clean_mol)
    qc_profile.min_energy_kcal = min_energy

    return _write_ligand(
        mol_3d,
        source="smiles",
        name=name,
        out_dir=out_dir,
        qc_profile=qc_profile,
        min_energy_kcal=min_energy,
        derivative_type=derivative_type,
        parent_id=parent_id,
        ligand_id=ligand_id,
    )


def prepare_ligand_from_file(
    path: Path,
    out_dir: Path,
    derivative_type: Literal["original", "aglycone", "stereoisomer"] = "original",
    parent_id: str | None = None,
    ligand_id: str | None = None,
    min_heavy: int = 3,
    max_heavy: int = 120,
) -> PreparedLigand:
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
                mol.UpdatePropertyCache(strict=False)
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

    # Desalt & filter
    clean_mol, fail_reason = desalt_and_filter(mol, min_heavy=min_heavy, max_heavy=max_heavy)
    if clean_mol is None:
        raise PreparationError(fail_reason or f"配体脱盐与元素/尺寸过滤失败: {path.name}")

    qc_profile = scan_pains_and_tcm_alerts(clean_mol)
    qc_profile.derivative_type = derivative_type
    qc_profile.parent_glycoside_id = parent_id

    # 3D Coordinates & Conformer Baseline
    has_3d = clean_mol.GetNumConformers() > 0 and any(
        abs(clean_mol.GetConformer().GetAtomPosition(i).z) > 1e-4
        for i in range(clean_mol.GetNumAtoms())
    )
    if not has_3d or clean_mol.GetNumAtoms() != Chem.AddHs(clean_mol).GetNumAtoms():
        mol_3d, min_energy = generate_conformers_and_min_energy(clean_mol)
    else:
        _, min_energy = generate_conformers_and_min_energy(clean_mol)
        mol_3d = clean_mol

    qc_profile.min_energy_kcal = min_energy
    name = clean_mol.GetProp("_Name") if clean_mol.HasProp("_Name") else path.stem
    return _write_ligand(
        mol_3d,
        source="file",
        name=name or None,
        out_dir=out_dir,
        qc_profile=qc_profile,
        min_energy_kcal=min_energy,
        derivative_type=derivative_type,
        parent_id=parent_id,
        ligand_id=ligand_id,
    )
