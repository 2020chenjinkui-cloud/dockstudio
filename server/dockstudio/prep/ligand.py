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
    standardize_tautomer_and_charge,
)

log = get_logger("prep.ligand")


def _embed_3d(mol: Chem.Mol) -> tuple[Chem.Mol, float]:
    """Generates ETKDGv3 conformers, MMFF94s minimization, and returns (best_mol, min_energy)."""
    return generate_conformers_and_min_energy(mol)


def _relax_existing_3d(mol: Chem.Mol, max_its: int = 500) -> tuple[Chem.Mol, float]:
    """P2-1: relax an already-3D molecule in place (on a copy) instead of re-embedding.

    When an SDF/MOL already carries explicit hydrogens and true 3D coordinates there is
    nothing to sample: a single MMFF94s (UFF fallback) minimisation on the copy gives a
    physically consistent energy baseline in a fraction of the time, and the relaxation
    is exactly the same force field used for ``min_energy_kcal`` downstream.
    """
    work = Chem.Mol(mol)
    if work.GetNumConformers() == 0:
        raise PreparationError("分子不含 3D 构象，无法直接弛豫")
    force_field = None
    if AllChem.MMFFHasAllMoleculeParams(work):
        props = AllChem.MMFFGetMoleculeProperties(work, mmffVariant="MMFF94s")
        if props is not None:
            force_field = AllChem.MMFFGetMoleculeForceField(work, props, confId=0)
    if force_field is None:
        force_field = AllChem.UFFGetMoleculeForceField(work, confId=0)
    if force_field is None:
        raise PreparationError("无法为该分子建立 MMFF94s/UFF 力场")
    force_field.Minimize(maxIts=max_its)
    return work, float(force_field.CalcEnergy())


def _prepare_core_mol(raw_mol: Chem.Mol, min_heavy: int, max_heavy: int, source_label: str) -> Chem.Mol:
    """Desalt + tautomer/charge standardisation shared by the SMILES and file paths."""
    clean_mol, fail_reason = desalt_and_filter(raw_mol, min_heavy=min_heavy, max_heavy=max_heavy)
    if clean_mol is None:
        if fail_reason == "multi_component_mixture":
            n_components = len(Chem.GetMolFrags(raw_mol, asMols=True)) if raw_mol is not None else 0
            raise PreparationError(
                f"检测到多个独立组分混合物（{n_components} 个不同组分），"
                "请先拆分为单分子后再准备。"
            )
        raise PreparationError(fail_reason or f"配体脱盐与元素/尺寸过滤失败{source_label}")
    try:
        standardized = standardize_tautomer_and_charge(clean_mol)
    except Exception as exc:  # pragma: no cover - defensive
        log.warning("standardisation failed for %s: %s", source_label or "SMILES", exc)
        return clean_mol
    # Tautomer canonicalisation may rebuild the molecule without conformers; carry the
    # input 3D coordinates over so P2-1 (relax instead of re-embed) still applies.
    if (
        standardized.GetNumConformers() == 0
        and clean_mol.GetNumConformers() > 0
        and standardized.GetNumAtoms() == clean_mol.GetNumAtoms()
    ):
        try:
            standardized.AddConformer(Chem.Conformer(clean_mol.GetConformer(0)), assignId=True)
        except Exception as exc:  # pragma: no cover - defensive
            log.debug("could not carry over 3D conformer: %s", exc)
    return standardized


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

    # Intelligent desalting (organic salt counter-ions are stripped; genuine
    # multi-component active mixtures are rejected by desalt_and_filter with
    # reason == "multi_component_mixture") + tautomer/charge standardisation.
    clean_mol = _prepare_core_mol(mol, min_heavy=min_heavy, max_heavy=max_heavy, source_label="")

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

    # Desalt (+ organic counter-ion stripping) & tautomer/charge standardisation
    clean_mol = _prepare_core_mol(
        mol, min_heavy=min_heavy, max_heavy=max_heavy, source_label=f": {path.name}"
    )

    qc_profile = scan_pains_and_tcm_alerts(clean_mol)
    qc_profile.derivative_type = derivative_type
    qc_profile.parent_glycoside_id = parent_id

    # 3D Coordinates & Conformer Baseline (P2-1)
    has_3d = clean_mol.GetNumConformers() > 0 and any(
        abs(clean_mol.GetConformer().GetAtomPosition(i).z) > 1e-4
        for i in range(clean_mol.GetNumAtoms())
    )
    has_all_hydrogens = clean_mol.GetNumAtoms() == Chem.AddHs(clean_mol).GetNumAtoms()
    if has_3d and has_all_hydrogens:
        # Already a true 3D structure: relax it once instead of embedding 10 new
        # conformers and throwing them away (P2-1 double work + energy mismatch).
        mol_3d, min_energy = _relax_existing_3d(clean_mol)
    else:
        mol_3d, min_energy = generate_conformers_and_min_energy(clean_mol)

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
