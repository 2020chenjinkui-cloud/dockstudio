"""TCM specific chemical operations: desalting, deglycosylation, conformer baseline, and PAINS/TCM alerts."""

from __future__ import annotations

import os
from typing import Literal

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, SaltRemover
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

try:
    from ..logging import get_logger
    log = get_logger("prep.tcm")
except Exception:
    import logging
    log = logging.getLogger("prep.tcm")

try:
    from ..models import LigandQCProfile
except ImportError:
    from pydantic import BaseModel, Field

    class LigandQCProfile(BaseModel):  # type: ignore[no-redef]
        """配体理化性质、天然产物特征与假阳性警示画像。"""
        heavy_atoms: int = Field(default=0, ge=0, description="重原子数")
        mw: float = Field(default=0.0, ge=0.0, description="分子量 (Da)")
        clogp: float = Field(default=0.0, description="计算脂水分配系数 cLogP")
        tpsa: float = Field(default=0.0, ge=0.0, description="拓扑极性表面积 (Å²)")
        hbd: int = Field(default=0, ge=0, description="氢键供体数")
        hba: int = Field(default=0, ge=0, description="氢键受体数")
        rotatable_bonds: int = Field(default=0, ge=0, description="可旋转键数")
        fsp3: float = Field(default=0.0, ge=0.0, le=1.0, description="sp3 碳杂化占比")
        qed: float = Field(default=0.0, ge=0.0, le=1.0, description="定量类药性评分 QED")
        inchikey: str = Field(default="", description="标准化 InChIKey")
        murcko_scaffold: str = Field(default="", description="Bemis-Murcko 骨架 SMILES")
        is_glycoside: bool = Field(default=False, description="是否为糖苷类化合物")
        parent_glycoside_id: str | None = Field(default=None, description="衍生苷元对应的原形糖苷 ID")
        derivative_type: Literal["original", "aglycone", "stereoisomer"] = Field(
            default="original", description="分子形态类型"
        )
        min_energy_kcal: float | None = Field(default=None, description="自由态 MMFF94s 全局最低能量基准 (kcal/mol)")
        pains_alerts: list[str] = Field(default_factory=list, description="PAINS 结构警示列表")
        pains_names: list[str] = Field(default_factory=list, description="PAINS 结构警示列表 (别名)")
        brenk_names: list[str] = Field(default_factory=list, description="Brenk 结构警示列表")
        tcm_warnings: list[str] = Field(
            default_factory=list,
            description="中药专属假阳性警示标签，如 Catechol_Polyphenol, Quinone_Redox, Michael_Acceptor, Tannin_Polyphenol"
        )
        tcm_alerts: list[str] = Field(default_factory=list, description="中药专属假阳性警示标签 (别名)")


ALLOWED_ORGANIC_ELEMENTS = frozenset(
    {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"}
)
ALLOWED_COUNTERIONS = frozenset({"Na", "K", "Li", "Ca", "Mg"})

_PAT_CATECHOL = Chem.MolFromSmarts("c1([OH])c([OH])cccc1")
_PAT_QUINONE_P = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")
_PAT_QUINONE_O = Chem.MolFromSmarts("O=C1C(=O)C=CC=C1")
_PAT_QUINONE_P_GEN = Chem.MolFromSmarts("O=C1~[#6]~[#6]~C(=O)~[#6]~[#6]~1")
_PAT_QUINONE_O_GEN = Chem.MolFromSmarts("O=C1~[#6](=O)~[#6]~[#6]~[#6]~[#6]~1")
_PAT_MICHAEL = Chem.MolFromSmarts("[C;!R]=[C;!R]-[C,S](=[O])")
_PAT_PHENOL_OH = Chem.MolFromSmarts("[c][OH]")

_FILTER_CATALOG: FilterCatalog | None = None


def _get_filter_catalog() -> FilterCatalog:
    global _FILTER_CATALOG
    if _FILTER_CATALOG is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        _FILTER_CATALOG = FilterCatalog(params)
    return _FILTER_CATALOG


def desalt_and_filter(
    mol: Chem.Mol,
    min_heavy: int = 5,
    max_heavy: int = 120,
) -> tuple[Chem.Mol | None, str | None]:
    """Strips non-covalent salts, waters, and solvates.

    Rejects molecules containing non-organic or heavy metal elements outside the whitelist
    {H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I}.
    Rejects molecules with heavy atom counts outside [min_heavy, max_heavy].

    Returns (clean_mol, failure_reason).
    """
    if mol is None:
        return None, "无效的分子对象 (None)"

    # 1. Reject transition/heavy metals and disallowed elements upfront (before desalting)
    raw_disallowed = sorted(
        {atom.GetSymbol() for atom in mol.GetAtoms()}
        - (ALLOWED_ORGANIC_ELEMENTS | ALLOWED_COUNTERIONS)
    )
    if raw_disallowed:
        elem_str = ", ".join(raw_disallowed)
        return None, f"含非有机/金属元素，不适用标准分子对接力场 (含: {elem_str})"

    # 2. Strip salts and crystal waters
    try:
        remover = SaltRemover.SaltRemover()
        stripped = remover.StripMol(mol, dontRemoveEverything=True)
    except Exception as exc:
        log.debug("SaltRemover exception: %s", exc)
        stripped = mol

    # 3. Choose largest organic fragment
    try:
        chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        chosen = chooser.choose(stripped)
    except Exception as exc:
        return None, f"脱盐与母核选择失败: {exc}"

    if chosen is None or chosen.GetNumAtoms() == 0:
        return None, "脱盐后无有效分子结构"

    # 4. Organic element whitelist verification
    disallowed = sorted({atom.GetSymbol() for atom in chosen.GetAtoms()} - ALLOWED_ORGANIC_ELEMENTS)
    if disallowed:
        elem_str = ", ".join(disallowed)
        return None, f"含非有机/金属元素，不适用标准分子对接力场 (含: {elem_str})"

    # 5. Heavy atom range verification
    n_heavy = chosen.GetNumHeavyAtoms()
    if n_heavy < min_heavy or n_heavy > max_heavy:
        return None, f"重原子数超出范围 [{min_heavy}, {max_heavy}] (当前: {n_heavy})"

    return chosen, None


def identify_sugar_rings(mol: Chem.Mol) -> list[set[int]]:
    """Identifies pyranose and furanose sugar rings via ring topology.

    A valid sugar ring must:
    - Be a 5-membered (furanose) or 6-membered (pyranose) non-aromatic saturated ring.
    - Contain exactly one ring oxygen atom and carbon for the remainder.
    - Bear at least 3 exocyclic heteroatom/carbon-oxygen substituents (or 2 for pentoses).
    """
    sugar_rings: list[set[int]] = []
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(a.GetIsAromatic() for a in atoms):
            continue
        # Ring bonds must all be single bonds
        bonds = [
            mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)])
            for i in range(len(ring))
        ]
        if any(b.GetBondType() != Chem.BondType.SINGLE for b in bonds):
            continue
        symbols = [a.GetSymbol() for a in atoms]
        if symbols.count("O") != 1 or symbols.count("C") != len(ring) - 1:
            continue

        hetero_subst = 0
        for a_idx in ring:
            atom = mol.GetAtomWithIdx(a_idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring:
                    if nbr.GetSymbol() in ("O", "N"):
                        hetero_subst += 1
                        break
                    elif nbr.GetSymbol() == "C" and any(
                        nn.GetSymbol() in ("O", "N") for nn in nbr.GetNeighbors()
                    ):
                        hetero_subst += 1
                        break

        min_subst = 3 if len(ring) == 6 else 2
        if hetero_subst >= min_subst:
            sugar_rings.append(set(ring))

    return sugar_rings


def deglycosylate_tcm_mol(
    mol: Chem.Mol,
    min_aglycone_heavy: int = 10,
) -> tuple[Chem.Mol | None, list[str]]:
    """Hydrolyzes O-/N-/ester-glycosidic linkages to recover the aglycone.

    - Protects direct C-glycosides (e.g. Puerarin, Vitexin) with direct C-C bonds without cleavage.
    - Cleaves bonds between sugar anomeric carbon and exocyclic heteroatoms (-O-, -N-).
    - Caps aglycone side with -OH, -NH2, or -COOH.
    - Requires aglycone heavy atoms >= min_aglycone_heavy.

    Returns (aglycone_mol, list_of_sugar_names).
    """
    if mol is None:
        return None, []

    sugar_rings = identify_sugar_rings(mol)
    if not sugar_rings:
        return None, []

    sugar_names = ["pyranose" if len(r) == 6 else "furanose" for r in sugar_rings]

    bonds_to_cleave = []
    for s_ring in sugar_rings:
        for a_idx in s_ring:
            atom = mol.GetAtomWithIdx(a_idx)
            if atom.GetSymbol() != "C":
                continue
            # Anomeric carbon is directly bonded to the ring oxygen
            is_anomeric = any(
                nbr.GetIdx() in s_ring and nbr.GetSymbol() == "O"
                for nbr in atom.GetNeighbors()
            )
            if not is_anomeric:
                continue

            for nbr in atom.GetNeighbors():
                # Glycosidic linkage: exocyclic heteroatom bridge (O or N)
                if nbr.GetIdx() not in s_ring and nbr.GetSymbol() in ("O", "N"):
                    # Must be connected to another carbon outside this sugar ring
                    other_carbons = [
                        nn
                        for nn in nbr.GetNeighbors()
                        if nn.GetIdx() != a_idx and nn.GetSymbol() == "C"
                    ]
                    if other_carbons:
                        b = mol.GetBondBetweenAtoms(a_idx, nbr.GetIdx())
                        bonds_to_cleave.append(b.GetIdx())

    if not bonds_to_cleave:
        # C-glycoside or unattached sugars
        return None, []

    bonds_to_cleave = list(set(bonds_to_cleave))
    frag_mol = Chem.FragmentOnBonds(
        mol, bonds_to_cleave, dummyLabels=[(0, 0)] * len(bonds_to_cleave)
    )
    frags = Chem.GetMolFrags(frag_mol, asMols=True)

    true_aglycones = []
    for f in frags:
        # Identify attachment atoms that are aromatic nitrogen in a 5-membered ring
        # (e.g. N9 of purine/adenosine/cordycepin, indole N, imidazole N)
        for a in f.GetAtoms():
            if a.GetAtomicNum() == 0:
                for nbr in a.GetNeighbors():
                    if nbr.GetSymbol() == "N" and nbr.GetIsAromatic() and nbr.IsInRingSize(5):
                        nbr.SetBoolProp("_AromaticN_Capping", True)

        capped = Chem.DeleteSubstructs(f, Chem.MolFromSmiles("*"))

        for a in capped.GetAtoms():
            if a.HasProp("_AromaticN_Capping"):
                a.SetNumExplicitHs(1)
                a.ClearProp("_AromaticN_Capping")

        try:
            Chem.SanitizeMol(capped)
        except Exception:
            continue
        if capped.GetNumHeavyAtoms() >= min_aglycone_heavy:
            cand_s_rings = []
            for r in capped.GetRingInfo().AtomRings():
                if len(r) in (5, 6):
                    atoms = [capped.GetAtomWithIdx(i) for i in r]
                    if not any(a.GetIsAromatic() for a in atoms) and [
                        a.GetSymbol() for a in atoms
                    ].count("O") == 1:
                        cand_s_rings.append(r)
            has_aromatic = any(a.GetIsAromatic() for a in capped.GetAtoms())
            cand_total_rings = capped.GetRingInfo().NumRings()
            if (
                has_aromatic
                or cand_total_rings > len(cand_s_rings)
                or (cand_total_rings == 0 and capped.GetNumHeavyAtoms() >= min_aglycone_heavy)
            ):
                true_aglycones.append(capped)

    if not true_aglycones:
        return None, sugar_names

    true_aglycones.sort(key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
    return true_aglycones[0], sugar_names


def generate_conformers_and_min_energy(
    mol: Chem.Mol,
    num_confs: int = 10,
    seed: int = 42,
) -> tuple[Chem.Mol, float]:
    """Generates an ETKDGv3 conformer ensemble and calculates the MMFF94s global minimum energy.

    - Uses useMacrocycleTorsions=True and useSmallRingTorsions=True.
    - Minimizes all conformers using MMFF94s (fallback to UFF).
    - Returns (best_single_conformer_mol, global_min_energy_kcal).
    """
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useMacrocycleTorsions = True
    params.useSmallRingTorsions = True
    params.randomSeed = seed

    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params))
    if not cids:
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params))
        if not cids:
            if AllChem.EmbedMolecule(mol, params) == 0:
                cids = [0]
            else:
                raise ValueError("3D 构象嵌入失败")

    has_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
    energies: dict[int, float] = {}

    if has_mmff:
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        for cid in cids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            if ff is not None:
                ff.Minimize(maxIts=500)
                energies[cid] = ff.CalcEnergy()
            else:
                energies[cid] = float("inf")
    else:
        for cid in cids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            if ff is not None:
                ff.Minimize(maxIts=500)
                energies[cid] = ff.CalcEnergy()
            else:
                energies[cid] = float("inf")

    best_cid = min(energies, key=energies.get)
    best_energy = energies[best_cid]

    # Return Mol containing strictly the global minimum conformer
    best_mol = Chem.Mol(mol)
    best_mol.RemoveAllConformers()
    best_mol.AddConformer(mol.GetConformer(best_cid), assignId=True)
    return best_mol, float(best_energy)


def scan_pains_and_tcm_alerts(mol: Chem.Mol) -> LigandQCProfile:
    """Scans FilterCatalog (PAINS_A/B/C, Brenk) and 4 TCM false-positive warning alerts.

    TCM Alerts:
    - Catechol_Polyphenol: c1([OH])c([OH])cccc1
    - Quinone_Redox: O=C1C=CC(=O)C=C1 / ortho-quinone / fused quinone
    - Michael_Acceptor: [C;!R]=[C;!R]-[C,S](=[O])
    - Tannin_Polyphenol: MW > 500 and phenolic OHs count >= 5

    Returns populated LigandQCProfile.
    """
    catalog = _get_filter_catalog()
    pains_matches = catalog.GetMatches(mol)
    pains_alerts = sorted({m.GetDescription() for m in pains_matches})
    brenk_names = [d for d in pains_alerts if "brenk" in d.lower()]
    pains_names = [d for d in pains_alerts if d not in brenk_names]

    tcm_warnings: list[str] = []
    if mol.HasSubstructMatch(_PAT_CATECHOL):
        tcm_warnings.append("Catechol_Polyphenol")
    if (
        mol.HasSubstructMatch(_PAT_QUINONE_P)
        or mol.HasSubstructMatch(_PAT_QUINONE_O)
        or mol.HasSubstructMatch(_PAT_QUINONE_P_GEN)
        or mol.HasSubstructMatch(_PAT_QUINONE_O_GEN)
    ):
        tcm_warnings.append("Quinone_Redox")
    if mol.HasSubstructMatch(_PAT_MICHAEL):
        tcm_warnings.append("Michael_Acceptor")
    if Descriptors.MolWt(mol) > 500 and len(mol.GetSubstructMatches(_PAT_PHENOL_OH)) >= 5:
        tcm_warnings.append("Tannin_Polyphenol")

    sugar_rings = identify_sugar_rings(mol)
    is_glyco = len(sugar_rings) > 0

    scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    inchikey = Chem.MolToInchiKey(mol)

    return LigandQCProfile(
        heavy_atoms=mol.GetNumHeavyAtoms(),
        mw=round(Descriptors.MolWt(mol), 3),
        clogp=round(Descriptors.MolLogP(mol), 3),
        tpsa=round(Descriptors.TPSA(mol), 3),
        hbd=Descriptors.NumHDonors(mol),
        hba=Descriptors.NumHAcceptors(mol),
        rotatable_bonds=Descriptors.NumRotatableBonds(mol),
        fsp3=round(Descriptors.FractionCSP3(mol), 3),
        qed=round(QED.qed(mol), 4),
        inchikey=inchikey,
        murcko_scaffold=scaffold,
        is_glycoside=is_glyco,
        pains_alerts=pains_alerts,
        pains_names=pains_names,
        brenk_names=brenk_names,
        tcm_warnings=tcm_warnings,
        tcm_alerts=tcm_warnings,
    )


def standardize_tautomer_and_charge(mol: Chem.Mol) -> Chem.Mol:
    """Normalizes tautomers and standardizes protonation state."""
    if mol is None:
        return mol
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(mol)
    te = rdMolStandardize.TautomerEnumerator()
    try:
        canonical = te.Canonicalize(uncharged)
        return canonical
    except Exception as exc:
        log.warning("Tautomer canonicalization failed: %s", exc)
        return uncharged


def enumerate_stereoisomers_tcm(
    mol: Chem.Mol,
    max_isomers: int = 4,
) -> tuple[list[Chem.Mol], list[str]]:
    """Enumerates stereoisomers for unassigned chiral centers.

    If 1-2 centers are undefined, enumerates up to max_isomers.
    If >2 centers are undefined, avoids combinatorial explosion and returns [mol] with warning.
    """
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    unassigned = [c for c in centers if c[1] == "?"]
    warnings: list[str] = []

    if not unassigned:
        return [mol], []

    if len(unassigned) > 2:
        warnings.append(
            f"含 {len(unassigned)} 个未定义手性中心，已保留默认构型避免组合爆炸，建议核查天然绝对构型"
        )
        return [mol], warnings

    opts = StereoEnumerationOptions(maxIsomers=max_isomers, onlyUnassigned=True)
    try:
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        return (isomers if isomers else [mol]), []
    except Exception as exc:
        warnings.append(f"手性中心枚举失败: {exc}")
        return [mol], warnings
