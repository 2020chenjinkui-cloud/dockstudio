"""TCM specific chemical operations: desalting, deglycosylation, conformer baseline, and PAINS/TCM alerts."""

from __future__ import annotations

import os
from typing import Literal

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
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

# ---------------------------------------------------------------------------
# Pharmaceutical counter-ions / inorganic ions / crystallization solvents (P1-2)
# ---------------------------------------------------------------------------
# The pharmaceutical salt/solvate forms in real TCM & drug libraries are built from
# a small, well-known family of acids and solvents.  All entries below are stored as
# *neutral, stereo-free canonical SMILES keys* so that the charged salt form
# (e.g. tartrate dianion, mesylate anion) and the free acid map onto one key.

#: Mineral acids / organic acids used as pharmaceutical counter-ions.
_ACID_COUNTERION_SOURCE: tuple[str, ...] = (
    # inorganic mineral acids: sulfate, phosphate, nitrate, carbonate
    "O=S(=O)(O)O",
    "O=P(O)(O)O",
    "O=[N+]([O-])O",
    "O=C(O)O",
    # aliphatic / aromatic organic acids
    "CC(=O)O",                     # acetic
    "O=C(O)C(F)(F)F",              # trifluoroacetic
    "CS(=O)(=O)O",                 # methanesulfonic (mesylate)
    "CCS(=O)(=O)O",                # ethanesulfonic (esylate)
    "O=S(=O)(O)c1ccccc1",          # benzenesulfonic (besylate)
    "Cc1ccc(S(=O)(=O)O)cc1",       # p-toluenesulfonic (tosylate)
    "O=C(O)C(=O)O",                # oxalic
    "O=C(O)CC(=O)O",               # malonic
    "O=C(O)CCC(=O)O",              # succinic
    "O=C(O)C=CC(=O)O",             # maleic / fumaric
    "O=C(O)CC(O)C(=O)O",           # malic
    "O=C(O)C(O)C(O)C(=O)O",        # tartaric
    "O=C(O)CC(O)(CC(=O)O)C(=O)O",  # citric
    "CC(O)C(=O)O",                 # lactic
    "O=C(O)C(O)C(O)C(O)C(O)CO",    # gluconic
)

#: Crystallization solvents frequently found as solvates / hydrates.
_SOLVENT_SOURCE: tuple[str, ...] = (
    "O",           # water
    "CO",          # methanol
    "CCO",         # ethanol
    "CC(C)O",      # isopropanol
    "CC(=O)C",     # acetone
    "CC#N",        # acetonitrile
    "CS(=O)C",     # DMSO
    "CN(C)C=O",    # DMF
    "ClCCl",       # dichloromethane
    "ClC(Cl)Cl",   # chloroform
)

#: Inorganic ions & halides.
_INORGANIC_ION_SOURCE: tuple[str, ...] = (
    "[Na+]", "[K+]", "[Li+]", "[Ca+2]", "[Mg+2]", "[Zn+2]", "[NH4+]",
    "[F-]", "[Cl-]", "[Br-]", "[I-]", "Cl", "Br", "I",
)

#: Fragments at or above this heavy-atom count are "substantial" organic
#: components when deciding whether a multi-fragment input is a real mixture.
MIN_SUBSTANTIAL_FRAGMENT_HEAVY = 3

_STANDARDIZER_UNCHARGER = rdMolStandardize.Uncharger()


def _neutral_smiles_key(mol: Chem.Mol | None) -> str | None:
    """Canonical, stereo-free, uncharged SMILES key used for whitelist matching."""
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    work = Chem.Mol(mol)
    if any(atom.GetAtomicNum() == 1 for atom in work.GetAtoms()):
        try:
            work = Chem.RemoveHs(work)
        except Exception as exc:  # pragma: no cover - defensive
            log.debug("RemoveHs failed on fragment: %s", exc)
    try:
        work = _STANDARDIZER_UNCHARGER.uncharge(work)
    except Exception as exc:  # pragma: no cover - defensive
        log.debug("Uncharger failed on fragment: %s", exc)
    try:
        return Chem.MolToSmiles(work, isomericSmiles=False)
    except Exception:  # pragma: no cover - defensive
        return None


def _keys_of(smiles_source: tuple[str, ...]) -> frozenset[str]:
    keys: set[str] = set()
    for smiles in smiles_source:
        parsed = Chem.MolFromSmiles(smiles)
        if parsed is None:  # pragma: no cover - static table is valid
            log.debug("counter-ion/solvent SMILES unparsable: %s", smiles)
            continue
        key = _neutral_smiles_key(parsed)
        if key:
            keys.add(key)
    return frozenset(keys)


_ACID_COUNTERION_KEYS: frozenset[str] = _keys_of(_ACID_COUNTERION_SOURCE)
_SOLVENT_KEYS: frozenset[str] = _keys_of(_SOLVENT_SOURCE)
_INORGANIC_ION_KEYS: frozenset[str] = _keys_of(_INORGANIC_ION_SOURCE)

#: Public contract (HVTS_FIX_CONTRACT.md §3.1): neutral, stereo-free canonical
#: SMILES of all pharmaceutical counter-ions, inorganic ions and crystal solvents.
PHARMACEUTICAL_COUNTERIONS_AND_SOLVENTS: frozenset[str] = frozenset(
    _ACID_COUNTERION_KEYS | _SOLVENT_KEYS | _INORGANIC_ION_KEYS
)


def _is_inorganic_fragment(frag_mol: Chem.Mol) -> bool:
    """True for inorganic ions / mineral acids (no carbon skeleton) or <= 2 heavy atoms."""
    if frag_mol is None:
        return True
    if frag_mol.GetNumHeavyAtoms() <= 2:
        return True
    return not any(atom.GetSymbol() == "C" for atom in frag_mol.GetAtoms())


def is_counterion_or_solvent(frag_mol: Chem.Mol) -> bool:
    """Whether a single fragment is a pharmaceutical counter-ion, an inorganic
    ion / mineral acid, a common crystallization solvent, or has <= 2 heavy atoms."""
    if frag_mol is None or frag_mol.GetNumAtoms() == 0:
        return True
    if _is_inorganic_fragment(frag_mol):
        return True
    key = _neutral_smiles_key(frag_mol)
    if key is None:  # pragma: no cover - defensive
        return False
    return key in PHARMACEUTICAL_COUNTERIONS_AND_SOLVENTS

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


def _select_core_fragment(mol: Chem.Mol) -> tuple[Chem.Mol | None, str | None]:
    """Split a (possibly multi-fragment) molecule into its single bioactive core.

    Policy (pharmaceutical salt/solvate aware):
    1. inorganic ions / mineral acids / halides / water are always dropped;
    2. pharmaceutical acid counter-ions (tartrate, citrate, maleate, mesylate, ...)
       are dropped whenever a parent organic core survives;
    3. neutral crystallization solvents (EtOH, iPrOH, acetone, MeCN, DMSO, DMF, ...)
       are dropped only when they are clearly smaller than the main component;
    4. if >= 2 distinct substantial organic fragments survive, the input is a genuine
       multi-component active mixture -> ``(None, "multi_component_mixture")``.
    """
    try:
        frags = [f for f in Chem.GetMolFrags(mol, asMols=True) if f.GetNumAtoms() > 0]
    except Exception as exc:  # pragma: no cover - defensive
        log.debug("fragment split failed: %s", exc)
        return mol, None
    if len(frags) <= 1:
        return (frags[0] if frags else mol), None

    # (1) inorganic ions / mineral acids / water / halides
    organics = [f for f in frags if not _is_inorganic_fragment(f)]
    if not organics:
        # Degenerate input made only of inorganic/solvate pieces: keep them so the
        # heavy-atom range check below can report an informative failure.
        organics = list(frags)

    # (2)+(3) pharmaceutical acid counter-ions and small neutral solvents
    if len(organics) > 1:
        main_heavy = max(f.GetNumHeavyAtoms() for f in organics)
        kept: list[Chem.Mol] = []
        for frag in organics:
            heavy = frag.GetNumHeavyAtoms()
            key = _neutral_smiles_key(frag)
            if key is not None and key in _ACID_COUNTERION_KEYS and heavy < main_heavy:
                continue
            if key is not None and key in _SOLVENT_KEYS and 2 * heavy < main_heavy:
                continue
            kept.append(frag)
        if kept:
            organics = kept

    # (4) genuine multi-component active mixture detection (stereo-aware dedup)
    unique: dict[str, Chem.Mol] = {}
    for frag in organics:
        try:
            species_key = Chem.MolToSmiles(frag)
        except Exception:  # pragma: no cover - defensive
            species_key = f"frag_{len(unique)}_{frag.GetNumAtoms()}_{frag.GetNumHeavyAtoms()}"
        if (
            species_key not in unique
            or frag.GetNumHeavyAtoms() > unique[species_key].GetNumHeavyAtoms()
        ):
            unique[species_key] = frag
    core_frags = list(unique.values())
    substantial = [
        f for f in core_frags if f.GetNumHeavyAtoms() >= MIN_SUBSTANTIAL_FRAGMENT_HEAVY
    ]
    if len(substantial) >= 2:
        log.info(
            "rejecting genuine multi-component mixture: %d substantial organic fragments",
            len(substantial),
        )
        return None, "multi_component_mixture"
    if not core_frags:
        return None, "脱盐后无有效分子结构"

    combined = core_frags[0]
    try:
        for extra in core_frags[1:]:
            combined = Chem.CombineMols(combined, extra)
        chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        chosen = chooser.choose(combined)
    except Exception as exc:
        return None, f"脱盐与母核选择失败: {exc}"
    if chosen is None or chosen.GetNumAtoms() == 0:
        return None, "脱盐后无有效分子结构"
    return chosen, None


def desalt_and_filter(
    mol: Chem.Mol,
    min_heavy: int = 5,
    max_heavy: int = 120,
) -> tuple[Chem.Mol | None, str | None]:
    """Strips non-covalent salts, crystal waters, solvates and organic counter-ions.

    Rejects molecules containing non-organic or heavy metal elements outside the whitelist
    {H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I}.
    Rejects molecules with heavy atom counts outside [min_heavy, max_heavy].
    Rejects genuine multi-component active mixtures with the failure reason literal
    ``"multi_component_mixture"`` (pharmaceutical organic salts are *not* mixtures:
    tartrate / citrate / maleate / mesylate / oxalate / acetate are stripped instead).

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

    # 2. Counter-ion / solvate stripping + genuine mixture rejection
    chosen, reason = _select_core_fragment(mol)
    if chosen is None:
        return None, reason or "脱盐后无有效分子结构"
    if chosen.GetNumAtoms() == 0:
        return None, "脱盐后无有效分子结构"

    # 3. Organic element whitelist verification
    disallowed = sorted({atom.GetSymbol() for atom in chosen.GetAtoms()} - ALLOWED_ORGANIC_ELEMENTS)
    if disallowed:
        elem_str = ", ".join(disallowed)
        return None, f"含非有机/金属元素，不适用标准分子对接力场 (含: {elem_str})"

    # 4. Heavy atom range verification
    n_heavy = chosen.GetNumHeavyAtoms()
    if n_heavy < min_heavy or n_heavy > max_heavy:
        return None, f"重原子数超出范围 [{min_heavy}, {max_heavy}] (当前: {n_heavy})"

    return chosen, None


def _ring_is_isolated(mol: Chem.Mol, ring: tuple[int, ...]) -> bool:
    """Condition A: the ring must be an isolated monocycle (no fusion/bridge/spiro).

    Any ring bond shared with another ring (fused/bridged systems) or any ring atom
    that also belongs to another ring (spiro) disqualifies the candidate.  This alone
    rules out flavonoid C-rings, catechin, podophyllotoxin, artemisinin, ginkgolides,
    tetrahydrofuran lignans and other caged polycycles.
    """
    ring_info = mol.GetRingInfo()
    n_ring_atoms = len(ring)
    for i in range(n_ring_atoms):
        bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % n_ring_atoms])
        if bond is None:
            return False
        if ring_info.NumBondRings(bond.GetIdx()) > 1:
            return False
    for atom_idx in ring:
        memberships = 0
        for candidate in ring_info.AtomRings():
            if atom_idx in candidate:
                memberships += 1
                if memberships > 1:
                    return False
    return True


def _count_sugar_polar_groups(mol: Chem.Mol, ring: set[int]) -> int:
    """Condition B: free -OH / -CH2OH / -C(=O)OH(-C(=O)O-) on ring carbons."""
    n_polar = 0
    for atom_idx in ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "C":
            continue
        polar = False
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in ring:
                continue
            if nbr.GetSymbol() == "O":
                heavy_nbrs = [x for x in nbr.GetNeighbors() if x.GetAtomicNum() > 1]
                if len(heavy_nbrs) == 1 and (
                    nbr.GetTotalNumHs() >= 1 or nbr.GetFormalCharge() < 0
                ):
                    polar = True  # free hydroxyl / alkoxide
                continue
            if nbr.GetSymbol() != "C":
                continue
            substituents = [x for x in nbr.GetNeighbors() if x.GetIdx() != atom_idx]
            if not substituents:
                continue
            double_oxy = []
            hydroxyl_oxy = []
            for sub in substituents:
                if sub.GetSymbol() != "O":
                    continue
                bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_oxy.append(sub)
                    continue
                heavy_sub = [x for x in sub.GetNeighbors() if x.GetAtomicNum() > 1]
                if len(heavy_sub) == 1 and (
                    sub.GetTotalNumHs() >= 1 or sub.GetFormalCharge() < 0
                ):
                    hydroxyl_oxy.append(sub)
            if double_oxy and hydroxyl_oxy and len(substituents) <= 2:
                polar = True  # carboxylic acid / carboxylate (uronic acid)
            elif len(substituents) == 1 and hydroxyl_oxy:
                polar = True  # hydroxymethyl
        if polar:
            n_polar += 1
    return n_polar


def _ring_has_anomeric_carbon(mol: Chem.Mol, ring: set[int], ring_oxygen: int) -> bool:
    """Ring carbon bonded to the ring oxygen *and* to an exocyclic O/N (acetal C1)."""
    for atom_idx in ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "C":
            continue
        if not any(nbr.GetIdx() == ring_oxygen for nbr in atom.GetNeighbors()):
            continue
        if any(
            nbr.GetIdx() not in ring and nbr.GetSymbol() in ("O", "N")
            for nbr in atom.GetNeighbors()
        ):
            return True
    return False


def _count_exocyclic_aryl_substituents(mol: Chem.Mol, ring: set[int]) -> int:
    """Condition C: ring carbons directly bearing an aromatic carbon substituent."""
    n_aryl = 0
    for atom_idx in ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "C":
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in ring:
                continue
            if nbr.GetSymbol() == "C" and nbr.GetIsAromatic():
                n_aryl += 1
                break
    return n_aryl


def identify_sugar_rings(mol: Chem.Mol) -> list[set[int]]:
    """Identifies pyranose and furanose sugar rings via topology + substituent pattern.

    A valid sugar ring must simultaneously satisfy:
    - be a 5-membered (furanose) or 6-membered (pyranose) saturated ring carrying
      exactly one ring oxygen;
    - **Condition A** be an isolated monocycle - no fused/bridged/spiro ring bond
      (excludes flavonoid C-rings, catechin, podophyllotoxin, artemisinin,
      ginkgolides and every caged lactone core);
    - **Condition B** carry >= 3 free -OH/-CH2OH/-COOH groups on ring carbons, or
      >= 2 such groups plus a genuine anomeric C1 (deoxy sugars such as rhamnose,
      and C-glycosides with >= 3 hydroxyls);
    - **Condition C** not bear >= 2 aromatic ring substituents (dibenzylbutyrolactone
      / tetrahydrofuran lignan cores are never sugars).
    """
    if mol is None or mol.GetNumAtoms() == 0:
        return []
    try:
        Chem.FastFindRings(mol)
    except Exception as exc:  # pragma: no cover - defensive
        log.debug("FastFindRings failed: %s", exc)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return []

    sugar_rings: list[set[int]] = []
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(a.GetIsAromatic() for a in atoms):
            continue
        # Ring bonds must all be single bonds (no enol/lactone/aromatic ring bonds)
        bonds = [
            mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)])
            for i in range(len(ring))
        ]
        if any(b is None or b.GetBondType() != Chem.BondType.SINGLE for b in bonds):
            continue
        symbols = [a.GetSymbol() for a in atoms]
        if symbols.count("O") != 1 or symbols.count("C") != len(ring) - 1:
            continue

        # Condition A: isolated monocycle only
        if not _ring_is_isolated(mol, ring):
            continue
        # Condition C: lignan core guard
        if _count_exocyclic_aryl_substituents(mol, ring) >= 2:
            continue

        ring_oxygen = ring[symbols.index("O")]
        ring_set = set(ring)
        n_polar = _count_sugar_polar_groups(mol, ring_set)
        has_anomeric = _ring_has_anomeric_carbon(mol, ring_set, ring_oxygen)
        # Condition B: pyranose needs >= 3 polar groups, furanose >= 2 + anomeric C1;
        # deoxy sugars / C-glycosides may trade the third hydroxyl for a C1 acetal.
        if n_polar >= 3 or (n_polar >= 2 and has_anomeric):
            sugar_rings.append(ring_set)

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


_PAT_CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def _apply_physiological_ph(mol: Chem.Mol, ph: float) -> Chem.Mol:
    """Opt-in carboxylate deprotonation correction for a physiological pH model.

    Kept deliberately minimal: only explicit carboxylic acids are deprotonated above
    their typical pKa (~4.8), which is enough to make pH-conditional runs reproducible
    without perturbing the neutral standard state used by default.
    """
    if ph < 5.5:
        return mol
    try:
        editable = Chem.RWMol(mol)
        for match in mol.GetSubstructMatches(_PAT_CARBOXYLIC_ACID):
            oxy_idx = match[2]
            oxygen = editable.GetAtomWithIdx(oxy_idx)
            if oxygen.GetFormalCharge() != 0:
                continue
            oxygen.SetFormalCharge(-1)
            oxygen.SetNumExplicitHs(0)
            oxygen.SetNoImplicit(True)
        corrected = editable.GetMol()
        Chem.SanitizeMol(corrected)
        return corrected
    except Exception as exc:  # pragma: no cover - defensive
        log.debug("physiological pH correction skipped: %s", exc)
        return mol


def standardize_tautomer_and_charge(
    mol: Chem.Mol,
    *,
    physiological_ph: float | None = None,
    max_tautomers: int = 50,
    max_transforms: int = 200,
) -> Chem.Mol:
    """Neutralises formal charges and picks the canonical tautomer (P0-4).

    1. ``rdMolStandardize.Uncharger`` removes ionisable charges while keeping
       permanent ones (nitro groups, quaternary ammonium / berberine-type alkaloids);
    2. ``TautomerEnumerator`` with bounded effort (``maxTautomers=50``,
       ``maxTransforms=200``) collapses keto/enol and lactam/lactim pairs onto one
       canonical tautomer, with sp3 stereochemistry preserved;
    3. ``physiological_ph=None`` (default) keeps the neutral standard state so the
       resulting InChIKey stays consistent with downstream Meeko PDBQT hydrogen
       typing; pass e.g. ``physiological_ph=7.4`` for an explicit pH 7.4 correction.

    Must be applied *before* InChIKey deduplication and *before* 3D embedding.
    """
    if mol is None:
        return mol
    work = Chem.Mol(mol)
    # Tautomer enumeration moves hydrogens between atoms and is only well defined for
    # implicit-H graphs; explicit-H inputs (3D SDF/MOL) must be hydrogen-suppressed
    # first, otherwise a migrated H can over-valence the acceptor atom.
    explicit_hydrogens = any(atom.GetAtomicNum() == 1 for atom in work.GetAtoms())
    if explicit_hydrogens:
        try:
            work = Chem.RemoveHs(work)
        except Exception as exc:  # pragma: no cover - defensive
            log.debug("RemoveHs failed before standardisation: %s", exc)
    try:
        work = rdMolStandardize.Uncharger().uncharge(work)
    except Exception as exc:
        log.warning("charge neutralisation failed: %s", exc)

    try:
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(int(max_tautomers))
        enumerator.SetMaxTransforms(int(max_transforms))
        enumerator.SetRemoveSp3Stereo(False)  # keep natural absolute configuration
        canonical = enumerator.Canonicalize(work)
        if canonical is not None and canonical.GetNumAtoms() > 0:
            work = canonical
    except Exception as exc:
        log.warning("Tautomer canonicalization failed: %s", exc)

    if physiological_ph is not None:
        work = _apply_physiological_ph(work, float(physiological_ph))

    if explicit_hydrogens:
        try:
            work = Chem.AddHs(work, addCoords=work.GetNumConformers() > 0)
        except Exception as exc:  # pragma: no cover - defensive
            log.debug("AddHs failed after standardisation: %s", exc)

    if mol.HasProp("_Name") and not work.HasProp("_Name"):
        try:
            work.SetProp("_Name", mol.GetProp("_Name"))
        except Exception:  # pragma: no cover - defensive
            pass
    return work


def enumerate_stereoisomers_tcm(
    mol: Chem.Mol,
    max_unspecified_centers: int = 3,
    max_isomers: int = 4,
) -> tuple[list[Chem.Mol], bool]:
    """Enumerates stereoisomers for centros with unassigned configuration.

    Returns ``(isomers, chiral_ok)``:

    * ``chiral_ok is False`` when the molecule carries more than
      ``max_unspecified_centers`` undefined stereocentres - the input is then
      unusable (the returned list is ``[mol]`` unchanged) and callers must reject
      the molecule (``library.py`` counts it as ``n_chiral_fail``);
    * otherwise the enumerated isomers (up to ``max_isomers``) are returned, with
      ``[mol]`` when nothing needed enumerating.
    """
    if mol is None:
        return [], False
    try:
        centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False
        )
    except Exception:  # pragma: no cover - defensive
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    unassigned = [idx for idx, tag in centers if tag == "?"]
    if len(unassigned) > max(0, int(max_unspecified_centers)):
        log.info(
            "%d undefined stereocentres (> %d) - chiral QC failed",
            len(unassigned),
            max_unspecified_centers,
        )
        return [mol], False
    if not unassigned:
        return [mol], True

    options = StereoEnumerationOptions(maxIsomers=int(max_isomers), onlyUnassigned=True)
    try:
        isomers = list(EnumerateStereoisomers(mol, options=options))
    except Exception as exc:
        log.warning("stereoisomer enumeration failed: %s", exc)
        return [mol], True
    return (isomers or [mol]), True
