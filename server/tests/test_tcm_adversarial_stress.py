"""Adversarial stress-test suite for Milestone 1: TCM Ligand Prep & PAINS Pipeline.

Independent empirical verification harness created by Challenger 1.
Covers 5 critical focus areas:
1. Multi-salt and complex ionic compounds (multi-counterions, zwitterions, ammonium salts).
2. Unusual glycosides & C-glycosides (Vitexin, Orientin, Chlorogenic acid, Ester-glycosides, N-glycosides).
3. Heavy metal & exotic element traps (Organo-As, Organo-Hg, Ferrocene, Gadolinium, Pt, Sn, Au).
4. Macrocycle torsion & conformer edge cases (14/15/16-membered macrolides, caged polycycles).
5. TCM warning flags (Ortho- vs para-quinones, disguised Michael acceptors, Tannin MW/OH boundary).
"""

from __future__ import annotations

import math
import pytest
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import Descriptors

from dockstudio.prep import tcm


# =============================================================================
# 1. Focus Area 1: Multi-Salt and Complex Ionic Compounds
# =============================================================================

@pytest.mark.parametrize(
    "name,smiles,expected_core_heavy",
    [
        (
            "Disodium phthalate",
            "[Na+].[Na+].[O-]C(=O)c1ccccc1C(=O)[O-]",
            12,
        ),
        (
            "Dipotassium tartrate",
            "[K+].[K+].[O-]C(=O)[C@@H](O)[C@H](O)C(=O)[O-]",
            10,
        ),
        (
            "Berberine sulfate (2:1 cation/sulfate)",
            "[O-]S(=O)(=O)[O-].COc1ccc2cc3[n+](cc2c1OC)CCc1c3ccc2c1OCO2.COc1ccc2cc3[n+](cc2c1OC)CCc1c3ccc2c1OCO2",
            25,
        ),
        (
            "Trisodium citrate",
            "[Na+].[Na+].[Na+].[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]",
            13,
        ),
    ],
)
def test_multisalt_counterion_stripping_and_core_retention(
    name: str, smiles: str, expected_core_heavy: int
) -> None:
    """Verifies that multi-counterion salts strip all inorganic counterions completely
    and retain the single largest bioactive organic scaffold."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse SMILES for {name}"

    clean_mol, err = tcm.desalt_and_filter(mol)
    assert err is None, f"{name} failed desalt: {err}"
    assert clean_mol is not None

    # Check heavy atom count matches the desalted organic core
    assert clean_mol.GetNumHeavyAtoms() == expected_core_heavy, (
        f"{name}: expected {expected_core_heavy} heavy atoms, got {clean_mol.GetNumHeavyAtoms()}"
    )

    # Check that inorganic metal/sulfate counterions are eliminated
    symbols = {a.GetSymbol() for a in clean_mol.GetAtoms()}
    assert "Na" not in symbols
    assert "K" not in symbols
    if "Berberine" in name:
        assert "S" not in symbols


@pytest.mark.parametrize(
    "name,smiles,expected_heavy",
    [
        ("L-Phenylalanine zwitterion", "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]", 12),
        ("L-Tryptophan zwitterion", "[NH3+][C@@H](Cc1c[nH]c2ccccc12)C(=O)[O-]", 15),
        ("Betaine", "C[N+](C)(C)CC(=O)[O-]", 8),
        ("Choline chloride", "[Cl-].C[N+](C)(C)CCO", 7),
        ("Tetrabutylammonium bromide", "[Br-].CCCC[N+](CCCC)(CCCC)CCCC", 17),
        ("Cetylpyridinium chloride", "[Cl-].CCCCCCCCCCCCCCCC[n+]1ccccc1", 22),
    ],
)
def test_zwitterions_and_quaternary_ammonium_salts(
    name: str, smiles: str, expected_heavy: int
) -> None:
    """Verifies that zwitterions and quaternary ammonium salts are properly handled:
    counterions are stripped, ionic cores pass filtering, 3D conformers generate with
    finite E_0, and valid Meeko PDBQT is produced."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse SMILES for {name}"

    clean_mol, err = tcm.desalt_and_filter(mol)
    assert err is None, f"{name} failed desalt: {err}"
    assert clean_mol is not None
    assert clean_mol.GetNumHeavyAtoms() == expected_heavy

    # Conformer baseline
    mol_3d, min_energy = tcm.generate_conformers_and_min_energy(clean_mol, num_confs=3, seed=42)
    assert math.isfinite(min_energy), f"{name}: E_0 is non-finite: {min_energy}"

    # Meeko PDBQT preparation
    setups = MoleculePreparation().prepare(mol_3d)
    assert setups, f"{name}: Meeko produced no setups"
    pdbqt, ok, perr = PDBQTWriterLegacy.write_string(setups[0])
    assert ok and pdbqt, f"{name}: Meeko PDBQT writing failed: {perr}"
    assert "TORSDOF" in pdbqt


# =============================================================================
# 2. Focus Area 2: Unusual Glycosides & C-Glycosides
# =============================================================================

@pytest.mark.parametrize(
    "name,smiles",
    [
        (
            "Vitexin (8-C-glucosyl-apigenin)",
            "O=c1cc(-c2ccc(O)cc2)oc2c(O)c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc12",
        ),
        (
            "Orientin (8-C-glucosyl-luteolin)",
            "O=c1cc(-c2ccc(O)c(O)c2)oc2c(O)c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc12",
        ),
    ],
)
def test_c_glycosides_vitexin_and_orientin_never_cleaved(name: str, smiles: str) -> None:
    """Verifies that C-glycosides (Vitexin, Orientin) possessing direct aromatic C-C
    bonds to the sugar unit are NEVER cleaved, preserving the intact natural compound."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse SMILES for {name}"

    sugar_rings = tcm.identify_sugar_rings(mol)
    assert len(sugar_rings) >= 1, f"{name}: expected sugar ring to be topologically detected"

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol)
    # Direct C-C bond must never be cleaved: aglycone must be None and sugars empty
    assert aglycone is None, f"{name}: C-glycoside was falsely cleaved!"
    assert sugars == [], f"{name}: sugars list should be empty when C-glycoside is not cleaved"


def test_chlorogenic_acid_non_sugar_retention() -> None:
    """Verifies that Chlorogenic acid (caffeoyl quinic acid ester) is recognized as a
    non-sugar carbocycle (quinic acid is a cyclohexane ring with 0 ring oxygens) and
    is NOT cleaved or shredded as a pseudo-glycoside."""
    chlorogenic_smi = "O=C(O)[C@@]1(O)C[C@@H](O)[C@H](OC(=O)/C=C/c2ccc(O)c(O)c2)[C@H](O)C1"
    mol = Chem.MolFromSmiles(chlorogenic_smi)
    assert mol is not None

    sugar_rings = tcm.identify_sugar_rings(mol)
    assert len(sugar_rings) == 0, "Quinic acid carbocycle must not be identified as a sugar ring"

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol)
    assert aglycone is None
    assert sugars == []


def test_true_ester_glycoside_cleavage() -> None:
    """Verifies that true ester-glycosides (1-O-acyl glycosides, e.g. 1-O-galloylglucose)
    are cleaved at the anomeric ester linkage into carboxylic acid aglycone (Gallic acid)."""
    galloyl_glucose_smi = "O=C(OC1OC(CO)C(O)C(O)C1O)c1cc(O)c(O)c(O)c1"
    mol = Chem.MolFromSmiles(galloyl_glucose_smi)
    assert mol is not None

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol, min_aglycone_heavy=10)
    assert aglycone is not None
    assert "pyranose" in sugars
    # Aglycone is gallic acid (12 heavy atoms)
    assert aglycone.GetNumHeavyAtoms() == 12
    agl_smiles = Chem.MolToSmiles(aglycone)
    assert "C(=O)O" in agl_smiles or "C(=O)[O-]" in agl_smiles or "c1cc(O)c(O)c(O)c1" in agl_smiles


def test_ginsenoside_rg1_poly_glycoside_deglycosylation() -> None:
    """Verifies that multi-glycosylated natural saponins (Ginsenoside Rg1 with 2 pyranose units)
    are fully deglycosylated into the protopanaxatriol aglycone core."""
    rg1_smi = (
        "CC(=CCC[C@@](C)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H]2CC[C@]3(C)[C@@H]2"
        "[C@@H](O)C[C@@H]4[C@@]5(C)C[C@@H](O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)"
        "C(C)(C)[C@@H]5CC[C@]34C)C"
    )
    mol = Chem.MolFromSmiles(rg1_smi)
    assert mol is not None

    sugar_rings = tcm.identify_sugar_rings(mol)
    assert len(sugar_rings) == 2, f"Ginsenoside Rg1 should contain 2 sugar rings, found {len(sugar_rings)}"

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol)
    assert aglycone is not None
    assert len(sugars) == 2
    assert aglycone.GetNumHeavyAtoms() == 34


def test_aromatic_n_glycoside_adenosine_behavior() -> None:
    """Verifies the behavior of aromatic N-glycosides (Adenosine, Cordycepin).
    Identifies whether deglycosylation successfully yields the purine aglycone or preserves
    the intact nucleoside when aromatic N kekulization requires explicit hydrogen."""
    adenosine_smi = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    mol = Chem.MolFromSmiles(adenosine_smi)
    assert mol is not None

    sugar_rings = tcm.identify_sugar_rings(mol)
    assert len(sugar_rings) == 1

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol)
    # Document observed empirical behavior:
    # If aglycone is None, adenosine remains intact in screening pool without crashing library prep.
    if aglycone is not None:
        assert aglycone.GetNumHeavyAtoms() == 10  # Adenine
    else:
        # Intact nucleoside retained safely
        assert sugars == ["furanose"]


# =============================================================================
# 3. Focus Area 3: Heavy Metal & Exotic Element Traps
# =============================================================================

@pytest.mark.parametrize(
    "name,smiles,metal_symbol",
    [
        ("Roxarsone (Organo-As)", "O=[As](O)(O)c1ccc(O)c([N+](=O)[O-])c1", "As"),
        ("Melarsoprol (Organo-As)", "c1ccc(cc1)[As]1SC(CO)CS1", "As"),
        ("Arsanilic acid (Organo-As)", "Nc1ccc([As](=O)(O)O)cc1", "As"),
        ("Thiomersal (Organo-Hg)", "CC[Hg]Sc1ccccc1C(=O)[O-].[Na+]", "Hg"),
        ("Methylmercury chloride", "C[Hg]Cl", "Hg"),
        ("Auranofin (Organo-Au)", "CC(=O)OCC1OC(S[Au]P(CC)(CC)CC)C(OC(=O)C)C(OC(=O)C)C1OC(=O)C", "Au"),
        ("Tributyltin chloride", "CCCC[Sn](CCCC)(CCCC)Cl", "Sn"),
        ("Cisplatin", "N[Pt](N)(Cl)Cl", "Pt"),
    ],
)
def test_heavy_metal_covalent_rejection(name: str, smiles: str, metal_symbol: str) -> None:
    """Verifies that covalent heavy metal and exotic organometallic compounds are 100%
    intercepted with clear explanatory rejection messages specifying the forbidden element."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse SMILES for {name}"

    clean_mol, err = tcm.desalt_and_filter(mol)
    assert clean_mol is None, f"{name} should be rejected, but passed with SMILES: {Chem.MolToSmiles(clean_mol)}"
    assert err is not None
    assert "非有机/金属元素" in err
    assert metal_symbol in err


@pytest.mark.parametrize(
    "name,smiles,metal_symbol",
    [
        (
            "Ferrocene (ionic representation)",
            "[cH-]1cccc1.[cH-]1cccc1.[Fe+2]",
            "Fe",
        ),
        (
            "Gadopentetate (Gd-DTPA contrast agent)",
            "[Gd+3].[Na+].[Na+].[O-]C(=O)CN(CCN(CC(=O)[O-])CC(=O)[O-])CCN(CC(=O)[O-])CC(=O)[O-]",
            "Gd",
        ),
        (
            "Gadobutrol (Gadovist contrast agent)",
            "[Gd+3].CC(CO)(O)CN1CCN(CC(=O)[O-])CCN(CC(=O)[O-])CCN(CC(=O)[O-])1",
            "Gd",
        ),
    ],
)
def test_challenge_heavy_metal_ionic_complex_traps(
    name: str, smiles: str, metal_symbol: str
) -> None:
    """ADVERSARIAL CHALLENGE TEST:
    Verifies whether heavy metal complexes and contrast agents (Ferrocene, Gd-DTPA, Gadobutrol)
    are 100% rejected, OR whether SaltRemover/LargestFragmentChooser strips the metal ion
    and allows the organic chelator to slip through undetected.
    
    This test asserts strict rejection as mandated by Milestone 1 requirements."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse SMILES for {name}"

    clean_mol, err = tcm.desalt_and_filter(mol)
    # The requirement mandates: "Confirm 100% rejection with clear messages."
    assert clean_mol is None, (
        f"CRITICAL VULNERABILITY: {name} contains heavy metal {metal_symbol}, "
        f"but was NOT rejected! The metal was stripped and '{Chem.MolToSmiles(clean_mol)}' passed."
    )
    assert err is not None
    assert metal_symbol in err or "非有机/金属元素" in err


# =============================================================================
# 4. Focus Area 4: Macrocycle Torsion & Conformer Edge Cases
# =============================================================================

@pytest.mark.parametrize(
    "name,smiles",
    [
        ("14-lactone (Oxacyclotetradecan-2-one)", "O=C1CCCCCCCCCCCCCO1"),
        ("Erythronolide A (14-macrolide core)", "CC1CC(C)C(=O)C(C)C(O)C(C)(O)CC(C)C(=O)C(C)C(O)C(C)C(=O)O1"),
        ("Azalide core (15-macrolide core)", "CC1CC(C)C(=O)C(C)C(O)C(C)(O)CC(C)CN(C)C(C)C(O)C(C)C(=O)O1"),
        ("Tylactone (16-macrolide core)", "CCC1C(=O)OC(CC=CC(=O)C(CC(C=CC1=O)C)C)C"),
        ("Ginkgolide B (TCM caged trilactone)", "CC(C)(C)C12CC3C4C(C1=O)OC(=O)C4(O)C(O)C5(C3O2)C(=O)OC5=O"),
        ("Strychnine (heptacyclic caged alkaloid)", "O=C1CC2CC3N4CCc5c(cccc5)C45CC(C=C2)C3C(O1)C5"),
        ("Paclitaxel (Taxol, highly constrained tetracycle)", "CC(=O)OC1C(=O)C2(C)C(O)CC3OCC3(OC(=O)C)C2C(OC(=O)c2ccccc2)C4(O)C(OC(=O)C(O)C(NC(=O)c3ccccc3)c3ccccc3)CC=C14"),
        ("Camptothecin (pentacyclic planar alkaloid)", "CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3ccccc3nc21"),
        ("Cubane-1,4-dicarboxylate (caged cube)", "O=C(O)C12C3C4C1C5C2C3C45C(=O)O"),
    ],
)
def test_macrocycles_and_constrained_polycycles_conformer_stability(
    name: str, smiles: str
) -> None:
    """Verifies that ETKDGv3 macrocycle torsion sampling and MMFF94s minimization
    handle large rings (14-, 15-, 16-membered) and highly strained/caged polycycles
    without crashing or returning non-finite/inf energies."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse SMILES for {name}"

    best_mol, min_energy = tcm.generate_conformers_and_min_energy(mol, num_confs=5, seed=42)
    assert best_mol is not None
    assert best_mol.GetNumConformers() == 1
    assert isinstance(min_energy, float)
    assert math.isfinite(min_energy), f"{name}: calculated energy {min_energy} is not finite!"


# =============================================================================
# 5. Focus Area 5: TCM Warning Flags & Boundary Conditions
# =============================================================================

@pytest.mark.parametrize(
    "name,smiles,expected_flag",
    [
        ("1,4-Benzoquinone (para)", "O=C1C=CC(=O)C=C1", "Quinone_Redox"),
        ("1,4-Naphthoquinone (para)", "O=C1C=CC(=O)c2ccccc12", "Quinone_Redox"),
        ("Anthraquinone (para)", "O=C1c2ccccc2C(=O)c2ccccc12", "Quinone_Redox"),
        ("Juglone (para)", "O=C1C=CC(=O)c2c1cccc2O", "Quinone_Redox"),
        ("Emodin (para)", "Cc1cc2c(c(O)c1)C(=O)c1c(O)cc(O)cc1C2=O", "Quinone_Redox"),
        ("Rhein (para)", "O=C(O)c1cc2c(c(O)c1)C(=O)c1c(O)cccc1C2=O", "Quinone_Redox"),
        ("Shikonin (para)", "CC(=CCC1=C(O)C(=O)c2c(O)ccc(O)c2C1=O)C", "Quinone_Redox"),
        ("1,2-Benzoquinone (ortho)", "O=C1C(=O)C=CC=C1", "Quinone_Redox"),
        ("1,2-Naphthoquinone (ortho)", "O=C1C(=O)c2ccccc2C=C1", "Quinone_Redox"),
        ("Tanshinone IIA (ortho-quinone, Salvia miltiorrhiza)", "CC1(C)CCc2c(ccc3c2C(=O)C(=O)c2c(C)coc23)C1", "Quinone_Redox"),
        ("Cryptotanshinone (ortho-quinone)", "CC1CCc2c(C1)C(=O)C(=O)c1c(C)coc12", "Quinone_Redox"),
    ],
)
def test_ortho_and_para_quinones_flagged_consistently(
    name: str, smiles: str, expected_flag: str
) -> None:
    """Verifies that both ortho-quinones and para-quinones (including fused phenanthrenequinones
    such as Tanshinone IIA from Danshen) consistently trigger Quinone_Redox."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Failed to parse {name}"

    qc = tcm.scan_pains_and_tcm_alerts(mol)
    assert expected_flag in qc.tcm_warnings, (
        f"{name} failed to trigger {expected_flag}. Found: {qc.tcm_warnings}"
    )


@pytest.mark.parametrize(
    "name,smiles,expect_michael",
    [
        ("Chalcone", "O=C(/C=C/c1ccccc1)c1ccccc1", True),
        ("Curcumin", "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1", True),
        ("Methyl cinnamate (ester)", "COC(=O)/C=C/c1ccccc1", True),
        ("N-Phenylacrylamide (amide)", "C=CC(=O)Nc1ccccc1", True),
        ("Phenyl vinyl sulfone (sulfone)", "C=CS(=O)(=O)c1ccccc1", True),
    ],
)
def test_disguised_michael_acceptors(
    name: str, smiles: str, expect_michael: bool
) -> None:
    """Verifies that open-chain alpha,beta-unsaturated carbonyls, esters, amides, and sulfones
    are flagged as Michael_Acceptor."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None

    qc = tcm.scan_pains_and_tcm_alerts(mol)
    has_michael = "Michael_Acceptor" in qc.tcm_warnings
    assert has_michael == expect_michael, f"{name}: expected Michael_Acceptor={expect_michael}, got {has_michael}"


def test_tannin_polyphenol_mw_and_phenolic_oh_boundary() -> None:
    """Verifies strict boundary condition for Tannin_Polyphenol:
    Requires BOTH MW > 500 AND phenolic OH count >= 5."""
    # 1. High OH, but MW < 500: EGCG (MW=458.4, 8 phenolic OHs) -> Not a tannin
    egcg_smi = "Oc1cc(O)c2c(c1)O[C@@H](c1cc(O)c(O)c(O)c1)[C@H](OC(=O)c1cc(O)c(O)c(O)c1)C2"
    m_egcg = Chem.MolFromSmiles(egcg_smi)
    assert Descriptors.MolWt(m_egcg) < 500
    qc_egcg = tcm.scan_pains_and_tcm_alerts(m_egcg)
    assert "Tannin_Polyphenol" not in qc_egcg.tcm_warnings

    # 2. Both satisfied: Trigalloylglucose (MW=636.5, 9 phenolic OHs) -> Tannin
    trigalloyl_smi = "O=C(OCC1OC(OC(=O)c2cc(O)c(O)c(O)c2)C(O)C(OC(=O)c2cc(O)c(O)c(O)c2)C1O)c1cc(O)c(O)c(O)c1"
    m_tri = Chem.MolFromSmiles(trigalloyl_smi)
    assert Descriptors.MolWt(m_tri) > 500
    qc_tri = tcm.scan_pains_and_tcm_alerts(m_tri)
    assert "Tannin_Polyphenol" in qc_tri.tcm_warnings

    # 3. Both satisfied: Procyanidin B2 (MW=578.5, 8 phenolic OHs) -> Tannin
    procyanidin_smi = "Oc1cc(O)c2c(c1)O[C@@H](c1ccc(O)c(O)c1)[C@H](O)C2c1c(O)cc(O)c2c1O[C@@H](c1ccc(O)c(O)c1)[C@H](O)C2"
    m_pro = Chem.MolFromSmiles(procyanidin_smi)
    assert Descriptors.MolWt(m_pro) > 500
    qc_pro = tcm.scan_pains_and_tcm_alerts(m_pro)
    assert "Tannin_Polyphenol" in qc_pro.tcm_warnings

    # 4. Low MW, low OH: Gallic acid (MW=170.1, 3 OHs) -> Not a tannin
    m_gal = Chem.MolFromSmiles("O=C(O)c1cc(O)c(O)c(O)c1")
    qc_gal = tcm.scan_pains_and_tcm_alerts(m_gal)
    assert "Tannin_Polyphenol" not in qc_gal.tcm_warnings
