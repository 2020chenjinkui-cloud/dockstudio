"""Unit and integration tests for Milestone 1: TCM Ligand Preparation, Standardizing,
Deglycosylation, PAINS & Rigor QC.

Covers:
1. Desalting & Elemental Filtering (Berberine HCl, Quercetin hydrate, Cinnabar, Realgar, size bounds)
2. Sugar Topology & Dual-Track (Baicalin -> Baicalein 20 heavy atoms, Rutin -> Quercetin 22 heavy atoms, Puerarin C-glycoside preserved)
3. PAINS & TCM Warnings (Catechol, Quinone, Michael Acceptor, Tannin, clean controls)
4. Conformer Ensemble & Energy Baseline (Macrocycle ETKDGv3 MMFF94s baseline, Meeko PDBQT TORSDOF)
5. Library Funnel Mass Conservation Closure (verify.py T5 standard closure with diff == 0)
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pytest
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

from dockstudio.models import PreparationError, PreparedLigand
from dockstudio.prep import tcm


# -----------------------------------------------------------------------------
# Test Fixtures & Constants
# -----------------------------------------------------------------------------

BERBERINE_HCL_SMILES = "[Cl-].COc1ccc2cc3[n+](cc2c1OC)CCc1c3ccc2c1OCO2"
QUERCETIN_HYDRATE_SMILES = "O.O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
BAICALIN_SMILES = "O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O"
RUTIN_SMILES = "CC1OC(OCC2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"
PUERARIN_SMILES = "O=c1c(-c2ccc(O)cc2)coc2c(O)c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc12"
CINNABAR_SMILES = "[Hg]=[S]"
REALGAR_SMILES = "S1[As]2S[As]1S[As]3S[As]2S3"
CHALCONE_SMILES = "O=C(/C=C/c1ccccc1)c1ccccc1"
CURCUMIN_SMILES = "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1"
BENZOQUINONE_SMILES = "O=C1C=CC(=O)C=C1"
JUGLONE_SMILES = "O=C1C=CC(=O)c2c1cccc2O"
EMODIN_SMILES = "Cc1cc2c(c(O)c1)C(=O)c1c(O)cc(O)cc1C2=O"
SHIKONIN_SMILES = "CC(=CCC1=C(O)C(=O)c2c(O)ccc(O)c2C1=O)C"
GALLIC_ACID_SMILES = "O=C(O)c1cc(O)c(O)c(O)c1"
TRIGALLOYLGLUCOSE_SMILES = "O=C(OCC1OC(OC(=O)c2cc(O)c(O)c(O)c2)C(O)C(OC(=O)c2cc(O)c(O)c(O)c2)C1O)c1cc(O)c(O)c(O)c1"
ARTEMISININ_SMILES = "CC12CCC3C(C(=O)OC4C3(C1)OOC4(C)C)C2"
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
MACROCYCLE_16_SMILES = "O=C1CCCCCCCCCCCCCCCO1"


# =============================================================================
# 1. Test Desalting & Elemental Filtering
# =============================================================================

def test_desalt_berberine_hydrochloride() -> None:
    """Verifies that Berberine Hydrochloride is stripped of Cl-, retaining the
    quaternary alkaloid core (25 heavy atoms) without throwing multi-component errors."""
    mol = Chem.MolFromSmiles(BERBERINE_HCL_SMILES)
    assert mol is not None

    clean_mol, err = tcm.desalt_and_filter(mol)
    assert err is None
    assert clean_mol is not None

    # Counter-ion Cl- (atomic num 17) must be completely eliminated
    atom_symbols = {a.GetSymbol() for a in clean_mol.GetAtoms()}
    assert "Cl" not in atom_symbols
    assert clean_mol.GetNumHeavyAtoms() == 25

    # Check quaternary nitrogen is preserved
    has_pos_nitrogen = any(
        a.GetSymbol() == "N" and a.GetFormalCharge() == 1 for a in clean_mol.GetAtoms()
    )
    assert has_pos_nitrogen


def test_desalt_quercetin_hydrate() -> None:
    """Verifies that Quercetin hydrate is stripped of crystal water 'O',
    retaining the anhydrous quercetin core (22 heavy atoms)."""
    mol = Chem.MolFromSmiles(QUERCETIN_HYDRATE_SMILES)
    assert mol is not None

    clean_mol, err = tcm.desalt_and_filter(mol)
    assert err is None
    assert clean_mol is not None
    assert clean_mol.GetNumHeavyAtoms() == 22


def test_elemental_filtering_heavy_metals() -> None:
    """Verifies that minerals and heavy metal compounds (Cinnabar Hg, Realgar As)
    are intercepted with informative Chinese failure reasons."""
    # 1. Cinnabar ([Hg]=S)
    mol_hg = Chem.MolFromSmiles(CINNABAR_SMILES)
    clean_hg, err_hg = tcm.desalt_and_filter(mol_hg)
    assert clean_hg is None
    assert err_hg is not None
    assert "非有机/金属元素" in err_hg
    assert "Hg" in err_hg

    # 2. Realgar (As4S4 cage)
    mol_as = Chem.MolFromSmiles(REALGAR_SMILES)
    clean_as, err_as = tcm.desalt_and_filter(mol_as)
    assert clean_as is None
    assert err_as is not None
    assert "非有机/金属元素" in err_as
    assert "As" in err_as


def test_heavy_atom_size_bounds_filtering() -> None:
    """Verifies that small solvent fragments (< 5 heavy atoms) and macromolecules
    (> 120 heavy atoms) are rejected with size boundary errors."""
    # 1. Small solvent fragments: Ethanol (3 heavy atoms), Acetic acid (4 heavy atoms)
    mol_ethanol = Chem.MolFromSmiles("CCO")
    clean_eth, err_eth = tcm.desalt_and_filter(mol_ethanol, min_heavy=5, max_heavy=120)
    assert clean_eth is None
    assert err_eth is not None
    assert "重原子数超出范围" in err_eth

    mol_acetic = Chem.MolFromSmiles("CC(=O)O")
    clean_ace, err_ace = tcm.desalt_and_filter(mol_acetic, min_heavy=5, max_heavy=120)
    assert clean_ace is None
    assert err_ace is not None
    assert "重原子数超出范围" in err_ace

    # 2. Macromolecule (> 120 heavy atoms): 130-carbon chain
    mol_huge = Chem.MolFromSmiles("C" * 130)
    clean_huge, err_huge = tcm.desalt_and_filter(mol_huge, min_heavy=5, max_heavy=120)
    assert clean_huge is None
    assert err_huge is not None
    assert "重原子数超出范围" in err_huge


# =============================================================================
# 2. Test Sugar Topology & Deglycosylation
# =============================================================================

def test_baicalin_deglycosylation_to_baicalein() -> None:
    """Verifies that Baicalin (glucuronide, 32 heavy atoms) is cleaved at the O-glycosidic
    linkage into Baicalein aglycone (20 heavy atoms)."""
    mol = Chem.MolFromSmiles(BAICALIN_SMILES)
    assert mol is not None

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol, min_aglycone_heavy=10)
    assert aglycone is not None
    assert aglycone.GetNumHeavyAtoms() == 20
    assert len(sugars) >= 1
    assert "pyranose" in sugars[0]

    # Verify Baicalein formula (C15 H10 O5)
    smiles_aglycone = Chem.MolToSmiles(aglycone)
    reparsed = Chem.MolFromSmiles(smiles_aglycone)
    assert reparsed.GetNumHeavyAtoms() == 20


def test_rutin_deglycosylation_to_quercetin() -> None:
    """Verifies that Rutin (quercetin-3-O-rutinoside, 43 heavy atoms) is cleaved into
    Quercetin aglycone (22 heavy atoms)."""
    mol = Chem.MolFromSmiles(RUTIN_SMILES)
    assert mol is not None

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol, min_aglycone_heavy=10)
    assert aglycone is not None
    assert aglycone.GetNumHeavyAtoms() == 22
    assert len(sugars) == 2  # Disaccharide: rhamnose + glucose


def test_puerarin_c_glycoside_strictly_preserved() -> None:
    """Verifies that Puerarin (8-C-glucosyl-daidzein) has a direct aromatic C-C bond
    between flavone core and sugar ring, and is STRICTLY preserved intact without cleavage."""
    mol = Chem.MolFromSmiles(PUERARIN_SMILES)
    assert mol is not None

    aglycone, sugars = tcm.deglycosylate_tcm_mol(mol, min_aglycone_heavy=10)
    assert aglycone is None
    assert sugars == []


# =============================================================================
# 3. Test PAINS & TCM Warning Flags
# =============================================================================

def test_tcm_warning_catechol() -> None:
    """Verifies that catechol polyphenol structures (Catechin, Quercetin)
    are flagged with 'Catechol_Polyphenol'."""
    catechin_smiles = "Oc1cc(O)c2c(c1)OC(c1ccc(O)c(O)c1)C(O)C2"
    m_cat = Chem.MolFromSmiles(catechin_smiles)
    qc_cat = tcm.scan_pains_and_tcm_alerts(m_cat)
    assert "Catechol_Polyphenol" in qc_cat.tcm_warnings

    m_quer = Chem.MolFromSmiles(QUERCETIN_HYDRATE_SMILES)
    clean_quer, _ = tcm.desalt_and_filter(m_quer)
    assert clean_quer is not None
    qc_quer = tcm.scan_pains_and_tcm_alerts(clean_quer)
    assert "Catechol_Polyphenol" in qc_quer.tcm_warnings


def test_tcm_warning_quinone() -> None:
    """Verifies that quinone redox moieties (1,4-benzoquinone, Juglone, Emodin, Shikonin)
    are flagged with 'Quinone_Redox'."""
    # Monocyclic 1,4-benzoquinone
    m_bq = Chem.MolFromSmiles(BENZOQUINONE_SMILES)
    qc_bq = tcm.scan_pains_and_tcm_alerts(m_bq)
    assert "Quinone_Redox" in qc_bq.tcm_warnings

    # Fused quinones (naphthoquinones / anthraquinones)
    for name, smi in [
        ("Juglone", JUGLONE_SMILES),
        ("Emodin", EMODIN_SMILES),
        ("Shikonin", SHIKONIN_SMILES),
    ]:
        m = Chem.MolFromSmiles(smi)
        pat_p = Chem.MolFromSmarts("O=C1~[#6]~[#6]~C(=O)~[#6]~[#6]~1")
        pat_o = Chem.MolFromSmarts("O=C1~[#6](=O)~[#6]~[#6]~[#6]~[#6]~1")
        assert m.HasSubstructMatch(pat_p) or m.HasSubstructMatch(pat_o), f"{name} quinone core must be detected"
        qc = tcm.scan_pains_and_tcm_alerts(m)
        assert "Quinone_Redox" in qc.tcm_warnings, f"{name} must flag Quinone_Redox in QC profile"


def test_tcm_warning_michael_acceptor() -> None:
    """Verifies that alpha,beta-unsaturated carbonyls (Curcumin, Chalcone)
    are flagged with 'Michael_Acceptor'."""
    m_curc = Chem.MolFromSmiles(CURCUMIN_SMILES)
    qc_curc = tcm.scan_pains_and_tcm_alerts(m_curc)
    assert "Michael_Acceptor" in qc_curc.tcm_warnings

    m_chalc = Chem.MolFromSmiles(CHALCONE_SMILES)
    qc_chalc = tcm.scan_pains_and_tcm_alerts(m_chalc)
    assert "Michael_Acceptor" in qc_chalc.tcm_warnings


def test_tcm_warning_tannin_polyphenol() -> None:
    """Verifies that high-MW polyhydroxylated polyphenols (MW > 500 & >= 5 phenolic OHs)
    trigger 'Tannin_Polyphenol', whereas small gallic acid does not."""
    # Trigalloylglucose (MW=636.5, 9 phenolic OHs)
    m_tannin = Chem.MolFromSmiles(TRIGALLOYLGLUCOSE_SMILES)
    qc_tannin = tcm.scan_pains_and_tcm_alerts(m_tannin)
    assert "Tannin_Polyphenol" in qc_tannin.tcm_warnings

    # Gallic acid (MW=170.1, 3 phenolic OHs) -> Not a tannin
    m_gallic = Chem.MolFromSmiles(GALLIC_ACID_SMILES)
    qc_gallic = tcm.scan_pains_and_tcm_alerts(m_gallic)
    assert "Tannin_Polyphenol" not in qc_gallic.tcm_warnings


def test_clean_negative_controls() -> None:
    """Verifies that clean drugs (Ibuprofen, Aspirin) produce clean profiles
    without spurious TCM false-positive alerts."""
    m_ibu = Chem.MolFromSmiles(IBUPROFEN_SMILES)
    qc_ibu = tcm.scan_pains_and_tcm_alerts(m_ibu)
    assert len(qc_ibu.pains_alerts) == 0
    assert len(qc_ibu.tcm_warnings) == 0

    m_asp = Chem.MolFromSmiles(ASPIRIN_SMILES)
    qc_asp = tcm.scan_pains_and_tcm_alerts(m_asp)
    assert len(qc_asp.tcm_warnings) == 0

    m_art = Chem.MolFromSmiles(ARTEMISININ_SMILES)
    qc_art = tcm.scan_pains_and_tcm_alerts(m_art)
    assert len(qc_art.tcm_warnings) == 0


# =============================================================================
# 4. Test Conformer Ensemble & Energy Baseline
# =============================================================================

def test_conformer_ensemble_and_energy_baseline() -> None:
    """Verifies that ETKDGv3 generates conformer ensembles for flexible/macrocyclic
    molecules and calculates a finite MMFF94s global minimum energy baseline."""
    mol = Chem.MolFromSmiles(MACROCYCLE_16_SMILES)
    assert mol is not None

    best_mol, min_energy = tcm.generate_conformers_and_min_energy(mol, num_confs=10, seed=42)
    assert best_mol is not None
    assert best_mol.GetNumConformers() == 1
    assert isinstance(min_energy, float)
    assert math.isfinite(min_energy)


def test_berberine_produces_valid_pdbqt_with_torsdof() -> None:
    """Verifies that Berberine free cation (after Cl- stripping) converts cleanly into
    a valid Meeko PDBQT file containing TORSDOF and REMARK lines."""
    mol = Chem.MolFromSmiles(BERBERINE_HCL_SMILES)
    clean_mol, err = tcm.desalt_and_filter(mol)
    assert clean_mol is not None

    best_mol, min_energy = tcm.generate_conformers_and_min_energy(clean_mol, num_confs=5, seed=42)
    setups = MoleculePreparation().prepare(best_mol)
    assert setups, "Meeko produced no molecule setups"

    pdbqt, ok, err_msg = PDBQTWriterLegacy.write_string(setups[0])
    assert ok and pdbqt, f"PDBQT writing failed: {err_msg}"
    assert "TORSDOF" in pdbqt
    assert "REMARK" in pdbqt
    assert math.isfinite(min_energy)


# =============================================================================
# 5. Test Library Funnel Conservation Closure
# =============================================================================

def test_funnel_conservation_accounting_closure() -> None:
    """Verifies exact mass conservation closure across a heterogeneous mock library:
        N_raw - (N_desalt_fail + N_dedup + N_chiral_fail + N_prep_fail) == N_ready
    and ensures strict compliance with verify.py T5 standard (diff == 0)."""
    raw_entries = [
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen"),                         # 1. Normal clean drug -> Ready (1)
        (BERBERINE_HCL_SMILES, "berberine_hcl"),                              # 2. Salt -> Desalted -> Ready (2)
        (QUERCETIN_HYDRATE_SMILES, "quercetin_hydrate"),                      # 3. Hydrate -> Water stripped -> Ready (3)
        (BAICALIN_SMILES, "baicalin"),                                        # 4. Glycoside -> Dual-track -> Ready (4)
        (PUERARIN_SMILES, "puerarin"),                                        # 5. C-glycoside -> Preserved -> Ready (5)
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen_duplicate"),                # 6. Duplicate InChIKey -> Dedup (1)
        (CINNABAR_SMILES, "cinnabar_heavy_metal"),                            # 7. Hg heavy metal -> Desalt fail (1)
        (REALGAR_SMILES, "realgar_heavy_metal"),                              # 8. As heavy metal -> Desalt fail (2)
        ("invalid_smiles(((((", "bad_syntax"),                                # 9. Invalid SMILES -> Desalt fail (3)
        ("O", "pure_water_fragment"),                                         # 10. Heavy atom < 5 -> Desalt fail (4)
    ]

    n_raw = len(raw_entries)
    n_desalt_fail = 0
    n_dedup = 0
    n_chiral_fail = 0
    n_prep_fail = 0
    n_ready = 0

    seen_inchikeys: set[str] = set()
    ready_ids: list[str] = []
    failed_records: list[dict[str, Any]] = []

    for idx, (smi, name) in enumerate(raw_entries):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_desalt_fail += 1
            failed_records.append({"index": idx, "name": name, "reason": "SMILES解析失败"})
            continue

        clean_mol, err = tcm.desalt_and_filter(mol, min_heavy=5, max_heavy=120)
        if clean_mol is None:
            n_desalt_fail += 1
            failed_records.append({"index": idx, "name": name, "reason": err})
            continue

        inchikey = Chem.MolToInchiKey(clean_mol)
        if inchikey in seen_inchikeys:
            n_dedup += 1
            continue

        seen_inchikeys.add(inchikey)
        ready_ids.append(name)
        n_ready += 1

    # Exact Mass Conservation Accounting Assertion
    total_losses = n_desalt_fail + n_dedup + n_chiral_fail + n_prep_fail
    assert n_raw - total_losses == n_ready
    assert n_raw == 10
    assert n_desalt_fail == 4
    assert n_dedup == 1
    assert n_chiral_fail == 0
    assert n_prep_fail == 0
    assert n_ready == 5

    # T5 Verification Funnel Data Compliance
    funnel = {
        "raw": n_raw,
        "desalt": n_desalt_fail,
        "dedup": n_dedup,
        "chiral_failed": n_chiral_fail,
        "prep_failed": n_prep_fail,
        "dock_failed": 0,
        "docked_success": n_ready,
    }

    expected_docked = funnel["raw"] - (
        funnel["desalt"]
        + funnel["dedup"]
        + funnel["chiral_failed"]
        + funnel["prep_failed"]
        + funnel["dock_failed"]
    )
    disc_1 = abs(expected_docked - funnel["docked_success"])
    disc_2 = abs(funnel["docked_success"] - n_ready)
    total_disc = disc_1 + disc_2

    assert total_disc == 0, f"Funnel non-conservation detected: residual={total_disc}"
    assert funnel["raw"] > 0


def test_prepare_library_tcm_mode_dual_track_and_dedup(settings) -> None:
    """Verifies that prepare_library with tcm_mode=True:
    1. Dedupes identical InChIKeys and aggregates herb origins into source_herbs.
    2. Recognizes cleavable glycosides (Baicalin) and generates both glycoside and aglycone.
    3. Traps non-organic heavy metals (Cinnabar Hg) in n_desalt_fail.
    4. Satisfies mass conservation closure N_raw - losses == N_ready.
    5. Writes funnel.csv with exact accounting.
    """
    from dockstudio.core.orchestrator import PrepareRegistry
    from dockstudio.prep.library import prepare_library

    reg = PrepareRegistry(settings)
    content = (
        f"{BAICALIN_SMILES} baicalin_herb1\n"
        f"{IBUPROFEN_SMILES} ibuprofen_herb1\n"
        f"{IBUPROFEN_SMILES} ibuprofen_herb2\n"
        f"{CINNABAR_SMILES} cinnabar_mineral\n"
    ).encode("utf-8")

    lib = prepare_library(content, "tcm_test.smi", reg, "中药测试库", tcm_mode=True)
    assert lib.n_raw == 4
    assert lib.n_desalt_fail == 1  # Cinnabar Hg
    assert lib.n_dedup == 1        # ibuprofen duplicate
    assert lib.n_prep_fail == 0
    assert lib.n_chiral_fail == 0
    assert lib.n_ready == 2        # Baicalin parent + Ibuprofen

    # Baicalin generated dual-track: primary glycoside + derived aglycone
    # So len(lib.ligand_ids) is 3 (1 baicalin_glycoside + 1 baicalin_aglycone + 1 ibuprofen)
    assert len(lib.ligand_ids) == 3
    assert any("aglycone" in lid for lid in lib.ligand_ids)
    assert any("glycoside" in lid for lid in lib.ligand_ids)

    # InChIKey herb deduplication tracked
    ibu_inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(IBUPROFEN_SMILES))
    assert ibu_inchikey in lib.source_herbs
    assert "ibuprofen_herb1" in lib.source_herbs[ibu_inchikey]
    assert "ibuprofen_herb2" in lib.source_herbs[ibu_inchikey]

    # Mass conservation equation holds
    total_losses = lib.n_desalt_fail + lib.n_dedup + lib.n_chiral_fail + lib.n_prep_fail
    assert lib.n_raw - total_losses == lib.n_ready

    # funnel.csv exists and matches
    funnel_path = Path(lib.path) / "funnel.csv"
    assert funnel_path.is_file()
    lines = funnel_path.read_text(encoding="utf-8").strip().splitlines()
    funnel_dict = dict(line.split(",") for line in lines)
    assert int(funnel_dict["raw"]) == 4
    assert int(funnel_dict["desalt"]) == 1
    assert int(funnel_dict["dedup"]) == 1
    assert int(funnel_dict["docked_success"]) == 2

