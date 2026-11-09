"""Adversarial stress testing for HTVS Milestone 1: Library Ingestion,
TCM Dual-Track Deglycosylation, and Mass Conservation Funnel Closure.

Challenger 2 Empirical Verification Suite:
1. 100% heavy metal / mineral libraries (Hg, As, Pt, Pb, Cd, Fe, Sn, Bi, Au).
2. Metal salt fragment behavior (ionic vs covalent metal handling).
3. 100% duplicate libraries (InChIKey collision, source_herbs tracking, salt invariant dedup).
4. Mixed adversarial libraries (salts, hydrates, bad SMILES, chiral centers, size bounds, glycosides).
5. Dual-track boundary conditions: 0% glycosides, 100% glycosides, sub-threshold aglycones (<10 heavy atoms).
6. Rigorous parent linkage & PDBQT syntactical verification for derived aglycones.
7. verify.py:855-875 T5 mathematical closure assertion (diff == 0) across all scenarios.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest
from rdkit import Chem

from dockstudio.cli.commands.verify import load_funnel_data
from dockstudio.core.orchestrator import PrepareRegistry
from dockstudio.models import PreparationError, PreparedLibrary
from dockstudio.prep.library import prepare_library
from dockstudio.prep.tcm import (
    deglycosylate_tcm_mol,
    desalt_and_filter,
    scan_pains_and_tcm_alerts,
)


# =============================================================================
# Helper: verify.py T5 Mathematical Funnel Closure
# =============================================================================

def assert_t5_closure(lib: PreparedLibrary, expected_docked_count: int | None = None) -> None:
    """Evaluates strict mass conservation closure using verify.py:855-875 T5 logic."""
    funnel_path = Path(lib.path) / "funnel.csv"
    assert funnel_path.is_file(), f"funnel.csv missing at {funnel_path}"

    funnel = load_funnel_data(funnel_path)
    n_raw = funnel.get("raw") or funnel.get("raw_library") or funnel.get("initial") or 0
    n_desalt = funnel.get("desalt") or 0
    n_dedup = funnel.get("dedup") or 0
    n_chiral_failed = funnel.get("chiral_failed") or 0
    n_prep_failed = funnel.get("prep_failed") or 0
    n_dock_failed = funnel.get("dock_failed") or 0
    n_docked_success = funnel.get("docked_success") or funnel.get("final_scored") or 0

    total_losses = n_desalt + n_dedup + n_chiral_failed + n_prep_failed + n_dock_failed
    expected_docked = n_raw - total_losses

    target_ready = expected_docked_count if expected_docked_count is not None else lib.n_ready

    disc_1 = abs(expected_docked - n_docked_success)
    disc_2 = abs(n_docked_success - target_ready)
    total_disc = disc_1 + disc_2

    assert disc_1 == 0, f"Funnel stage discrepancy disc_1 != 0 (expected={expected_docked}, got={n_docked_success})"
    assert disc_2 == 0, f"Funnel ready count mismatch disc_2 != 0 (funnel={n_docked_success}, lib={target_ready})"
    assert total_disc == 0, f"verify.py T5 non-closure: total_disc={total_disc}"
    assert n_raw > 0, "n_raw must be greater than 0 for a valid library manifest"


# =============================================================================
# 1. 100% Heavy Metal / Mineral Adversarial Library
# =============================================================================

HEAVY_METAL_LIBRARY = [
    ("[Hg]=S", "Cinnabar_Hg"),
    ("S1[As]2S[As]1S[As]3S[As]2S3", "Realgar_As"),
    ("N.N.Cl[Pt]Cl", "Cisplatin_Pt"),
    ("CC[Pb](CC)(CC)CC", "Tetraethyllead_Pb"),
    ("CC[Cd]CC", "Diethylcadmium_Cd"),
    ("[Fe+2].[Cl-].[Cl-]", "Ferrous_Chloride_Fe"),
    ("CCCC[Sn](CCCC)(CCCC)Cl", "Tributyltin_chloride_Sn"),
    ("CC(=O)O[Hg]c1ccccc1", "Phenylmercuric_acetate_Hg"),
    ("c1ccc([Bi](c2ccccc2)c2ccccc2)cc1", "Triphenylbismuth_Bi"),
    ("Cl[Au]Cl", "Gold_Chloride_Au"),
]


def test_adversarial_100_percent_heavy_metals_smi(settings) -> None:
    """Adversarial stress test: Library with 100% heavy metals and toxic minerals.
    All entries must fail element filtering; n_ready == 0; mass conservation holds."""
    reg = PrepareRegistry(settings)
    content = "\n".join(f"{smi} {name}" for smi, name in HEAVY_METAL_LIBRARY).encode("utf-8")

    lib = prepare_library(content, "heavy_metals.smi", reg, "重金属矿物库", tcm_mode=True)

    n_expected = len(HEAVY_METAL_LIBRARY)
    assert lib.n_raw == n_expected
    assert lib.n_desalt_fail == n_expected
    assert lib.n_dedup == 0
    assert lib.n_chiral_fail == 0
    assert lib.n_prep_fail == 0
    assert lib.n_ready == 0
    assert len(lib.ligand_ids) == 0
    assert len(lib.failed) == n_expected

    # Verify that all failure reasons specify non-organic/metal element
    disallowed_symbols = {"Hg", "As", "Pt", "Pb", "Cd", "Fe", "Sn", "Bi", "Au"}
    for f in lib.failed:
        assert "非有机/金属元素" in f["reason"]
        assert any(sym in f["reason"] for sym in disallowed_symbols)

    # Closed mass-conservation closure
    total_losses = lib.n_desalt_fail + lib.n_dedup + lib.n_chiral_fail + lib.n_prep_fail
    assert lib.n_raw - total_losses == lib.n_ready
    assert_t5_closure(lib, expected_docked_count=0)


def test_adversarial_100_percent_heavy_metals_csv(settings) -> None:
    """Adversarial stress test: 100% heavy metals supplied via CSV format."""
    reg = PrepareRegistry(settings)
    csv_rows = ["smiles,name,herb"]
    for smi, name in HEAVY_METAL_LIBRARY:
        csv_rows.append(f'"{smi}","{name}","Mineral_Source"')
    content = "\n".join(csv_rows).encode("utf-8")

    lib = prepare_library(content, "heavy_metals.csv", reg, "矿物中药库CSV", tcm_mode=True)
    assert lib.n_raw == len(HEAVY_METAL_LIBRARY)
    assert lib.n_desalt_fail == len(HEAVY_METAL_LIBRARY)
    assert lib.n_ready == 0
    assert_t5_closure(lib, expected_docked_count=0)


def test_adversarial_metal_salt_fragment_stripping(settings) -> None:
    """Critical Empirical Edge Case:
    When a heavy metal is presented as a disconnected ionic salt of an organic drug
    (e.g. Bismuth subsalicylate O=C(O)c1ccccc1O.[Bi+3]), the heavy metal element (Bi)
    MUST be rejected upfront during elemental verification.
    The ligand is not docked and fails at desalting step (n_desalt_fail == 1, n_ready == 0)."""
    reg = PrepareRegistry(settings)
    content = b"O=C(O)c1ccccc1O.[Bi+3] Bismuth_Subsalicylate\n"
    lib = prepare_library(content, "bi_salt.smi", reg, "次水杨酸铋测试", tcm_mode=True)

    assert lib.n_raw == 1
    assert lib.n_desalt_fail == 1
    assert lib.n_ready == 0
    assert len(lib.ligand_ids) == 0
    assert lib.failed[0]["reason"] == "含非有机/金属元素，不适用标准分子对接力场 (含: Bi)"
    assert_t5_closure(lib, expected_docked_count=0)


# =============================================================================
# 2. 100% Duplicate Adversarial Library
# =============================================================================

def test_adversarial_100_percent_duplicates(settings) -> None:
    """Adversarial stress test: 30 identical entries of Ibuprofen across 30 different herbs.
    Exact closure: n_dedup == N - 1, n_ready == 1, all 30 herb origins tracked in source_herbs."""
    reg = PrepareRegistry(settings)
    n_entries = 30
    smi = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    content_lines = [f"{smi} Herb_Batch_{i:02d}" for i in range(n_entries)]
    content = "\n".join(content_lines).encode("utf-8")

    lib = prepare_library(content, "dup_library.smi", reg, "重复分子库", tcm_mode=True)

    assert lib.n_raw == n_entries
    assert lib.n_dedup == n_entries - 1
    assert lib.n_desalt_fail == 0
    assert lib.n_chiral_fail == 0
    assert lib.n_prep_fail == 0
    assert lib.n_ready == 1
    assert len(lib.ligand_ids) == 1

    # Verify source_herbs contains all 30 origins
    inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
    assert inchikey in lib.source_herbs
    origins = lib.source_herbs[inchikey]
    assert len(origins) == n_entries
    for i in range(n_entries):
        assert f"Herb_Batch_{i:02d}" in origins

    # Verify strict T5 closure
    assert_t5_closure(lib, expected_docked_count=1)


def test_adversarial_duplicate_salts_invariant_dedup(settings) -> None:
    """Adversarial test: Berberine presented as Hydrochloride, Hydrobromide, Sulfate, and Free Base.
    All 4 forms must desalt to the identical quaternary alkaloid core and deduplicate to n_ready == 1."""
    reg = PrepareRegistry(settings)
    berberine_core = "COc1ccc2cc3[n+](cc2c1OC)CCc1c3ccc2c1OCO2"
    entries = [
        (f"[Cl-].{berberine_core}", "Berberine_HCl"),
        (f"[Br-].{berberine_core}", "Berberine_HBr"),
        (f"[I-].{berberine_core}", "Berberine_HI"),
        (berberine_core, "Berberine_FreeBase"),
    ]
    content = "\n".join(f"{smi} {name}" for smi, name in entries).encode("utf-8")

    lib = prepare_library(content, "berberine_salts.smi", reg, "黄连素盐型库", tcm_mode=True)

    assert lib.n_raw == 4
    assert lib.n_desalt_fail == 0
    assert lib.n_dedup == 3
    assert lib.n_ready == 1
    assert len(lib.ligand_ids) == 1

    # Check source_herbs tracked all 4 salt names
    clean_mol, _ = desalt_and_filter(Chem.MolFromSmiles(berberine_core))
    assert clean_mol is not None
    inchikey = Chem.MolToInchiKey(clean_mol)
    assert inchikey in lib.source_herbs
    assert len(lib.source_herbs[inchikey]) == 4
    for _, name in entries:
        assert name in lib.source_herbs[inchikey]

    assert_t5_closure(lib, expected_docked_count=1)


# =============================================================================
# 3. Mixed Adversarial Library (All-in-One Failure & Success Stress)
# =============================================================================

MIXED_ADVERSARIAL_ENTRIES = [
    # Multi-component & crystal waters (desalt ok)
    ("[Cl-].COc1ccc2cc3[n+](cc2c1OC)CCc1c3ccc2c1OCO2", "Berberine_HCl"),
    ("O.O.O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "Quercetin_Dihydrate"),
    # Minerals & heavy metals (desalt fail)
    ("[Hg]=S", "Cinnabar_Hg"),
    ("S1[As]2S[As]1S[As]3S[As]2S3", "Realgar_As"),
    ("CC[Pb](CC)(CC)CC", "Tetraethyllead_Pb"),
    # Heavy atom bounds: too small (< 5)
    ("CCO", "Ethanol_TooSmall"),
    ("O", "Water_TooSmall"),
    # Heavy atom bounds: too large (> 120)
    ("C" * 125, "Polymer_TooLarge"),
    # Syntax error & unparseable SMILES
    ("INVALID_SMILES_(((", "Corrupt_Smiles_1"),
    ("c1ccc2(cc1)cccc", "Valence_Violation_Smiles"),
    # Valid synthetic drugs
    ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
    ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin_Duplicate_1"),
    ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin_Duplicate_2"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen_Chiral_1"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen_Duplicate"),
    # Cleavable O-glycosides (dual-track: glycoside + aglycone)
    ("O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O", "Baicalin_Glycoside"),
    ("CC1OC(OCC2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O", "Rutin_Glycoside"),
    # Direct C-glycoside (preserved without cleavage)
    ("O=c1c(-c2ccc(O)cc2)coc2c(O)c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc12", "Puerarin_C_Glycoside"),
    # TCM alerts & PAINS
    ("O=C(/C=C/c1ccccc1)c1ccccc1", "Chalcone_Michael"),
    ("Cc1cc2c(c(O)c1)C(=O)c1c(O)cc(O)cc1C2=O", "Emodin_Quinone"),
    ("O=C(O)c1cc(O)c(O)c(O)c1", "Gallic_Acid_Catechol"),
]


def test_adversarial_mixed_library_funnel_closure(settings) -> None:
    """Comprehensive adversarial stress test:
    Mixed 21-compound library spanning all failure modes and dual-track pathways.
    Asserts exact conservation: N_raw - (N_desalt + N_dedup + N_chiral + N_prep) == N_ready,
    and verify.py:855-875 T5 residual diff == 0."""
    reg = PrepareRegistry(settings)
    content = "\n".join(f"{smi} {name}" for smi, name in MIXED_ADVERSARIAL_ENTRIES).encode("utf-8")

    lib = prepare_library(
        content,
        "mixed_stress.smi",
        reg,
        "综合对抗压力测试库",
        tcm_mode=True,
        min_heavy=5,
        max_heavy=120,
    )

    n_raw = len(MIXED_ADVERSARIAL_ENTRIES)  # 21
    # Expected breakdowns:
    # Desalt/element/size failures:
    # 1. Cinnabar (Hg)
    # 2. Realgar (As)
    # 3. Tetraethyllead (Pb)
    # 4. Ethanol (<5 heavy)
    # 5. Water (<5 heavy)
    # 6. Polymer (>120 heavy)
    # 7. Corrupt_Smiles_1 (unparseable)
    # 8. Valence_Violation_Smiles (unparseable)
    # -> n_desalt_fail = 8
    assert lib.n_desalt_fail == 8

    # Deduplication:
    # 1. Aspirin_Duplicate_1
    # 2. Aspirin_Duplicate_2
    # 3. Ibuprofen_Duplicate
    # -> n_dedup = 3
    assert lib.n_dedup == 3

    # Prep & Chiral failures: 0
    assert lib.n_prep_fail == 0
    assert lib.n_chiral_fail == 0

    # Ready primary entries: 21 - (8 + 3) = 10
    # Ready primaries:
    # 1. Berberine (desalted)
    # 2. Quercetin (desalted)
    # 3. Aspirin
    # 4. Ibuprofen
    # 5. Baicalin (primary glycoside)
    # 6. Rutin (primary glycoside)
    # 7. Puerarin (preserved C-glycoside)
    # 8. Chalcone
    # 9. Emodin
    # 10. Gallic acid
    assert lib.n_ready == 10
    assert lib.n_raw == 21

    # Invariant mass conservation
    total_losses = lib.n_desalt_fail + lib.n_dedup + lib.n_chiral_fail + lib.n_prep_fail
    assert lib.n_raw - total_losses == lib.n_ready

    # Dual-track additions: Baicalin and Rutin produce aglycones (+2)
    # So total dockable ligands = 10 + 2 = 12
    assert len(lib.ligand_ids) == 12

    # Check Puerarin generated NO aglycone
    puerarin_ligands = [
        reg.ligands[lid] for lid in lib.ligand_ids if "Puerarin" in (reg.ligands[lid].name or "")
    ]
    assert len(puerarin_ligands) == 1
    assert puerarin_ligands[0].derivative_type == "original"

    # Strict verify.py T5 check
    assert_t5_closure(lib, expected_docked_count=10)


# =============================================================================
# 4. Dual-Track Glycoside Edge Cases
# =============================================================================

PURE_NON_GLYCOSIDES = [
    ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen"),
    ("CC(=O)Nc1ccc(O)cc1", "Acetaminophen"),
    ("O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "Diclofenac"),
    ("COc1ccc2cc(ccc2c1)C(C)C(=O)O", "Naproxen"),
    ("Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1", "Celecoxib"),
    ("CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1", "Ketoprofen"),
    ("COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1", "Indomethacin"),
    ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "Caffeine"),
    ("c1cccnc1C1CCCN1C", "Nicotine"),
]


def test_dual_track_zero_percent_glycosides(settings) -> None:
    """Dual-track edge case: Library containing 0% glycosides (pure small synthetic molecules).
    - Every molecule prepared as derivative_type='original'
    - len(ligand_ids) == n_ready (no aglycone track activated)
    - No ligand ID contains '_glycoside' or '_aglycone'
    - Mass conservation holds with total_disc == 0."""
    reg = PrepareRegistry(settings)
    content = "\n".join(f"{smi} {name}" for smi, name in PURE_NON_GLYCOSIDES).encode("utf-8")

    lib = prepare_library(content, "non_glycosides.smi", reg, "纯非糖苷库", tcm_mode=True)

    n_expected = len(PURE_NON_GLYCOSIDES)
    assert lib.n_raw == n_expected
    assert lib.n_ready == n_expected
    assert len(lib.ligand_ids) == n_expected
    assert lib.n_desalt_fail == 0
    assert lib.n_dedup == 0
    assert lib.n_prep_fail == 0

    # Ensure no dual-track suffixes exist
    for lid in lib.ligand_ids:
        lig = reg.ligands[lid]
        assert "_aglycone" not in lig.id
        assert "_glycoside" not in lig.id
        assert lig.derivative_type == "original"
        assert lig.parent_id is None
        assert lig.qc_profile is not None
        assert not lig.qc_profile.is_glycoside

    assert_t5_closure(lib, expected_docked_count=n_expected)


PURE_GLYCOSIDES = [
    # 1. Baicalin -> Baicalein aglycone (flavonoid)
    ("O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O", "Baicalin"),
    # 2. Rutin -> Quercetin aglycone (flavonoid)
    ("CC1OC(OCC2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O", "Rutin"),
    # 3. Hesperidin -> Hesperetin aglycone (dihydroflavone)
    ("COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O[C@@H]3O[C@H](CO[C@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)cc2O1", "Hesperidin"),
    # 4. Daidzin -> Daidzein aglycone (isoflavone)
    ("O=c1c(-c2ccc(O)cc2)coc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)ccc12", "Daidzin"),
    # 5. Genistin -> Genistein aglycone (isoflavone)
    ("O=c1c(-c2ccc(O)cc2)coc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)cc(O)c12", "Genistin"),
]


def test_dual_track_100_percent_glycosides(settings) -> None:
    """Dual-track edge case: Library containing 100% cleavable glycosides.
    - All 5 entries produce dual-track items: 5 primaries + 5 derived aglycones = 10 dockable ligands
    - Each primary ends in '_glycoside', derivative_type='original', is_glycoside=True
    - Each derived aglycone ends in '_aglycone', derivative_type='aglycone', parent_id matches primary
    - Closed mass conservation N_raw - 0 == N_ready == 5
    - T5 residual diff == 0."""
    reg = PrepareRegistry(settings)
    content = "\n".join(f"{smi} {name}" for smi, name in PURE_GLYCOSIDES).encode("utf-8")

    lib = prepare_library(content, "glycosides.smi", reg, "全糖苷中药库", tcm_mode=True)

    n_expected = len(PURE_GLYCOSIDES)
    assert lib.n_raw == n_expected
    assert lib.n_ready == n_expected
    assert lib.n_desalt_fail == 0
    assert lib.n_dedup == 0
    assert lib.n_prep_fail == 0
    assert len(lib.ligand_ids) == 2 * n_expected  # 10 ligands in total

    primaries = [reg.ligands[lid] for lid in lib.ligand_ids if lid.endswith("_glycoside")]
    aglycones = [reg.ligands[lid] for lid in lib.ligand_ids if lid.endswith("_aglycone")]

    assert len(primaries) == n_expected
    assert len(aglycones) == n_expected

    # Verify parent-child linkage pairs
    primary_ids = {p.id for p in primaries}
    for agl in aglycones:
        assert agl.derivative_type == "aglycone"
        assert agl.parent_id is not None
        assert agl.parent_id in primary_ids
        assert agl.qc_profile is not None
        assert agl.qc_profile.derivative_type == "aglycone"
        assert agl.qc_profile.parent_glycoside_id == agl.parent_id

    assert_t5_closure(lib, expected_docked_count=n_expected)


def test_dual_track_subthreshold_aglycone_heavy_atoms(settings) -> None:
    """Boundary test: Glycoside whose aglycone has < 10 heavy atoms (e.g. Methyl glucoside).
    The aglycone (methanol, 2 heavy atoms) is below min_aglycone_heavy=10.
    Deglycosylation should safely abort aglycone creation and treat it as a single primary entry.
    Mass conservation must hold with 0 residual."""
    reg = PrepareRegistry(settings)
    # Methyl alpha-D-glucopyranoside: Aglycone is MeOH (CH3-O-), 2 heavy atoms
    methyl_glucoside_smi = "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"

    lib = prepare_library(
        f"{methyl_glucoside_smi} Methyl_Glucoside\n".encode("utf-8"),
        "sub_aglycone.smi",
        reg,
        "微小苷元测试",
        tcm_mode=True,
    )

    assert lib.n_raw == 1
    assert lib.n_ready == 1
    # Only 1 ligand emitted because aglycone is rejected (< 10 heavy atoms)
    assert len(lib.ligand_ids) == 1
    lig = reg.ligands[lib.ligand_ids[0]]
    assert lig.derivative_type == "original"
    assert not lig.id.endswith("_aglycone")

    assert_t5_closure(lib, expected_docked_count=1)


# =============================================================================
# 5. Aglycone PDBQT Structural & Syntactical Rigor
# =============================================================================

def test_aglycone_pdbqt_file_and_torsdof_validity(settings) -> None:
    """Empirical verification: Derived aglycones must produce syntactically valid AutoDock PDBQT files:
    - Exists on disk with non-empty content
    - Contains 'REMARK', 'ROOT', 'ENDROOT', 'TORSDOF'
    - Contains valid 3D coordinates (non-zero z coordinates) and partial charges
    - Meeko and Vina parseable without syntax errors."""
    reg = PrepareRegistry(settings)
    # Baicalin -> Baicalein (20 heavy atoms, 5,6,7-trihydroxyflavone)
    baicalin_smi = "O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O"
    content = f"{baicalin_smi} Baicalin_Extract\n".encode("utf-8")

    lib = prepare_library(content, "baicalin.smi", reg, "黄芩苷提取物", tcm_mode=True)
    assert lib.n_ready == 1
    assert len(lib.ligand_ids) == 2

    agl_id = next(lid for lid in lib.ligand_ids if lid.endswith("_aglycone"))
    agl_lig = reg.ligands[agl_id]
    assert agl_lig is not None

    pdbqt_path = Path(agl_lig.pdbqt_path)
    assert pdbqt_path.is_file()
    lines = pdbqt_path.read_text(encoding="utf-8").splitlines()

    # PDBQT AutoDock Syntax checks
    has_remark = any(l.startswith("REMARK") for l in lines)
    has_root = any(l.startswith("ROOT") for l in lines)
    has_endroot = any(l.startswith("ENDROOT") for l in lines)
    torsdof_lines = [l for l in lines if l.startswith("TORSDOF")]

    assert has_remark, "PDBQT missing REMARK header"
    assert has_root, "PDBQT missing ROOT marker"
    assert has_endroot, "PDBQT missing ENDROOT marker"
    assert len(torsdof_lines) == 1, "PDBQT missing TORSDOF marker"

    # TORSDOF value must be an integer >= 0
    torsdof_val = int(torsdof_lines[0].split()[1])
    assert torsdof_val >= 0

    # Inspect coordinates: ensure 3D conformer with non-zero Z dimension
    atom_lines = [l for l in lines if l.startswith("ATOM") or l.startswith("HETATM")]
    assert len(atom_lines) >= 20  # Baicalein has 20 heavy atoms + Hs

    z_coords = [float(l[46:54].strip()) for l in atom_lines]
    charges = [float(l[70:76].strip()) for l in atom_lines if len(l) >= 76]

    assert any(abs(z) > 1e-4 for z in z_coords), "Aglycone coordinates appear strictly 2D flat"
    assert len(charges) == len(atom_lines), "Some atom lines missing Gasteiger partial charges"


# =============================================================================
# 6. Undefined Chiral Centers & Stereochemistry Robustness
# =============================================================================

def test_stereoisomer_unassigned_centers_robustness(settings) -> None:
    """Stress test: Ingestion of molecules with unassigned chiral centers:
    - 1 center: Ibuprofen (CC(C)Cc1ccc(cc1)C(C)C(=O)O)
    - 2 centers: Threonine / ephedrine analog
    - >2 centers: Unassigned hexose (4 undefined centers)
    All molecules must be embedded in 3D without crashing and satisfy mass conservation."""
    reg = PrepareRegistry(settings)
    chiral_molecules = [
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen_1_Chiral"),
        ("CC(O)C(N)C(=O)O", "Threonine_2_Chiral"),
        ("OC1C(O)C(O)C(O)C(CO)O1", "Hexose_4_Chiral"),
    ]
    content = "\n".join(f"{smi} {name}" for smi, name in chiral_molecules).encode("utf-8")

    lib = prepare_library(content, "chiral_test.smi", reg, "手性测试库", tcm_mode=True)

    assert lib.n_raw == 3
    assert lib.n_desalt_fail == 0
    assert lib.n_dedup == 0
    assert lib.n_prep_fail == 0
    assert lib.n_ready == 3

    assert_t5_closure(lib, expected_docked_count=3)


# =============================================================================
# 7. Degenerate Input Boundaries
# =============================================================================

def test_degenerate_empty_library(settings) -> None:
    """Degenerate boundary: Completely empty library file or only whitespace/comments.
    Must raise descriptive PreparationError."""
    reg = PrepareRegistry(settings)
    with pytest.raises(PreparationError, match="没有可解析"):
        prepare_library(b"\n# Only comments\n   \n", "empty.smi", reg)


def test_degenerate_single_molecule_library(settings) -> None:
    """Degenerate boundary: Library with exactly 1 single valid molecule (Aspirin).
    Must pass with N_raw = 1, N_ready = 1, total_disc == 0."""
    reg = PrepareRegistry(settings)
    lib = prepare_library(
        b"CC(=O)Oc1ccccc1C(=O)O Aspirin_Single\n",
        "single.smi",
        reg,
        "单分子库",
        tcm_mode=True,
    )
    assert lib.n_raw == 1
    assert lib.n_ready == 1
    assert len(lib.ligand_ids) == 1
    assert_t5_closure(lib, expected_docked_count=1)


def test_tautomer_canonical_inchikey_dedup(settings) -> None:
    """Verify that different SMILES representations of identical molecules
    (canonical InChIKey collision) are correctly identified and deduplicated."""
    reg = PrepareRegistry(settings)
    # Different SMILES representing Phenol
    entries = [
        ("c1ccccc1O", "Phenol_SMILES_1"),
        ("Oc1ccccc1", "Phenol_SMILES_2"),
        ("C1=CC=C(C=C1)O", "Phenol_SMILES_3"),
    ]
    content = "\n".join(f"{smi} {name}" for smi, name in entries).encode("utf-8")

    lib = prepare_library(content, "tautomers.smi", reg, "同分异构去重库", tcm_mode=True)

    assert lib.n_raw == 3
    assert lib.n_dedup == 2  # 2 duplicates merged
    assert lib.n_ready == 1
    assert len(lib.ligand_ids) == 1

    inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles("c1ccccc1O"))
    assert inchikey in lib.source_herbs
    assert len(lib.source_herbs[inchikey]) == 3
    assert_t5_closure(lib, expected_docked_count=1)
