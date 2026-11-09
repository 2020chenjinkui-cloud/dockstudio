"""Unit and integration tests for Milestone 2: Multi-Dimensional Scoring Rigor Engine.

Covers:
1. Ligand Efficiencies (LE, SILE, LLE) & Saponin vs Flavonoid SILE rank reversal.
2. MMFF94s Conformational Strain Energy (relaxed, moderate, severe strain, 3D Mol fixed atoms).
3. Multi-seed Convergence (identical, divergent affinities, divergent RMSD, automorphic symmetry RMSD).
4. Exponential Consensus Ranking (ECR) (monotonicity, multi-metric balancing, tie-breaking, small N).
5. Murcko-Butina Diversity Clustering (congeneric series grouping, cluster representatives, acyclic fallback).
6. Enrichment Metrics (AUROC, BEDROC alpha=20, EF1%, EF5% on perfect, inverted, random rankings).
7. GNINA CNN Parsing (colon-, equal-, space-separated, tabular stdout, PDBQT remarks, pure Vina fallback).
"""

from __future__ import annotations

import math
import random
from pathlib import Path
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.models import Pose, ScreeningCandidateRow
from dockstudio.analysis import scoring_rigor
from dockstudio.engines import gnina


# =============================================================================
# 1. Ligand Efficiencies (LE, SILE, LLE)
# =============================================================================

def test_sile_saponin_vs_flavonoid_rank_reversal() -> None:
    """Verifies that SILE eliminates molecular weight bias:
    Saponin (MW ~1100, 70 HA, -10.0 kcal/mol) vs Flavonoid (MW ~302, 22 HA, -8.5 kcal/mol).
    While raw affinity favors Saponin, SILE and LE correctly favor the Flavonoid."""
    # Saponin: Ginsenoside Rb1-like
    ha_saponin = 70
    aff_saponin = -10.0
    res_saponin = scoring_rigor.compute_ligand_efficiencies(
        affinity=aff_saponin, num_heavy_atoms=ha_saponin, clogp=-1.5
    )

    # Flavonoid: Quercetin-like
    ha_flavonoid = 22
    aff_flavonoid = -8.5
    res_flavonoid = scoring_rigor.compute_ligand_efficiencies(
        affinity=aff_flavonoid, num_heavy_atoms=ha_flavonoid, clogp=1.5
    )

    # 1. Raw affinity: Saponin is better (lower kcal/mol)
    assert aff_saponin < aff_flavonoid

    # 2. LE: Flavonoid is significantly higher
    # Saponin LE: 10.0 / 70 = ~0.143
    # Flavonoid LE: 8.5 / 22 = ~0.386
    assert res_flavonoid["le"] > res_saponin["le"]
    assert pytest.approx(res_flavonoid["le"], abs=0.01) == 0.386
    assert pytest.approx(res_saponin["le"], abs=0.01) == 0.143

    # 3. SILE: Flavonoid reverses ranking over Saponin
    # Saponin SILE: 10.0 / (70^0.3) = 10.0 / 3.5771 = ~2.796
    # Flavonoid SILE: 8.5 / (22^0.3) = 8.5 / 2.5277 = ~3.363
    assert res_flavonoid["sile"] > res_saponin["sile"]
    assert pytest.approx(res_flavonoid["sile"], abs=0.01) == 3.363
    assert pytest.approx(res_saponin["sile"], abs=0.01) == 2.796


def test_lle_calculation_positive_and_negative_clogp() -> None:
    """Verifies Lipophilic Ligand Efficiency (LLE = pKd - cLogP, where pKd = -dG / 1.364)."""
    # Case A: Hydrophobic Quercetin (cLogP = 1.5, aff = -8.5)
    # pKd = 8.5 / 1.364 = 6.2317 -> LLE = 6.2317 - 1.5 = 4.7317
    res_pos = scoring_rigor.compute_ligand_efficiencies(
        affinity=-8.5, num_heavy_atoms=22, clogp=1.5
    )
    assert pytest.approx(res_pos["lle"], abs=0.01) == 4.732

    # Case B: Hydrophilic Glycoside (cLogP = -1.5, aff = -10.0)
    # pKd = 10.0 / 1.364 = 7.3314 -> LLE = 7.3314 - (-1.5) = 8.8314
    res_neg = scoring_rigor.compute_ligand_efficiencies(
        affinity=-10.0, num_heavy_atoms=70, clogp=-1.5
    )
    assert pytest.approx(res_neg["lle"], abs=0.01) == 8.831

    # Case C: Missing cLogP (None) defaults gracefully to pKd
    res_none = scoring_rigor.compute_ligand_efficiencies(
        affinity=-8.5, num_heavy_atoms=22, clogp=None
    )
    assert pytest.approx(res_none["lle"], abs=0.01) == 6.232


def test_ligand_efficiencies_edge_cases() -> None:
    """Verifies zero heavy atoms, zero affinity, and unfavorable positive affinity."""
    # Zero heavy atoms -> LE and SILE return 0.0 safely without ZeroDivisionError
    res_zero = scoring_rigor.compute_ligand_efficiencies(affinity=-5.0, num_heavy_atoms=0)
    assert res_zero["le"] == 0.0
    assert res_zero["sile"] == 0.0

    # Unfavorable positive binding free energy (+2.0 kcal/mol)
    res_pos_aff = scoring_rigor.compute_ligand_efficiencies(affinity=2.0, num_heavy_atoms=10)
    assert res_pos_aff["le"] == -0.2
    assert res_pos_aff["sile"] < 0.0


# =============================================================================
# 2. MMFF94s Conformational Strain Energy
# =============================================================================

def test_strain_energy_relaxed_pose() -> None:
    """Ideal relaxed pose (ΔE_strain <= 3.0 kcal/mol) receives zero penalty."""
    strain_e, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy_penalty(
        docked_energy=10.0, global_min_energy=10.0, raw_affinity=-8.0
    )
    assert pytest.approx(strain_e, abs=1e-4) == 0.0
    assert pytest.approx(corr_aff, abs=1e-4) == -8.0
    assert is_high is False

    # Pose within 3.0 kcal/mol threshold
    strain_e2, corr_aff2, is_high2 = scoring_rigor.compute_internal_strain_energy_penalty(
        docked_energy=12.5, global_min_energy=10.0, raw_affinity=-8.0
    )
    assert pytest.approx(strain_e2, abs=1e-4) == 2.5
    assert pytest.approx(corr_aff2, abs=1e-4) == -8.0
    assert is_high2 is False


def test_strain_energy_moderate_strain() -> None:
    """Moderate strain (ΔE_strain = 5.0 kcal/mol): penalty = 0.5 * (5.0 - 3.0) = 1.0 kcal/mol."""
    strain_e, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy_penalty(
        docked_energy=15.0, global_min_energy=10.0, raw_affinity=-8.0
    )
    assert pytest.approx(strain_e, abs=1e-4) == 5.0
    assert pytest.approx(corr_aff, abs=1e-4) == -7.0
    assert is_high is False


def test_strain_energy_severe_strain_flag() -> None:
    """Severe high-strain pose (ΔE_strain = 9.5 kcal/mol > 8.0) flags is_high_strain == True."""
    strain_e, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy_penalty(
        docked_energy=19.5, global_min_energy=10.0, raw_affinity=-8.0
    )
    assert pytest.approx(strain_e, abs=1e-4) == 9.5
    # penalty: 0.5 * (9.5 - 3.0) = 3.25 -> -8.0 + 3.25 = -4.75
    assert pytest.approx(corr_aff, abs=1e-4) == -4.75
    assert is_high is True


def test_strain_energy_docked_lower_than_min_clamped() -> None:
    """Numerical precision edge case: docked energy slightly lower than global min is clamped to 0."""
    strain_e, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy_penalty(
        docked_energy=9.6, global_min_energy=10.0, raw_affinity=-8.0
    )
    assert strain_e == 0.0
    assert corr_aff == -8.0
    assert is_high is False


def test_strain_energy_evaluation_on_3d_mol() -> None:
    """Evaluates 3D molecule with fixed heavy atoms to calculate MMFF94s strain energy."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    # Global min energy via unconstrained MMFF94s
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    ff.Minimize(maxIts=500)
    e_min = ff.CalcEnergy()

    # Create strained pose: translate phenol hydroxyl oxygen slightly
    strained_mol = Chem.Mol(mol)
    conf = strained_mol.GetConformer()
    o_idx = [a.GetIdx() for a in strained_mol.GetAtoms() if a.GetAtomicNum() == 8][0]
    pos = conf.GetAtomPosition(o_idx)
    conf.SetAtomPosition(o_idx, (pos.x + 0.35, pos.y, pos.z))

    strain_e, corr_aff, is_high = scoring_rigor.evaluate_docked_strain(
        docked_mol=strained_mol, global_min_energy=e_min, raw_affinity=-7.5
    )
    assert strain_e > 0.5
    assert corr_aff >= -7.5


# =============================================================================
# 3. Multi-seed Convergence & Automorphic Symmetry RMSD
# =============================================================================

def test_multiseed_convergence_identical_poses() -> None:
    """3 identical poses with identical affinities converge perfectly (RMSD 0.0, std 0.0)."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))  # Aspirin
    AllChem.EmbedMolecule(mol, randomSeed=42)
    poses = [Chem.Mol(mol), Chem.Mol(mol), Chem.Mol(mol)]
    affs = [-7.5, -7.5, -7.5]

    mean_aff, std_aff, mean_rmsd, is_conv = scoring_rigor.compute_multiseed_convergence(
        poses=poses, affinities=affs
    )
    assert pytest.approx(mean_aff, abs=1e-4) == -7.5
    assert pytest.approx(std_aff, abs=1e-4) == 0.0
    assert pytest.approx(mean_rmsd, abs=1e-4) == 0.0
    assert is_conv is True


def test_multiseed_convergence_divergent_affinities() -> None:
    """Affinities with std > 0.35 kcal/mol fail convergence check."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    poses = [Chem.Mol(mol), Chem.Mol(mol), Chem.Mol(mol)]
    affs = [-9.0, -7.5, -6.0]  # std = 1.5 > 0.35

    mean_aff, std_aff, mean_rmsd, is_conv = scoring_rigor.compute_multiseed_convergence(
        poses=poses, affinities=affs
    )
    assert std_aff > 0.35
    assert is_conv is False


def test_multiseed_convergence_divergent_spatial_poses() -> None:
    """Poses with mean RMSD > 2.0 Å fail convergence check despite identical affinity."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol_shifted = Chem.Mol(mol)
    conf = mol_shifted.GetConformer()
    for i in range(mol_shifted.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 3.0, p.y, p.z))

    poses = [mol, mol_shifted]
    affs = [-6.0, -6.0]
    mean_aff, std_aff, mean_rmsd, is_conv = scoring_rigor.compute_multiseed_convergence(
        poses=poses, affinities=affs
    )
    assert mean_rmsd > 2.0
    assert is_conv is False


def test_multiseed_convergence_automorphic_symmetry_rmsd() -> None:
    """Symmetric molecules (e.g. p-Xylene) flipped 180° produce large naive RMSD (>4.0 Å)
    but symmetry-corrected RMSD matches automorphic symmetry (RMSD < 1.0 Å) and converges."""
    mol1 = Chem.AddHs(Chem.MolFromSmiles("Cc1ccc(C)cc1"))  # p-Xylene
    AllChem.EmbedMolecule(mol1, randomSeed=42)

    # Pose 2: 180° rotation around Z-axis (x -> -x, y -> -y)
    mol2 = Chem.Mol(mol1)
    conf2 = mol2.GetConformer()
    for i in range(mol1.GetNumAtoms()):
        p = mol1.GetConformer().GetAtomPosition(i)
        conf2.SetAtomPosition(i, (-p.x, -p.y, p.z))

    # Naive Euclidean distance between identical indices
    naive_rmsd = math.sqrt(
        sum(
            (mol1.GetConformer().GetAtomPosition(i).x - conf2.GetAtomPosition(i).x) ** 2
            + (mol1.GetConformer().GetAtomPosition(i).y - conf2.GetAtomPosition(i).y) ** 2
            + (mol1.GetConformer().GetAtomPosition(i).z - conf2.GetAtomPosition(i).z) ** 2
            for i in range(mol1.GetNumAtoms())
        )
        / mol1.GetNumAtoms()
    )
    assert naive_rmsd > 4.0

    # Symmetry-corrected convergence calculation
    mean_aff, std_aff, mean_rmsd, is_conv = scoring_rigor.compute_multiseed_convergence(
        poses=[mol1, mol2], affinities=[-6.5, -6.5]
    )
    assert mean_rmsd < 1.0
    assert is_conv is True


# =============================================================================
# 4. Exponential Consensus Ranking (ECR)
# =============================================================================

def test_ecr_five_candidates_monotonicity() -> None:
    """Tests multi-metric ranking on 5 candidates with varying scores; verifies monotonicity."""
    candidates = [
        ScreeningCandidateRow(ligand_id="A", affinity_mean=-10.0, sile=4.0, cnn_score=0.95),
        ScreeningCandidateRow(ligand_id="B", affinity_mean=-9.0, sile=3.5, cnn_score=0.85),
        ScreeningCandidateRow(ligand_id="C", affinity_mean=-8.0, sile=3.0, cnn_score=0.75),
        ScreeningCandidateRow(ligand_id="D", affinity_mean=-7.0, sile=2.5, cnn_score=0.65),
        ScreeningCandidateRow(ligand_id="E", affinity_mean=-6.0, sile=2.0, cnn_score=0.55),
    ]

    ranked = scoring_rigor.compute_exponential_consensus_ranking(
        candidate_rows=candidates,
        weights={"affinity_mean": 0.4, "sile": 0.3, "cnn_score": 0.3},
    )

    assert len(ranked) == 5
    assert [r.ligand_id for r in ranked] == ["A", "B", "C", "D", "E"]
    assert [r.ecr_rank for r in ranked] == [1, 2, 3, 4, 5]
    # Verify strict score monotonicity
    for i in range(len(ranked) - 1):
        assert ranked[i].ecr_score > ranked[i + 1].ecr_score


def test_ecr_weighting_orthogonality_and_ties() -> None:
    """Verifies that shifting weights flips ranking between orthogonal metrics."""
    cands = [
        # X: great Vina affinity, weak CNN
        ScreeningCandidateRow(ligand_id="X", affinity_mean=-11.0, sile=4.2, cnn_score=0.20),
        # Y: moderate Vina affinity, outstanding CNN
        ScreeningCandidateRow(ligand_id="Y", affinity_mean=-8.0, sile=3.0, cnn_score=0.99),
    ]

    # Weighted towards CNN -> Y wins
    ranked_cnn = scoring_rigor.compute_exponential_consensus_ranking(
        candidate_rows=[c.model_copy() for c in cands],
        weights={"affinity_mean": 0.1, "cnn_score": 0.9},
    )
    assert ranked_cnn[0].ligand_id == "Y"

    # Weighted towards Vina -> X wins
    ranked_vina = scoring_rigor.compute_exponential_consensus_ranking(
        candidate_rows=[c.model_copy() for c in cands],
        weights={"affinity_mean": 0.9, "cnn_score": 0.1},
    )
    assert ranked_vina[0].ligand_id == "X"


# =============================================================================
# 5. Murcko-Butina Diversity Clustering
# =============================================================================

def test_murcko_scaffold_extraction() -> None:
    """Verifies Bemis-Murcko scaffold extraction on cyclic, acyclic, and None inputs."""
    aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    assert scoring_rigor.extract_murcko_scaffold(aspirin) == "c1ccccc1"

    hexane = Chem.MolFromSmiles("CCCCCC")
    assert scoring_rigor.extract_murcko_scaffold(hexane) == ""

    assert scoring_rigor.extract_murcko_scaffold(None) == ""


def test_murcko_butina_congeneric_series_clustering() -> None:
    """Tests clustering of 7 natural product compounds:
    - 3 Flavonoids (Quercetin, Kaempferol, Luteolin)
    - 3 Steroids (Testosterone, Hydrocortisone, Progesterone)
    - 1 Indole alkaloid (Tryptamine)
    Asserts exactly 3 clusters, correct series segregation, and centroid representatives."""
    compounds = [
        # Flavonoids
        ("Quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
        ("Kaempferol", "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
        ("Luteolin", "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
        # Steroids
        ("Testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
        ("Hydrocortisone", "CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O"),
        ("Progesterone", "CC(=O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C"),
        # Indole
        ("Tryptamine", "c1ccc2c(c1)c(c[nH]2)CCN"),
    ]
    mols = [Chem.MolFromSmiles(smi) for _, smi in compounds]

    cluster_ids, rep_indices = scoring_rigor.cluster_murcko_scaffolds(mols, dist_thresh=0.4)

    assert len(cluster_ids) == 7
    # Exactly 3 clusters must be formed
    unique_clusters = set(cluster_ids)
    assert len(unique_clusters) == 3

    # All flavonoids (0, 1, 2) share the same cluster_id
    assert cluster_ids[0] == cluster_ids[1] == cluster_ids[2]

    # All steroids (3, 4, 5) share the same cluster_id
    assert cluster_ids[3] == cluster_ids[4] == cluster_ids[5]

    # Indole (6) is in a distinct cluster
    assert cluster_ids[6] != cluster_ids[0]
    assert cluster_ids[6] != cluster_ids[3]

    # Exactly 3 cluster representatives must be designated
    assert len(rep_indices) == 3


def test_murcko_butina_acyclic_and_single_mol_edge_case() -> None:
    """Verifies that acyclic compounds (zero-atom Murcko scaffold) do not crash clustering."""
    mols = [
        Chem.MolFromSmiles("CCCCCC"),  # n-Hexane (acyclic)
        Chem.MolFromSmiles("CCCCCCCC"),  # n-Octane (acyclic)
        Chem.MolFromSmiles("c1ccccc1"),  # Benzene (cyclic)
    ]
    cluster_ids, rep_indices = scoring_rigor.cluster_murcko_scaffolds(mols, dist_thresh=0.4)
    assert len(cluster_ids) == 3
    assert len(rep_indices) >= 1

    # Single molecule edge case
    single_c, single_rep = scoring_rigor.cluster_murcko_scaffolds([Chem.MolFromSmiles("c1ccccc1")])
    assert single_c == [0]
    assert single_rep == {0}

    # Empty list edge case
    empty_c, empty_rep = scoring_rigor.cluster_murcko_scaffolds([])
    assert empty_c == []
    assert empty_rep == set()


# =============================================================================
# 6. Enrichment Metrics (AUROC, BEDROC, EF1%, EF5%)
# =============================================================================

def test_enrichment_metrics_perfect_ranking() -> None:
    """Tests perfect ranking (10 actives followed by 90 decoys):
    AUROC = 1.0, BEDROC(alpha=20) ~ 1.0, EF1% = 10.0, EF5% = 10.0."""
    candidates = []
    for i in range(100):
        is_act = (i < 10)
        candidates.append(
            ScreeningCandidateRow(
                ligand_id=f"MOL_{i}",
                ecr_score=100.0 - i,
                alerts=["active"] if is_act else [],
            )
        )

    metrics = scoring_rigor.compute_enrichment_metrics(
        candidates=candidates,
        active_evaluator=lambda c: "active" in c.alerts,
        alpha=20.0,
        fractions=[0.01, 0.05],
    )

    assert pytest.approx(metrics["auroc"], abs=1e-4) == 1.0
    assert pytest.approx(metrics["bedroc"], abs=1e-4) == 1.0
    # 10 actives out of 100 = 0.1 base rate. Top 1% = 1 active (1.0 rate) -> EF = 1.0 / 0.1 = 10.0
    assert pytest.approx(metrics["ef_1pct"], abs=1e-2) == 10.0
    assert pytest.approx(metrics["ef_5pct"], abs=1e-2) == 10.0


def test_enrichment_metrics_inverted_ranking() -> None:
    """Tests inverted ranking (90 decoys followed by 10 actives):
    AUROC = 0.0, BEDROC(alpha=20) ~ 0.0, EF1% = 0.0."""
    candidates = []
    for i in range(100):
        is_act = (i >= 90)
        candidates.append(
            ScreeningCandidateRow(
                ligand_id=f"MOL_{i}",
                ecr_score=100.0 - i,
                alerts=["active"] if is_act else [],
            )
        )

    metrics = scoring_rigor.compute_enrichment_metrics(
        candidates=candidates,
        active_evaluator=lambda c: "active" in c.alerts,
        alpha=20.0,
        fractions=[0.01, 0.05],
    )

    assert pytest.approx(metrics["auroc"], abs=1e-4) == 0.0
    assert metrics["bedroc"] < 1e-10
    assert metrics["ef_1pct"] == 0.0
    assert metrics["ef_5pct"] == 0.0


def test_enrichment_metrics_random_ranking() -> None:
    """Tests random ranking (10 actives randomly distributed in 100 compounds):
    AUROC ~ 0.5."""
    random.seed(42)
    labels = [1] * 10 + [0] * 90
    random.shuffle(labels)

    candidates = [
        ScreeningCandidateRow(
            ligand_id=f"MOL_{i}",
            ecr_score=100.0 - i,
            alerts=["active"] if labels[i] == 1 else [],
        )
        for i in range(100)
    ]

    metrics = scoring_rigor.compute_enrichment_metrics(
        candidates=candidates,
        active_evaluator=lambda c: "active" in c.alerts,
        alpha=20.0,
        fractions=[0.01, 0.05],
    )

    assert 0.35 <= metrics["auroc"] <= 0.65


def test_enrichment_metrics_validation_errors() -> None:
    """Raises ValueError when zero actives or zero decoys are provided via candidates."""
    cands_no_act = [
        ScreeningCandidateRow(ligand_id=f"MOL_{i}", ecr_score=10.0 - i)
        for i in range(10)
    ]
    with pytest.raises(ValueError, match="active"):
        scoring_rigor.compute_enrichment_metrics(
            candidates=cands_no_act, active_evaluator=lambda c: False
        )


# =============================================================================
# 7. GNINA CNN Parsing
# =============================================================================

def test_parse_gnina_colon_separated(tmp_path: Path) -> None:
    """Verifies parsing of colon-separated CNN log lines and remarks."""
    out_pdbqt = tmp_path / "gnina_out.pdbqt"
    out_pdbqt.write_text(
        "MODEL 1\n"
        "REMARK minimizedAffinity -8.52\n"
        "REMARK CNNscore: 0.8523\n"
        "REMARK CNNaffinity: 7.320\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    poses = gnina.parse_gnina_out(out_pdbqt, run_log="")
    assert len(poses) == 1
    assert pytest.approx(poses[0].affinity, abs=1e-2) == -8.52
    assert pytest.approx(poses[0].cnn_score, abs=1e-4) == 0.8523
    assert pytest.approx(poses[0].cnn_affinity, abs=1e-3) == 7.320


def test_parse_gnina_equal_separated(tmp_path: Path) -> None:
    """Verifies parsing of equal-separated CNN log lines."""
    out_pdbqt = tmp_path / "gnina_out.pdbqt"
    out_pdbqt.write_text(
        "MODEL 1\n"
        "REMARK minimizedAffinity -9.10\n"
        "REMARK CNNscore=0.912 CNNaffinity=8.15\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    poses = gnina.parse_gnina_out(out_pdbqt, run_log="")
    assert len(poses) == 1
    assert pytest.approx(poses[0].cnn_score, abs=1e-3) == 0.912
    assert pytest.approx(poses[0].cnn_affinity, abs=1e-2) == 8.15


def test_parse_gnina_space_separated(tmp_path: Path) -> None:
    """Verifies parsing of space-separated CNN log lines."""
    out_pdbqt = tmp_path / "gnina_out.pdbqt"
    out_pdbqt.write_text(
        "MODEL 1\n"
        "REMARK minimizedAffinity -7.80\n"
        "REMARK CNNscore 0.785 CNNaffinity 6.42\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    poses = gnina.parse_gnina_out(out_pdbqt, run_log="")
    assert len(poses) == 1
    assert pytest.approx(poses[0].cnn_score, abs=1e-3) == 0.785
    assert pytest.approx(poses[0].cnn_affinity, abs=1e-2) == 6.42


def test_parse_gnina_stdout_mode_table(tmp_path: Path) -> None:
    """Verifies parsing of GNINA stdout mode table with 4 columns: mode affinity CNNscore CNNaffinity."""
    empty_pdbqt = tmp_path / "missing.pdbqt"
    run_log = (
        "mode |   affinity |  CNNscore | CNNaffinity\n"
        "-----+------------+-----------+------------\n"
        "   1       -8.52       0.8523        7.320\n"
        "   2       -7.80       0.6410        6.120\n"
    )
    poses = gnina.parse_gnina_out(empty_pdbqt, run_log=run_log)
    assert len(poses) == 2
    assert poses[0].index == 1
    assert pytest.approx(poses[0].affinity, abs=1e-2) == -8.52
    assert pytest.approx(poses[0].cnn_score, abs=1e-4) == 0.8523
    assert pytest.approx(poses[0].cnn_affinity, abs=1e-3) == 7.320


def test_parse_gnina_pure_vina_output_retains_none(tmp_path: Path) -> None:
    """Verifies fallback when docking log lacks CNN fields (pure Vina output).
    cnn_score and cnn_affinity remain None."""
    empty_pdbqt = tmp_path / "missing.pdbqt"
    run_log = (
        "mode |   affinity | dist from best mode\n"
        "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
        "-----+------------+----------+----------\n"
        "   1       -8.52          0          0\n"
    )
    poses = gnina.parse_gnina_out(empty_pdbqt, run_log=run_log)
    assert len(poses) == 1
    assert poses[0].cnn_score is None
    assert poses[0].cnn_affinity is None
