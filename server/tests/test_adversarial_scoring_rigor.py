"""Adversarial stress-testing suite for Milestone 2: Scoring Rigor Engine.

Empirically challenges:
1. LE, SILE, LLE edge cases (0, 1, 500 heavy atoms; negative/zero/extreme affinities; NaN/extreme cLogP).
2. MMFF94s strain energy & UFF fallback (Boron, Silicon, Selenium; severely clashed conformers; planar vs 3D; unsupported elements).
3. Multi-seed convergence with automorphic symmetry (benzene, cyclohexane, p-xylene, mesitylene; empty, single, NaN inputs).
4. ECR, Murcko clustering, and enrichment analytics extreme conditions.
"""

from __future__ import annotations

import math
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from dockstudio.analysis import scoring_rigor
from dockstudio.models import ScreeningCandidateRow


# =============================================================================
# 1. Ligand Efficiency Edge Cases (LE, SILE, LLE)
# =============================================================================

class TestLigandEfficiencyEdgeCases:
    """Rigorous boundary and stress tests for LE, SILE, and LLE formulas."""

    @pytest.mark.parametrize("n_heavy", [0, -1, -5, -100])
    def test_zero_and_negative_heavy_atoms(self, n_heavy: int) -> None:
        """Verifies no ZeroDivisionError or negative exponent error when heavy atoms <= 0."""
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=-8.5, num_heavy_atoms=n_heavy, clogp=2.0
        )
        assert res["le"] == 0.0
        assert res["sile"] == 0.0
        # LLE does not depend on heavy atom count: pkd = 8.5 / 1.364 = 6.2317 -> 6.2317 - 2.0 = 4.2317
        assert pytest.approx(res["lle"], abs=0.01) == 4.232

    def test_single_heavy_atom(self) -> None:
        """With 1 heavy atom (e.g. methane, formaldehyde), LE and SILE must be mathematically identical."""
        aff = -4.5
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=aff, num_heavy_atoms=1, clogp=0.5
        )
        # -(-4.5) / 1 = 4.5
        # -(-4.5) / (1^0.3) = 4.5
        assert res["le"] == 4.5
        assert res["sile"] == 4.5
        assert res["le"] == res["sile"]

    @pytest.mark.parametrize("n_heavy", [100, 500, 2000])
    def test_massive_heavy_atoms(self, n_heavy: int) -> None:
        """Verifies large macromolecules (peptides, macrolides) do not overflow or produce non-finite values."""
        aff = -12.0
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=aff, num_heavy_atoms=n_heavy, clogp=3.0
        )
        assert res["le"] is not None and res["le"] > 0.0
        assert res["sile"] is not None and res["sile"] > 0.0
        assert math.isfinite(res["le"])
        assert math.isfinite(res["sile"])
        # SILE must always be strictly greater than LE for N > 1
        assert res["sile"] > res["le"]

    def test_zero_affinity(self) -> None:
        """Affinity = 0.0 kcal/mol (no net binding energy)."""
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=0.0, num_heavy_atoms=25, clogp=1.2
        )
        assert res["le"] == 0.0
        assert res["sile"] == 0.0
        # pkd = 0.0 -> LLE = 0.0 - 1.2 = -1.2
        assert res["lle"] == -1.2

    def test_repulsive_positive_affinity(self) -> None:
        """Affinity = +15.0 kcal/mol (strongly unfavorable/repulsive interaction)."""
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=15.0, num_heavy_atoms=20, clogp=2.5
        )
        assert res["le"] == -0.75
        assert res["sile"] < 0.0
        assert res["lle"] < 0.0

    def test_ultra_tight_picomolar_affinity(self) -> None:
        """Affinity = -35.0 kcal/mol (ultra-tight binding, e.g. Streptavidin-Biotin)."""
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=-35.0, num_heavy_atoms=30, clogp=0.8
        )
        assert pytest.approx(res["le"], abs=0.01) == 1.167
        assert pytest.approx(res["sile"], abs=0.01) == 12.607
        assert pytest.approx(res["lle"], abs=0.01) == 24.859

    def test_clogp_nan_handling(self) -> None:
        """Verifies that NaN cLogP is handled gracefully with or without default fallback."""
        # Case A: default_none_clogp_to_zero = True (default) -> fallback to pkd
        res_default = scoring_rigor.compute_ligand_efficiencies(
            affinity=-8.5, num_heavy_atoms=20, clogp=float("nan"), default_none_clogp_to_zero=True
        )
        assert res_default["lle"] is not None
        assert not math.isnan(res_default["lle"])
        assert pytest.approx(res_default["lle"], abs=0.01) == 6.232

        # Case B: default_none_clogp_to_zero = False -> return None
        res_no_default = scoring_rigor.compute_ligand_efficiencies(
            affinity=-8.5, num_heavy_atoms=20, clogp=float("nan"), default_none_clogp_to_zero=False
        )
        assert res_no_default["lle"] is None

    def test_clogp_extreme_values(self) -> None:
        """Verifies extreme hydrophilic (cLogP = -15.0) and extreme lipophilic (cLogP = +15.0)."""
        res_hydro = scoring_rigor.compute_ligand_efficiencies(
            affinity=-8.5, num_heavy_atoms=20, clogp=-15.0
        )
        # pkd = 6.2317 - (-15.0) = 21.2317
        assert pytest.approx(res_hydro["lle"], abs=0.01) == 21.232

        res_lipo = scoring_rigor.compute_ligand_efficiencies(
            affinity=-8.5, num_heavy_atoms=20, clogp=15.0
        )
        # pkd = 6.2317 - 15.0 = -8.7683
        assert pytest.approx(res_lipo["lle"], abs=0.01) == -8.768

    def test_nan_affinity_safely_propagates(self) -> None:
        """When affinity is NaN, functions should not raise unhandled exceptions."""
        res = scoring_rigor.compute_ligand_efficiencies(
            affinity=float("nan"), num_heavy_atoms=15, clogp=1.0
        )
        assert math.isnan(res["le"])
        assert math.isnan(res["sile"])
        assert math.isnan(res["lle"])


# =============================================================================
# 2. MMFF94s Conformational Strain Energy & UFF Fallback
# =============================================================================

class TestMMFF94sStrainEnergyStress:
    """Stress tests MMFF94s parameter boundaries, UFF fallback, and physical distortions."""

    def test_boron_molecule_graceful_uff_fallback(self) -> None:
        """Boron compounds lack MMFF94s parameters in RDKit.
        Verifies that AllChem.MMFFHasAllMoleculeParams is False, but UFF fallback succeeds cleanly."""
        # Phenylboronic acid: OB(O)c1ccccc1
        mol = Chem.AddHs(Chem.MolFromSmiles("OB(O)c1ccccc1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        assert not AllChem.MMFFHasAllMoleculeParams(mol)

        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-7.5
        )
        assert math.isfinite(strain)
        assert strain >= 0.0
        assert math.isfinite(corr_aff)
        assert isinstance(is_high, bool)

    def test_bortezomib_core_uff_fallback(self) -> None:
        """Bortezomib boronic acid core lacks MMFF parameters; confirms UFF handles peptide boronate."""
        smi = "CC(C)CC(NC(=O)C(Cc1ccccc1)NC(=O)c1cnccn1)B(O)O"
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        assert not AllChem.MMFFHasAllMoleculeParams(mol)

        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-9.0
        )
        assert math.isfinite(strain)
        assert strain >= 0.0
        assert corr_aff >= -9.0

    def test_selenium_molecule_ebselen(self) -> None:
        """Ebselen (Se-containing drug) lacks MMFF parameters; verifies UFF handles selenium."""
        mol = Chem.AddHs(Chem.MolFromSmiles("O=C1Nc2ccccc2[Se]1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        assert not AllChem.MMFFHasAllMoleculeParams(mol)

        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-8.0
        )
        assert math.isfinite(strain)
        assert math.isfinite(corr_aff)

    def test_silane_molecules(self) -> None:
        """Silanes (Si) parameters check."""
        mol = Chem.AddHs(Chem.MolFromSmiles("C[Si](C)(C)c1ccccc1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-6.5
        )
        assert math.isfinite(strain)
        assert strain >= 0.0

    def test_severely_distorted_impossible_conformer(self) -> None:
        """Constructs an unphysical conformer where two bonded carbon atoms overlap at 0.05 Å.
        Verifies extreme strain (> 1000 kcal/mol) and triggers is_high_strain = True."""
        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
        conf.SetAtomPosition(1, Point3D(0.05, 0.0, 0.0))  # Extreme overlap clash

        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-6.0, high_strain_threshold=8.0
        )
        assert strain > 1000.0
        assert is_high is True
        assert corr_aff > 500.0

    def test_eclipsed_vs_anti_butane_torsional_strain(self) -> None:
        """Eclipsed butane (dihedral = 0°) has ~5.5-6.0 kcal/mol torsional strain compared to anti (180°)."""
        mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        # Force dihedral 0-1-2-3 to 0.0 deg (syn-periplanar / eclipsed)
        AllChem.SetDihedralDeg(mol.GetConformer(), 0, 1, 2, 3, 0.0)

        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-5.0, high_strain_threshold=8.0
        )
        # Expected torsional strain is between 4.5 and 7.5 kcal/mol
        assert 4.5 <= strain <= 7.5
        # High strain cutoff is 8.0 kcal/mol, so eclipsed butane should NOT be marked high strain
        assert is_high is False
        # Penalty should be 0.5 * (strain - 3.0)
        expected_penalty = 0.5 * (strain - 3.0)
        assert pytest.approx(corr_aff, abs=0.01) == -5.0 + expected_penalty

    def test_planar_cyclohexane_strain_against_3d_chair_baseline(self) -> None:
        """Planar cyclohexane (all atoms in z=0 plane) evaluated against 3D chair baseline.
        Must show enormous angle and eclipsing strain (> 100 kcal/mol) and trigger is_high_strain."""
        # 1. Relaxed 3D chair baseline
        m_chair = Chem.AddHs(Chem.MolFromSmiles("C1CCCCC1"))
        AllChem.EmbedMolecule(m_chair, randomSeed=42)
        ff_chair = AllChem.MMFFGetMoleculeForceField(
            m_chair, AllChem.MMFFGetMoleculeProperties(m_chair)
        )
        ff_chair.Minimize()
        e_chair = ff_chair.CalcEnergy()

        # 2. Flat planar cyclohexane
        m_flat = Chem.AddHs(Chem.MolFromSmiles("C1CCCCC1"))
        AllChem.Compute2DCoords(m_flat)

        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            m_flat, global_min_energy=e_chair, affinity=-5.0, high_strain_threshold=8.0
        )
        assert strain > 100.0
        assert is_high is True
        assert corr_aff > 50.0

    def test_docked_lower_than_min_clamped_to_zero(self) -> None:
        """If docked energy is numerically lower than baseline min, delta_e is clamped to 0.0."""
        strain, corr, is_high = scoring_rigor.compute_internal_strain_energy_penalty(
            docked_energy=-25.0, global_min_energy=-20.0, raw_affinity=-8.0
        )
        assert strain == 0.0
        assert corr == -8.0
        assert is_high is False

    def test_exotic_unsupported_atom_fallback(self) -> None:
        """An unsupported exotic atom (e.g. Uranium [U]) where force fields fail gracefully."""
        mol = Chem.MolFromSmiles("[U]")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        strain, corr_aff, is_high = scoring_rigor.compute_internal_strain_energy(
            mol, affinity=-4.0
        )
        assert strain == 0.0
        assert corr_aff == -4.0
        assert is_high is False

    def test_empty_and_none_mol_safety(self) -> None:
        """None mol or mol with 0 conformers returns safely without error."""
        s1, c1, h1 = scoring_rigor.compute_internal_strain_energy(None, affinity=-6.0)
        assert (s1, c1, h1) == (0.0, -6.0, False)

        empty_mol = Chem.MolFromSmiles("CC")
        assert empty_mol.GetNumConformers() == 0
        s2, c2, h2 = scoring_rigor.compute_internal_strain_energy(empty_mol, affinity=-6.0)
        assert (s2, c2, h2) == (0.0, -6.0, False)


# =============================================================================
# 3. Multi-Seed Convergence & Automorphic Symmetry RMSD
# =============================================================================

class TestMultiSeedConvergenceAutomorphicStress:
    """Stress tests multi-seed convergence with ring automorphisms and edge cases."""

    def test_benzene_rotational_symmetry_convergence(self) -> None:
        """Benzene (D6h symmetry) rotated by 60°, 120°, 180° in the pocket.
        Naive Euclidean RMSDs are 1.40 Å, 2.42 Å, 2.79 Å, but automorphic RMSD must be < 0.05 Å."""
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        m_heavy = Chem.RemoveHs(mol)
        pts = m_heavy.GetConformer().GetPositions()
        center = pts.mean(axis=0)

        poses = []
        for deg in [0, 60, 120, 180]:
            theta = np.radians(deg)
            R = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1]
            ])
            m_rot = Chem.Mol(m_heavy)
            new_pts = (pts - center) @ R.T + center
            for i, pt in enumerate(new_pts):
                m_rot.GetConformer().SetAtomPosition(i, pt)
            poses.append(m_rot)

        affs = [-5.50, -5.52, -5.49, -5.51]
        mean_aff, std_aff, mean_rmsd, is_converged = scoring_rigor.compute_multiseed_convergence(
            poses, affs, max_sigma_dG=0.35, max_rmsd_seed=2.0, in_place=True
        )

        assert pytest.approx(mean_aff, abs=0.01) == -5.505
        assert std_aff < 0.05
        assert mean_rmsd < 0.10  # Automorphic symmetry resolves 60/120/180 deg rotations
        assert is_converged is True

    def test_pxylene_twofold_symmetry_flip(self) -> None:
        """p-Xylene (Cc1ccc(C)cc1) flipped 180° around the ring normal.
        Naive RMSD is 3.78 Å, but automorphic RMSD must be < 0.15 Å and pass convergence."""
        mol = Chem.AddHs(Chem.MolFromSmiles("Cc1ccc(C)cc1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        m_heavy = Chem.RemoveHs(mol)
        pts = m_heavy.GetConformer().GetPositions()
        center = pts.mean(axis=0)

        # Ring normal vector
        v1 = pts[1] - pts[0]
        v2 = pts[2] - pts[0]
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)

        # 180 degree rotation around ring normal
        theta = np.pi
        K = np.array([
            [0, -normal[2], normal[1]],
            [normal[2], 0, -normal[0]],
            [-normal[1], normal[0], 0]
        ])
        R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)

        m_rot = Chem.Mol(m_heavy)
        new_pts = (pts - center) @ R.T + center
        for i, pt in enumerate(new_pts):
            m_rot.GetConformer().SetAtomPosition(i, pt)

        affs = [-6.20, -6.22]
        mean_aff, std_aff, mean_rmsd, is_converged = scoring_rigor.compute_multiseed_convergence(
            [m_heavy, m_rot], affs, max_sigma_dG=0.35, max_rmsd_seed=2.0, in_place=True
        )

        assert mean_rmsd < 0.15
        assert is_converged is True

    def test_mesitylene_threefold_symmetry_rotation(self) -> None:
        """Mesitylene (1,3,5-trimethylbenzene) rotated 120° around ring normal.
        Naive RMSD is 3.52 Å, but automorphic RMSD must be < 0.15 Å and pass convergence."""
        mol = Chem.AddHs(Chem.MolFromSmiles("Cc1cc(C)cc(C)c1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        m_heavy = Chem.RemoveHs(mol)
        pts = m_heavy.GetConformer().GetPositions()
        center = pts.mean(axis=0)

        v1 = pts[1] - pts[0]
        v2 = pts[2] - pts[0]
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)

        theta = 2 * np.pi / 3  # 120 deg
        K = np.array([
            [0, -normal[2], normal[1]],
            [normal[2], 0, -normal[0]],
            [-normal[1], normal[0], 0]
        ])
        R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)

        m_rot = Chem.Mol(m_heavy)
        new_pts = (pts - center) @ R.T + center
        for i, pt in enumerate(new_pts):
            m_rot.GetConformer().SetAtomPosition(i, pt)

        affs = [-6.80, -6.83]
        mean_aff, std_aff, mean_rmsd, is_converged = scoring_rigor.compute_multiseed_convergence(
            [m_heavy, m_rot], affs, max_sigma_dG=0.35, max_rmsd_seed=2.0, in_place=True
        )
        assert mean_rmsd < 0.15
        assert is_converged is True

    def test_cyclohexane_chair_automorphic_permutation(self) -> None:
        """Cyclohexane chair with permuted atom indices matching symmetry operation.
        Verifies automorphic RMSD = 0.0 Å and passes convergence."""
        mol = Chem.AddHs(Chem.MolFromSmiles("C1CCCCC1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        m_heavy = Chem.RemoveHs(mol)
        pts = m_heavy.GetConformer().GetPositions()

        # Permute atoms by chair symmetry (0->2, 1->3, 2->4, 3->5, 4->0, 5->1)
        m_perm = Chem.Mol(m_heavy)
        perm = [2, 3, 4, 5, 0, 1]
        for i in range(6):
            m_perm.GetConformer().SetAtomPosition(perm[i], pts[i])

        affs = [-4.50, -4.50]
        mean_aff, std_aff, mean_rmsd, is_converged = scoring_rigor.compute_multiseed_convergence(
            [m_heavy, m_perm], affs, max_sigma_dG=0.35, max_rmsd_seed=2.0, in_place=True
        )
        assert mean_rmsd == 0.0
        assert is_converged is True

    def test_divergent_spatial_poses_fail_convergence(self) -> None:
        """Poses translated by 3.5 Å in space must fail convergence due to high RMSD (> 2.0 Å)."""
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        m_heavy = Chem.RemoveHs(mol)
        pts = m_heavy.GetConformer().GetPositions()

        m_shifted = Chem.Mol(m_heavy)
        for i, pt in enumerate(pts + np.array([3.5, 0.0, 0.0])):
            m_shifted.GetConformer().SetAtomPosition(i, pt)

        affs = [-5.50, -5.50]  # Identical affinity, but different binding sites
        mean_aff, std_aff, mean_rmsd, is_converged = scoring_rigor.compute_multiseed_convergence(
            [m_heavy, m_shifted], affs, max_sigma_dG=0.35, max_rmsd_seed=2.0, in_place=True
        )
        assert mean_rmsd > 3.0
        assert is_converged is False

    def test_empty_and_nan_affinities_handling(self) -> None:
        """Tests empty list, single element, and NaN values in affinities list."""
        # 1. Empty lists
        assert scoring_rigor.compute_multiseed_convergence([], []) == (0.0, 0.0, 0.0, False)

        # 2. Single pose
        res_single = scoring_rigor.compute_multiseed_convergence([None], [-7.5])
        assert res_single[0] == -7.5
        assert res_single[1] == 0.0
        assert res_single[3] is True

        # 3. NaN in affinities
        res_nan = scoring_rigor.compute_multiseed_convergence([], [float("nan"), -7.5])
        assert math.isnan(res_nan[0])
        assert math.isnan(res_nan[1])
        assert res_nan[3] is False

    def test_in_place_false_uses_superimposed_get_best_rms(self) -> None:
        """When in_place=False, poses are superimposed using GetBestRMS."""
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        m_heavy = Chem.RemoveHs(mol)
        pts = m_heavy.GetConformer().GetPositions()

        # Shift one pose by 5.0 Å
        m_shifted = Chem.Mol(m_heavy)
        for i, pt in enumerate(pts + np.array([5.0, 0.0, 0.0])):
            m_shifted.GetConformer().SetAtomPosition(i, pt)

        affs = [-6.0, -6.0]
        # in_place=False aligns them, yielding RMSD ~0.0
        _, _, rmsd_super, conv_super = scoring_rigor.compute_multiseed_convergence(
            [m_heavy, m_shifted], affs, in_place=False
        )
        assert rmsd_super < 0.01
        assert conv_super is True

        # in_place=True does not align, yielding RMSD 5.0 Å
        _, _, rmsd_inplace, conv_inplace = scoring_rigor.compute_multiseed_convergence(
            [m_heavy, m_shifted], affs, in_place=True
        )
        assert rmsd_inplace >= 5.0
        assert conv_inplace is False


# =============================================================================
# 4. Exponential Consensus Ranking (ECR) & Diversity Stress
# =============================================================================

class TestScoringRigorECRAndDiversityStress:
    """Stress tests ECR consensus and Murcko-Butina clustering under edge conditions."""

    def test_ecr_empty_input(self) -> None:
        """Verifies empty list returns empty list without error."""
        assert scoring_rigor.compute_exponential_consensus_ranking([]) == []

    def test_ecr_missing_metrics_handled_gracefully(self) -> None:
        """Candidate rows with missing / None values for all or some metrics are penalized appropriately."""
        c1 = {"name": "c1", "affinity": -9.0, "sile": 3.5, "cnn_score": 0.8}
        c2 = {"name": "c2", "affinity": None, "sile": None, "cnn_score": None}
        c3 = {"name": "c3", "affinity": -7.0, "sile": 2.5, "cnn_score": 0.5}

        res = scoring_rigor.compute_exponential_consensus_ranking([c1, c2, c3])
        assert len(res) == 3
        # c1 has best metrics across board, c2 is missing all, c3 is moderate
        assert res[0]["name"] == "c1"
        assert res[1]["name"] == "c3"
        assert res[2]["name"] == "c2"
        assert res[0]["ecr_rank"] == 1
        assert res[2]["ecr_rank"] == 3

    def test_murcko_scaffold_inorganics_and_acyclic(self) -> None:
        """Inorganic salts, metals, and acyclic molecules have empty Murcko scaffold SMILES."""
        assert scoring_rigor.extract_murcko_scaffold(None) == ""
        assert scoring_rigor.extract_murcko_scaffold(Chem.MolFromSmiles("[Na+].[Cl-]")) == ""
        assert scoring_rigor.extract_murcko_scaffold(Chem.MolFromSmiles("CCCC")) == ""
        # Benzene Murcko scaffold is benzene
        assert scoring_rigor.extract_murcko_scaffold(Chem.MolFromSmiles("c1ccccc1")) == "c1ccccc1"

    def test_murcko_butina_all_acyclic_fallback(self) -> None:
        """When all molecules are acyclic, clustering falls back to whole-molecule Morgan fingerprints."""
        mols = [
            Chem.MolFromSmiles("CCCC"),
            Chem.MolFromSmiles("CCCCC"),
            Chem.MolFromSmiles("O=C(O)CC"),
        ]
        cluster_ids, rep_indices = scoring_rigor.cluster_murcko_scaffolds(mols, cutoff=0.5)
        assert len(cluster_ids) == 3
        assert len(rep_indices) >= 1

    def test_enrichment_metrics_empty_and_single_label(self) -> None:
        """Enrichment metrics handle empty lists and all-active / zero-active gracefully."""
        res_empty = scoring_rigor.compute_enrichment_metrics([], [])
        assert res_empty["auroc"] == 0.0

        # All inactive
        res_no_act = scoring_rigor.compute_enrichment_metrics([1.0, 2.0], [False, False])
        assert res_no_act["auroc"] == 0.0

        # All active
        res_all_act = scoring_rigor.compute_enrichment_metrics([1.0, 2.0], [True, True])
        assert res_all_act["auroc"] == 1.0
