"""Multi-dimensional scoring rigor engine: Ligand Efficiency (LE, SILE, LLE),
MMFF94s conformational strain energy, multi-seed convergence evaluation,
Exponential Consensus Ranking (ECR), Bemis-Murcko Butina diversity clustering,
and virtual screening enrichment analytics (BEDROC, AUROC, EF1%, EF5%)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from rdkit.ML.Scoring import Scoring

try:
    from scipy.stats import rankdata
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

if TYPE_CHECKING:
    from ..models import ScreeningCandidateRow


# ---------------------------------------------------------------------------
# 1. Ligand Efficiency Metrics (LE, SILE, LLE)
# ---------------------------------------------------------------------------

def compute_ligand_efficiencies(
    affinity: float,
    num_heavy_atoms: int,
    clogp: float | None = None,
    default_none_clogp_to_zero: bool = True,
) -> dict[str, float | None]:
    """Computes Ligand Efficiency (LE), Size-Independent LE (SILE), and Lipophilic LE (LLE).

    Formulas:
        LE   = -affinity / N_heavy
        SILE = -affinity / (N_heavy ** 0.3)
        LLE  = (-affinity / 1.364) - cLogP

    Args:
        affinity: Docking binding free energy ΔG in kcal/mol (negative is favorable).
        num_heavy_atoms: Number of non-hydrogen heavy atoms.
        clogp: Calculated octanol-water partition coefficient.
        default_none_clogp_to_zero: If True, missing cLogP defaults to 0.0 (LLE = pKd).

    Returns:
        dict with keys "le", "sile", "lle".
    """
    if num_heavy_atoms <= 0:
        le = 0.0
        sile = 0.0
    else:
        le = round(-affinity / num_heavy_atoms, 4)
        sile = round(-affinity / (num_heavy_atoms ** 0.3), 4)

    pkd = -affinity / 1.364
    if clogp is not None and not math.isnan(clogp):
        lle: float | None = round(pkd - clogp, 4)
    elif default_none_clogp_to_zero:
        lle = round(pkd, 4)
    else:
        lle = None

    return {"le": le, "sile": sile, "lle": lle}


# ---------------------------------------------------------------------------
# 2. MMFF94s Conformational Strain Energy
# ---------------------------------------------------------------------------

def compute_internal_strain_energy_penalty(
    docked_energy: float,
    global_min_energy: float,
    raw_affinity: float = 0.0,
    high_strain_threshold: float = 8.0,
) -> tuple[float, float, bool]:
    """Calculates internal conformational strain energy ΔE_strain and corrected affinity.

    Formulas:
        ΔE_strain = max(0.0, docked_energy - global_min_energy)
        ΔG_corr   = raw_affinity + 0.5 * max(0.0, ΔE_strain - 3.0)
        is_high_strain = bool(ΔE_strain > high_strain_threshold)

    Args:
        docked_energy: Single-point potential energy of the docked pose (kcal/mol).
        global_min_energy: Baseline relaxed minimum potential energy (kcal/mol).
        raw_affinity: Docking binding free energy ΔG in kcal/mol.
        high_strain_threshold: Cutoff above which poses are tagged as high strain.

    Returns:
        tuple (strain_energy, corrected_affinity, is_high_strain_flag).
    """
    delta_e = max(0.0, docked_energy - global_min_energy)
    penalty = 0.5 * max(0.0, delta_e - 3.0)
    corr_aff = raw_affinity + penalty
    is_high = bool(delta_e > high_strain_threshold)
    return round(delta_e, 4), round(corr_aff, 4), is_high


def compute_internal_strain_energy(
    docked_mol: Chem.Mol,
    global_min_energy: float | None = None,
    affinity: float = 0.0,
    high_strain_threshold: float = 8.0,
    force_constant: float = 100.0,
    **kwargs: Any,
) -> tuple[float, float, bool]:
    """Calculates internal conformational strain energy ΔE_strain and corrected affinity for a 3D Mol.

    Freezes all heavy atoms in their docked 3D coordinates, relaxes hydrogens with
    MMFF94s (fallback to UFF), and computes strain relative to the global/local
    relaxed minimum.

    Args:
        docked_mol: RDKit Mol with 3D conformer coordinates.
        global_min_energy: Baseline global minimum energy in kcal/mol. If None,
            computed via unconstrained force field relaxation.
        affinity: Raw docking binding free energy ΔG in kcal/mol.
        high_strain_threshold: Cutoff for marking artifact poses (default 8.0 kcal/mol).
        force_constant: Restraint force constant (kept for API compatibility).
        kwargs: Accepts optional raw_affinity alias for affinity.

    Returns:
        tuple (strain_energy, corrected_affinity, is_high_strain_flag).
    """
    raw_affinity = kwargs.get("raw_affinity", affinity)
    if docked_mol is None or docked_mol.GetNumConformers() == 0:
        return 0.0, raw_affinity, False

    mol = Chem.Mol(docked_mol)
    has_hs = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    if not has_hs:
        mol = Chem.AddHs(mol, addCoords=True)

    has_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
    if has_mmff:
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    else:
        props = None
        ff = AllChem.UFFGetMoleculeForceField(mol)

    if ff is None:
        return 0.0, raw_affinity, False

    # Freeze heavy atoms so only hydrogens relax
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            ff.AddFixedPoint(a.GetIdx())
    ff.Minimize(maxIts=200)
    e_docked = ff.CalcEnergy()

    # Determine baseline relaxed energy
    if global_min_energy is None or math.isnan(global_min_energy):
        if has_mmff and props is not None:
            ff_free = AllChem.MMFFGetMoleculeForceField(mol, props)
        else:
            ff_free = AllChem.UFFGetMoleculeForceField(mol)
        if ff_free is not None:
            ff_free.Minimize(maxIts=500)
            global_min_energy = ff_free.CalcEnergy()
        else:
            global_min_energy = e_docked

    return compute_internal_strain_energy_penalty(
        docked_energy=e_docked,
        global_min_energy=global_min_energy,
        raw_affinity=raw_affinity,
        high_strain_threshold=high_strain_threshold,
    )


# Alias for evaluate_docked_strain
evaluate_docked_strain = compute_internal_strain_energy


# ---------------------------------------------------------------------------
# 3. Multi-Seed Convergence Evaluation
# ---------------------------------------------------------------------------

def compute_multiseed_convergence(
    poses: list[Chem.Mol],
    affinities: list[float],
    max_sigma_dG: float = 0.35,
    max_rmsd_seed: float = 2.0,
    in_place: bool = True,
) -> tuple[float, float, float, bool]:
    """Evaluates multi-seed docking convergence across independent random seeds.

    Calculates:
        1. Mean affinity μ_ΔG
        2. Standard deviation σ_ΔG (must be <= max_sigma_dG, default 0.35 kcal/mol)
        3. Mean pairwise symmetry-corrected automorphic RMSD (must be <= max_rmsd_seed, default 2.0 Å)

    Args:
        poses: List of 3D Mol poses from independent seeds.
        affinities: List of docking affinities corresponding to each seed.
        max_sigma_dG: Maximum allowed affinity standard deviation (kcal/mol).
        max_rmsd_seed: Maximum allowed mean pairwise RMSD (Å).
        in_place: If True, uses in-pocket CalcRMS; if False, uses superimposed GetBestRMS.

    Returns:
        tuple (mean_affinity, std_affinity, mean_rmsd, is_converged).
    """
    if not affinities:
        return 0.0, 0.0, 0.0, False

    mean_aff = float(np.mean(affinities))
    std_aff = float(np.std(affinities, ddof=1)) if len(affinities) > 1 else 0.0

    valid_poses = [p for p in poses if p is not None and p.GetNumConformers() > 0]
    if len(valid_poses) < 2:
        is_converged = bool(std_aff <= max_sigma_dG)
        return round(mean_aff, 4), round(std_aff, 4), 0.0, is_converged

    # Remove hydrogens to prevent automorphism mapping combinatorial explosion
    heavy_poses = [Chem.RemoveHs(p) for p in valid_poses]
    rmsds: list[float] = []
    n = len(heavy_poses)
    for i in range(n):
        for j in range(i + 1, n):
            p1, p2 = heavy_poses[i], heavy_poses[j]
            try:
                if in_place:
                    r = AllChem.CalcRMS(p1, p2)
                else:
                    r = AllChem.GetBestRMS(Chem.Mol(p1), Chem.Mol(p2))
            except Exception:
                r = 999.0
            rmsds.append(r)

    mean_rmsd = float(np.mean(rmsds)) if rmsds else 0.0
    is_converged = bool((std_aff <= max_sigma_dG) and (mean_rmsd <= max_rmsd_seed))

    return round(mean_aff, 4), round(std_aff, 4), round(mean_rmsd, 4), is_converged


# ---------------------------------------------------------------------------
# 4. Exponential Consensus Ranking (ECR)
# ---------------------------------------------------------------------------

def _rank_array(arr: np.ndarray, higher_is_better: bool) -> np.ndarray:
    """Computes 1-based fractional ranks for an array using average ties."""
    target = -arr if higher_is_better else arr
    if _HAS_SCIPY:
        return rankdata(target, method="average")
    order = np.argsort(target)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(arr) + 1, dtype=float)
    return ranks


def compute_exponential_consensus_ranking(
    candidate_rows: list[Any],
    weights: dict[str, float] | None = None,
    sigma_r: float | None = None,
) -> list[Any]:
    """Computes Exponential Consensus Ranking (ECR) across heterogeneous scoring metrics.

    Formula:
        ECR(i) = sum_m w_m * exp(-r_{i,m} / sigma_r)
        sigma_r = max(1.0, 0.05 * N)

    Args:
        candidate_rows: List of objects (e.g. ScreeningCandidateRow) or dicts.
        weights: Optional mapping of metric name to weight.
        sigma_r: Optional decay scale parameter. Defaults to max(1.0, 0.05 * N).

    Returns:
        List of candidate rows sorted descending by ecr_score, with ecr_score
        and ecr_rank populated on each row.
    """
    if not candidate_rows:
        return []

    n_candidates = len(candidate_rows)
    if sigma_r is None:
        sigma_r = max(1.0, 0.05 * n_candidates)

    metric_directions = {
        "affinity_mean": False,
        "affinity": False,
        "strain_energy": False,
        "sile": True,
        "le": True,
        "lle": True,
        "cnn_score": True,
        "cnn_affinity": True,
        "anchor_match_ratio": True,
        "ecr_score": True,
    }

    if weights is None:
        weights = {
            "affinity_mean": 0.3,
            "sile": 0.3,
            "cnn_score": 0.2,
            "anchor_match_ratio": 0.2,
        }

    def _get_val(row: Any, key: str) -> Any:
        if isinstance(row, dict):
            return row.get(key)
        return getattr(row, key, None)

    def _set_val(row: Any, key: str, val: Any) -> None:
        if isinstance(row, dict):
            row[key] = val
        else:
            setattr(row, key, val)

    # Filter weights to metrics present on at least some candidates
    active_weights: dict[str, float] = {}
    for m, w in weights.items():
        if any(_get_val(r, m) is not None for r in candidate_rows):
            active_weights[m] = float(w)

    if not active_weights:
        for m in ("affinity_mean", "affinity", "sile"):
            if any(_get_val(r, m) is not None for r in candidate_rows):
                active_weights[m] = 1.0
                break

    total_w = sum(active_weights.values())
    if total_w <= 0:
        norm_weights = {m: 1.0 / len(active_weights) for m in active_weights}
    else:
        norm_weights = {m: w / total_w for m, w in active_weights.items()}

    # Compute ranks for each active metric
    ranks_by_metric: dict[str, np.ndarray] = {}
    for m in norm_weights:
        higher_is_better = metric_directions.get(m, True)
        vals = []
        for r in candidate_rows:
            v = _get_val(r, m)
            if v is None:
                vals.append(-1e9 if higher_is_better else 1e9)
            else:
                vals.append(float(v))
        ranks_by_metric[m] = _rank_array(np.array(vals), higher_is_better)

    # Calculate ECR score
    for idx, r in enumerate(candidate_rows):
        ecr = 0.0
        for m, w in norm_weights.items():
            rank_val = ranks_by_metric[m][idx]
            ecr += w * math.exp(-rank_val / sigma_r)
        _set_val(r, "ecr_score", round(ecr, 5))

    # Sort descending by ecr_score
    candidate_rows.sort(key=lambda x: _get_val(x, "ecr_score"), reverse=True)
    for rank_idx, r in enumerate(candidate_rows, start=1):
        _set_val(r, "ecr_rank", rank_idx)
        if hasattr(r, "rank"):
            _set_val(r, "rank", rank_idx)

    return candidate_rows


# Alias for compute_ecr_consensus
compute_ecr_consensus = compute_exponential_consensus_ranking


# ---------------------------------------------------------------------------
# 5. Bemis-Murcko Scaffold Butina Diversity Clustering
# ---------------------------------------------------------------------------

def extract_murcko_scaffold(mol: Chem.Mol | None, generic: bool = False) -> str:
    """Extracts the Bemis-Murcko scaffold from an RDKit Mol and returns its canonical SMILES.

    Args:
        mol: Input molecule.
        generic: If True, simplifies scaffold to generic carbon skeleton.

    Returns:
        Canonical SMILES of the scaffold, or empty string "" if acyclic or None.
    """
    if mol is None:
        return ""
    try:
        scaff = MurckoScaffold.GetScaffoldForMol(mol)
        if scaff is None or scaff.GetNumAtoms() == 0:
            return ""
        if generic:
            scaff = MurckoScaffold.MakeScaffoldGeneric(scaff)
        return Chem.MolToSmiles(scaff)
    except Exception:
        return ""


def cluster_murcko_scaffolds(
    mols: list[Chem.Mol | None],
    cutoff: float = 0.4,
    dist_thresh: float | None = None,
    generic: bool = False,
    reordering: bool = False,
) -> tuple[list[int], set[int]]:
    """Clusters molecules based on Bemis-Murcko scaffold ECFP4 fingerprints using Butina.

    Args:
        mols: List of RDKit molecules.
        cutoff: Distance cutoff (1 - Tanimoto). Default 0.4.
        dist_thresh: Optional alias for cutoff.
        generic: If True, uses generic carbon skeleton scaffolds.
        reordering: Whether to reorder clusters dynamically.

    Returns:
        tuple (cluster_ids, representative_indices).
    """
    if not mols:
        return [], set()

    n = len(mols)
    if n == 1:
        return [0], {0}

    threshold = dist_thresh if dist_thresh is not None else cutoff

    scaffolds: list[Chem.Mol] = []
    for m in mols:
        if m is None or m.GetNumAtoms() == 0:
            scaffolds.append(Chem.Mol())
            continue
        try:
            sc = MurckoScaffold.GetScaffoldForMol(m)
            if sc is None or sc.GetNumAtoms() == 0:
                scaffolds.append(m)
            else:
                if generic:
                    sc = MurckoScaffold.MakeScaffoldGeneric(sc)
                scaffolds.append(sc)
        except Exception:
            scaffolds.append(m)

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps: list[Any] = []
    for s in scaffolds:
        if s is not None and s.GetNumAtoms() > 0:
            try:
                fps.append(gen.GetFingerprint(s))
            except Exception:
                fps.append(None)
        else:
            fps.append(None)

    valid_indices = [i for i, fp in enumerate(fps) if fp is not None]
    if not valid_indices:
        return [0] * n, {0}

    valid_fps = [fps[i] for i in valid_indices]
    num_valid = len(valid_fps)

    dists: list[float] = []
    for i in range(num_valid):
        for j in range(i):
            sim = DataStructs.TanimotoSimilarity(valid_fps[i], valid_fps[j])
            dists.append(1.0 - sim)

    clusters = Butina.ClusterData(
        dists, num_valid, distThresh=threshold, isDistData=True, reordering=reordering
    )

    cluster_ids = [-1] * n
    representative_indices: set[int] = set()

    for c_id, members in enumerate(clusters):
        rep_orig_idx = valid_indices[members[0]]
        representative_indices.add(rep_orig_idx)
        for m_idx in members:
            cluster_ids[valid_indices[m_idx]] = c_id

    next_cid = len(clusters)
    for i in range(n):
        if cluster_ids[i] == -1:
            cluster_ids[i] = next_cid
            representative_indices.add(i)
            next_cid += 1

    return cluster_ids, representative_indices


def assign_clusters_to_candidates(
    candidates: list[Any],
    mols: list[Chem.Mol | None],
    cutoff: float = 0.4,
    dist_thresh: float | None = None,
    generic: bool = False,
) -> tuple[list[int], set[int]]:
    """Assigns cluster_id and is_cluster_rep onto candidates in-place."""
    cluster_ids, rep_indices = cluster_murcko_scaffolds(
        mols, cutoff=cutoff, dist_thresh=dist_thresh, generic=generic
    )
    for i, cand in enumerate(candidates):
        if i < len(cluster_ids):
            c_id = cluster_ids[i]
            is_rep = bool(i in rep_indices)
            if isinstance(cand, dict):
                cand["cluster_id"] = c_id
                cand["is_cluster_rep"] = is_rep
            else:
                setattr(cand, "cluster_id", c_id)
                setattr(cand, "is_cluster_rep", is_rep)
    return cluster_ids, rep_indices


# ---------------------------------------------------------------------------
# 6. Virtual Screening Enrichment Analytics
# ---------------------------------------------------------------------------

def compute_enrichment_metrics(
    scores: list[float] | None = None,
    is_active: list[bool] | list[int] | None = None,
    labels: list[bool] | list[int] | None = None,
    higher_is_better: bool = False,
    alpha: float = 20.0,
    fractions: list[float] | None = None,
    *,
    candidates: list[Any] | None = None,
    active_evaluator: Any | None = None,
) -> dict[str, float]:
    """Calculates virtual screening benchmark metrics: BEDROC, AUROC, and Enrichment Factors.

    Supports both raw score/label lists and candidates list with active_evaluator predicate.

    Args:
        scores: Docking or consensus scores.
        is_active: Boolean or 0/1 indicator of true active status.
        labels: Optional alias for is_active.
        higher_is_better: False for binding energy ΔG, True for ECR / SILE / CNNscore.
        alpha: Exponential weighting factor for BEDROC (default 20.0).
        fractions: List of fraction cutoffs for EF (default [0.01, 0.05]).
        candidates: Optional list of ScreeningCandidateRow objects.
        active_evaluator: Optional predicate function for candidate active status.

    Returns:
        dict with keys "bedroc", "bedroc_20", "auroc", "ef_1%", "ef_1pct", "ef_5%", "ef_5pct".
    """
    if fractions is None:
        fractions = [0.01, 0.05]

    if candidates is not None:
        if active_evaluator is None:
            active_evaluator = lambda c: bool(getattr(c, "is_active", False))
        scores = []
        is_active = []
        for c in candidates:
            if hasattr(c, "ecr_score"):
                scores.append(float(c.ecr_score))
            elif hasattr(c, "affinity_mean"):
                scores.append(float(c.affinity_mean))
            elif isinstance(c, dict):
                scores.append(float(c.get("ecr_score", c.get("affinity", 0.0))))
            else:
                scores.append(0.0)
            is_active.append(bool(active_evaluator(c)))
        higher_is_better = True

    if is_active is None and labels is not None:
        is_active = labels

    if scores is None or is_active is None or len(scores) == 0:
        return {
            "bedroc": 0.0,
            "bedroc_20": 0.0,
            "auroc": 0.0,
            "ef_1%": 0.0,
            "ef_1pct": 0.0,
            "ef_5%": 0.0,
            "ef_5pct": 0.0,
        }

    if len(scores) != len(is_active):
        raise ValueError(
            f"Length mismatch: {len(scores)} scores does not match {len(is_active)} labels"
        )

    n_actives = sum(1 for a in is_active if a)
    n_total = len(scores)

    if candidates is not None and (n_actives == 0 or n_actives == n_total):
        raise ValueError("At least one active and one inactive candidate required for enrichment metrics")

    if n_actives == 0:
        return {
            "bedroc": 0.0,
            "bedroc_20": 0.0,
            "auroc": 0.0,
            "ef_1%": 0.0,
            "ef_1pct": 0.0,
            "ef_5%": 0.0,
            "ef_5pct": 0.0,
        }
    if n_actives == n_total:
        return {
            "bedroc": 1.0,
            "bedroc_20": 1.0,
            "auroc": 1.0,
            "ef_1%": 1.0,
            "ef_1pct": 1.0,
            "ef_5%": 1.0,
            "ef_5pct": 1.0,
        }

    paired = list(zip(scores, [1 if a else 0 for a in is_active]))
    paired.sort(key=lambda x: x[0], reverse=higher_is_better)
    data = [[s, a] for s, a in paired]

    bedroc = float(Scoring.CalcBEDROC(data, 1, alpha))
    auroc = float(Scoring.CalcAUC(data, 1))
    efs = Scoring.CalcEnrichment(data, 1, fractions)

    ef_1 = float(efs[0]) if len(efs) > 0 else 0.0
    ef_5 = float(efs[1]) if len(efs) > 1 else 0.0

    return {
        "bedroc": round(bedroc, 4),
        "bedroc_20": round(bedroc, 4),
        "auroc": round(auroc, 4),
        "ef_1%": round(ef_1, 2),
        "ef_1pct": round(ef_1, 2),
        "ef_5%": round(ef_5, 2),
        "ef_5pct": round(ef_5, 2),
    }


def compute_enrichment_from_ranked_labels(
    ranked_labels: list[bool | int],
    alpha: float = 20.0,
    fractions: list[float] | None = None,
) -> dict[str, float]:
    """Computes enrichment metrics directly from an already rank-ordered list of active labels."""
    scores = list(range(len(ranked_labels)))
    return compute_enrichment_metrics(
        scores=scores,
        is_active=ranked_labels,
        higher_is_better=False,
        alpha=alpha,
        fractions=fractions,
    )
