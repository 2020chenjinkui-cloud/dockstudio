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

#: Number of ETKDGv3 conformers used to locate the global minimum when the
#: caller does not supply a trustworthy ``global_min_energy`` baseline.
_STRAIN_EMBED_NUM_CONFS = 10
#: Deterministic embedding seed (contract Task Fix 3.1: ``randomSeed=42``).
_STRAIN_EMBED_RANDOM_SEED = 42
#: Unconstrained minimisation steps for each ensemble member / local relaxation.
_STRAIN_MINIMIZE_MAX_ITS = 500

# ---- Polar-hydrogen (H-bond network) rotamer search -----------------------
# Task Fix 3.1b: after a PDBQT round trip the pose's hydrogens are re-generated
# with ``AddHs(addCoords=True)`` and can sit in a hydrogen-bond network that is
# many kcal/mol above the prep baseline. Relaxing hydrogens *alone* cannot cross
# O-H/N-H rotation barriers, which used to fake 4-26 kcal/mol of "strain" for
# polyphenols (quercetin / rutin / baicalin / arbutin) and vetoed legitimate
# Tier-2 poses. The docked heavy-atom geometry is authoritative, so only the
# polar hydrogens are re-oriented, deterministically:
_STRAIN_H_ROTAMER_STEPS = 12  # 30 deg grid around the donor-heavy/upstream axis
_STRAIN_MAX_POLAR_H_DONORS = 10  # hard cap (rutin needs 10; see task report)
_STRAIN_SCAN_MIN_ITS = 80  # H relaxation per rotamer candidate
#: Total minimisation budget for one call (Task Fix 3.1c raised it from 80 to 240
#: so the worst-case coupled H-bond networks still have head-room; the cheap path
#: only uses ~one minimisation per donor).
_STRAIN_MAX_MINIMIZATIONS = 240
_STRAIN_FALLBACK_RESIDUAL = 1.5  # kcal/mol above which multi-start is attempted
#: Coordinated all-donor rotations (lock-step patterns) for coupled networks.
#: These are applied to *every* donor with the same offset, so the sweep is
#: independent of donor ordering (Task Fix 3.1c).
_STRAIN_MULTISTART_PATTERNS = 16
#: Additional minimisations allowed for the multi-start fallback (bounded).
_STRAIN_MAX_MULTISTART_MINIMIZATIONS = 24
#: Tolerance (kcal/mol) used to decide whether a fully relaxed pose lands back
#: near the caller's baseline (=> genuine strain) or nowhere near it
#: (=> the prep baseline was evaluated in a different force-field state).
_STRAIN_STATE_TOLERANCE = 5.0
_POLAR_DONOR_ELEMENTS = (7, 8, 16)  # N, O, S


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


def _is_finite_number(value: Any) -> bool:
    """True when ``value`` is a finite real number (None/NaN/inf -> False)."""
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


#: Coordinate resolution used to break canonical-ranking ties (0.0001 Å).
_CANONICAL_COORD_DECIMALS = 4


def _canonical_atom_order(mol: Chem.Mol) -> list[int]:
    """Total atom order derived only from the molecular graph and 3D coordinates.

    ``Chem.CanonicalRankAtoms(breakTies=False)`` leaves symmetry-equivalent atoms
    tied; the tie is broken here by the atom's own 3D position (which travels
    with the atom through ``Chem.RenumberAtoms``). The resulting order is
    therefore identical for every input numbering of the same conformer, which is
    what makes the heuristic searches below atom-order independent.
    """
    n_atoms = mol.GetNumAtoms()
    try:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        return list(range(n_atoms))
    if len(ranks) != n_atoms:
        return list(range(n_atoms))
    try:
        positions = np.array(mol.GetConformer().GetPositions(), dtype=float)
    except Exception:
        positions = np.zeros((n_atoms, 3), dtype=float)
    quantized = np.round(positions, _CANONICAL_COORD_DECIMALS)
    return sorted(
        range(n_atoms),
        key=lambda idx: (
            int(ranks[idx]),
            float(quantized[idx][0]),
            float(quantized[idx][1]),
            float(quantized[idx][2]),
        ),
    )


def _canonicalize_atom_order(mol: Chem.Mol) -> Chem.Mol:
    """Renumbers ``mol`` canonically (conformer coordinates move with the atoms)."""
    try:
        order = _canonical_atom_order(mol)
        if len(order) != mol.GetNumAtoms():
            return mol
        renumbered = Chem.RenumberAtoms(mol, order)
        return renumbered if renumbered is not None else mol
    except Exception:
        return mol


def _sanitized_copy(mol: Chem.Mol) -> Chem.Mol | None:
    """Full RDKit sanitisation of a copy (only used as a force-field fallback)."""
    try:
        work = Chem.Mol(mol)
        work.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(work)
        return work
    except Exception:
        try:
            work = Chem.Mol(mol)
            work.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(
                work,
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
            )
            return work
        except Exception:
            return None


def _force_field_available(mol: Chem.Mol) -> bool:
    """True when an MMFF94s (UFF fallback) force field can be built for ``mol``.

    The probe runs on a *copy*: ``MMFFGetMoleculeProperties`` mutates the
    molecule's aromaticity perception, and that must never leak into the
    molecule whose force-field state we are trying to preserve.
    """
    try:
        probe = Chem.Mol(mol)
    except Exception:
        return False
    try:
        use_mmff = bool(AllChem.MMFFHasAllMoleculeParams(probe))
    except Exception:
        use_mmff = False
    return _build_force_field(probe, use_mmff) is not None


def _prepare_mol_for_force_field(docked_mol: Chem.Mol) -> Chem.Mol | None:
    """Returns a canonically ordered copy with explicit hydrogens and kept coords.

    Two properties are essential for a *comparable* strain measurement
    (Task Fix 3.1c):

    1. **Atom-order invariance** - the molecule is renumbered canonically before
       any search runs, so the rotamer/ensemble trajectories (and therefore the
       energies) no longer depend on the atom ordering produced by the prep or
       the library pipeline.
    2. **Force-field state preservation** - a full ``Chem.SanitizeMol`` is *not*
       run unconditionally. Sanitisation re-perceives RDKit aromaticity and can
       flip a Kekulé ring (e.g. the chromone pyranone ring of baicalin) back to
       "aromatic", which shifts MMFF94s energies by ~18 kcal/mol relative to the
       prep ``min_energy_kcal`` baseline and produced the residual library-path
       false-strain veto. The caller's bond representation is kept; sanitisation
       is only applied as a fallback when no force field can be built otherwise.
    """
    try:
        mol = Chem.Mol(docked_mol)
    except Exception:
        return None
    if mol is None:
        return None

    mol = _canonicalize_atom_order(mol)

    # Primary path: use the molecule exactly as delivered. Only if no force field
    # can be built do we touch its state (property cache / ring info, then full
    # sanitisation) - both of those can re-perceive aromaticity, so they are
    # strictly fallbacks. NOTE: even ``Chem.FastFindRings`` can rewrite bond
    # aromaticity on some prep molecules, so it is never applied in-place here.
    if not _force_field_available(mol):
        candidate = Chem.Mol(mol)
        try:
            candidate.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(candidate)
        except Exception:
            pass
        if _force_field_available(candidate):
            mol = candidate
        else:
            sanitized = _sanitized_copy(mol)
            if sanitized is not None and _force_field_available(sanitized):
                mol = sanitized

    # ``Chem.AddHs`` re-perceives aromaticity as a side effect, so it is only
    # called when hydrogens are actually missing (a molecule that already carries
    # its full hydrogen set must be left untouched).
    try:
        has_explicit_h = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
        has_implicit_h = sum(a.GetNumImplicitHs() for a in mol.GetAtoms()) > 0
    except Exception:
        has_explicit_h = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
        has_implicit_h = False
    if has_implicit_h or not has_explicit_h:
        try:
            mol = Chem.AddHs(mol, addCoords=True)
        except Exception:
            try:
                mol = Chem.AddHs(mol)
            except Exception:
                pass
    return mol


def _build_force_field(
    mol: Chem.Mol,
    use_mmff: bool,
    conf_id: int = -1,
) -> Any | None:
    """Builds an MMFF94s (preferred) or UFF force field for the same topology."""
    try:
        if use_mmff:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            if props is None:
                return None
            return AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
        return AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    except Exception:
        return None


def _relaxed_energy(mol: Chem.Mol, use_mmff: bool, conf_id: int = -1) -> float | None:
    """Unconstrained minimisation energy of one conformer with the shared force field."""
    ff = _build_force_field(mol, use_mmff, conf_id=conf_id)
    if ff is None:
        return None
    try:
        ff.Minimize(maxIts=_STRAIN_MINIMIZE_MAX_ITS)
        energy = float(ff.CalcEnergy())
    except Exception:
        return None
    return energy if math.isfinite(energy) else None


def _search_conformational_global_minimum(mol: Chem.Mol, use_mmff: bool) -> float | None:
    """ETKDGv3 conformer ensemble search for the global relaxed minimum.

    Contract (Task Fix 3.1): when no ``global_min_energy`` is supplied, relaxing
    only the *current* pose can get trapped in a nearby local minimum and thus
    severely underestimate the strain of e.g. a high-energy rotamer. A 10-member
    deterministic ETKDGv3 ensemble is therefore embedded on a copy of the
    molecule and every member is relaxed with the *same* force field; the lowest
    energy is returned (``None`` when embedding/FF fails for all members).
    """
    mol_copy = Chem.Mol(mol)
    cids: list[int] = []
    for use_random_coords in (False, True):
        try:
            mol_copy.RemoveAllConformers()
        except Exception:
            return None
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = _STRAIN_EMBED_RANDOM_SEED
            if use_random_coords:
                # Fallback prescribed by the contract for rigid/unusual skeletons
                params.useRandomCoords = True
            cids = list(
                AllChem.EmbedMultipleConfs(
                    mol_copy,
                    numConfs=_STRAIN_EMBED_NUM_CONFS,
                    params=params,
                )
            )
        except Exception:
            cids = []
        if cids:
            break

    energies: list[float] = []
    for cid in cids:
        energy = _relaxed_energy(mol_copy, use_mmff, conf_id=cid)
        if energy is not None:
            energies.append(energy)
    return min(energies) if energies else None


def _polar_hydrogen_rotamers(mol: Chem.Mol) -> list[tuple[int, int, int]]:
    """Rotatable polar hydrogens as ``(h_idx, donor_heavy_idx, upstream_heavy_idx)``.

    A rotatable polar hydrogen is an explicit H bonded to O/N/S whose donor heavy
    atom has at least one heavy neighbour (so an O-H / N-H / S-H torsion axis
    exists). The list is ordered by hydrogen atom index, which keeps the search
    fully deterministic.
    """
    rotamers: list[tuple[int, int, int]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        neighbours = atom.GetNeighbors()
        if len(neighbours) != 1:
            continue
        donor = neighbours[0]
        if donor.GetAtomicNum() not in _POLAR_DONOR_ELEMENTS:
            continue
        heavy_neighbours = sorted(
            n.GetIdx() for n in donor.GetNeighbors() if n.GetAtomicNum() > 1
        )
        if not heavy_neighbours:
            continue
        rotamers.append((atom.GetIdx(), donor.GetIdx(), heavy_neighbours[0]))
    return rotamers


def _rotate_point_about_axis(
    position: Any,
    origin: Any,
    axis: Any,
    angle: float,
) -> np.ndarray:
    """Rodrigues rotation of ``position`` about the (origin, axis) line."""
    vector = np.asarray(position, dtype=float) - np.asarray(origin, dtype=float)
    direction = np.asarray(axis, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm < 1e-9:
        return np.asarray(position, dtype=float)
    direction = direction / norm
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    rotated = (
        vector * cos_a
        + np.cross(direction, vector) * sin_a
        + direction * float(np.dot(direction, vector)) * (1.0 - cos_a)
    )
    return np.asarray(origin, dtype=float) + rotated


def _hydrogen_relaxation_force_field(mol: Chem.Mol, use_mmff: bool) -> Any | None:
    """Force field with every heavy atom frozen (only hydrogens may move)."""
    ff = _build_force_field(mol, use_mmff)
    if ff is None:
        return None
    try:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 1:
                ff.AddFixedPoint(atom.GetIdx())
        ff.Initialize()
    except Exception:
        return None
    return ff


def _relax_and_energy(ff: Any, max_its: int) -> float | None:
    """Minimises ``ff`` and returns the (finite) energy, or ``None``."""
    try:
        ff.Minimize(maxIts=max_its)
        energy = float(ff.CalcEnergy())
    except Exception:
        return None
    return energy if math.isfinite(energy) else None


def _optimise_docked_hydrogen_network(
    mol: Chem.Mol,
    use_mmff: bool,
    baseline: float | None,
) -> float | None:
    """Best heavy-atom-frozen energy over polar-hydrogen rotamer states.

    The docked pose's heavy-atom coordinates stay authoritative (frozen); only
    the polar hydrogens (O-H / N-H / S-H) are re-oriented so that the pose is
    not punished for a hydrogen-bond network that ``AddHs(addCoords=True)``
    happened to place in an arbitrary orientation (Task Fix 3.1b).

    Three deterministic stages, all sharing one force field and one topology:

    1. Local relaxation of all hydrogens from the pose as given.
    2. Coordinate descent over the polar-H rotamers: each donor hydrogen is spun
       around its donor-heavy/upstream-heavy axis on a fixed 30 deg grid, the
       best single-point rotamer is relaxed with all hydrogens free, and the
       move is kept only when the relaxed energy improves.
    3. Multi-start fallback when the residual above ``baseline`` is still large
       (strongly coupled networks, e.g. adjacent catechol OH groups): every
       donor is rotated by a deterministic offset pattern and relaxed, and the
       lowest energy is kept.

    Returns ``None`` when no force field can be built. No random numbers are
    involved, so identical input always produces the identical output.
    """
    work = Chem.Mol(mol)
    ff = _hydrogen_relaxation_force_field(work, use_mmff)
    if ff is None:
        return None
    best_energy = _relax_and_energy(ff, _STRAIN_MINIMIZE_MAX_ITS)
    if best_energy is None:
        return None
    minimizations = 1

    try:
        conf = work.GetConformer()
        start_positions = np.array(conf.GetPositions(), dtype=float)
    except Exception:
        return best_energy
    snapshot = start_positions.copy()

    def _restore_snapshot() -> None:
        for idx in range(snapshot.shape[0]):
            conf.SetAtomPosition(idx, snapshot[idx].tolist())

    donors = _polar_hydrogen_rotamers(work)[:_STRAIN_MAX_POLAR_H_DONORS]
    if not donors:
        return best_energy

    step_angle = 2.0 * math.pi / _STRAIN_H_ROTAMER_STEPS

    # ---- Stage 2: deterministic per-donor rotamer coordinate descent --------
    for h_idx, donor_idx, axis_idx in donors:
        if minimizations >= _STRAIN_MAX_MINIMIZATIONS:
            break
        origin = conf.GetAtomPosition(donor_idx)
        axis = np.asarray(conf.GetAtomPosition(axis_idx)) - np.asarray(origin)
        start_pos = np.asarray(conf.GetAtomPosition(h_idx))
        best_scan_energy = None
        best_scan_pos = start_pos
        for step in range(1, _STRAIN_H_ROTAMER_STEPS + 1):
            position = _rotate_point_about_axis(
                start_pos, origin, axis, step * step_angle
            )
            conf.SetAtomPosition(h_idx, position.tolist())
            try:
                scan_energy = float(ff.CalcEnergy())
            except Exception:
                scan_energy = math.inf
            if best_scan_energy is None or scan_energy < best_scan_energy:
                best_scan_energy = scan_energy
                best_scan_pos = position
        conf.SetAtomPosition(h_idx, np.asarray(best_scan_pos).tolist())
        relaxed = _relax_and_energy(ff, _STRAIN_SCAN_MIN_ITS)
        minimizations += 1
        if relaxed is not None and relaxed < best_energy - 1e-4:
            best_energy = relaxed
            snapshot = np.array(conf.GetPositions(), dtype=float)
        else:
            _restore_snapshot()

    # ---- Stage 3: coordinated multi-start for coupled H-bond networks -------
    residual = best_energy - baseline if _is_finite_number(baseline) else None
    if residual is not None and residual > _STRAIN_FALLBACK_RESIDUAL:
        multistart_budget = min(
            _STRAIN_MAX_MULTISTART_MINIMIZATIONS,
            max(1, _STRAIN_MAX_MINIMIZATIONS - minimizations),
        )
        minimizations_used = 0
        lock_step_energies: list[float] = []

        def _reset_to_start() -> None:
            for idx in range(start_positions.shape[0]):
                conf.SetAtomPosition(idx, start_positions[idx].tolist())

        def _rotate_donors(offset_for_donor: Any) -> None:
            for donor_pos, (h_idx, donor_idx, axis_idx) in enumerate(donors):
                origin = conf.GetAtomPosition(donor_idx)
                axis = np.asarray(conf.GetAtomPosition(axis_idx)) - np.asarray(origin)
                start_pos = np.asarray(conf.GetAtomPosition(h_idx))
                offset = offset_for_donor(donor_pos)
                position = _rotate_point_about_axis(start_pos, origin, axis, offset)
                conf.SetAtomPosition(h_idx, position.tolist())

        def _try_state() -> float | None:
            nonlocal best_energy, snapshot
            energy = _relax_and_energy(ff, _STRAIN_MINIMIZE_MAX_ITS)
            if energy is not None and energy < best_energy - 1e-4:
                best_energy = energy
                snapshot = np.array(conf.GetPositions(), dtype=float)
            return energy

        # Sweep A - lock-step patterns: every donor is rotated by the *same*
        # offset, so the sweep is completely independent of the donor ordering.
        for pattern in range(1, _STRAIN_MULTISTART_PATTERNS):
            if minimizations_used >= multistart_budget:
                break
            _reset_to_start()
            offset = pattern * (2.0 * math.pi / _STRAIN_MULTISTART_PATTERNS)
            _rotate_donors(lambda _donor_pos, _offset=offset: _offset)
            energy = _try_state()
            minimizations += 1
            minimizations_used += 1
            if energy is not None:
                lock_step_energies.append(energy)

        # Convergence criterion (Task Fix 3.1c): if the two best lock-step states
        # still differ by more than the fallback threshold, the coupled network is
        # not resolved - add an independent staggered sweep.
        if len(lock_step_energies) >= 2:
            ordered = sorted(lock_step_energies)
            lock_step_gap = ordered[1] - ordered[0]
        else:
            lock_step_gap = math.inf if len(lock_step_energies) < 2 else 0.0
        if lock_step_gap > _STRAIN_FALLBACK_RESIDUAL:
            for pattern in range(1, _STRAIN_MULTISTART_PATTERNS):
                if minimizations_used >= multistart_budget:
                    break
                _reset_to_start()

                def _staggered(donor_pos: int, _pattern: int = pattern) -> float:
                    step_index = (_pattern * (donor_pos + 1) * 5) % _STRAIN_H_ROTAMER_STEPS
                    return step_index * step_angle

                _rotate_donors(_staggered)
                _try_state()
                minimizations += 1
                minimizations_used += 1

    return best_energy


def _ring_synchronised_copy(mol: Chem.Mol) -> Chem.Mol | None:
    """Copy whose property cache and ring info are (re)synchronised.

    ``Chem.FastFindRings`` may rewrite bond aromaticity flags on molecules whose
    ring bookkeeping came from a different pipeline. This copy is used only as a
    force-field-state fallback (see :func:`_state_consistent_docked_energy`).
    """
    try:
        work = Chem.Mol(mol)
        work.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(work)
        return work
    except Exception:
        return None


def _state_consistent_docked_energy(
    mol: Chem.Mol,
    use_mmff: bool,
    baseline: float | None,
) -> float | None:
    """Heavy-atom-frozen energy that is comparable to the caller's baseline.

    The prep ``min_energy_kcal`` baseline is produced by the ligand preparation
    pipeline, whose repeated ``AddHs`` / ETKDG / MMFF calls progressively rewrite
    the aromaticity bookkeeping of the molecule (RDKit kekulises the chromone
    pyranone ring on the way). A pose restored from that ligand's PDBQT is
    usually in a *different* bookkeeping state, which shifts MMFF94s energies by
    ~18 kcal/mol even at identical heavy-atom coordinates (Task Fix 3.1c).

    The delivered state is always evaluated first. Only when the residual above
    the baseline is suspiciously large *and* fully relaxing the pose does not
    land back near the baseline (the signature of a state mismatch rather than
    real strain) is the molecule re-evaluated in a ring-synchronised state; the
    lower of the two energies is kept so a genuine mismatch cannot fake strain.
    """
    best = _optimise_docked_hydrogen_network(mol, use_mmff, baseline)
    if best is None:
        return None
    if not _is_finite_number(baseline):
        return best
    if best - baseline <= _STRAIN_FALLBACK_RESIDUAL:
        return best

    # Distinguish genuine strain from a force-field state mismatch: relaxing the
    # pose back to its own minimum in the same state must land near the baseline
    # when the state is consistent with the prep baseline.
    relaxed = _relaxed_energy(Chem.Mol(mol), use_mmff)
    if relaxed is not None and abs(relaxed - baseline) <= _STRAIN_STATE_TOLERANCE:
        return best

    alternative = _ring_synchronised_copy(mol)
    if alternative is None:
        return best
    alternative_energy = _optimise_docked_hydrogen_network(alternative, use_mmff, baseline)
    if alternative_energy is not None and alternative_energy < best - 1e-4:
        return alternative_energy
    return best


def compute_internal_strain_energy(
    docked_mol: Chem.Mol,
    global_min_energy: float | None = None,
    affinity: float = 0.0,
    high_strain_threshold: float = 8.0,
    force_constant: float = 100.0,
    **kwargs: Any,
) -> tuple[float, float, bool]:
    """Calculates internal conformational strain energy ΔE_strain and corrected affinity for a 3D Mol.

    Freezes all heavy atoms in their docked 3D coordinates, optimises the
    **polar-hydrogen (H-bond) network** deterministically, evaluates MMFF94s
    (UFF fallback) and computes strain relative to the relaxed conformational
    minimum evaluated with the very same force field.

    Both ``e_docked`` and ``global_min_energy`` are therefore guaranteed to share
    one force field (MMFF94s preferred, UFF fallback) and one atom topology
    (canonically ordered molecule with explicit hydrogens; the caller's bond /
    aromaticity representation is preserved so the prep baseline stays comparable
    - see :func:`_prepare_mol_for_force_field`):

    * ``e_docked``: energy of the docked pose with the heavy atoms frozen and the
      polar-hydrogen rotamers optimised (see
      :func:`_optimise_docked_hydrogen_network`). This removes the spurious
      4-26 kcal/mol "strain" that a PDBQT round trip used to inject into
      polyphenols by re-generating hydrogens outside their H-bond network
      (Task Fix 3.1b) while genuine heavy-atom torsional strain is preserved.
    * ``global_min_energy``: caller-supplied baseline when finite; otherwise the
      minimum of (a) the fully relaxed *current* conformer and (b) the lowest
      member of a deterministic 10-conformer ETKDGv3 ensemble relaxed with the
      same force field. This prevents the classic local-minimum underestimation.

    Args:
        docked_mol: RDKit Mol with 3D conformer coordinates.
        global_min_energy: Baseline global minimum energy in kcal/mol. If None
            (or NaN/non-finite), searched via conformer ensemble relaxation.
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

    mol = _prepare_mol_for_force_field(docked_mol)
    if mol is None or mol.GetNumConformers() == 0:
        return 0.0, raw_affinity, False

    # Re-evaluate MMFF availability *after* hydrogen completion: the hydrogen
    # topology changes atom typing, so a pre-hydrogen decision is not reliable.
    # Both probes run on copies - MMFF property construction mutates the
    # molecule's aromaticity bookkeeping and must not disturb the pose state.
    use_mmff = bool(AllChem.MMFFHasAllMoleculeParams(Chem.Mol(mol)))
    if not _force_field_available(mol):
        return 0.0, raw_affinity, False

    baseline: float
    if not _is_finite_number(global_min_energy):
        # No trustworthy baseline: search the conformational space instead of
        # only locally relaxing the current pose (which traps local minima).
        # NOTE: work on copies - the docked conformer must stay untouched.
        candidates: list[float] = []
        local_relaxed = _relaxed_energy(Chem.Mol(mol), use_mmff)
        if local_relaxed is not None:
            candidates.append(local_relaxed)
        ensemble_min = _search_conformational_global_minimum(mol, use_mmff)
        if ensemble_min is not None:
            candidates.append(ensemble_min)
        baseline = min(candidates) if candidates else math.nan
    else:
        baseline = float(global_min_energy)

    e_docked = _state_consistent_docked_energy(mol, use_mmff, baseline)
    if e_docked is None or not math.isfinite(e_docked):
        return 0.0, raw_affinity, False

    if not math.isfinite(baseline):
        baseline = e_docked
    # Cross-baseline guard: a baseline above the docked single-point energy
    # cannot come from the same force field/topology (e.g. the caller mixed
    # an MMFF baseline with a UFF pose). Clamping keeps ΔE >= 0 and makes
    # the max(0, ...) semantics explicit.
    if baseline > e_docked:
        baseline = e_docked

    return compute_internal_strain_energy_penalty(
        docked_energy=e_docked,
        global_min_energy=baseline,
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
