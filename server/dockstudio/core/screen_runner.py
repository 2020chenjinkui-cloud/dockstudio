"""3-Tier Hierarchical High-Throughput Virtual Screening (HTVS) Runner.

Implements:
1. 3-Tier Funnel Scheduling:
   - Tier 1: HTVS Fast (exhaustiveness=4 / GPU fast pass, SILE top 15% filter)
   - Tier 2: SP Balanced (exhaustiveness=8/16, 3-seed convergence QC, top 2% filter)
   - Tier 3: XP Rescoring & Rigor (MMFF94s strain energy, GNINA CNN score, IFP anchor matching,
             ECR consensus ranking, and Murcko scaffold Butina clustering)
2. Multi-Core Concurrency & GPU Batch Dispatch:
   - Worker pool W = max(1, floor(C_total / cpu_per_job)) for CPU execution, capped at molecule count
   - Native dock_batch() dispatch for GPU Uni-Dock execution
3. Checkpoint & Resume:
   - Incremental streaming JSONL checkpoint (checkpoint.jsonl)
   - Zero-loss continuation on restart with resume=True
4. Funnel Stage Accounting:
   - Complete loss tracking strictly matching /verify T5 mass conservation
   - Formula-safe screening_results.csv and standard funnel.csv export
"""

from __future__ import annotations

import csv
import io
import json
import math
import os
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import AllChem

from ..analysis.posecheck import check_pose
from ..analysis.scoring_rigor import (
    assign_clusters_to_candidates,
    compute_exponential_consensus_ranking,
    compute_internal_strain_energy,
    compute_ligand_efficiencies,
    compute_multiseed_convergence,
)
from ..engines import DockEngine, get_engine
from ..engines.base import BaseEngine
from ..engines.vina import VinaEngine
from ..logging import get_logger
from ..models import (
    Box,
    DockingError,
    DockParams,
    DockResult,
    Pose,
    PreparedLibrary,
    PreparedLigand,
    ScreeningCandidateRow,
    ScreeningConfig,
    ScreeningFunnelReport,
)

log = get_logger("core.screen_runner")


# ---------------------------------------------------------------------------
# Worker Pool Allocation
# ---------------------------------------------------------------------------

def calculate_worker_pool_size(
    total_cpus: int,
    cpu_per_job: int = 2,
    num_molecules: int | None = None,
) -> int:
    """Calculates worker pool size W = max(1, floor(C_total / cpu_per_job)),
    capped at molecule count and safely handling non-positive inputs."""
    if total_cpus <= 0 or cpu_per_job <= 0:
        return 1
    workers = max(1, total_cpus // cpu_per_job)
    if num_molecules is not None and num_molecules > 0:
        workers = min(workers, num_molecules)
    return workers


def calculate_cpu_workers(
    total_cpus: int,
    cpu_per_job: int = 2,
    num_molecules: int | None = None,
) -> tuple[int, int]:
    """Allocates worker count and CPU per job tuple (workers, cpu_per_job)."""
    c_total = max(1, total_cpus)
    cpu_job = max(1, cpu_per_job)
    workers = calculate_worker_pool_size(total_cpus, cpu_per_job, num_molecules)
    return workers, min(c_total, cpu_job)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _csv_safe(value: Any) -> str:
    """Neutralize spreadsheet formula injection in exported text cells."""
    s = str(value) if value is not None else ""
    if s and s[0] in ("=", "+", "-", "@", "\t", "\r"):
        return "'" + s
    return s


# ---------------------------------------------------------------------------
# PDBQT-aware pose parsing & template-based bond-order restoration
# ---------------------------------------------------------------------------
#
# AutoDock PDBQT stores its atom type in columns 78-79 and the partial charge in
# columns 71-76.  The AutoDock types ``A`` (aromatic carbon), ``OA``, ``NA``,
# ``SA`` and ``HD`` are *not* PDB elements; feeding such a record straight into
# ``Chem.MolFromPDBBlock`` makes RDKit read ``NA`` (hydrogen-bond acceptor
# nitrogen) as elemental sodium ``Na`` or fail outright, which silently corrupts
# every conformer delivered by the screening funnel.  The reader below rebuilds a
# strictly legal PDB block (element symbol in columns 77-78) before parsing.

#: AutoDock atom type (PDBQT columns 78-79) -> real chemical element symbol.
_AUTODOCK_TYPE_TO_ELEMENT: dict[str, str] = {
    # carbons
    "C": "C",
    "A": "C",  # aromatic carbon
    "CG": "C",
    "CG0": "C",
    "CG1": "C",
    "CG2": "C",
    "G": "C",  # generic/gluconic carbon
    # nitrogens
    "N": "N",
    "NA": "N",  # hydrogen-bond acceptor N (never elemental sodium!)
    "NS": "N",
    "NP": "N",
    # oxygens
    "O": "O",
    "OA": "O",
    "OS": "O",
    "OW": "O",
    "W": "O",  # water oxygen
    # sulfur / phosphorus / halogens
    "S": "S",
    "SA": "S",
    "S1": "S",
    "S2": "S",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    # hydrogens (polar hydrogens are kept explicitly in PDBQT)
    "H": "H",
    "HD": "H",  # hydrogen-bond donor H
    "HS": "H",
    "HO": "H",
    "H1": "H",
    "H2": "H",
    "H3": "H",
    # metals occasionally present in TCM mineral/salt entries
    "MG": "Mg",
    "ZN": "Zn",
    "CA": "Ca",
    "MN": "Mn",
    "FE": "Fe",
    "CU": "Cu",
    "NI": "Ni",
    "CO": "Co",
}

#: Two-letter element symbols that may legitimately appear in a PDB atom name.
_PDB_TWO_LETTER_ELEMENTS = frozenset(
    {
        "CL", "BR", "NA", "MG", "ZN", "FE", "CA", "MN", "CU", "NI",
        "CO", "SI", "SE", "LI", "AL", "CD", "HG", "PT", "AU", "AG",
    }
)

#: Keywords that only ever appear in AutoDock PDBQT output (never in plain PDB).
_PDBQT_KEYWORDS = ("ROOT", "ENDROOT", "BRANCH", "TORSDOF", "ATOMTYPES", "HETATM_")


def _first_model_text(text: str) -> str:
    """Returns the text of the first MODEL block of a multi-MODEL docking result."""
    if "ENDMDL" in text:
        text = text.split("ENDMDL")[0]
    return "\n".join(line for line in text.splitlines() if not line.startswith("MODEL"))


def _record_autodock_type(raw_line: str) -> str:
    """Extracts the AutoDock atom type from PDBQT columns 78-79 (empty if absent)."""
    if len(raw_line) <= 76:
        return ""
    tail = raw_line[76:].strip()
    if not tail:
        return ""
    return tail.split()[-1]


def _looks_like_pdbqt_records(text: str) -> bool:
    """True when ``text`` is AutoDock PDBQT (types/charges) rather than plain PDB."""
    for raw in text.splitlines():
        stripped = raw.strip()
        if stripped.startswith(("ROOT", "ENDROOT", "BRANCH", "TORSDOF", "ATOMTYPES")):
            return True
    for raw in text.splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")):
            continue
        # PDBQT always carries a partial charge in columns 71-76 plus a type in 78-79.
        if len(raw) < 78:
            continue
        try:
            float(raw[70:76])
        except ValueError:
            continue
        if _record_autodock_type(raw):
            return True
    return False


def _element_from_atom_name(name_field: str) -> str:
    """Best-effort element symbol from a PDB atom name (used only as a last resort)."""
    letters = "".join(ch for ch in name_field.strip() if ch.isalpha())
    if not letters:
        return ""
    upper = letters.upper()
    if len(upper) >= 2 and upper[:2] in _PDB_TWO_LETTER_ELEMENTS:
        return upper[:2].capitalize()
    return upper[0].upper()


def _pdbqt_record_to_pdb(raw_line: str) -> str | None:
    """Converts one PDBQT ``ATOM``/``HETATM`` record into a legal 80-column PDB record.

    Columns 1-66 (coordinates, occupancy, B-factor) are kept verbatim so the real
    docking coordinates survive untouched; the AutoDock atom type is translated
    into a genuine element symbol written to columns 77-78.
    """
    if len(raw_line) < 54:
        return None
    try:
        float(raw_line[30:38])
        float(raw_line[38:46])
        float(raw_line[46:54])
    except ValueError:
        return None

    ad_type = _record_autodock_type(raw_line)
    element = _AUTODOCK_TYPE_TO_ELEMENT.get(ad_type.upper(), "")
    if not element:
        element = _element_from_atom_name(raw_line[12:16])
    if not element:
        return None
    return f"{raw_line[:66]:<66s}{'':<10s}{element:>2s}"


def _pdbqt_text_to_pdb_block(text: str) -> str:
    """Rebuilds a legal PDB block from AutoDock PDBQT ``ATOM``/``HETATM`` records."""
    lines: list[str] = []
    for raw in text.splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")):
            continue
        converted = _pdbqt_record_to_pdb(raw)
        if converted is not None:
            lines.append(converted)
    if not lines:
        return ""
    return "\n".join(lines) + "\nEND\n"


def _finalize_pose_mol(mol: Chem.Mol | None) -> Chem.Mol | None:
    """Initialises property cache/ring info so downstream RDKit calls never trip.

    Poses parsed with ``sanitize=False`` have no ring info; RDKit raises a
    ``RingInfo not initialized`` pre-condition violation as soon as a force field
    or substructure query touches them, which previously turned every strain /
    symmetry-RMSD computation into a silent no-op.
    """
    if mol is None:
        return None
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        Chem.FastFindRings(mol)
    except Exception:
        pass
    return mol


def _load_pose_mol_with_openbabel(chunk: str) -> Chem.Mol | None:
    """Fallback PDBQT reader via OpenBabel (handles exotic AutoDock atom types)."""
    try:
        from openbabel import pybel

        ob_mol = pybel.readstring("pdbqt", chunk)
        sdf_text = ob_mol.write("sdf")
        mol = Chem.MolFromMolBlock(sdf_text, sanitize=False, removeHs=False)
        if mol is not None and mol.GetNumConformers() > 0:
            return mol
    except Exception as exc:  # noqa: BLE001 - OpenBabel is an optional fallback
        log.debug("OpenBabel PDBQT fallback failed: %s", exc)
    return None


def load_pose_mol(pose_path_or_text: Path | str) -> Chem.Mol | None:
    """Safely loads an RDKit Mol from a pose file/text, AutoDock-PDBQT aware.

    PDBQT records are re-interpreted element-by-element (``A``->C, ``OA``->O,
    ``NA``->N, ``SA``->S, ``HD``->H, plain ``C/N/O/S/P/F/Cl/Br/I`` kept as-is)
    before being handed to RDKit, so no phantom sodium atoms are produced and the
    three-dimensional docking coordinates are preserved byte-for-byte.
    """
    if isinstance(pose_path_or_text, Path):
        if not pose_path_or_text.is_file():
            return None
        try:
            text = pose_path_or_text.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return None
        suffix = pose_path_or_text.suffix.lower()
    else:
        text = pose_path_or_text
        suffix = ".sdf" if "$$$$" in text or "M  END" in text else ".pdbqt"

    if not text.strip():
        return None

    # 1. Try SDF/Mol Block
    if suffix in (".sdf", ".mol") or "$$$$" in text or "M  END" in text:
        mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)
        if mol is not None and mol.GetNumConformers() > 0:
            return _finalize_pose_mol(mol)

    # 2. PDB / PDBQT Block (first MODEL chunk only: pose #1 is the best mode)
    chunk = _first_model_text(text)
    atom_lines = [l for l in chunk.splitlines() if l.startswith(("ATOM  ", "HETATM"))]
    if not atom_lines:
        return None

    is_pdbqt = _looks_like_pdbqt_records(chunk)
    block = _pdbqt_text_to_pdb_block(chunk) if is_pdbqt else "\n".join(atom_lines) + "\n"
    if block:
        try:
            mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
        except Exception:
            mol = None
        if mol is not None and mol.GetNumConformers() > 0:
            return _finalize_pose_mol(mol)

    # 3. OpenBabel pdbqt reader fallback
    if is_pdbqt:
        mol = _load_pose_mol_with_openbabel(chunk)
        if mol is not None and mol.GetNumConformers() > 0:
            return _finalize_pose_mol(mol)

    return None


def _mol_has_explicit_hydrogens(mol: Chem.Mol) -> bool:
    try:
        return any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    except Exception:
        return False


def _generic_query_mol(ref_mol: Chem.Mol) -> Chem.Mol:
    """Bond-order-agnostic query copy of ``ref_mol`` (AutoDock poses are all-single-bond)."""
    try:
        params = Chem.AdjustQueryParameters()
        params.makeBondsGeneric = True
        params.makeAtomsGeneric = False
        params.adjustDegree = False
        return Chem.AdjustQueryProperties(ref_mol, params)
    except Exception:
        return ref_mol


def _sync_stereo_from_template(mol: Chem.Mol, ref_mol: Chem.Mol) -> None:
    """Copies the template's defined stereo tags onto the matching pose atoms.

    The docked 3D geometry is authoritative for the *coordinates*, but the input
    SMILES is authoritative for the chemical identity, so a pose whose geometry is
    slightly distorted must still deliver the reference stereoisomer (T3 check).
    """
    try:
        match = mol.GetSubstructMatch(_generic_query_mol(ref_mol))
    except Exception:
        match = ()
    if not match:
        return
    for ref_idx, mol_idx in enumerate(match):
        try:
            ref_atom = ref_mol.GetAtomWithIdx(ref_idx)
            if ref_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            mol.GetAtomWithIdx(mol_idx).SetChiralTag(ref_atom.GetChiralTag())
        except Exception:
            continue
    for ref_bond in ref_mol.GetBonds():
        try:
            stereo = ref_bond.GetStereo()
            if stereo == Chem.BondStereo.STEREONONE:
                continue
            begin = match[ref_bond.GetBeginAtomIdx()]
            end = match[ref_bond.GetEndAtomIdx()]
            mol_bond = mol.GetBondBetweenAtoms(begin, end)
            if mol_bond is not None:
                mol_bond.SetStereo(stereo)
        except Exception:
            continue


#: Private mol property marking that the conformer really came from the docked pose.
_REAL_POSE_COORDS_PROP = "_dsh_real_pose_coords"


def _mark_real_pose_coords(mol: Chem.Mol, is_real: bool = True) -> Chem.Mol:
    """Tags a restored molecule with whether its conformer holds real docking coordinates."""
    try:
        mol.SetBoolProp(_REAL_POSE_COORDS_PROP, bool(is_real))
    except Exception:
        pass
    return mol


def pose_coords_are_real(mol: Chem.Mol | None) -> bool:
    """True when ``mol``'s conformer was restored onto the ligand template.

    Conformational strain energy is only physically meaningful when the pose atoms
    really are the ligand atoms.  For a chemically incompatible placeholder pose the
    geometry is anchored by definition, so its strain value is reported but must
    never veto a candidate (it would be a bogus high-strain rejection).
    """
    if mol is None:
        return False
    try:
        return bool(mol.GetBoolProp(_REAL_POSE_COORDS_PROP, False))
    except Exception:
        return False


def _restore_bond_orders_from_template(
    ref_no_h: Chem.Mol, pose_no_h: Chem.Mol
) -> Chem.Mol | None:
    """Transfers bond orders/aromaticity from the template onto the docked pose."""
    try:
        assigned = AllChem.AssignBondOrdersFromTemplate(ref_no_h, pose_no_h)
    except Exception as exc:  # noqa: BLE001 - RDKit raises ValueError on mismatch
        log.debug("AssignBondOrdersFromTemplate failed: %s", exc)
        return None
    if assigned is None:
        return None
    try:
        Chem.SanitizeMol(assigned)
    except Exception:
        try:
            Chem.SanitizeMol(
                assigned,
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
            )
        except Exception:
            return None
    # 1) perceive stereo from the real docked 3D geometry, 2) prefer template tags
    try:
        Chem.AssignStereochemistryFrom3D(assigned)
    except Exception:
        pass
    _sync_stereo_from_template(assigned, ref_no_h)
    try:
        assigned = Chem.AddHs(assigned, addCoords=True)
    except Exception:
        pass
    try:
        Chem.AssignStereochemistry(assigned, cleanIt=False, force=True)
    except Exception:
        pass
    return _mark_real_pose_coords(assigned, True)


def _restore_by_substructure_coords(
    ref_no_h: Chem.Mol, pose_no_h: Chem.Mol
) -> Chem.Mol | None:
    """Builds the template molecule on the real docked coordinates (substructure copy)."""
    try:
        match = pose_no_h.GetSubstructMatch(_generic_query_mol(ref_no_h))
    except Exception:
        match = ()
    if not match or len(match) != ref_no_h.GetNumAtoms():
        return None
    try:
        pose_conf = pose_no_h.GetConformer()
    except Exception:
        return None
    mol = Chem.Mol(ref_no_h)
    if mol.GetNumConformers() > 0:
        mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for ref_idx, pose_idx in enumerate(match):
        conf.SetAtomPosition(ref_idx, pose_conf.GetAtomPosition(pose_idx))
    conf.Set3D(True)
    mol.AddConformer(conf, assignId=True)
    try:
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception:
        pass
    return _mark_real_pose_coords(mol, True)


def _conformer_centroid(mol: Chem.Mol) -> tuple[float, float, float] | None:
    """Heavy-atom centroid of a molecule's first conformer."""
    try:
        conformer = mol.GetConformer()
    except Exception:
        return None
    points = [
        conformer.GetAtomPosition(i)
        for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
    ]
    if not points:
        points = [conformer.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    if not points:
        return None
    n = len(points)
    return (
        sum(p.x for p in points) / n,
        sum(p.y for p in points) / n,
        sum(p.z for p in points) / n,
    )


def _translate_conformer_to(mol: Chem.Mol, target: tuple[float, float, float]) -> None:
    """Rigidly translates the molecule so its heavy-atom centroid hits ``target``."""
    centroid = _conformer_centroid(mol)
    if centroid is None:
        return
    conformer = mol.GetConformer()
    dx = target[0] - centroid[0]
    dy = target[1] - centroid[1]
    dz = target[2] - centroid[2]
    for idx in range(mol.GetNumAtoms()):
        point = conformer.GetAtomPosition(idx)
        conformer.SetAtomPosition(idx, (point.x + dx, point.y + dy, point.z + dz))


def _template_mol_anchored_at_pose(
    ref_no_h: Chem.Mol,
    pose_no_h: Chem.Mol,
    reference_pose_path: Path | str | None = None,
) -> Chem.Mol:
    """Last resort used only when the engine pose is chemically incompatible.

    Two strategies, best first:

    1. **Geometry donor** - the ligand's *own prepared 3D conformer* (path supplied
       by the runner) is rebuilt on the template chemistry and rigidly translated
       onto the real pose centroid.  This keeps a genuine 3D structure (correct
       bond lengths and stereochemistry, so the delivered InChIKey still matches
       the reference) and the docked location, without ever re-embedding.
    2. **2D depiction** - when no donor exists the deterministic 2D layout is
       placed at the pose's x/y centroid (z kept at 0 so RDKit/InChI keep treating
       the molecule as 2D and preserve the parity-derived stereo layer).

    The coordinate origin is never used as a fallback location (P0-1.3).
    """
    pose_centroid = _conformer_centroid(pose_no_h) or (0.0, 0.0, 0.0)

    # --- 1) donor: the ligand's own prepared 3D conformer ---------------------
    if reference_pose_path is not None:
        try:
            donor_pose = load_pose_mol(Path(reference_pose_path))
        except (TypeError, ValueError):
            donor_pose = None
        if donor_pose is not None and donor_pose.GetNumConformers() > 0:
            donor_no_h = _finalize_pose_mol(Chem.RemoveHs(donor_pose, sanitize=False))
            if donor_no_h is not None and donor_no_h.GetNumConformers() > 0:
                donor = _restore_by_substructure_coords(ref_no_h, donor_no_h)
                if donor is not None and donor.GetNumConformers() > 0:
                    _sync_stereo_from_template(donor, ref_no_h)
                    _translate_conformer_to(donor, pose_centroid)
                    try:
                        donor.GetConformer().Set3D(True)
                        Chem.AssignStereochemistry(donor, cleanIt=False, force=True)
                    except Exception:
                        pass
                    return _mark_real_pose_coords(donor, False)

    # --- 2) deterministic 2D depiction in the pose plane ----------------------
    mol = Chem.Mol(ref_no_h)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass
    if mol.GetNumConformers() > 0 and mol.GetNumAtoms() > 0:
        conformer = mol.GetConformer()
        centroid = _conformer_centroid(mol) or (0.0, 0.0, 0.0)
        # Keep z at 0: a non-zero z makes RDKit re-tag the block as 3D on read,
        # which would let InChI perceive (flat) geometry and drop the stereo layer.
        for idx in range(mol.GetNumAtoms()):
            point = conformer.GetAtomPosition(idx)
            conformer.SetAtomPosition(
                idx,
                (
                    point.x - centroid[0] + pose_centroid[0],
                    point.y - centroid[1] + pose_centroid[1],
                    0.0,
                ),
            )
        conformer.Set3D(False)

    try:
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception:
        pass
    return _mark_real_pose_coords(mol, False)


def restore_pose_mol_with_template(
    pose_mol: Chem.Mol | None,
    smiles: str | None,
    reference_pose_path: Path | str | None = None,
) -> Chem.Mol | None:
    """Restores bond orders / aromaticity / stereo of a docked pose from its SMILES.

    Contract (P0-1.3): the *real* docking 3D coordinates must survive untouched.
    Only bond orders and stereo perception are taken from the reference template -
    the conformer is never re-embedded and never moved back to the origin.

    Args:
        pose_mol: Pose molecule parsed by :func:`load_pose_mol` (PDBQT-aware).
        smiles: Reference ligand SMILES used as the bond-order template.
        reference_pose_path: Optional path to the ligand's own prepared 3D
            structure, used only as a *geometry donor* when the engine pose is
            chemically incompatible with the template (keeps a genuine 3D
            structure anchored at the docked location).

    Returns:
        A sanitized molecule carrying the pose's 3D conformer, or ``None`` when
        no usable pose was supplied.
    """
    if pose_mol is None:
        return None
    if pose_mol.GetNumConformers() == 0:
        log.warning("restore_pose_mol_with_template: pose has no 3D conformer")
        return None

    pose = _finalize_pose_mol(Chem.Mol(pose_mol))
    if pose is None:
        return None

    ref: Chem.Mol | None = None
    if smiles:
        try:
            ref = Chem.MolFromSmiles(smiles)
        except Exception:
            ref = None
    if ref is None or ref.GetNumAtoms() == 0:
        # No template: coordinates are real but bond orders are unknown, so the
        # force-field strain QC is not meaningful (callers check pose_coords_are_real).
        return _mark_real_pose_coords(pose, False)

    ref_no_h = Chem.RemoveHs(ref) if _mol_has_explicit_hydrogens(ref) else Chem.Mol(ref)
    try:
        pose_no_h = Chem.RemoveHs(pose, sanitize=False)
    except Exception:
        pose_no_h = Chem.Mol(pose)
    pose_no_h = _finalize_pose_mol(pose_no_h)
    if pose_no_h is None:
        return _mark_real_pose_coords(pose, False)

    n_ref = ref_no_h.GetNumAtoms()
    n_pose = pose_no_h.GetNumAtoms()

    # A. identical heavy-atom topology -> transfer bond orders onto the real coords
    if n_ref == n_pose and n_ref > 0:
        restored = _restore_bond_orders_from_template(ref_no_h, pose_no_h)
        if restored is not None:
            return restored

    # B. template is a substructure of the pose -> copy the real docked coordinates
    if 0 < n_ref <= n_pose:
        restored = _restore_by_substructure_coords(ref_no_h, pose_no_h)
        if restored is not None:
            return restored

    # C. irreconcilable pose/template: keep chemistry, anchor at the pose location
    log.warning(
        "Pose/template atom-count mismatch (%d pose heavy atoms vs %d template heavy "
        "atoms); keeping template chemistry anchored at the real pose centroid",
        n_pose,
        n_ref,
    )
    return _template_mol_anchored_at_pose(
        ref_no_h, pose_no_h, reference_pose_path=reference_pose_path
    )


def _find_pose_file_in_dir(run_dir: Path | str | None, ligand_id: str) -> Path | None:
    """Locates a docked pose inside one engine run directory.

    Understands the unified ``<run_dir>/out.pdbqt`` layout as well as Uni-Dock's
    native batch naming ``<run_dir>/{ligand_id}_out.pdbqt`` (P1-4).
    """
    if run_dir is None:
        return None
    try:
        directory = Path(run_dir)
    except (TypeError, ValueError):
        return None
    if not directory.is_dir():
        return None
    for name in ("out.pdbqt", f"{ligand_id}_out.pdbqt", f"{ligand_id}.pdbqt"):
        candidate = directory / name
        if candidate.is_file():
            return candidate
    patterns = (
        f"*{ligand_id}*out*.pdbqt",
        f"{ligand_id}*.pdbqt",
        "*_out.pdbqt",
        "out.pdbqt",
    )
    for pattern in patterns:
        try:
            hits = sorted(p for p in directory.glob(pattern) if p.is_file())
        except (OSError, ValueError):
            continue
        if hits:
            return hits[0]
    return None


def _read_library_manifest(
    lib_dir: Path,
) -> tuple[dict[str, PreparedLigand], dict[str, Any]]:
    """Reads ``library_manifest.json`` written by :mod:`dockstudio.prep.library`.

    Contract §3.2: the persisted top-level object is
    ``{"ligands": [...], "n_raw": …, "n_ready": …, "n_aglycone_generated": …}``.
    A bare ``ligands`` list is accepted as well for forward compatibility.
    """
    path = lib_dir / "library_manifest.json"
    if not path.is_file():
        return {}, {}
    try:
        data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError) as exc:
        log.warning("Could not read library manifest %s: %s", path, exc)
        return {}, {}

    meta: dict[str, Any] = data if isinstance(data, dict) else {}
    entries = data.get("ligands") if isinstance(data, dict) else data
    if not isinstance(entries, list):
        return {}, meta

    ligands: dict[str, PreparedLigand] = {}
    for item in entries:
        if not isinstance(item, dict):
            continue
        try:
            ligand = PreparedLigand.model_validate(item)
        except Exception:
            continue
        if not ligand.id:
            continue
        pdbqt_path = Path(ligand.pdbqt_path) if ligand.pdbqt_path else None
        if pdbqt_path is not None and not pdbqt_path.is_absolute():
            pdbqt_path = lib_dir / pdbqt_path
        if pdbqt_path is None or not pdbqt_path.is_file():
            for candidate in (
                lib_dir / ligand.id / "ligand.pdbqt",
                lib_dir / f"{ligand.id}.pdbqt",
                lib_dir.parent / ligand.id / "ligand.pdbqt",
            ):
                if candidate.is_file():
                    pdbqt_path = candidate
                    break
        if pdbqt_path is None or not pdbqt_path.is_file():
            continue
        ligands[ligand.id] = ligand.model_copy(update={"pdbqt_path": str(pdbqt_path)})
    return ligands, meta


def evaluate_anchor_match(
    receptor_pdb_or_pdbqt: Path,
    pose_mol: Chem.Mol | None,
    anchor_residues: list[str],
    cutoff_dist: float = 4.0,
) -> float:
    """Evaluates the fraction of user-designated pocket anchor residues contacted by pose."""
    if not anchor_residues:
        return 1.0
    if pose_mol is None or pose_mol.GetNumConformers() == 0 or not receptor_pdb_or_pdbqt.is_file():
        return 0.0

    try:
        rec_text = receptor_pdb_or_pdbqt.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 0.0

    rec_atoms: list[tuple[str, str, int, tuple[float, float, float]]] = []
    for line in rec_text.splitlines():
        if (line.startswith("ATOM  ") or line.startswith("HETATM")) and len(line) >= 54:
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip().upper() or "_"
            try:
                resnr = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                rec_atoms.append((resname, chain, resnr, (x, y, z)))
            except ValueError:
                continue

    if not rec_atoms:
        return 0.0

    conf = pose_mol.GetConformer()
    lig_coords: list[tuple[float, float, float]] = []
    for i in range(pose_mol.GetNumAtoms()):
        atom = pose_mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() > 1:
            p = conf.GetAtomPosition(i)
            lig_coords.append((p.x, p.y, p.z))

    if not lig_coords:
        return 0.0

    cutoff_sq = cutoff_dist * cutoff_dist
    matched_anchors = 0

    for anchor in anchor_residues:
        norm_anchor = anchor.strip().upper()
        matching_coords: list[tuple[float, float, float]] = []
        for r_name, r_chain, r_num, coords in rec_atoms:
            cand_keys = [
                f"{r_name}{r_num}:{r_chain}",
                f"{r_name}{r_num}",
                f"{r_chain}:{r_num}",
                f"{r_num}:{r_chain}",
                f"{r_name}:{r_num}",
                str(r_num),
            ]
            if norm_anchor in cand_keys:
                matching_coords.append(coords)

        is_hit = False
        for lx, ly, lz in lig_coords:
            for rx, ry, rz in matching_coords:
                dist_sq = (lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2
                if dist_sq <= cutoff_sq:
                    is_hit = True
                    break
            if is_hit:
                break

        if is_hit:
            matched_anchors += 1

    return round(matched_anchors / len(anchor_residues), 3)


# ---------------------------------------------------------------------------
# Progress Event & Checkpoint Manager
# ---------------------------------------------------------------------------

class ScreeningProgressEvent(BaseModel):
    """Event emitted periodically to notify UI or CLI of screening progress."""

    tier: int
    stage_name: str
    current: int
    total: int
    current_ligand_id: str = ""
    speed_mol_per_sec: float = 0.0
    eta_seconds: float = 0.0
    message: str = ""


class CheckpointState(dict):
    """Container holding checkpoint data by tier and flat by ligand ID."""

    def __init__(self, by_tier: dict[int, dict[str, dict[str, Any]]], by_id: dict[str, dict[str, Any]]):
        super().__init__(by_tier)
        self.by_id = by_id


class CheckpointManager:
    """Manages streaming JSONL checkpoint and zero-loss restart reconstruction."""

    def __init__(self, checkpoint_path: Path):
        self.path = checkpoint_path
        self._lock = threading.Lock()
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def record_single(
        self,
        ligand_id: str,
        affinity: float,
        sile: float = 0.0,
        le: float = 0.0,
        lle: float = 0.0,
        run_dir: str = "",
        n_modes: int = 1,
    ) -> None:
        """Records a single-pass docking result for single-tier runs."""
        row = {
            "ligand_id": ligand_id,
            "affinity": affinity,
            "sile": sile,
            "le": le,
            "lle": lle,
            "run_dir": run_dir,
            "n_modes": n_modes,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        self._write_line(row)

    def record_tier1(
        self,
        ligand_id: str,
        affinity: float,
        sile: float,
        le: float,
        lle: float,
        run_dir: str,
        n_modes: int = 1,
    ) -> None:
        row = {
            "tier": 1,
            "ligand_id": ligand_id,
            "affinity": affinity,
            "sile": sile,
            "le": le,
            "lle": lle,
            "run_dir": run_dir,
            "n_modes": n_modes,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        self._write_line(row)

    def record_tier2(
        self,
        ligand_id: str,
        affinity_mean: float,
        affinity_std: float,
        pose_rmsd_seed: float,
        is_converged: bool,
        seed_affinities: list[float],
        seed_run_dirs: list[str],
    ) -> None:
        row = {
            "tier": 2,
            "ligand_id": ligand_id,
            "affinity_mean": affinity_mean,
            "affinity_std": affinity_std,
            "pose_rmsd_seed": pose_rmsd_seed,
            "is_converged": is_converged,
            "seed_affinities": seed_affinities,
            "seed_run_dirs": seed_run_dirs,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        self._write_line(row)

    def record_tier3(self, ligand_id: str, candidate_dict: dict[str, Any]) -> None:
        row = {
            "tier": 3,
            "ligand_id": ligand_id,
            "candidate": candidate_dict,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        self._write_line(row)

    def _write_line(self, data: dict[str, Any]) -> None:
        line = json.dumps(data, ensure_ascii=False)
        with self._lock:
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(line + "\n")
                f.flush()

    def load_completed(self) -> CheckpointState:
        """Loads completed records grouped by tier and by ligand_id."""
        by_tier: dict[int, dict[str, dict[str, Any]]] = {1: {}, 2: {}, 3: {}}
        by_id: dict[str, dict[str, Any]] = {}
        if not self.path.is_file():
            return CheckpointState(by_tier, by_id)

        with open(self.path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    if not isinstance(data, dict):
                        continue
                    lig_id = str(data.get("ligand_id", ""))
                    if not lig_id:
                        continue
                    by_id[lig_id] = data
                    tier = int(data.get("tier", 1))
                    if tier in by_tier:
                        by_tier[tier][lig_id] = data
                except (json.JSONDecodeError, ValueError):
                    continue

        return CheckpointState(by_tier, by_id)


# ---------------------------------------------------------------------------
# CSV Exporters
# ---------------------------------------------------------------------------

def export_funnel_csv(report: ScreeningFunnelReport, csv_path: Path) -> Path:
    """Exports stage loss accounting CSV strictly matching /verify T5 expectations."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        for stage, count in report.to_funnel_rows():
            writer.writerow([stage, count])
    return csv_path


def export_screening_results_csv(
    candidates: list[ScreeningCandidateRow], csv_path: Path
) -> Path:
    """Exports final candidate ranking with formula-injection neutralization,
    including backward-compatible Affinity and InChIKey headers for /verify compatibility.

    The ``Affinity`` compatibility column carries the **primary-seed** single-run
    docking score (``affinity_seed0``) so that /verify T1 - which re-docks with
    ``seed = tier2_seeds[0]`` from conf.txt - reproduces it bit-for-bit; the
    multi-seed mean stays available in ``Affinity_Mean_kcal_mol``.  The strain
    corrected score and the real best-pose path are appended as new trailing
    columns (existing columns are never removed).
    """
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank",
            "Ligand_ID",
            "Ligand_Name",
            "Derivative_Type",
            "Parent_ID",
            "Affinity_Mean_kcal_mol",
            "Affinity_Std",
            "Affinity",  # Backward-compatible alias for /verify and legacy tools
            "Pose_RMSD_Seed_A",
            "SILE",
            "LE",
            "LLE",
            "Strain_Energy_kcal_mol",
            "CNN_Score",
            "CNN_Affinity",
            "Anchor_Match_Ratio",
            "PoseBusters_OK",
            "ECR_Score",
            "ECR_Rank",
            "Cluster_ID",
            "Is_Cluster_Rep",
            "Alerts",
            "InChIKey",  # For /verify T3 structural identity verification
            "SMILES",
            "Affinity_Corrected",  # dG_corr = raw + 0.5*max(0, dE_strain - 3.0)
            "Best_Pose_Path",  # Real docked out.pdbqt behind the delivered conformer
        ])
        for rank, c in enumerate(candidates, start=1):
            inchikey = getattr(c, "inchikey", "") or ""
            if not inchikey and c.smiles:
                try:
                    m = Chem.MolFromSmiles(c.smiles)
                    if m is not None:
                        inchikey = Chem.MolToInchiKey(m)
                except Exception:
                    pass
            affinity_compat = (
                c.affinity_seed0 if c.affinity_seed0 is not None else c.affinity_mean
            )
            writer.writerow([
                rank,
                _csv_safe(c.ligand_id),
                _csv_safe(c.ligand_name or ""),
                _csv_safe(c.derivative_type),
                _csv_safe(c.parent_id or ""),
                f"{c.affinity_mean:.3f}",
                f"{c.affinity_std:.3f}",
                f"{affinity_compat:.3f}",  # Affinity alias mirrors the primary seed
                f"{c.pose_rmsd_seed:.3f}" if c.pose_rmsd_seed is not None else "",
                f"{c.sile:.4f}",
                f"{c.le:.4f}",
                f"{c.lle:.4f}",
                f"{c.strain_energy:.3f}",
                f"{c.cnn_score:.4f}" if c.cnn_score is not None else "",
                f"{c.cnn_affinity:.4f}" if c.cnn_affinity is not None else "",
                f"{c.anchor_match_ratio:.2f}",
                1 if c.posebusters_ok else 0,
                f"{c.ecr_score:.5f}",
                c.ecr_rank,
                c.cluster_id,
                1 if c.is_cluster_rep else 0,
                _csv_safe(";".join(c.alerts)),
                _csv_safe(inchikey),
                _csv_safe(c.smiles or ""),
                f"{c.affinity_corrected:.3f}" if c.affinity_corrected is not None else "",
                _csv_safe(c.best_pose_path or ""),
            ])
    return csv_path


# ---------------------------------------------------------------------------
# Main ScreeningRunner
# ---------------------------------------------------------------------------

class ScreeningRunner:
    """3-Tier Hierarchical High-Throughput Virtual Screening Runner."""

    def __init__(
        self,
        config: ScreeningConfig | None = None,
        receptor_pdbqt: Path | str | None = None,
        ligands: list[PreparedLigand] | dict[str, PreparedLigand] | PreparedLibrary | list[Path] | None = None,
        out_dir: Path | str | None = None,
        engine: DockEngine | str | None = None,
        protein_pdb: Path | str | None = None,
        resume: bool = False,
        cpu_per_job: int = 2,
        batch_size: int = 500,
        callback: Callable[[ScreeningProgressEvent], None] | None = None,
        prepared_library: PreparedLibrary | None = None,
        **kwargs: Any,
    ):
        self.config = config or ScreeningConfig()
        rec_path = receptor_pdbqt or kwargs.get("receptor") or ""
        self.receptor_pdbqt = Path(rec_path)
        out_path = out_dir or kwargs.get("out_dir") or Path.cwd() / "dock_screen_out"
        self.out_dir = Path(out_path)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.protein_pdb = Path(protein_pdb) if protein_pdb else None
        self.resume = resume
        self.cpu_per_job = max(1, cpu_per_job)
        self.batch_size = max(1, batch_size)
        self.callback = callback
        self._cancel_requested = False

        # Tracking state
        self.tier1_passed_ids: list[str] = []
        self.tier2_passed_ids: list[str] = []
        # Real (cross-tier) docking failure accounting for the funnel report: engine
        # crash, missing pose file or a >=900 sentinel affinity.  A set keeps the
        # count idempotent when a ligand is retried in a later tier.
        self._dock_failed_ids: set[str] = set()
        self._library_manifest_meta: dict[str, Any] = {}

        # Resolve engine
        if engine is None or isinstance(engine, str):
            engine_name = engine or self.config.engine or "vina"
            found_engine = get_engine(engine_name)
            self.engine = found_engine if found_engine is not None else VinaEngine()
        else:
            self.engine = engine

        # Checkpoint manager
        self.checkpoint = CheckpointManager(self.out_dir / "checkpoint.jsonl")

        # Catalogue ligands and loss tracking baseline
        self.prepared_lib: PreparedLibrary | None = prepared_library
        self.ligand_dict: dict[str, PreparedLigand] = {}

        raw_ligands = ligands if ligands is not None else prepared_library
        # An empty in-memory dict means "no ligands supplied": fall back to the
        # persistent PreparedLibrary so the runner never silently screens 0 molecules.
        if isinstance(raw_ligands, dict) and not raw_ligands and prepared_library is not None:
            raw_ligands = prepared_library
        if isinstance(raw_ligands, PreparedLibrary):
            self.prepared_lib = raw_ligands
            lib_path = Path(raw_ligands.path)
            # Restore the full PreparedLigand metadata (SMILES, QC profile, global
            # minimum energy) persisted by prep.library (contract §3.2 / P1-6.3)
            # instead of degrading every heavy-atom count to a blind constant.
            manifest_ligands, manifest_meta = _read_library_manifest(lib_path)
            self._library_manifest_meta = manifest_meta
            if manifest_ligands:
                wanted = set(raw_ligands.ligand_ids) if raw_ligands.ligand_ids else None
                for lid, ligand in manifest_ligands.items():
                    if wanted is None or lid in wanted:
                        self.ligand_dict[lid] = ligand
                log.info(
                    "Recovered %d prepared ligand(s) from %s",
                    len(self.ligand_dict),
                    lib_path / "library_manifest.json",
                )
            for lid in raw_ligands.ligand_ids:
                if lid in self.ligand_dict:
                    continue
                p_file = lib_path / lid / "ligand.pdbqt"
                if not p_file.is_file() and (lib_path.parent / lid / "ligand.pdbqt").is_file():
                    p_file = lib_path.parent / lid / "ligand.pdbqt"
                dtype = "aglycone" if lid.endswith("_aglycone") else "original"
                self.ligand_dict[lid] = PreparedLigand(
                    id=lid,
                    source="file",
                    pdbqt_path=str(p_file),
                    name=lid,
                    derivative_type=dtype,
                )
        elif isinstance(raw_ligands, dict):
            self.ligand_dict = dict(raw_ligands)
        elif isinstance(raw_ligands, list):
            for item in raw_ligands:
                if isinstance(item, PreparedLigand):
                    self.ligand_dict[item.id] = item
                elif isinstance(item, Path):
                    self.ligand_dict[item.stem] = PreparedLigand(
                        id=item.stem,
                        source="file",
                        pdbqt_path=str(item),
                        name=item.stem,
                    )
                elif isinstance(item, tuple) and len(item) == 2:
                    lid, p = item
                    self.ligand_dict[str(lid)] = PreparedLigand(
                        id=str(lid),
                        source="file",
                        pdbqt_path=str(p),
                        name=str(lid),
                    )

    def cancel(self) -> None:
        """Signal runner to terminate early gracefully."""
        self._cancel_requested = True

    # -----------------------------------------------------------------------
    # Execution helpers: failure accounting, batch capability, pose QC
    # -----------------------------------------------------------------------

    def _mark_dock_failed(self, ligand_id: str) -> None:
        """Records a genuine docking failure (engine crash / no pose / 999 sentinel).

        Deliberately separate from funnel ratio cut-offs and QC rejections: those
        are counted as ``n_funnel_filtered`` so "正常淘汰" never masquerades as
        "对接程序失败" (report P1-6.1).
        """
        self._dock_failed_ids.add(str(ligand_id))

    def _engine_has_native_batch(self) -> bool:
        """True only for engines with a *native* ``dock_batch`` implementation.

        ``BaseEngine.dock_batch`` is merely a sequential fallback that also breaks
        incremental checkpointing (a crash would discard every completed ligand of
        the chunk), so CPU engines keep using the thread pool.  Uni-Dock (GPU) and
        any adapter that really implements batching are dispatched through
        ``dock_batch`` in ``self.batch_size`` chunks (P1-4).
        """
        if not getattr(self.engine.manifest, "supports_batch_ligands", False):
            return False
        if not hasattr(self.engine, "dock_batch"):
            return False
        instance_impl = getattr(self.engine, "__dict__", {}).get("dock_batch")
        if instance_impl is not None:
            return True
        type_impl = getattr(type(self.engine), "dock_batch", None)
        return type_impl is not None and type_impl is not BaseEngine.dock_batch

    def _restore_pose_for_ligand(
        self, lig: PreparedLigand, pose_path: Path | str | None
    ) -> tuple[Chem.Mol | None, Path | None]:
        """Loads a docked pose and restores bond orders from the ligand template."""
        if pose_path is None:
            return None, None
        pose_file = Path(pose_path)
        if not pose_file.is_file():
            return None, None
        pose_mol = load_pose_mol(pose_file)
        restored = restore_pose_mol_with_template(
            pose_mol,
            lig.smiles,
            reference_pose_path=lig.pdbqt_path,
        )
        return restored, pose_file

    def _pose_geometry_qc(
        self,
        ligand_id: str,
        pose_mol: Chem.Mol | None,
        work_dir: Path | str | None,
    ) -> tuple[bool, list[str]]:
        """Exports the pose to SDF and runs the real ``check_pose`` QC (P1-5).

        Returns ``(posebusters_ok, alerts)``.  Failures never abort the run: the
        reason is surfaced through the candidate ``Alerts`` column instead of the
        previous hard-coded "always pass" verdict.
        """
        if pose_mol is None or pose_mol.GetNumConformers() == 0:
            return True, []
        try:
            base_dir = Path(work_dir) if work_dir is not None else self.out_dir
        except (TypeError, ValueError):
            base_dir = self.out_dir
        pose_sdf = base_dir / f"pose_{ligand_id}.sdf"
        try:
            pose_sdf.parent.mkdir(parents=True, exist_ok=True)
            pose_sdf.write_text(Chem.MolToMolBlock(pose_mol), encoding="utf-8")
        except Exception as exc:  # noqa: BLE001 - QC must never break screening
            log.warning("Could not export pose SDF for %s: %s", ligand_id, exc)
            return True, []
        try:
            pose_chk = check_pose(
                pose_sdf,
                box=self.config.box,
                protein_pdb=self.protein_pdb,
            )
        except Exception as exc:  # noqa: BLE001 - QC must never break screening
            log.warning("check_pose failed for %s: %s", ligand_id, exc)
            return True, []
        if pose_chk.ok:
            return True, []
        alerts = [
            f"POSECHECK_{issue.rule}: {issue.message}" for issue in pose_chk.issues
        ]
        return False, alerts

    def _emit_progress(
        self,
        tier: int,
        stage_name: str,
        current: int,
        total: int,
        current_ligand_id: str,
        start_time: float,
        message: str = "",
    ) -> None:
        if self.callback is None:
            return
        elapsed = max(0.001, time.time() - start_time)
        speed = current / elapsed
        rem = max(0, total - current)
        eta = rem / max(0.001, speed)
        event = ScreeningProgressEvent(
            tier=tier,
            stage_name=stage_name,
            current=current,
            total=total,
            current_ligand_id=current_ligand_id,
            speed_mol_per_sec=round(speed, 2),
            eta_seconds=round(eta, 1),
            message=message,
        )
        try:
            self.callback(event)
        except Exception:
            pass

    def run(self, resume: bool | None = None) -> ScreeningFunnelReport:
        """Executes the screening workflow (hierarchical 3-tier or direct single-tier)."""
        if resume is not None:
            self.resume = resume
        self._cancel_requested = False

        log.info(
            "Starting screening run in %s (engine=%s, resume=%s, tiered=%s)",
            self.out_dir,
            self.engine.manifest.name,
            self.resume,
            self.config.enable_tiered_funnel,
        )

        # Baseline loss accounting from prepared library if available
        if self.prepared_lib is not None:
            n_raw = self.prepared_lib.n_raw
            n_desalt = self.prepared_lib.n_desalt_fail
            n_dedup = self.prepared_lib.n_dedup
            n_chiral = self.prepared_lib.n_chiral_fail
            n_prep = self.prepared_lib.n_prep_fail
            n_ready = self.prepared_lib.n_ready
        else:
            n_ready = len(self.ligand_dict)
            n_raw = n_ready
            n_desalt = 0
            n_dedup = 0
            n_chiral = 0
            n_prep = 0

        # Load completed checkpoint states if resuming
        checkpoint_state = self.checkpoint.load_completed() if self.resume else CheckpointState({1: {}, 2: {}, 3: {}}, {})

        # Branch: Single-tier standard docking
        if not self.config.enable_tiered_funnel:
            return self._run_single_tier(
                n_raw, n_desalt, n_dedup, n_chiral, n_prep, n_ready, checkpoint_state
            )

        # -------------------------------------------------------------------
        # TIER 1: HTVS Fast Docking & SILE Cutoff
        # -------------------------------------------------------------------
        tier1_candidates = self._run_tier1(checkpoint_state.get(1, {}))
        self.tier1_passed_ids = [r["ligand"].id for r in tier1_candidates]
        n_tier1_in = len(self.ligand_dict)
        n_tier1_out = len(tier1_candidates)

        if not tier1_candidates or self._cancel_requested:
            return self._finalize_report(
                n_raw, n_desalt, n_dedup, n_chiral, n_prep, n_ready,
                n_tier1_in, n_tier1_out, 0, 0, 0, 0, []
            )

        # -------------------------------------------------------------------
        # TIER 2: SP Balanced Multi-Seed Convergence
        # -------------------------------------------------------------------
        tier2_candidates = self._run_tier2(tier1_candidates, checkpoint_state.get(2, {}))
        self.tier2_passed_ids = [r["ligand"].id for r in tier2_candidates]
        n_tier2_in = n_tier1_out
        n_tier2_out = len(tier2_candidates)

        if not tier2_candidates or self._cancel_requested:
            return self._finalize_report(
                n_raw, n_desalt, n_dedup, n_chiral, n_prep, n_ready,
                n_tier1_in, n_tier1_out, n_tier2_in, n_tier2_out, 0, 0, []
            )

        # -------------------------------------------------------------------
        # TIER 3: XP Rescoring, ECR Consensus & Butina Clustering
        # -------------------------------------------------------------------
        final_candidates = self._run_tier3(tier2_candidates, checkpoint_state.get(3, {}))
        n_tier3_in = n_tier2_out
        n_tier3_out = len(final_candidates)

        return self._finalize_report(
            n_raw, n_desalt, n_dedup, n_chiral, n_prep, n_ready,
            n_tier1_in, n_tier1_out, n_tier2_in, n_tier2_out, n_tier3_in, n_tier3_out,
            final_candidates,
        )

    # -----------------------------------------------------------------------
    # Single-Tier Fallback
    # -----------------------------------------------------------------------

    def _run_single_tier(
        self,
        n_raw: int,
        n_desalt: int,
        n_dedup: int,
        n_chiral: int,
        n_prep: int,
        n_ready: int,
        checkpoint_state: CheckpointState,
    ) -> ScreeningFunnelReport:
        """Executes a direct single-pass docking without hierarchical funnel reduction."""
        dock_dir = self.out_dir / "docking"
        dock_dir.mkdir(parents=True, exist_ok=True)
        items = list(self.ligand_dict.values())
        total = len(items)
        t_start = time.time()

        completed_by_id = checkpoint_state.by_id
        results_by_id: dict[str, dict[str, Any]] = {}
        pending_items: list[PreparedLigand] = []

        for lig in items:
            if lig.id in completed_by_id:
                rec = completed_by_id[lig.id]
                run_dir = rec.get("run_dir", str(dock_dir / f"dock_{lig.id}"))
                pose_file = _find_pose_file_in_dir(run_dir, lig.id)
                affinity = float(rec.get("affinity", rec.get("affinity_mean", 0.0)))
                if affinity >= 900.0:
                    self._mark_dock_failed(lig.id)
                results_by_id[lig.id] = {
                    "ligand": lig,
                    "affinity": affinity,
                    "sile": float(rec.get("sile", 0.0)),
                    "le": float(rec.get("le", 0.0)),
                    "lle": float(rec.get("lle", 0.0)),
                    "run_dir": run_dir,
                    "pose_file": str(pose_file) if pose_file is not None else None,
                }
            else:
                pending_items.append(lig)

        log.info("Single-tier: %d completed from checkpoint, %d pending", len(results_by_id), len(pending_items))

        if pending_items and not self._cancel_requested:
            workers, cpu_per_job = calculate_cpu_workers(self.config.cpu, self.cpu_per_job, len(pending_items))
            box = self.config.box or Box(center=(0, 0, 0), size=(20, 20, 20))
            params = DockParams(exhaustiveness=8, num_modes=9, cpu=cpu_per_job, scoring="vina", seed=42)

            if self._engine_has_native_batch():
                # GPU engines: chunked native dock_batch() dispatch (P1-4), never a
                # per-molecule serial dock() loop on the GPU adapter.
                self._run_single_tier_batch(pending_items, dock_dir, results_by_id, total, t_start, box, params)
            else:
                def _do_single_dock(lig: PreparedLigand) -> dict[str, Any]:
                    if self._cancel_requested:
                        raise KeyboardInterrupt("Docking cancelled")
                    sub_dir = dock_dir / f"dock_{lig.id}"
                    sub_dir.mkdir(parents=True, exist_ok=True)
                    try:
                        res = self._dock_single_safe(lig, box, params, sub_dir)
                    except BaseException:
                        self._cancel_requested = True
                        raise
                    return self._collect_single_tier_result(lig, res, sub_dir)

                with ThreadPoolExecutor(max_workers=workers) as executor:
                    future_to_lig = {executor.submit(_do_single_dock, lig): lig for lig in pending_items}
                    for fut in as_completed(future_to_lig):
                        lig = future_to_lig[fut]
                        try:
                            res_dict = fut.result()
                            results_by_id[lig.id] = res_dict
                        except BaseException as exc:
                            self._cancel_requested = True
                            self._mark_dock_failed(lig.id)
                            if isinstance(exc, KeyboardInterrupt):
                                raise
                            log.warning("Single-tier worker error on %s: %s", lig.id, exc)
                        if self._cancel_requested:
                            break
                        self._emit_progress(
                            1, "Single-tier Docking", len(results_by_id), total, lig.id, t_start
                        )

        # Build candidate rows for successful docks
        candidates: list[ScreeningCandidateRow] = []
        for lig in items:
            if lig.id not in results_by_id:
                continue
            r = results_by_id[lig.id]
            if r["affinity"] >= 900.0:
                continue
            n_heavy = self._get_heavy_atom_count(lig)
            clogp = getattr(lig.qc_profile, "clogp", None) if lig.qc_profile else None
            effs = compute_ligand_efficiencies(r["affinity"], n_heavy, clogp)
            run_dir = r.get("run_dir") or str(dock_dir / f"dock_{lig.id}")

            alerts: list[str] = []
            if lig.qc_profile:
                alerts.extend(lig.qc_profile.pains_alerts)
                alerts.extend(lig.qc_profile.tcm_warnings)

            # Restore the real docked conformer (bond orders from SMILES, 3D kept)
            pose_path = r.get("pose_file")
            if pose_path is None:
                pose_path = _find_pose_file_in_dir(run_dir, lig.id)
            restored_mol, pose_file = self._restore_pose_for_ligand(lig, pose_path)

            smi = lig.smiles or ""
            if not smi and restored_mol is not None:
                try:
                    smi = Chem.MolToSmiles(Chem.RemoveHs(restored_mol))
                except Exception:
                    smi = ""

            # Real PoseBusters / builtin geometry QC instead of a hard-coded True
            posebusters_ok, qc_alerts = self._pose_geometry_qc(lig.id, restored_mol, run_dir)
            alerts.extend(qc_alerts)

            cand = ScreeningCandidateRow(
                rank=0,
                ligand_id=lig.id,
                ligand_name=lig.name or lig.id,
                derivative_type=getattr(lig, "derivative_type", "original"),
                parent_id=getattr(lig, "parent_id", None),
                smiles=smi,
                affinity_mean=round(r["affinity"], 3),
                affinity_std=0.0,
                affinity_seed0=round(r["affinity"], 3),
                pose_rmsd_seed=None,
                le=float(effs["le"] or 0.0),
                sile=float(effs["sile"] or 0.0),
                lle=float(effs["lle"] or 0.0),
                strain_energy=0.0,
                affinity_corrected=None,
                cnn_score=None,
                cnn_affinity=None,
                anchor_match_ratio=1.0,
                posebusters_ok=posebusters_ok,
                ecr_score=0.0,
                ecr_rank=0,
                cluster_id=0,
                is_cluster_rep=False,
                best_pose_path=str(pose_file) if pose_file is not None else None,
                alerts=alerts,
            )
            candidates.append(cand)

        # Apply ECR and Murcko-Butina clustering on final candidates
        if candidates:
            compute_exponential_consensus_ranking(candidates)
            mols = [Chem.MolFromSmiles(c.smiles) if c.smiles else None for c in candidates]
            assign_clusters_to_candidates(candidates, mols, cutoff=0.4)

            # Designate highest ECR ranked candidate per cluster as cluster representative
            cluster_map: dict[int, list[ScreeningCandidateRow]] = {}
            for c in candidates:
                c.is_cluster_rep = False
                cluster_map.setdefault(c.cluster_id, []).append(c)

            for cid, cands in cluster_map.items():
                best = min(
                    cands,
                    key=lambda x: (x.ecr_rank if x.ecr_rank > 0 else 999999, -x.ecr_score, x.affinity_mean),
                )
                best.is_cluster_rep = True

        n_success = len(candidates)
        return self._finalize_report(
            n_raw, n_desalt, n_dedup, n_chiral, n_prep, n_ready,
            n_ready, n_ready, n_ready, n_ready, n_ready, n_success,
            candidates,
        )

    def _collect_single_tier_result(
        self, lig: PreparedLigand, dock_res: DockResult, fallback_run_dir: Path
    ) -> dict[str, Any]:
        """Normalises one single-tier docking result (effs + checkpoint + pose path)."""
        top_aff = float(dock_res.poses[0].affinity) if dock_res.poses else 999.0
        run_dir = str(dock_res.run_dir) if dock_res.run_dir else str(fallback_run_dir)
        pose_file = _find_pose_file_in_dir(run_dir, lig.id)
        if pose_file is None:
            pose_file = _find_pose_file_in_dir(fallback_run_dir, lig.id)
        if not dock_res.poses or top_aff >= 900.0:
            self._mark_dock_failed(lig.id)
        n_heavy = self._get_heavy_atom_count(lig)
        clogp = getattr(lig.qc_profile, "clogp", None) if lig.qc_profile else None
        effs = compute_ligand_efficiencies(top_aff, n_heavy, clogp)
        self.checkpoint.record_single(
            ligand_id=lig.id,
            affinity=top_aff,
            sile=float(effs["sile"] or 0.0),
            le=float(effs["le"] or 0.0),
            lle=float(effs["lle"] or 0.0),
            run_dir=run_dir,
            n_modes=len(dock_res.poses),
        )
        return {
            "ligand": lig,
            "affinity": top_aff,
            "sile": float(effs["sile"] or 0.0),
            "le": float(effs["le"] or 0.0),
            "lle": float(effs["lle"] or 0.0),
            "run_dir": run_dir,
            "pose_file": str(pose_file) if pose_file is not None else None,
        }

    def _run_single_tier_batch(
        self,
        pending: list[PreparedLigand],
        work_dir: Path,
        results: dict[str, dict[str, Any]],
        total: int,
        t_start: float,
        box: Box,
        params: DockParams,
    ) -> None:
        """Chunked native ``dock_batch()`` dispatch for single-tier GPU screening (P1-4)."""
        for chunk_idx in range(0, len(pending), self.batch_size):
            if self._cancel_requested:
                break
            chunk = pending[chunk_idx : chunk_idx + self.batch_size]
            batch_dir = work_dir / f"batch_{chunk_idx // self.batch_size}"
            batch_dir.mkdir(parents=True, exist_ok=True)
            batch_tuples = [(lig.id, Path(lig.pdbqt_path)) for lig in chunk]
            try:
                batch_res = self.engine.dock_batch(
                    self.receptor_pdbqt, batch_tuples, box, params, batch_dir
                )
            except Exception as exc:  # noqa: BLE001 - engine failure must not abort it all
                log.warning("Single-tier dock_batch failed (%s); isolating per ligand", exc)
                batch_res = {}
            for lig in chunk:
                if self._cancel_requested:
                    break
                dock_res = batch_res.get(lig.id)
                if dock_res is None:
                    self._mark_dock_failed(lig.id)
                    continue
                results[lig.id] = self._collect_single_tier_result(
                    lig, dock_res, batch_dir / f"dock_{lig.id}"
                )
                self._emit_progress(
                    1, "Single-tier GPU Batch", len(results), total, lig.id, t_start
                )

    # -----------------------------------------------------------------------
    # Tier 1 Implementation
    # -----------------------------------------------------------------------

    def _run_tier1(self, completed: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
        """Tier 1: Fast HTVS docking pass with SILE ranking."""
        tier1_dir = self.out_dir / "tier1"
        tier1_dir.mkdir(parents=True, exist_ok=True)
        items = list(self.ligand_dict.values())
        total = len(items)
        t_start = time.time()

        tier1_results: list[dict[str, Any]] = []
        pending_items: list[PreparedLigand] = []

        # 1. Recover completed from checkpoint
        for lig in items:
            if lig.id in completed:
                rec = completed[lig.id]
                run_dir = rec.get("run_dir", str(tier1_dir / f"dock_{lig.id}"))
                affinity = float(rec.get("affinity", 0.0))
                if affinity >= 900.0 or _find_pose_file_in_dir(run_dir, lig.id) is None:
                    self._mark_dock_failed(lig.id)
                tier1_results.append({
                    "ligand": lig,
                    "affinity": affinity,
                    "sile": float(rec.get("sile", 0.0)),
                    "le": float(rec.get("le", 0.0)),
                    "lle": float(rec.get("lle", 0.0)),
                    "run_dir": run_dir,
                })
            else:
                pending_items.append(lig)

        log.info("Tier 1: %d completed from checkpoint, %d pending", len(tier1_results), len(pending_items))

        # 2. Dock pending items
        if pending_items and not self._cancel_requested:
            if self._engine_has_native_batch():
                self._run_tier1_batch(pending_items, tier1_dir, tier1_results, total, t_start)
            else:
                self._run_tier1_workers(pending_items, tier1_dir, tier1_results, total, t_start)

        # 3. Filter and sort by SILE
        valid = [r for r in tier1_results if r["affinity"] < 900.0]
        valid.sort(key=lambda x: x["sile"], reverse=True)

        if not self.config.enable_tiered_funnel:
            return valid

        # Apply Tier 1 keep ratio (default 15%)
        keep_ratio = self.config.tier1_keep_ratio
        k = max(1, min(len(valid), int(math.ceil(len(valid) * keep_ratio))))
        k = min(k, self.config.tier1_max_candidates)
        passed = valid[:k]
        log.info("Tier 1 complete: %d/%d passed cutoff (k=%d)", len(passed), total, k)
        return passed

    def _run_tier1_workers(
        self,
        pending: list[PreparedLigand],
        work_dir: Path,
        results: list[dict[str, Any]],
        total: int,
        t_start: float,
    ) -> None:
        """Executes Tier 1 single docks across multi-core worker thread pool."""
        workers, cpu_per_job = calculate_cpu_workers(self.config.cpu, self.cpu_per_job, len(pending))
        box = self.config.box or Box(center=(0, 0, 0), size=(20, 20, 20))
        params = DockParams(exhaustiveness=4, num_modes=1, cpu=cpu_per_job, scoring="vina", seed=42)

        def _do_dock(lig: PreparedLigand) -> dict[str, Any]:
            sub_dir = work_dir / f"dock_{lig.id}"
            sub_dir.mkdir(parents=True, exist_ok=True)
            res = self._dock_single_safe(lig, box, params, sub_dir)
            top_aff = res.poses[0].affinity if res.poses else 999.0
            run_dir = str(res.run_dir) if res.run_dir else str(sub_dir)
            if not res.poses or top_aff >= 900.0:
                self._mark_dock_failed(lig.id)
            n_heavy = self._get_heavy_atom_count(lig)
            clogp = getattr(lig.qc_profile, "clogp", None) if lig.qc_profile else None
            effs = compute_ligand_efficiencies(top_aff, n_heavy, clogp)
            self.checkpoint.record_tier1(
                ligand_id=lig.id,
                affinity=top_aff,
                sile=float(effs["sile"] or 0.0),
                le=float(effs["le"] or 0.0),
                lle=float(effs["lle"] or 0.0),
                run_dir=run_dir,
                n_modes=len(res.poses),
            )
            return {
                "ligand": lig,
                "affinity": top_aff,
                "sile": float(effs["sile"] or 0.0),
                "le": float(effs["le"] or 0.0),
                "lle": float(effs["lle"] or 0.0),
                "run_dir": run_dir,
            }

        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_lig = {executor.submit(_do_dock, lig): lig for lig in pending}
            for fut in as_completed(future_to_lig):
                if self._cancel_requested:
                    break
                lig = future_to_lig[fut]
                try:
                    res_dict = fut.result()
                    results.append(res_dict)
                except Exception as exc:
                    self._mark_dock_failed(lig.id)
                    log.warning("Tier 1 worker error on %s: %s", lig.id, exc)
                self._emit_progress(
                    1, "Tier 1: HTVS Fast Docking", len(results), total,
                    lig.id, t_start
                )

    def _run_tier1_batch(
        self,
        pending: list[PreparedLigand],
        work_dir: Path,
        results: list[dict[str, Any]],
        total: int,
        t_start: float,
    ) -> None:
        """Executes Tier 1 GPU batch docking via engine.dock_batch()."""
        box = self.config.box or Box(center=(0, 0, 0), size=(20, 20, 20))
        params = DockParams(exhaustiveness=4, num_modes=1, scoring="vina", seed=42)

        for chunk_idx in range(0, len(pending), self.batch_size):
            if self._cancel_requested:
                break
            chunk = pending[chunk_idx : chunk_idx + self.batch_size]
            batch_tuples = [(lig.id, Path(lig.pdbqt_path)) for lig in chunk]
            batch_work_dir = work_dir / f"batch_{chunk_idx // self.batch_size}"
            batch_work_dir.mkdir(parents=True, exist_ok=True)

            try:
                batch_res = self.engine.dock_batch(
                    self.receptor_pdbqt,
                    batch_tuples,
                    box,
                    params,
                    batch_work_dir,
                )
            except Exception as exc:
                log.warning("Tier 1 dock_batch failed: %s; falling back to worker pool", exc)
                self._run_tier1_workers(chunk, work_dir, results, total, t_start)
                continue

            for lig in chunk:
                dock_res = batch_res.get(lig.id)
                top_aff = dock_res.poses[0].affinity if dock_res and dock_res.poses else 999.0
                # Unified per-ligand run directory: DockResult.run_dir (Uni-Dock now
                # guarantees <run_dir>/out.pdbqt) with the documented layout fallback.
                run_dir = (
                    str(dock_res.run_dir)
                    if dock_res is not None and dock_res.run_dir
                    else str(batch_work_dir / f"dock_{lig.id}")
                )
                pose_file = _find_pose_file_in_dir(run_dir, lig.id)
                if pose_file is None:
                    pose_file = _find_pose_file_in_dir(batch_work_dir / f"dock_{lig.id}", lig.id)
                if not dock_res or not dock_res.poses or top_aff >= 900.0:
                    self._mark_dock_failed(lig.id)
                n_heavy = self._get_heavy_atom_count(lig)
                clogp = getattr(lig.qc_profile, "clogp", None) if lig.qc_profile else None
                effs = compute_ligand_efficiencies(top_aff, n_heavy, clogp)
                self.checkpoint.record_tier1(
                    ligand_id=lig.id,
                    affinity=top_aff,
                    sile=float(effs["sile"] or 0.0),
                    le=float(effs["le"] or 0.0),
                    lle=float(effs["lle"] or 0.0),
                    run_dir=run_dir,
                    n_modes=len(dock_res.poses) if dock_res else 0,
                )
                results.append({
                    "ligand": lig,
                    "affinity": top_aff,
                    "sile": float(effs["sile"] or 0.0),
                    "le": float(effs["le"] or 0.0),
                    "lle": float(effs["lle"] or 0.0),
                    "run_dir": run_dir,
                })
                self._emit_progress(
                    1, "Tier 1: HTVS Fast GPU Batch", len(results), total,
                    lig.id, t_start
                )

    # -----------------------------------------------------------------------
    # Tier 2 Implementation
    # -----------------------------------------------------------------------

    def _run_tier2(
        self,
        tier1_passed: list[dict[str, Any]],
        completed: dict[str, dict[str, Any]],
    ) -> list[dict[str, Any]]:
        """Tier 2: SP Balanced multi-seed docking & physical convergence filtering."""
        tier2_dir = self.out_dir / "tier2"
        tier2_dir.mkdir(parents=True, exist_ok=True)
        total = len(tier1_passed)
        t_start = time.time()
        seeds = self.config.tier2_seeds or [42, 1024, 2026]
        box = self.config.box or Box(center=(0, 0, 0), size=(20, 20, 20))

        tier2_results: list[dict[str, Any]] = []
        pending_items: list[dict[str, Any]] = []

        # 1. Recover completed from checkpoint
        for r in tier1_passed:
            lig = r["ligand"]
            if lig.id in completed:
                rec = completed[lig.id]
                seed_dirs = [
                    str(d) for d in (rec.get("seed_run_dirs") or [r.get("run_dir", "")]) if d
                ]
                seed_affs = [
                    float(a) for a in (rec.get("seed_affinities") or []) if a is not None
                ]
                best_pose: Path | None = None
                for seed_dir in seed_dirs:
                    best_pose = _find_pose_file_in_dir(seed_dir, lig.id)
                    if best_pose is not None:
                        break
                tier2_results.append({
                    "ligand": lig,
                    "affinity_mean": float(rec.get("affinity_mean", r["affinity"])),
                    "affinity_std": float(rec.get("affinity_std", 0.0)),
                    "pose_rmsd_seed": rec.get("pose_rmsd_seed"),
                    "is_converged": bool(rec.get("is_converged", True)),
                    "sile": r["sile"],
                    "run_dir": seed_dirs[0] if seed_dirs else r.get("run_dir", str(tier2_dir / f"dock_{lig.id}")),
                    "seed_run_dirs": seed_dirs,
                    "seed_affinities": seed_affs,
                    "affinity_seed0": seed_affs[0] if seed_affs else None,
                    "affinity_corrected": None,
                    "strain_energy": 0.0,
                    "best_pose_path": str(best_pose) if best_pose is not None else None,
                    "is_strain_ok": True,
                })
            else:
                pending_items.append(r)

        log.info("Tier 2: %d completed from checkpoint, %d pending", len(tier2_results), len(pending_items))

        # 2. Multi-seed docking for pending items (native GPU batching or CPU pool)
        if pending_items and not self._cancel_requested:
            if self._engine_has_native_batch():
                self._run_tier2_batch(pending_items, tier2_dir, tier2_results, total, t_start, box, seeds)
            else:
                _, cpu_per_job = calculate_cpu_workers(
                    self.config.cpu, self.cpu_per_job, len(pending_items)
                )
                self._run_tier2_workers(
                    pending_items, tier2_dir, tier2_results, total, t_start, box, seeds, cpu_per_job
                )

        # 3. Convergence & Strain QC filtering (strict: no "pad up to 3" back door)
        converged = [
            r for r in tier2_results
            if r["affinity_mean"] < 900.0 and r["is_converged"] and r.get("is_strain_ok", True)
        ]
        if not converged and tier2_results:
            # Only when *nothing* passed the strict QC: keep the single best penalised
            # candidate and flag it loudly instead of silently waiving the thresholds
            # (report P0-2.2 / Task Fix 3.2.2).
            fallback_pool = [
                r for r in tier2_results if r["affinity_mean"] < 900.0
            ] or list(tier2_results)
            fallback_pool.sort(key=lambda x: x.get("sile", 0.0), reverse=True)
            fallback = fallback_pool[0]
            fallback["qc_fallback_warning"] = "UNCONVERGED_OR_HIGH_STRAIN"
            converged = [fallback]
            log.warning(
                "Tier 2: no candidate passed convergence/strain QC; degrading to best SILE "
                "candidate %s with HIGH_STRAIN_UNCONVERGED_FALLBACK warning",
                fallback["ligand"].id,
            )

        converged.sort(key=lambda x: x["sile"], reverse=True)

        if not self.config.enable_tiered_funnel:
            return converged

        # Apply Tier 2 keep ratio (default 2% of original ready library)
        n_base = len(self.ligand_dict)
        keep_ratio = self.config.tier2_keep_ratio
        k = max(1, min(len(converged), int(math.ceil(n_base * keep_ratio))))
        k = min(k, self.config.tier2_max_candidates)
        passed = converged[:k]
        log.info("Tier 2 complete: %d/%d passed convergence & cutoff (k=%d)", len(passed), total, k)
        return passed

    # --- Tier 2 helpers -----------------------------------------------------

    def _tier2_seed_from_dock_result(
        self,
        lig: PreparedLigand,
        dock_res: DockResult | None,
        fallback_dir: Path,
    ) -> tuple[float, Chem.Mol | None, Path | None]:
        """Converts one seed's DockResult into (affinity, restored pose, pose path)."""
        if dock_res is None or not dock_res.poses:
            return 999.0, None, None
        affinity = float(dock_res.poses[0].affinity)
        run_dir = Path(dock_res.run_dir) if dock_res.run_dir else fallback_dir
        pose_file = _find_pose_file_in_dir(run_dir, lig.id)
        if pose_file is None:
            pose_file = _find_pose_file_in_dir(fallback_dir, lig.id)
        pose_mol = load_pose_mol(pose_file) if pose_file is not None else None
        # Restore bond orders/aromaticity *before* RMSD & strain evaluation so both
        # the docked energy and the global-minimum baseline share one force field
        # (Task Fix 3.2.1 / P1-3).
        restored = restore_pose_mol_with_template(pose_mol, lig.smiles)
        return affinity, restored, pose_file

    def _evaluate_tier2_ligand(
        self,
        lig: PreparedLigand,
        tier1_row: dict[str, Any],
        lig_t2_dir: Path,
        seed_dirs: list[str],
        affinities: list[float],
        poses_mols: list[Chem.Mol | None],
        pose_files: list[Path | None],
    ) -> dict[str, Any]:
        """Computes convergence, strain-corrected affinity and efficiencies for Tier 2."""
        valid_affs = [a for a in affinities if a < 900.0]
        valid_mols = [m for m in poses_mols if m is not None]

        if valid_affs:
            mean_aff, std_aff, mean_rmsd, is_conv = compute_multiseed_convergence(
                valid_mols,
                valid_affs,
                max_sigma_dG=self.config.max_sigma_dG,
                max_rmsd_seed=self.config.max_rmsd_seed,
            )
        else:
            mean_aff, std_aff, mean_rmsd, is_conv = 999.0, 0.0, 999.0, False
            self._mark_dock_failed(lig.id)

        # Conformational strain & strain-corrected affinity (Task Fix 3.2.1)
        is_strain_ok = True
        strain_e = 0.0
        corr_aff: float | None = None
        if valid_mols:
            top_mol = valid_mols[0]
            min_e = lig.min_energy_kcal or (
                lig.qc_profile.min_energy_kcal if lig.qc_profile else None
            )
            try:
                strain_e, corr_aff, is_high_strain = compute_internal_strain_energy(
                    top_mol,
                    global_min_energy=min_e,
                    affinity=mean_aff,
                    high_strain_threshold=self.config.max_strain_kcal,
                )
                # Only a pose whose atoms really are the ligand atoms may *reject*
                # a candidate on strain grounds: a chemically incompatible
                # placeholder pose would otherwise produce a bogus high-strain veto.
                if pose_coords_are_real(top_mol) and (
                    is_high_strain or strain_e > self.config.max_strain_kcal
                ):
                    is_strain_ok = False
            except Exception as exc:  # noqa: BLE001 - strain QC is best effort
                log.warning("Strain energy calculation failed for %s: %s", lig.id, exc)

        # The strain penalty must actually drive the SILE ranking (P0-2.1).
        eff_affinity = corr_aff if corr_aff is not None else mean_aff
        n_heavy = self._get_heavy_atom_count(lig)
        clogp = getattr(lig.qc_profile, "clogp", None) if lig.qc_profile else None
        effs = compute_ligand_efficiencies(eff_affinity, n_heavy, clogp)
        new_sile = float(effs["sile"] or tier1_row.get("sile", 0.0))

        best_pose_path: str | None = None
        for pose_file in pose_files:  # seed[0] is the primary seed: preferred delivery
            if pose_file is not None:
                best_pose_path = str(pose_file)
                break
        affinity_seed0 = (
            affinities[0] if affinities and affinities[0] < 900.0 else None
        )

        self.checkpoint.record_tier2(
            ligand_id=lig.id,
            affinity_mean=mean_aff,
            affinity_std=std_aff,
            pose_rmsd_seed=mean_rmsd,
            is_converged=is_conv,
            seed_affinities=affinities,
            seed_run_dirs=seed_dirs,
        )

        return {
            "ligand": lig,
            "affinity_mean": mean_aff,
            "affinity_std": std_aff,
            "pose_rmsd_seed": mean_rmsd,
            "is_converged": is_conv,
            "is_strain_ok": is_strain_ok,
            "strain_energy": strain_e,
            "affinity_corrected": corr_aff,
            "affinity_seed0": affinity_seed0,
            "seed_affinities": affinities,
            "seed_run_dirs": seed_dirs,
            "sile": new_sile,
            "le": float(effs["le"] or 0.0),
            "lle": float(effs["lle"] or 0.0),
            "run_dir": seed_dirs[0] if seed_dirs else str(lig_t2_dir),
            "best_pose_path": best_pose_path,
        }

    def _run_tier2_workers(
        self,
        pending: list[dict[str, Any]],
        tier2_dir: Path,
        results: list[dict[str, Any]],
        total: int,
        t_start: float,
        box: Box,
        seeds: list[int],
        cpu_per_job: int,
    ) -> None:
        """CPU path: per-ligand multi-seed docking across the worker thread pool."""
        workers = calculate_cpu_workers(self.config.cpu, self.cpu_per_job, len(pending))[0]

        def _do_multiseed(row: dict[str, Any]) -> dict[str, Any]:
            lig = row["ligand"]
            lig_t2_dir = tier2_dir / f"dock_{lig.id}"
            lig_t2_dir.mkdir(parents=True, exist_ok=True)
            affinities: list[float] = []
            poses_mols: list[Chem.Mol | None] = []
            seed_dirs: list[str] = []
            pose_files: list[Path | None] = []
            for seed in seeds:
                seed_dir = lig_t2_dir / f"seed_{seed}"
                seed_dir.mkdir(parents=True, exist_ok=True)
                seed_dirs.append(str(seed_dir))
                params = DockParams(
                    exhaustiveness=8, num_modes=9, cpu=cpu_per_job, scoring="vina", seed=seed
                )
                dock_res = self._dock_single_safe(lig, box, params, seed_dir)
                affinity, pose_mol, pose_file = self._tier2_seed_from_dock_result(
                    lig, dock_res, seed_dir
                )
                affinities.append(affinity)
                poses_mols.append(pose_mol)
                pose_files.append(pose_file)
            return self._evaluate_tier2_ligand(
                lig, row, lig_t2_dir, seed_dirs, affinities, poses_mols, pose_files
            )

        with ThreadPoolExecutor(max_workers=workers) as executor:
            fut_map = {executor.submit(_do_multiseed, r): r for r in pending}
            for fut in as_completed(fut_map):
                if self._cancel_requested:
                    break
                row = fut_map[fut]
                try:
                    results.append(fut.result())
                except Exception as exc:  # noqa: BLE001 - isolate one ligand's failure
                    self._mark_dock_failed(row["ligand"].id)
                    log.warning("Tier 2 error on %s: %s", row["ligand"].id, exc)
                self._emit_progress(
                    2, "Tier 2: SP Balanced Multi-Seed", len(results), total,
                    row["ligand"].id, t_start
                )

    def _run_tier2_batch(
        self,
        pending: list[dict[str, Any]],
        tier2_dir: Path,
        results: list[dict[str, Any]],
        total: int,
        t_start: float,
        box: Box,
        seeds: list[int],
    ) -> None:
        """GPU path: one chunked native ``dock_batch()`` call per random seed (P1-4)."""
        accumulated: dict[str, dict[str, Any]] = {
            row["ligand"].id: {
                "row": row,
                "lig_t2_dir": tier2_dir / f"dock_{row['ligand'].id}",
                "affinities": [],
                "poses_mols": [],
                "seed_dirs": [],
                "pose_files": [],
            }
            for row in pending
        }

        for seed in seeds:
            for chunk_idx in range(0, len(pending), self.batch_size):
                if self._cancel_requested:
                    break
                chunk = pending[chunk_idx : chunk_idx + self.batch_size]
                batch_dir = tier2_dir / f"seed_{seed}" / f"batch_{chunk_idx // self.batch_size}"
                batch_dir.mkdir(parents=True, exist_ok=True)
                batch_tuples = [
                    (row["ligand"].id, Path(row["ligand"].pdbqt_path)) for row in chunk
                ]
                params = DockParams(
                    exhaustiveness=8, num_modes=9, scoring="vina", seed=seed
                )
                try:
                    batch_res = self.engine.dock_batch(
                        self.receptor_pdbqt, batch_tuples, box, params, batch_dir
                    )
                except Exception as exc:  # noqa: BLE001 - isolate batch failure per seed
                    log.warning("Tier 2 dock_batch (seed=%s) failed: %s", seed, exc)
                    batch_res = {}
                for row in chunk:
                    lig = row["ligand"]
                    bucket = accumulated[lig.id]
                    seed_dir = batch_dir / f"dock_{lig.id}"
                    seed_dir.mkdir(parents=True, exist_ok=True)
                    bucket["seed_dirs"].append(str(seed_dir))
                    affinity, pose_mol, pose_file = self._tier2_seed_from_dock_result(
                        lig, batch_res.get(lig.id), batch_dir
                    )
                    bucket["affinities"].append(affinity)
                    bucket["poses_mols"].append(pose_mol)
                    bucket["pose_files"].append(pose_file)

        for bucket in accumulated.values():
            if self._cancel_requested:
                break
            lig = bucket["row"]["ligand"]
            try:
                results.append(
                    self._evaluate_tier2_ligand(
                        lig,
                        bucket["row"],
                        bucket["lig_t2_dir"],
                        bucket["seed_dirs"],
                        bucket["affinities"],
                        bucket["poses_mols"],
                        bucket["pose_files"],
                    )
                )
            except Exception as exc:  # noqa: BLE001 - isolate one ligand's failure
                self._mark_dock_failed(lig.id)
                log.warning("Tier 2 (batch) error on %s: %s", lig.id, exc)
            self._emit_progress(
                2, "Tier 2: SP Balanced GPU Batch", len(results), total, lig.id, t_start
            )

    # -----------------------------------------------------------------------
    # Tier 3 Implementation
    # -----------------------------------------------------------------------

    def _run_tier3(
        self,
        tier2_passed: list[dict[str, Any]],
        completed: dict[str, dict[str, Any]],
    ) -> list[ScreeningCandidateRow]:
        """Tier 3: XP Rescoring, MMFF94s strain, GNINA CNN score, ECR ranking & clustering."""
        tier3_dir = self.out_dir / "tier3"
        tier3_dir.mkdir(parents=True, exist_ok=True)
        total = len(tier2_passed)
        t_start = time.time()

        candidates: list[ScreeningCandidateRow] = []
        pending_items: list[dict[str, Any]] = []

        # 1. Recover completed from checkpoint
        for r in tier2_passed:
            lig = r["ligand"]
            if lig.id in completed:
                c_data = completed[lig.id].get("candidate", {})
                try:
                    candidates.append(ScreeningCandidateRow(**c_data))
                except Exception:
                    pending_items.append(r)
            else:
                pending_items.append(r)

        log.info("Tier 3: %d completed from checkpoint, %d pending", len(candidates), len(pending_items))

        # Check GNINA availability once
        gnina_eng: DockEngine | None = None
        if self.config.enable_gnina_rescore:
            try:
                g = get_engine("gnina")
                if g is not None:
                    ok, _ = g.is_available()
                    if ok:
                        gnina_eng = g
            except Exception:
                pass

        # 2. Evaluate rigor properties on pending items
        for r in pending_items:
            if self._cancel_requested:
                break
            lig = r["ligand"]
            run_dir = Path(r.get("run_dir") or self.out_dir / "tier3" / f"dock_{lig.id}")
            top_pose_file = (
                Path(r["best_pose_path"])
                if r.get("best_pose_path")
                else _find_pose_file_in_dir(run_dir, lig.id)
            )
            pose_mol, top_pose_file = self._restore_pose_for_ligand(lig, top_pose_file)

            # A. Strain energy (bond orders restored above: one force field for both
            #    the docked pose and the relaxed global-minimum baseline)
            min_energy = lig.min_energy_kcal
            if min_energy is None and lig.qc_profile:
                min_energy = lig.qc_profile.min_energy_kcal

            strain_e = 0.0
            corr_aff: float | None = None
            if pose_mol is not None:
                try:
                    strain_e, corr_aff, _ = compute_internal_strain_energy(
                        pose_mol,
                        global_min_energy=min_energy,
                        affinity=r["affinity_mean"],
                        high_strain_threshold=self.config.max_strain_kcal,
                    )
                except Exception as exc:
                    log.warning("Strain energy calculation failed for %s: %s", lig.id, exc)
                if not pose_coords_are_real(pose_mol):
                    log.debug(
                        "Tier 3 strain for %s computed on a template-anchored conformer "
                        "(pose was not template-restorable); reporting only, no veto",
                        lig.id,
                    )

            # B. Ligand efficiencies on the strain-corrected affinity (P0-2.1)
            eff_affinity = corr_aff if corr_aff is not None else r["affinity_mean"]
            n_heavy = self._get_heavy_atom_count(lig)
            clogp = getattr(lig.qc_profile, "clogp", None) if lig.qc_profile else None
            effs = compute_ligand_efficiencies(eff_affinity, n_heavy, clogp)

            # C. Anchor residues match
            anchor_ratio = evaluate_anchor_match(
                self.protein_pdb or self.receptor_pdbqt,
                pose_mol,
                self.config.anchor_residues,
            )

            # D. GNINA CNN score rescore: in-place refinement of the *verified* Tier 2
            #    pose (--minimize) instead of a fresh low-exhaustiveness global search.
            cnn_score: float | None = None
            cnn_affinity: float | None = None
            if gnina_eng is not None and top_pose_file is not None and top_pose_file.is_file():
                try:
                    gnina_res = gnina_eng.dock(
                        self.receptor_pdbqt,
                        top_pose_file,
                        self.config.box or Box(center=(0, 0, 0), size=(20, 20, 20)),
                        DockParams(
                            exhaustiveness=1,
                            num_modes=1,
                            minimize=True,
                            score_only=False,
                        ),
                        run_dir / "gnina_rescore",
                    )
                    if gnina_res.poses:
                        cnn_score = gnina_res.poses[0].cnn_score
                        cnn_affinity = gnina_res.poses[0].cnn_affinity
                except Exception:
                    pass

            # E. Assemble Candidate Row
            alerts: list[str] = []
            if lig.qc_profile:
                alerts.extend(lig.qc_profile.pains_alerts)
                alerts.extend(lig.qc_profile.tcm_warnings)
            if r.get("qc_fallback_warning"):
                alerts.append(
                    "HIGH_STRAIN_UNCONVERGED_FALLBACK: 无任何分子通过多 seed 收敛/张力能 "
                    "QC，已按惩罚后 SILE 降级保留该分子，请人工复核构象真实性"
                )

            smi = lig.smiles or ""
            if not smi and pose_mol is not None:
                try:
                    smi = Chem.MolToSmiles(Chem.RemoveHs(pose_mol))
                except Exception:
                    smi = ""

            # F. Real PoseBusters / builtin geometry QC (P1-5)
            posebusters_ok, qc_alerts = self._pose_geometry_qc(lig.id, pose_mol, run_dir)
            alerts.extend(qc_alerts)

            cand = ScreeningCandidateRow(
                rank=0,
                ligand_id=lig.id,
                ligand_name=lig.name or lig.id,
                derivative_type=getattr(lig, "derivative_type", "original"),
                parent_id=getattr(lig, "parent_id", None),
                smiles=smi,
                affinity_mean=round(r["affinity_mean"], 3),
                affinity_std=round(r["affinity_std"], 3),
                affinity_seed0=(
                    round(r["affinity_seed0"], 3)
                    if r.get("affinity_seed0") is not None
                    else None
                ),
                pose_rmsd_seed=round(r["pose_rmsd_seed"], 3) if r["pose_rmsd_seed"] is not None else None,
                le=float(effs["le"] or 0.0),
                sile=float(effs["sile"] or 0.0),
                lle=float(effs["lle"] or 0.0),
                strain_energy=round(strain_e, 3),
                affinity_corrected=round(corr_aff, 3) if corr_aff is not None else None,
                cnn_score=cnn_score,
                cnn_affinity=cnn_affinity,
                anchor_match_ratio=anchor_ratio,
                posebusters_ok=posebusters_ok,
                ecr_score=0.0,
                ecr_rank=0,
                cluster_id=0,
                is_cluster_rep=False,
                best_pose_path=str(top_pose_file) if top_pose_file is not None else None,
                alerts=alerts,
            )
            candidates.append(cand)
            self.checkpoint.record_tier3(lig.id, cand.model_dump())
            self._emit_progress(
                3, "Tier 3: XP Rescoring & Rigor", len(candidates), total,
                lig.id, t_start
            )

        if not candidates:
            return []

        # 3. Exponential Consensus Ranking (ECR)
        compute_exponential_consensus_ranking(candidates)

        # 4. Murcko Scaffold Butina Diversity Clustering
        mols = [Chem.MolFromSmiles(c.smiles) if c.smiles else None for c in candidates]
        assign_clusters_to_candidates(candidates, mols, cutoff=0.4)

        # Designate highest ECR ranked candidate per cluster as cluster representative
        cluster_map: dict[int, list[ScreeningCandidateRow]] = {}
        for c in candidates:
            c.is_cluster_rep = False
            cluster_map.setdefault(c.cluster_id, []).append(c)

        for cid, cands in cluster_map.items():
            best = min(
                cands,
                key=lambda x: (x.ecr_rank if x.ecr_rank > 0 else 999999, -x.ecr_score, x.affinity_mean),
            )
            best.is_cluster_rep = True

        log.info("Tier 3 complete: %d candidates ranked and clustered", len(candidates))
        return candidates

    # -----------------------------------------------------------------------
    # Execution Helpers
    # -----------------------------------------------------------------------

    def _dock_single_safe(
        self,
        lig: PreparedLigand,
        box: Box,
        params: DockParams,
        work_dir: Path,
    ) -> DockResult:
        """Safely invokes engine.dock() catching all exceptions without crashing."""
        lig_path = Path(lig.pdbqt_path)
        if not lig_path.is_file():
            return DockResult(
                engine=self.engine.manifest.name,
                engine_version="unknown",
                params=params,
                box=box,
                poses=[],
                log=f"Ligand file not found: {lig_path}",
                run_dir=str(work_dir),
            )
        try:
            return self.engine.dock(
                self.receptor_pdbqt,
                lig_path,
                box,
                params,
                work_dir,
            )
        except Exception as exc:
            return DockResult(
                engine=self.engine.manifest.name,
                engine_version="unknown",
                params=params,
                box=box,
                poses=[],
                log=f"Docking error: {exc}",
                run_dir=str(work_dir),
            )

    @staticmethod
    def _count_heavy_atoms_in_pdbqt(pdbqt_path: str | None) -> int:
        """Counts non-hydrogen atoms in a prepared ligand PDBQT (AutoDock-type aware)."""
        if not pdbqt_path:
            return 0
        path = Path(pdbqt_path)
        if not path.is_file():
            return 0
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return 0
        count = 0
        for line in text.splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            ad_type = _record_autodock_type(line).upper()
            element = _AUTODOCK_TYPE_TO_ELEMENT.get(ad_type, "") or _element_from_atom_name(
                line[12:16]
            )
            if not element:
                count += 1
                continue
            if element.upper() == "H":
                continue
            count += 1
        return count

    def _get_heavy_atom_count(self, lig: PreparedLigand) -> int:
        """Real heavy-atom count: QC profile -> SMILES -> prepared PDBQT.

        Never returns the historical blind constant ``20`` that silently corrupted
        every SILE/LE ranking for ligands without an in-memory QC profile (P1-6.3).
        """
        if lig.qc_profile and lig.qc_profile.heavy_atoms > 0:
            return lig.qc_profile.heavy_atoms
        if lig.smiles:
            try:
                mol = Chem.MolFromSmiles(lig.smiles)
            except Exception:
                mol = None
            if mol is not None and mol.GetNumAtoms() > 0:
                return mol.GetNumHeavyAtoms()
        counted = self._count_heavy_atoms_in_pdbqt(lig.pdbqt_path)
        return counted if counted > 0 else 1

    def _finalize_report(
        self,
        n_raw: int,
        n_desalt: int,
        n_dedup: int,
        n_chiral: int,
        n_prep: int,
        n_ready: int,
        n_t1_in: int,
        n_t1_out: int,
        n_t2_in: int,
        n_t2_out: int,
        n_t3_in: int,
        n_t3_out: int,
        candidates: list[ScreeningCandidateRow],
    ) -> ScreeningFunnelReport:
        """Exports funnel.csv, screening_results.csv, conf.txt, receptor.pdbqt,
        and poses/ directory, ensuring strict /verify T1~T5 compliance."""
        n_docked_success = len(candidates)
        # Funnel mass-conservation (report P1-6.1): aglycone derivatives actually
        # added to the docking pool, genuinely failed docks and funnel threshold
        # rejections are three *distinct* quantities - never a single fudge factor.
        n_aglycone_generated = max(0, len(self.ligand_dict) - n_ready)
        n_dock_failed = len(self._dock_failed_ids)
        n_docking_pool = n_ready + n_aglycone_generated
        tiered = bool(self.config.enable_tiered_funnel)
        if tiered:
            n_funnel_filtered = max(
                0, n_docking_pool - n_docked_success - n_dock_failed
            )
            n_tier1_filtered = max(0, n_t1_in - n_t1_out)
            n_tier2_filtered = max(0, n_t2_in - n_t2_out)
            # Engine failures already have their own bucket above; keep the tier
            # detail rows free of them so that (tier1 + tier2) == funnel_filtered.
            excess = (n_tier1_filtered + n_tier2_filtered) - n_funnel_filtered
            if excess > 0:
                taken = min(excess, n_tier2_filtered)
                n_tier2_filtered -= taken
                excess -= taken
                if excess > 0:
                    n_tier1_filtered = max(0, n_tier1_filtered - excess)
        else:
            # Single-tier mode performs no ratio cut-off: nothing may be reported as
            # funnel-filtered (keeps dock_failed=1 / docked_success=9 semantics).
            n_funnel_filtered = 0
            n_tier1_filtered = 0
            n_tier2_filtered = 0

        funnel_csv = self.out_dir / "funnel.csv"
        results_csv = self.out_dir / "screening_results.csv"

        # 1. Export funnel.csv & screening_results.csv
        report = ScreeningFunnelReport(
            n_raw=n_raw,
            n_aglycone_generated=n_aglycone_generated,
            n_desalt_fail=n_desalt,
            n_dedup=n_dedup,
            n_chiral_fail=n_chiral,
            n_prep_fail=n_prep,
            n_ready=n_ready,
            n_tier1_in=n_t1_in,
            n_tier1_out=n_t1_out,
            n_tier2_in=n_t2_in,
            n_tier2_out=n_t2_out,
            n_tier3_in=n_t3_in,
            n_tier3_out=n_t3_out,
            n_tier1_filtered=n_tier1_filtered,
            n_tier2_filtered=n_tier2_filtered,
            n_funnel_filtered=n_funnel_filtered,
            n_dock_failed=n_dock_failed,
            n_docked_success=n_docked_success,
            funnel_csv_path=str(funnel_csv),
            checkpoint_path=str(self.checkpoint.path),
            results_csv_path=str(results_csv),
            candidates=candidates,
        )

        export_funnel_csv(report, funnel_csv)
        export_screening_results_csv(candidates, results_csv)

        # 2. Export conf.txt for /verify T1, T2, T4, D5
        conf_path = self.out_dir / "conf.txt"
        box = self.config.box or Box(center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0))
        seed_val = self.config.tier2_seeds[0] if (self.config and self.config.tier2_seeds) else 42
        workers_val = getattr(self, "workers", None) or self.config.cpu or max(1, os.cpu_count() or 4)
        conf_content = (
            f"center_x = {box.center[0]:.3f}\n"
            f"center_y = {box.center[1]:.3f}\n"
            f"center_z = {box.center[2]:.3f}\n"
            f"size_x = {box.size[0]:.3f}\n"
            f"size_y = {box.size[1]:.3f}\n"
            f"size_z = {box.size[2]:.3f}\n"
            f"seed = {seed_val}\n"
            f"exhaustiveness = 8\n"
            f"num_modes = 9\n"
            f"energy_range = 3.0\n"
            f"cpu = {workers_val}\n"
            f"scoring = vina\n"
        )
        conf_path.write_text(conf_content, encoding="utf-8")

        # 3. Copy receptor.pdbqt into out_dir for /verify T1, D1
        rec_dest = self.out_dir / "receptor.pdbqt"
        if self.receptor_pdbqt and Path(self.receptor_pdbqt).is_file():
            try:
                if rec_dest.resolve() != Path(self.receptor_pdbqt).resolve():
                    shutil.copyfile(self.receptor_pdbqt, rec_dest)
            except (OSError, ValueError):
                pass

        # 4. Assemble poses/ directory for /verify T1, T3, D7
        poses_dir = self.out_dir / "poses"
        poses_dir.mkdir(parents=True, exist_ok=True)
        for c in candidates:
            cid = c.ligand_id
            pose_dest_pdbqt = poses_dir / f"pose_{cid}.pdbqt"
            pose_dest_sdf = poses_dir / f"pose_{cid}.sdf"

            # Locate the *real* docked out.pdbqt: candidate best_pose_path first,
            # then the documented Tier-2/Tier-1 scratch layouts, then a per-ligand
            # search (also understands Uni-Dock's "<cid>_out.pdbqt" naming).
            found_pdbqt: Path | None = None
            if c.best_pose_path:
                candidate_path = Path(c.best_pose_path)
                if candidate_path.is_file():
                    found_pdbqt = candidate_path
            if found_pdbqt is None:
                scratch_candidates: list[Path] = []
                for seed in (self.config.tier2_seeds or [42]):
                    scratch_candidates.append(
                        self.out_dir / "tier2" / f"dock_{cid}" / f"seed_{seed}" / "out.pdbqt"
                    )
                    scratch_candidates.append(
                        self.out_dir / "tier2" / f"dock_{cid}" / f"seed_{seed}" / f"{cid}_out.pdbqt"
                    )
                scratch_candidates.append(self.out_dir / "tier2" / f"dock_{cid}" / "out.pdbqt")
                scratch_candidates.append(self.out_dir / "tier1" / f"dock_{cid}" / "out.pdbqt")
                scratch_candidates.append(self.out_dir / "docking" / f"dock_{cid}" / "out.pdbqt")
                scratch_candidates.append(self.out_dir / f"dock_{cid}" / "out.pdbqt")
                for p in scratch_candidates:
                    if p.is_file():
                        found_pdbqt = p
                        break
            if found_pdbqt is None:
                try:
                    candidate_dirs = sorted(self.out_dir.glob(f"**/*{cid}*/"))
                except (OSError, ValueError):
                    candidate_dirs = []
                for run_dir in candidate_dirs:
                    if not run_dir.is_dir():
                        continue
                    found_pdbqt = _find_pose_file_in_dir(run_dir, cid)
                    if found_pdbqt is not None:
                        break
            if found_pdbqt is None:
                log.warning(
                    "No docked out.pdbqt found for %s; refusing to substitute the "
                    "un-docked input ligand as the delivered conformer (P0-1)",
                    cid,
                )

            if found_pdbqt is not None and found_pdbqt.is_file():
                try:
                    shutil.copyfile(found_pdbqt, pose_dest_pdbqt)
                except OSError:
                    pass
                c.best_pose_path = str(found_pdbqt)

            # Write pose_<cid>.sdf from the REAL docked coordinates.  Bond orders and
            # stereo come from the reference SMILES; the conformer itself is never
            # re-embedded (the historical "embed a fresh random conformer at the
            # coordinate origin" fallback is gone - P0-1.3).
            try:
                pose_mol = load_pose_mol(found_pdbqt) if found_pdbqt is not None else None
                lig_item = self.ligand_dict.get(cid)
                restored = restore_pose_mol_with_template(
                    pose_mol,
                    c.smiles or None,
                    reference_pose_path=(
                        lig_item.pdbqt_path if lig_item is not None else None
                    ),
                )
                if restored is not None and restored.GetNumConformers() > 0:
                    pose_dest_sdf.write_text(Chem.MolToMolBlock(restored), encoding="utf-8")
                else:
                    log.warning("Could not write a conformer for %s (no usable pose)", cid)
            except Exception as exc:  # noqa: BLE001 - delivery is best effort
                log.warning("Failed to export pose SDF for %s: %s", cid, exc)

        log.info(
            "Screening report finalized: %d final candidates. Conserved: %s. Output: %s",
            len(candidates), report.is_conserved, self.out_dir
        )
        return report
