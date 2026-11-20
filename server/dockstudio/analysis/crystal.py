"""Co-crystal reference comparison: heavy-atom RMSD between crystal and docked pose.

Two different questions can be asked of a docked pose, and they give different
numbers for the same geometry:

* **in-place / 原位** — "is the pose still where the co-crystal ligand sits?"
  Computed with :func:`rdkit.Chem.rdMolAlign.CalcRMS`, which does *not* move the
  probe.  This is the question the D9 acceptance criterion asks
  (「共晶配体在原位点重新对接的 RMSD」), so it is the number the platform
  reports as the verdict (:data:`CRITERION_METRIC`).
* **superposed / 叠合** — "does the pose reproduce the crystal *conformation*,
  disregarding where it was placed?"  Computed with
  :func:`rdkit.Chem.rdMolAlign.GetBestRMS`, which superimposes first.  Useful as
  a physics/degeneracy diagnostic (:data:`SUPERPOSED_METRIC`) but it is
  mathematically blind to rigid-body drift: a perfectly copied pose translated
  10 Å away still scores 0.0 Å.

RDKit's own documentation makes the same distinction and recommends the
in-place form for exactly this use case::

    CalcRMS: "In contrast to getBestRMS, the RMS is computed 'in place', i.e.
    probe molecules are not aligned to the reference ahead of the RMS
    calculation. This is useful, for example, to compute the RMSD between
    docking poses and the co-crystallized ligand."

Both metrics are symmetry-corrected (RDKit enumerates the reference graph's
automorphisms and takes the minimum), so comparing them isolates the
superposition, not the atom matching.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from ..logging import get_logger

log = get_logger("analysis.crystal")

# --- metric identifiers: stable strings, safe to put in JSON/report payloads ---
INPLACE_METRIC = "inplace_calc_rms"
SUPERPOSED_METRIC = "superposed_best_rms"
#: The metric that answers D9.  Named so callers never have to guess.
CRITERION_METRIC = INPLACE_METRIC
#: The auxiliary metric, kept for diagnostics only.
AUXILIARY_METRIC = SUPERPOSED_METRIC

#: D9 acceptance threshold (Å): co-crystal ligand re-docked *at its original site*.
D9_THRESHOLD_A = 2.0

METRIC_LABELS_ZH = {
    INPLACE_METRIC: "原位 RMSD（不叠合，D9 判据）",
    SUPERPOSED_METRIC: "叠合 RMSD（先最优叠合，仅辅助）",
}


def _read_sdf(path: Path) -> Chem.Mol:
    mol = Chem.MolFromMolBlock(path.read_text(encoding="utf-8"), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not read SDF {path}")
    return mol


def _strip_bond_types(mol: Chem.Mol) -> Chem.Mol:
    """Copy the mol with all bonds as SINGLE — RMSD is about geometry/topology,
    not bond order; mixed typing would break substructure matching."""
    mol = Chem.Mol(mol)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    return mol


def _load_pair(crystal_sdf: Path, pose_sdf: Path) -> tuple[Chem.Mol, Chem.Mol]:
    """Read both SDFs and reduce them to comparable heavy-atom graphs.

    ``RemoveHs`` is a no-op while the molecule is unsanitised, and a co-crystal SDF
    produced by :func:`dockstudio.prep.receptor.restore_cocrystal_bond_orders` has
    already had its explicit hydrogens stripped, so a raw ``RemoveHs`` here would
    leave the two sides with different atom counts (1PPC: 37 vs 38) and the
    substructure match would fail with "No sub-structure match found".  Sanitise a
    copy first so both sides really are heavy-atom graphs.
    """
    crystal = _read_sdf(crystal_sdf)
    pose = _read_sdf(pose_sdf)
    ref = _strip_bond_types(_heavy_atom_graph(crystal))
    probe = _strip_bond_types(_heavy_atom_graph(pose))
    return ref, probe


def _heavy_atom_graph(mol: Chem.Mol) -> Chem.Mol:
    """Drop hydrogens deterministically, regardless of the file's sanitisation state."""
    if not any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        return Chem.RemoveHs(mol)
    work = Chem.RWMol(mol)
    for idx in sorted((a.GetIdx() for a in work.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        work.RemoveAtom(idx)
    return work.GetMol()


@dataclass(frozen=True)
class CrystalPoseRmsd:
    """Both RMSD flavours for one pose, with the metric names attached.

    Attributes
    ----------
    inplace_rmsd:
        Symmetry-corrected RMSD computed *without* moving the probe
        (``rdMolAlign.CalcRMS``).  This is the D9 verdict.
    superposed_rmsd:
        Symmetry-corrected RMSD after optimal rigid-body superposition
        (``rdMolAlign.GetBestRMS``).  Auxiliary diagnostic only.
    metric / auxiliary_metric:
        Which function produced each number.
    rmsd:
        Alias of :attr:`inplace_rmsd` — the single number a D9 consumer wants.
    """

    inplace_rmsd: float | None
    superposed_rmsd: float | None
    metric: str = CRITERION_METRIC
    auxiliary_metric: str = AUXILIARY_METRIC
    threshold_a: float = D9_THRESHOLD_A

    @property
    def rmsd(self) -> float | None:
        """D9 number (in-place). Never the superposed one."""
        return self.inplace_rmsd

    @property
    def d9_pass(self) -> bool | None:
        """``None`` when the metric could not be computed."""
        if self.inplace_rmsd is None:
            return None
        return bool(self.inplace_rmsd <= self.threshold_a)

    def d9_status(self) -> str:
        """``PASS`` / ``FAIL`` / ``SKIP`` for report and audit rendering."""
        verdict = self.d9_pass
        if verdict is None:
            return "SKIP"
        return "PASS" if verdict else "FAIL"

    def as_dict(self) -> dict:
        """JSON-friendly payload; field names carry the metric, not just the value."""
        return {
            "rmsd": self.rmsd,
            "metric": self.metric,
            "auxiliary_rmsd": self.superposed_rmsd,
            "auxiliary_metric": self.auxiliary_metric,
            "threshold_a": self.threshold_a,
            "d9_status": self.d9_status(),
        }


def crystal_pose_rmsd(crystal_sdf: Path, pose_sdf: Path) -> CrystalPoseRmsd:
    """Compute both heavy-atom RMSD flavours for one pose against the co-crystal.

    A failure in one flavour does not discard the other: each is caught
    separately so the caller still receives whatever could be measured, and the
    reason is logged (科学降级必须留痕).
    """
    ref, probe = _load_pair(crystal_sdf, pose_sdf)

    inplace: float | None = None
    try:
        inplace = round(float(rdMolAlign.CalcRMS(Chem.Mol(probe), Chem.Mol(ref))), 3)
    except Exception as exc:  # 降级兜底
        log.warning("in-place RMSD (CalcRMS) failed: %s", exc)

    superposed: float | None = None
    try:
        superposed = round(float(rdMolAlign.GetBestRMS(Chem.Mol(probe), Chem.Mol(ref))), 3)
    except Exception as exc:  # 降级兜底
        log.warning("superposed RMSD (GetBestRMS) failed: %s", exc)

    return CrystalPoseRmsd(inplace_rmsd=inplace, superposed_rmsd=superposed)


def heavy_atom_rmsd(crystal_sdf: Path, pose_sdf: Path) -> float | None:
    """Heavy-atom RMSD **after best-fit superimposition** (``GetBestRMS``).

    Deprecated as a D9 verdict: superimposing first makes the number blind to
    rigid-body drift — a crystal-identical pose translated 10 Å still reports
    ~0.0 Å.  It is kept because callers and stored results predate the split of
    the two metrics; new code must use :func:`crystal_pose_rmsd` and read
    :attr:`CrystalPoseRmsd.inplace_rmsd` (the D9 number) from it.
    """
    return crystal_pose_rmsd(crystal_sdf, pose_sdf).superposed_rmsd


def inplace_heavy_atom_rmsd(crystal_sdf: Path, pose_sdf: Path) -> float | None:
    """Heavy-atom RMSD **in place** (``rdMolAlign.CalcRMS``) — the D9 metric.

    The probe is not moved, so the number answers "is the pose at the crystal
    ligand's site?" rather than "does it have the crystal ligand's shape?".
    """
    return crystal_pose_rmsd(crystal_sdf, pose_sdf).inplace_rmsd


def superposed_heavy_atom_rmsd(crystal_sdf: Path, pose_sdf: Path) -> float | None:
    """Explicit alias of :func:`heavy_atom_rmsd` (superimposed; auxiliary only)."""
    return crystal_pose_rmsd(crystal_sdf, pose_sdf).superposed_rmsd


def cocrystal_reference(receptor_dir: Path) -> tuple[str | None, Path | None]:
    """Locate the first usable co-crystal ligand for a prepared receptor.

    The receptor preparation writes a sidecar ``cocrystals.json`` next to the
    receptor PDBQT (see :mod:`dockstudio.prep.receptor`); this reader is the
    single place where the API, the report and the audit bundle agree on where
    the D9 reference lives.
    """
    receptor_dir = Path(receptor_dir)
    meta_path = receptor_dir / "cocrystals.json"
    if not meta_path.is_file():
        return None, None
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:  # 降级兜底
        log.warning("could not read co-crystal sidecar %s: %s", meta_path, exc)
        return None, None
    if not meta:
        return None, None
    first = meta[0]
    reference = f"{first.get('resname', '?')}:{first.get('chain', '?')}"
    sdf_path = Path(str(first.get("sdf_path", "")))
    return reference, sdf_path if sdf_path.is_file() else None


def pose_rmsd_series(
    crystal_sdf: Path, pose_paths: list[Path], max_poses: int | None = None
) -> dict[int, CrystalPoseRmsd]:
    """Compute both RMSD flavours for the first ``max_poses`` docked poses.

    Keys are the 1-based pose indices that ``DockResult.poses[i].index`` uses, so
    callers can render the D9 verdict per pose without re-deriving the numbering.
    """
    limit = len(pose_paths) if max_poses is None else min(max_poses, len(pose_paths))
    out: dict[int, CrystalPoseRmsd] = {}
    for idx in range(1, limit + 1):
        try:
            out[idx] = crystal_pose_rmsd(crystal_sdf, pose_paths[idx - 1])
        except Exception as exc:  # 降级兜底
            log.warning("pose %d: co-crystal RMSD unavailable: %s", idx, exc)
    return out


def best_pose_rmsd(series: dict[int, CrystalPoseRmsd]) -> CrystalPoseRmsd | None:
    """Best-of-N: the pose whose **D9** metric is smallest.

    Ranking by the superposed metric would be wrong here — that number ignores
    placement, so "best" by it can be a pose sitting outside the pocket.
    """
    scored = [r for r in series.values() if r.inplace_rmsd is not None]
    return min(scored, key=lambda r: r.inplace_rmsd) if scored else None
