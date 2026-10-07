"""Protein-ligand interaction analysis: PLIP preferred, builtin geometry fallback."""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path

from openbabel import pybel
from rdkit import Chem

from ..logging import get_logger
from ..models import Box, Interaction, PoseAnalysis

log = get_logger("analysis.interactions")


# ---------------------------------------------------------------------------
# PLIP (optional)
# ---------------------------------------------------------------------------


def _plip_available() -> bool:
    try:
        import plip  # noqa: F401
    except ImportError:
        return False
    return True


def _patch_openbabel_for_plip() -> None:
    """openbabel-wheel ships without the inchikey format; PLIP calls it
    unconditionally. Swallow that one ValueError."""
    orig = pybel.Molecule.write

    def write(self, format="smi", filename=None, overwrite=False, opt=None):
        try:
            return orig(self, format, filename, overwrite, opt)
        except ValueError as exc:
            if "inchikey" in str(exc):
                return ""
            raise

    pybel.Molecule.write = write


def _build_complex_pdb(protein_pdb: Path, pose_sdf: Path, out_path: Path) -> None:
    """Concatenate protein PDB + pose (as LIG HETATM) into one PLIP-readable PDB."""
    mol = Chem.MolFromMolBlock(pose_sdf.read_text(encoding="utf-8"), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not read pose file {pose_sdf}")
    block = Chem.MolToPDBBlock(mol)
    lig_lines = [
        line[:17] + "LIG" + line[20:22] + "Z" + line[23:] if line.startswith("HETATM") else line
        for line in block.splitlines()
        if line.strip()
    ]
    protein_lines = [
        line
        for line in protein_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() not in ("END", "")
    ]
    out_path.write_text("\n".join(protein_lines + lig_lines + ["END", ""]), encoding="utf-8")


def _run_plip(protein_pdb: Path, pose_sdf: Path, work_dir: Path) -> list[Interaction]:
    _patch_openbabel_for_plip()
    import plip  # noqa: F401
    from plip.structure.preparation import PDBComplex

    logging.disable(logging.CRITICAL)
    try:
        complex_pdb = work_dir / "plip_complex.pdb"
        _build_complex_pdb(protein_pdb, pose_sdf, complex_pdb)
        complex_ = PDBComplex()
        complex_.load_pdb(str(complex_pdb))
        complex_.analyze()
    finally:
        logging.disable(logging.NOTSET)

    interactions: list[Interaction] = []
    for bs in complex_.interaction_sets.values():
        hbonds = list(bs.hbonds_ldon) + list(bs.hbonds_pdon)
        for hb in hbonds:
            interactions.append(
                Interaction(
                    type="hbond",
                    residue=f"{hb.restype}{hb.resnr}:{hb.reschain}",
                    distance=float(hb.distance_ad),
                    detail="蛋白供体" if hb.protisdon else "配体供体",
                )
            )
        for hc in bs.hydrophobic_contacts:
            interactions.append(
                Interaction(
                    type="hydrophobic",
                    residue=f"{hc.restype}{hc.resnr}:{hc.reschain}",
                    distance=float(hc.distance),
                )
            )
        for pi in bs.pistacking:
            interactions.append(
                Interaction(
                    type="pi_stacking",
                    residue=f"{pi.restype}{pi.resnr}:{pi.reschain}",
                    distance=float(pi.distance),
                )
            )
        for sb in list(bs.saltbridge_lneg) + list(bs.saltbridge_pneg):
            interactions.append(
                Interaction(
                    type="salt_bridge",
                    residue=f"{sb.restype}{sb.resnr}:{sb.reschain}",
                    distance=float(sb.distance),
                )
            )
        for xb in bs.halogen_bonds:
            interactions.append(
                Interaction(
                    type="halogen",
                    residue=f"{xb.restype}{xb.resnr}:{xb.reschain}",
                    distance=float(xb.distance),
                )
            )
    return interactions


# ---------------------------------------------------------------------------
# Builtin geometric detection
# ---------------------------------------------------------------------------

_DONOR_ACCEPTOR = {"N", "O", "S", "F"}
_POS_RES = {"ARG", "LYS", "HIP"}
_NEG_RES = {"ASP", "GLU"}
_AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS", "HIE", "HID", "HIP"}
_HYDROPHOBIC_RES = {
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "TYR",
    "HIS",
    "CYS",
}


class _Atom:
    __slots__ = ("chain", "element", "name", "pos", "resname", "resnr")

    def __init__(self, line: str):
        self.element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not self.element:
            self.element = line[12:16].strip().lstrip("0123456789")[:1].upper()
        self.pos = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.resname = line[17:20].strip()
        self.resnr = line[22:26].strip()
        self.chain = line[21:22].strip()
        self.name = line[12:16].strip()


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def _res_label(a: _Atom) -> str:
    return f"{a.resname}{a.resnr}:{a.chain or '_'}"


def _parse_protein_atoms(pdb_path: Path) -> list[_Atom]:
    atoms: list[_Atom] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM  ") and len(line) >= 54:
            atoms.append(_Atom(line))
    return atoms


def _ring_centroid(mof: Chem.Mol, ring: tuple[int, ...]) -> tuple[float, float, float]:
    conf = mof.GetConformer()
    xs = ys = zs = 0.0
    for i in ring:
        p = conf.GetAtomPosition(i)
        xs += p.x
        ys += p.y
        zs += p.z
    n = len(ring)
    return xs / n, ys / n, zs / n


def _run_builtin(protein_pdb: Path, pose_sdf: Path) -> list[Interaction]:
    pose = Chem.MolFromMolBlock(
        pose_sdf.read_text(encoding="utf-8"), sanitize=False, removeHs=False
    )
    if pose is None:
        raise ValueError(f"RDKit could not read pose {pose_sdf}")
    protein = _parse_protein_atoms(protein_pdb)
    conf = pose.GetConformer()

    lig_pos = [conf.GetAtomPosition(i) for i in range(pose.GetNumAtoms())]
    lig_elems = [pose.GetAtomWithIdx(i).GetSymbol().upper() for i in range(pose.GetNumAtoms())]
    ring_info = pose.GetRingInfo()
    aromatic_rings = [
        ring
        for ring in ring_info.AtomRings()
        if all(pose.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
    ]

    interactions: list[Interaction] = []
    seen_hb: set[str] = set()
    seen_hydrophobic: set[str] = set()
    seen_pi: set[str] = set()
    seen_sb: set[str] = set()

    for pa in protein:
        if pa.element == "H":
            continue
        for idx, (lp, le) in enumerate(zip(lig_pos, lig_elems)):
            if le == "H":
                continue
            d = _dist(pa.pos, (lp.x, lp.y, lp.z))
            label = _res_label(pa)
            # hydrogen bonds: both atoms must be donor/acceptor types
            if d < 3.5 and pa.element in _DONOR_ACCEPTOR and le in _DONOR_ACCEPTOR:
                if label not in seen_hb:
                    seen_hb.add(label)
                    interactions.append(
                        Interaction(
                            type="hbond",
                            residue=label,
                            distance=round(d, 2),
                            ligand_atom_ids=[idx],
                        )
                    )
            # hydrophobic: carbon-carbon contact with hydrophobic residue
            elif (
                d < 4.0
                and pa.element == "C"
                and le == "C"
                and pa.resname in _HYDROPHOBIC_RES
                and label not in seen_hydrophobic
            ):
                seen_hydrophobic.add(label)
                interactions.append(
                    Interaction(
                        type="hydrophobic",
                        residue=label,
                        distance=round(d, 2),
                        ligand_atom_ids=[idx],
                    )
                )
            # salt bridges: charged groups within 4 Å
            if (
                d < 4.0
                and pa.element in {"N", "O"}
                and le in {"N", "O"}
                and pa.resname in _POS_RES | _NEG_RES
                and label not in seen_sb
            ):
                lig_sym = pose.GetAtomWithIdx(idx).GetSymbol()
                res_charged = (pa.resname in _POS_RES and pa.element == "N") or (
                    pa.resname in _NEG_RES and pa.element == "O"
                )
                if res_charged and lig_sym in {"N", "O"}:
                    seen_sb.add(label)
                    interactions.append(
                        Interaction(
                            type="salt_bridge",
                            residue=label,
                            distance=round(d, 2),
                            ligand_atom_ids=[idx],
                        )
                    )

    # pi stacking: aromatic ring centroid distance < 5.5 Å against aromatic residues
    res_rings: dict[str, list[tuple[float, float, float]]] = {}
    if aromatic_rings:
        # group protein aromatic side-chain carbons by residue
        by_res: dict[str, list[tuple[float, float, float]]] = {}
        for pa in protein:
            if pa.resname in _AROMATIC_RES and pa.element == "C":
                by_res.setdefault(_res_label(pa), []).append(pa.pos)
        for label, coords in by_res.items():
            if len(coords) >= 5:
                n = len(coords)
                cx = sum(c[0] for c in coords) / n
                cy = sum(c[1] for c in coords) / n
                cz = sum(c[2] for c in coords) / n
                res_rings.setdefault(label, []).append((cx, cy, cz))
        for ring in aromatic_rings:
            c = _ring_centroid(pose, ring)
            for label, centroids in res_rings.items():
                for rc in centroids:
                    d = _dist(c, rc)
                    if d < 5.5 and label not in seen_pi:
                        seen_pi.add(label)
                        interactions.append(
                            Interaction(
                                type="pi_stacking",
                                residue=label,
                                distance=round(d, 2),
                                ligand_atom_ids=list(ring),
                            )
                        )

    interactions.sort(key=lambda i: (i.type, i.distance or 99.0))
    return interactions


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyze_pose(
    protein_pdb: Path,
    pose_sdf: Path,
    work_dir: Path,
    use_plip: bool = True,
    box: Box | None = None,
) -> PoseAnalysis:
    """Analyze one pose against the protein. Tries PLIP, falls back to builtin."""
    from .posecheck import check_pose

    pose_idx = int(pose_sdf.stem.rsplit("_", 1)[-1]) if "_" in pose_sdf.stem else 1
    warnings: list[str] = []
    pose_check = None
    try:
        pose_check = check_pose(pose_sdf, box=box, pose_index=pose_idx)
    except Exception as exc:  # noqa: BLE001
        log.warning("pose check failed (non-fatal): %s", exc)
        warnings.append(f"姿态质检失败（不影响分析）: {exc}")
    if use_plip and _plip_available():
        try:
            interactions = _run_plip(protein_pdb, pose_sdf, work_dir)
            return PoseAnalysis(
                pose_index=pose_idx,
                source="plip",
                interactions=interactions,
                pose_check=pose_check,
                warnings=warnings,
            )
        except Exception as exc:  # noqa: BLE001
            log.warning("PLIP failed, falling back to builtin: %s", exc)
            warnings.append(f"PLIP 失败，使用内置几何检测: {exc}")
    interactions = _run_builtin(protein_pdb, pose_sdf)
    return PoseAnalysis(
        pose_index=pose_idx,
        source="builtin",
        interactions=interactions,
        warnings=warnings,
        pose_check=pose_check,
    )


def analyze_job_poses(
    protein_pdb: Path,
    pose_paths: list[Path],
    job_dir: Path,
    use_plip: bool = True,
    pose_indices: list[int] | None = None,
    box: Box | None = None,
) -> list[PoseAnalysis]:
    """Analyze (and cache under job_dir/analysis/) the given poses."""
    cache_dir = job_dir / "analysis"
    cache_dir.mkdir(parents=True, exist_ok=True)
    analyses: list[PoseAnalysis] = []
    indices = pose_indices or list(range(1, len(pose_paths) + 1))
    for idx, pose_path in zip(indices, pose_paths):
        cache_path = cache_dir / f"pose_{idx}.json"
        if cache_path.is_file():
            try:
                analyses.append(PoseAnalysis(**json.loads(cache_path.read_text(encoding="utf-8"))))
                continue
            except (json.JSONDecodeError, ValueError):
                pass
        analysis = analyze_pose(protein_pdb, pose_path, job_dir, use_plip=use_plip, box=box)
        cache_path.write_text(analysis.model_dump_json(indent=2), encoding="utf-8")
        analyses.append(analysis)
    return analyses
