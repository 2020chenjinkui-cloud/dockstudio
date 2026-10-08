"""PLIP worker: runs in an ISOLATED process to escape rdkit<->openbabel ABI conflicts.

IMPORTANT: this file is executed directly by path (NOT as a package module),
so `dockstudio/__init__` (which preloads rdkit) must NOT run here — do not add
any `dockstudio.` imports. Only openbabel/plip/std imports allowed.

Usage: python plip_worker.py <protein.pdb> <pose.sdf> <work_dir>
Stdout: JSON array of interaction dicts {type, residue, distance, detail}.
Exit code: 0 on success (even with zero interactions); !=0 on failure.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path


def _patch_openbabel_for_plip() -> None:
    from openbabel import pybel

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
    """SDF pose -> PDB via pybel（worker 里不碰 rdkit）。string I/O 兜中文路径。"""
    from openbabel import pybel

    sd_text = pose_sdf.read_text(encoding="utf-8", errors="replace")
    mol = pybel.readstring("sdf", sd_text)
    if mol is None or mol.OBMol.NumAtoms() == 0:
        raise ValueError(f"OpenBabel could not read pose {pose_sdf}")
    lig_lines = []
    for line in mol.write("pdb").splitlines():
        if line.startswith("HETATM"):
            lig_lines.append(line[:17] + "LIG" + line[20:22] + "Z" + line[23:])
    protein_lines = [
        line
        for line in protein_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() not in ("END", "")
    ]
    out_path.write_text("\n".join(protein_lines + lig_lines + ["END", ""]), encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: plip_worker.py <protein.pdb> <pose.sdf> <work_dir>", file=sys.stderr)
        return 2
    protein_pdb, pose_sdf, work_dir = (Path(a) for a in sys.argv[1:4])
    work_dir.mkdir(parents=True, exist_ok=True)

    _patch_openbabel_for_plip()
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

    interactions: list[dict] = []
    for bs in complex_.interaction_sets.values():
        for hb in list(bs.hbonds_ldon) + list(bs.hbonds_pdon):
            interactions.append(
                {
                    "type": "hbond",
                    "residue": f"{hb.restype}{hb.resnr}:{hb.reschain}",
                    "distance": float(hb.distance_ad),
                    "detail": "蛋白供体" if hb.protisdon else "配体供体",
                }
            )
        for hc in bs.hydrophobic_contacts:
            interactions.append(
                {
                    "type": "hydrophobic",
                    "residue": f"{hc.restype}{hc.resnr}:{hc.reschain}",
                    "distance": float(hc.distance),
                }
            )
        for pi in bs.pistacking:
            interactions.append(
                {
                    "type": "pi_stacking",
                    "residue": f"{pi.restype}{pi.resnr}:{pi.reschain}",
                    "distance": float(pi.distance),
                }
            )
        for sb in list(bs.saltbridge_lneg) + list(bs.saltbridge_pneg):
            interactions.append(
                {
                    "type": "salt_bridge",
                    "residue": f"{sb.restype}{sb.resnr}:{sb.reschain}",
                    "distance": float(sb.distance),
                }
            )
        for xb in bs.halogen_bonds:
            interactions.append(
                {
                    "type": "halogen",
                    "residue": f"{xb.restype}{xb.resnr}:{xb.reschain}",
                    "distance": float(xb.distance),
                }
            )
    json.dump(interactions, sys.stdout, ensure_ascii=False)
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
