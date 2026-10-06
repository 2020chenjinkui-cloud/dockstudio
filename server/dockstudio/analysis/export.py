"""Pose export: split vina multi-MODEL PDBQT into per-pose SDF files."""

from __future__ import annotations

from pathlib import Path

from openbabel import pybel
from rdkit import Chem

from ..logging import get_logger
from ..models import DockingError

log = get_logger("analysis.export")


def _split_models(pdbqt_text: str) -> list[str]:
    """Split a vina multi-MODEL PDBQT into per-MODEL chunks (MODEL/ENDMDL)."""
    chunks: list[str] = []
    current: list[str] | None = None
    for line in pdbqt_text.splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            if current is not None:
                chunks.append("\n".join(current))
                current = None
        elif current is not None:
            current.append(line)
    return chunks


def export_poses(out_pdbqt: Path, run_dir: Path) -> list[Path]:
    """Convert each MODEL in a vina output PDBQT to pose_<n>.sdf in run_dir."""
    if not out_pdbqt.is_file():
        raise DockingError(f"vina output not found for export: {out_pdbqt}")
    text = out_pdbqt.read_text(encoding="utf-8", errors="replace")
    chunks = _split_models(text)
    if not chunks:
        raise DockingError(f"no MODEL records found in {out_pdbqt}")

    paths: list[Path] = []
    for idx, chunk in enumerate(chunks, start=1):
        try:
            mol = pybel.readstring("pdbqt", chunk)
        except OSError as exc:
            raise DockingError(f"OpenBabel could not parse MODEL {idx} of {out_pdbqt}") from exc
        pose_path = run_dir / f"pose_{idx}.sdf"
        pose_path.write_text(mol.write("sdf"), encoding="utf-8", newline="\n")
        # RDKit's SDMolSupplier cannot open non-ASCII paths on Windows; verify
        # by parsing the file content instead.
        if Chem.MolFromMolBlock(pose_path.read_text(encoding="utf-8"), sanitize=False) is None:
            raise DockingError(f"exported pose unreadable by RDKit: {pose_path}")
        paths.append(pose_path)
    log.info("exported %d poses from %s", len(paths), out_pdbqt)
    return paths
