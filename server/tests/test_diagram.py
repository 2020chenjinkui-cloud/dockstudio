"""2D interaction diagram: well-formed SVG with ligand drawing + residue labels."""

from __future__ import annotations

# module-scoped vina run reused by the svg-structure test above
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from dockstudio.analysis.diagram import generate_diagram_svg
from dockstudio.analysis.interactions import _run_builtin
from dockstudio.models import Box, DockParams, Interaction


def _engine_available() -> bool:
    from dockstudio.engines.vina import VinaEngine

    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


@pytest.fixture(scope="module")
def pose_files(tmp_path_factory, indinavir_smiles: str, hsg_pdb_bytes: bytes):
    if not _engine_available():
        pytest.skip("vina binary not available")
    from dockstudio.analysis.export import export_poses
    from dockstudio.engines.vina import VinaEngine
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    base = tmp_path_factory.mktemp("diagram")
    lig = prepare_ligand_from_smiles(indinavir_smiles, "indinavir", base / "lig")
    rec = prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", base / "rec")
    eng = VinaEngine()
    mk1 = rec.cocrystal_ligands[0]
    run = base / "run"
    eng.dock(
        Path(rec.pdbqt_path),
        Path(lig.pdbqt_path),
        Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=1, num_modes=3, cpu=max(1, (os.cpu_count() or 2) - 1), seed=42),
        run,
    )
    return rec, export_poses(run / "out.pdbqt", run), None


def test_svg_structure_with_interactions(pose_files) -> None:
    rec, poses, _ = pose_files
    protein = Path(rec.pdbqt_path).parent / "receptor_protein.pdb"
    interactions = _run_builtin(protein, poses[0])
    svg = generate_diagram_svg(poses[0], interactions, title="indinavir")
    root = ET.fromstring(svg)  # must be well-formed XML
    assert root.tag == "svg" or root.tag.endswith("}svg")
    # ligand drawing: RDKit emits path primitives
    assert "<path" in svg
    # >=1 residue pill drawn (vina+hash ordering can differ across machines)
    import re

    assert re.search(r"[A-Z]{3}\s?\d+:[A-Z_]", svg), "no residue label found in SVG"
    # legend uses interaction-type colors
    assert any(c in svg for c in ("#2563eb", "#059669", "#9333ea", "#dc2626"))


def test_svg_empty_interactions(tmp_path: Path) -> None:
    sdf = tmp_path / "pose_1.sdf"
    sdf.write_text(
        "e\n  DockStudio\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0\n"
        "    1.5000    0.0000    0.0000 O   0  0\n"
        "  1  2  1  0  0  0  0\nM  END\n$$$$\n",
        encoding="utf-8",
    )
    svg = generate_diagram_svg(sdf, [], title="empty")
    ET.fromstring(svg)
    assert "empty" in svg


def test_svg_interaction_label_color(tmp_path: Path) -> None:
    sdf = tmp_path / "pose_1.sdf"
    sdf.write_text(
        "e\n  DockStudio\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0\n"
        "    1.5000    0.0000    0.0000 O   0  0\n"
        "  1  2  1  0  0  0  0\nM  END\n$$$$\n",
        encoding="utf-8",
    )
    svg = generate_diagram_svg(sdf, [Interaction(type="hbond", residue="ASP25:A", distance=2.9)])
    root = ET.fromstring(svg)
    texts = [el.text or "" for el in root.iter() if el.tag.endswith("text")]
    assert any("ASP25:A" in t for t in texts)
    assert any("氢键" in t for t in texts)
    assert "#2563eb" in svg  # hbond blue
