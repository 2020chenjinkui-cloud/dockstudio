"""Golden redocking regression set (opt-in: -m golden or DOCKSTUDIO_GOLDEN=1).

Each target: prepare receptor (cocrystal ligand defines the box), prepare ligand
from the CCD SMILES, dock with cheap params, assert the top affinity lands in the
measured range in targets.json.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from dockstudio.engines.vina import VinaEngine
from dockstudio.models import Box, DockParams
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.prep.receptor import prepare_receptor

GOLDEN_DIR = Path(__file__).parent
FIXTURES = GOLDEN_DIR.parent / "fixtures" / "golden"

pytestmark = pytest.mark.golden


def _targets() -> list[dict]:
    data = json.loads((GOLDEN_DIR / "targets.json").read_text(encoding="utf-8"))
    return data["targets"]


def _params():
    return [pytest.param(t, id=t["pdb_id"]) for t in _targets()]


def _engine_available() -> bool:
    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


@pytest.mark.parametrize("target", _params())
def test_golden_target(target: dict, tmp_path: Path) -> None:
    pdb_id, ccd = target["pdb_id"], target["ccd"]
    pdb_path = FIXTURES / f"{pdb_id}.pdb"
    if not pdb_path.is_file():
        pytest.fail(f"missing golden fixture {pdb_path} (run the fetch script)")
    if not _engine_available():
        pytest.skip("vina binary not available")

    rec = prepare_receptor(pdb_path.read_bytes(), f"{pdb_id}.pdb", tmp_path / "rec")
    cocrystal = next(
        (c for c in rec.cocrystal_ligands if c.resname == ccd and c.n_atoms >= 5), None
    )
    if cocrystal is None:
        pytest.fail(f"no cocrystal ligand {ccd} detected in {pdb_id} (chains={rec.chains})")
    lig = prepare_ligand_from_smiles(target["smiles"], ccd, tmp_path / "lig")

    eng = VinaEngine()
    result = eng.dock(
        Path(rec.pdbqt_path),
        Path(lig.pdbqt_path),
        Box(center=cocrystal.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=2, num_modes=3, cpu=max(1, (os.cpu_count() or 2) - 1), seed=42),
        tmp_path / "run",
    )
    top = result.poses[0].affinity
    lo, hi = target["range"] or (None, None)
    assert lo is not None, f"targets.json has no measured range for {pdb_id}"
    assert lo <= top <= hi, f"{pdb_id} top affinity {top} outside [{lo}, {hi}]"
