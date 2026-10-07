"""Flexible-residue receptor variants: build, validate, dock with --flex."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from dockstudio.models import PreparationError
from dockstudio.prep.receptor import prepare_flex_variant, prepare_receptor


def _engine_available() -> bool:
    from dockstudio.engines.vina import VinaEngine

    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


@pytest.fixture(scope="module")
def base_receptor(tmp_path_factory, hsg_pdb_bytes: bytes):
    base = tmp_path_factory.mktemp("flex")
    return prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", base)


def test_flex_variant_build(base_receptor, tmp_path) -> None:
    variant = prepare_flex_variant(base_receptor, ["A:25"], tmp_path)
    assert variant.parent_id == base_receptor.id
    assert variant.flex_residues == ["A:25"]
    assert variant.flex_pdbqt_path and Path(variant.flex_pdbqt_path).is_file()
    assert variant.receptor_json_path and Path(variant.receptor_json_path).is_file()
    assert Path(variant.pdbqt_path).is_file()
    flex_text = Path(variant.flex_pdbqt_path).read_text(encoding="utf-8")
    assert "BEGIN_RES" in flex_text  # flexible residue branch marker


def test_flex_validation_errors(base_receptor, tmp_path) -> None:
    with pytest.raises(PreparationError, match="不存在"):
        prepare_flex_variant(base_receptor, ["A:999"], tmp_path)
    with pytest.raises(PreparationError, match="不存在"):
        prepare_flex_variant(base_receptor, ["Z:1"], tmp_path)
    with pytest.raises(PreparationError, match="格式"):
        prepare_flex_variant(base_receptor, ["A25"], tmp_path)
    with pytest.raises(PreparationError, match="最多"):
        prepare_flex_variant(base_receptor, [f"A:{i}" for i in range(1, 12)], tmp_path)


def test_flex_api(base_receptor, client, hsg_pdb_bytes: bytes) -> None:
    resp = client.post("/api/prepare/receptor", files={"file": ("1HSG.pdb", hsg_pdb_bytes)})
    rid = resp.json()["id"]
    resp = client.post(f"/api/prepare/receptor/{rid}/flex", json={"flex_residues": ["A:25"]})
    assert resp.status_code == 201, resp.text
    variant = resp.json()
    assert variant["parent_id"] == rid
    assert variant["flex_pdbqt_path"]
    # bad residue -> 400 with chain info
    resp = client.post(f"/api/prepare/receptor/{rid}/flex", json={"flex_residues": ["A:9999"]})
    assert resp.status_code == 400
    assert "A" in resp.json()["detail"]
    # variant appears in receptor listing
    listing = client.get("/api/prepare/receptors").json()
    assert any(r["id"] == variant["id"] for r in listing)


@pytest.mark.skipif(not _engine_available(), reason="vina binary not available")
def test_flex_dock_e2e(tmp_path, base_receptor, indinavir_smiles: str) -> None:
    from dockstudio.engines.vina import VinaEngine
    from dockstudio.models import Box, DockParams
    from dockstudio.prep.ligand import prepare_ligand_from_smiles

    variant = prepare_flex_variant(base_receptor, ["A:25"], tmp_path / "var")
    lig = prepare_ligand_from_smiles(indinavir_smiles, "indinavir", tmp_path / "lig")
    eng = VinaEngine()
    mk1 = base_receptor.cocrystal_ligands[0]
    result = eng.dock(
        Path(variant.pdbqt_path),
        Path(lig.pdbqt_path),
        Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=1, num_modes=3, cpu=max(1, (os.cpu_count() or 2) - 1), seed=42),
        tmp_path / "run",
        Path(variant.flex_pdbqt_path),
    )
    assert result.poses
    assert -14.0 <= result.poses[0].affinity <= -5.0
    assert (tmp_path / "run" / "flex.pdbqt").is_file()  # staged --flex input
