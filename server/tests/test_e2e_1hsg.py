"""End-to-end test: prepare 1HSG + indinavir, dock with real vina, verify poses."""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from rdkit import Chem

from dockstudio.engines.vina import VinaEngine


def _engine_available() -> bool:
    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


pytestmark = pytest.mark.skipif(not _engine_available(), reason="vina binary not available")


def test_e2e_1hsg(client, settings, indinavir_smiles: str, hsg_pdb_bytes: bytes) -> None:
    resp = client.get("/api/health")
    assert resp.status_code == 200 and resp.json()["status"] == "ok"

    engines = client.get("/api/engines").json()
    vina = next(e for e in engines if e["name"] == "vina")
    assert vina["available"], vina["unavailable_reason"]
    assert {"vina", "vinardo"} <= set(vina["scoring_functions"])

    resp = client.post(
        "/api/prepare/receptor",
        files={"file": ("1HSG.pdb", hsg_pdb_bytes, "application/octet-stream")},
    )
    assert resp.status_code == 200, resp.text
    receptor = resp.json()
    mk1 = next(c for c in receptor["cocrystal_ligands"] if c["resname"] == "MK1")

    resp = client.post(
        "/api/prepare/ligand", data={"smiles": indinavir_smiles, "name": "indinavir"}
    )
    assert resp.status_code == 200, resp.text
    ligand = resp.json()

    body = {
        "engine": "vina",
        "receptor_id": receptor["id"],
        "ligand_id": ligand["id"],
        "box": {"center": list(mk1["center"]), "size": [22.5, 22.5, 22.5]},
        "params": {
            "exhaustiveness": 4,
            "num_modes": 5,
            "cpu": max(1, (os.cpu_count() or 2) - 1),
            "seed": 42,
            "scoring": "vina",
        },
    }
    resp = client.post("/api/dock", json=body)
    assert resp.status_code == 200, resp.text
    result = resp.json()

    poses = result["poses"]
    assert len(poses) >= 3, result["log"][-2000:]
    assert -14.0 <= poses[0]["affinity"] <= -6.0
    # poses sorted best-first
    assert poses == sorted(poses, key=lambda p: p["affinity"])

    run_dir = Path(result["run_dir"])
    pose1 = run_dir / "pose_1.sdf"
    assert pose1.is_file()
    assert poses[0]["file"] and Path(poses[0]["file"]).is_file()
    assert Chem.MolFromMolBlock(pose1.read_text(encoding="utf-8"), sanitize=False) is not None

    # exported pose is fetchable through the sandboxed file endpoint
    resp = client.get("/api/file", params={"path": str(pose1)})
    assert resp.status_code == 200
    assert len(resp.text) > 0
