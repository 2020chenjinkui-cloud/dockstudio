"""End-to-end job lifecycle: projects -> jobs API -> analysis -> report (real vina)."""

from __future__ import annotations

import os
import time

import pytest

from dockstudio.engines.vina import VinaEngine


def _engine_available() -> bool:
    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


pytestmark = pytest.mark.skipif(not _engine_available(), reason="vina binary not available")


def _prepare(client, indinavir_smiles: str, hsg_pdb_bytes: bytes) -> tuple[dict, dict]:
    resp = client.post(
        "/api/prepare/receptor",
        files={"file": ("1HSG.pdb", hsg_pdb_bytes, "application/octet-stream")},
    )
    assert resp.status_code == 200, resp.text
    receptor = resp.json()
    resp = client.post("/api/prepare/ligand", data={"smiles": indinavir_smiles})
    assert resp.status_code == 200, resp.text
    return receptor, resp.json()


def test_job_lifecycle_e2e(client, indinavir_smiles: str, hsg_pdb_bytes: bytes) -> None:
    receptor, ligand = _prepare(client, indinavir_smiles, hsg_pdb_bytes)
    mk1 = next(c for c in receptor["cocrystal_ligands"] if c["resname"] == "MK1")

    # project CRUD
    resp = client.post("/api/projects", json={"name": "1HSG 印地那韦"})
    assert resp.status_code == 201
    project = resp.json()
    resp = client.get("/api/projects")
    assert project["id"] in [p["id"] for p in resp.json()]

    # prepare listings
    assert any(r["id"] == receptor["id"] for r in client.get("/api/prepare/receptors").json())
    assert any(l["id"] == ligand["id"] for l in client.get("/api/prepare/ligands").json())

    body = {
        "engine": "vina",
        "receptor_id": receptor["id"],
        "ligand_id": ligand["id"],
        "box": {"center": list(mk1["center"]), "size": [22.5, 22.5, 22.5]},
        "params": {
            "exhaustiveness": 2,
            "num_modes": 3,
            "cpu": max(1, (os.cpu_count() or 2) - 1),
            "seed": 42,
        },
    }
    resp = client.post(f"/api/projects/{project['id']}/jobs", json=body)
    assert resp.status_code == 201, resp.text
    job_id = resp.json()["job_id"]

    deadline = time.time() + 300
    while time.time() < deadline:
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["status"] in ("completed", "failed", "canceled", "interrupted"):
            break
        time.sleep(2)
    assert job["status"] == "completed", job.get("error")
    poses = job["result"]["poses"]
    # 姿态数受配体初始构象影响（rdkit 版本会漂移），生命周期考察不做硬性数量断言；
    # 亲和力区间才是真正的质量回归线
    assert len(poses) >= 1
    assert -14.0 <= poses[0]["affinity"] <= -6.0

    # analysis
    resp = client.get(f"/api/jobs/{job_id}/analysis", params={"pose_index": 1})
    assert resp.status_code == 200, resp.text
    analysis = resp.json()
    assert analysis["pose_index"] == 1
    assert analysis["source"] in ("plip", "builtin")
    assert analysis["interactions"], "no interactions found"

    # report
    resp = client.get(f"/api/jobs/{job_id}/report")
    assert resp.status_code == 200
    assert "attachment" in resp.headers.get("content-disposition", "")
    page = resp.text
    assert page.lstrip().startswith("<!DOCTYPE html>")
    assert "3Dmol" in page and "姿态" in page and "免责声明" in page

    # project detail aggregates job summary
    detail = client.get(f"/api/projects/{project['id']}").json()
    assert len(detail["jobs"]) == 1
    summary = detail["jobs"][0]
    assert summary["id"] == job_id and summary["status"] == "completed"
    assert summary["top_affinity"] == poses[0]["affinity"]

    # cancel of finished job is a no-op returning the job
    resp = client.post(f"/api/jobs/{job_id}/cancel")
    assert resp.status_code == 200
    assert resp.json()["status"] == "completed"


def test_projects_404_and_soft_delete(client) -> None:
    assert client.get("/api/projects/nope1234").status_code == 404
    pid = client.post("/api/projects", json={"name": "trash"}).json()["id"]
    assert client.delete(f"/api/projects/{pid}").status_code == 204
    assert client.get(f"/api/projects/{pid}").status_code == 404
    assert client.delete(f"/api/projects/{pid}").status_code == 404
