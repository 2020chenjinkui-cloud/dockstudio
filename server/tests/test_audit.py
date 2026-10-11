"""audit.zip + audit markdown + demo endpoints + CLI audit."""

from __future__ import annotations

import json
import time
import zipfile
from pathlib import Path

import pytest

import dockstudio.core.queue as queue_mod
from dockstudio.engines.base import EngineManifest
from dockstudio.models import DockResult, Pose

_FAKE_PDBQT = """MODEL 1
REMARK VINA RESULT: -8.5 0.0 0.0
ROOT
ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O
ENDROOT
ENDMDL
"""


class _FakeEngine:
    manifest = EngineManifest(name="fake", scoring_functions=["vina"])

    def is_available(self):
        return True, None

    def cancel_current(self):
        pass

    def dock(self, receptor_pdbqt, ligand_pdbqt, box, params, work_dir, flex_pdbqt=None):
        work_dir.mkdir(parents=True, exist_ok=True)
        time.sleep(0.2)
        (work_dir / "out.pdbqt").write_text(_FAKE_PDBQT, encoding="utf-8")
        return DockResult(
            engine="fake",
            engine_version="9.9",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-8.5)],
            log="fake engine log line",
            run_dir=str(work_dir),
        )


def _seed(client, monkeypatch) -> str:
    fake = _FakeEngine()
    original = queue_mod.get_engine
    monkeypatch.setattr(
        queue_mod, "get_engine", lambda name: fake if name == "fake" else original(name)
    )
    registry = client.app.state.registry
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    fixtures = Path(__file__).parent / "fixtures"
    rec = prepare_receptor((fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", registry.receptors_dir)
    registry.add_receptor(rec)
    smiles = (fixtures / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()
    lig = prepare_ligand_from_smiles(smiles, "indinavir", registry.ligands_dir)
    registry.add_ligand(lig)

    pid = client.post("/api/projects", json={"name": "audit"}).json()["id"]
    body = {
        "engine": "fake",
        "receptor_id": rec.id,
        "ligand_id": lig.id,
        "box": {"center": [13.0, 22.4, 5.5], "size": [22.5, 22.5, 22.5]},
        "params": {"num_modes": 1},
    }
    jid = client.post(f"/api/projects/{pid}/jobs", json=body).json()["job_id"]
    deadline = time.time() + 20
    while time.time() < deadline:
        j = client.get(f"/api/jobs/{jid}").json()
        if j["status"] == "completed":
            return jid, pid, rec, lig
        time.sleep(0.1)
    raise TimeoutError(jid)


def test_audit_zip_structure(client, monkeypatch) -> None:
    jid, _pid, _rec, _lig = _seed(client, monkeypatch)
    resp = client.get(f"/api/jobs/{jid}/audit.zip")
    assert resp.status_code == 200
    assert "attachment" in resp.headers.get("content-disposition", "")
    tmp = Path(client.app.state.settings.data_dir) / "audit_test.zip"
    tmp.write_bytes(resp.content)
    with zipfile.ZipFile(tmp) as zf:
        names = set(zf.namelist())
        # core entries always; staged input/output pdbqt only when the engine writes them
        assert {"snapshot.json", "events.jsonl", "engine.log", "report.md", "out.pdbqt"} <= names
        snapshot = json.loads(zf.read("snapshot.json"))
        assert snapshot["job"]["id"] == jid
        assert snapshot["job"]["env"]["rdkit"]
        assert snapshot["inputs"]["receptor"]
        report = zf.read("report.md").decode("utf-8")
        for needle in ("环境指纹", "SHA256", "姿态打分", "事件序列", "免责声明", "-8.50"):
            assert needle in report, f"missing {needle!r}"
        events_text = zf.read("events.jsonl").decode("utf-8")
        assert "job.created" in events_text and "job.completed" in events_text
    tmp.unlink()


def test_project_events_endpoint(client, monkeypatch) -> None:
    _jid, pid, _rec, _lig = _seed(client, monkeypatch)
    resp = client.get(f"/api/projects/{pid}/events.jsonl")
    assert resp.status_code == 200
    assert resp.headers["content-type"].startswith("text/plain")
    lines = [json.loads(line) for line in resp.text.strip().splitlines()]
    actions = [e["action"] for e in lines]
    assert "job.created" in actions
    assert "job.started" in actions
    assert "job.completed" in actions
    created = next(e for e in lines if e["action"] == "job.created")
    assert created["payload"]["input_sha256"]["receptor"]
    env = created["payload"]["env"]
    assert env["rdkit"] and env["meeko"]
    # since filter
    later = client.get(
        f"/api/projects/{pid}/events.jsonl", params={"since": lines[-1]["ts"]}
    ).text.strip()
    assert later == ""


def test_job_env_captured(client, monkeypatch) -> None:
    jid, *_ = _seed(client, monkeypatch)
    j = client.get(f"/api/jobs/{jid}").json()
    assert j["env"]["python"]
    assert j["env"]["rdkit"]


def test_demo_example_idempotent(client) -> None:
    r1 = client.get("/api/demo/example")
    assert r1.status_code == 200
    d1 = r1.json()
    r2 = client.get("/api/demo/example")
    d2 = r2.json()
    assert d1["receptor"]["id"] == d2["receptor"]["id"]
    assert d1["ligand"]["id"] == d2["ligand"]["id"]
    assert "MK1" in (c["resname"] for c in d1["receptor"]["cocrystal_ligands"])
    assert d1["box"]["size"] == [22.5, 22.5, 22.5]
    assert "1HSG" in d1["note"]


def test_demo_library(client) -> None:
    from rdkit import Chem

    data = client.get("/api/demo/library").json()
    assert data["name"] == "demo-library"
    assert len(data["items"]) == 10
    assert any(item["name"] == "indinavir" for item in data["items"])
    for item in data["items"]:
        assert Chem.MolFromSmiles(item["smiles"]) is not None, item["smiles"]


def test_crystal_ligand_and_rmsd(client, monkeypatch, tmp_path: Path) -> None:
    """1HSG indinavir redock: crystal SDF exists, real rmsd is finite and small."""
    from dockstudio.engines.vina import VinaEngine

    ok, _ = VinaEngine().is_available()
    if not ok:
        pytest.skip("vina not available")
    registry = client.app.state.registry
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    fixtures = Path(__file__).parent / "fixtures"
    rec = prepare_receptor((fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", registry.receptors_dir)
    registry.add_receptor(rec)
    smiles = (fixtures / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()
    lig = prepare_ligand_from_smiles(smiles, "indinavir", registry.ligands_dir)
    registry.add_ligand(lig)

    pid = client.post("/api/projects", json={"name": "crystal"}).json()["id"]
    body = {
        "engine": "vina",
        "receptor_id": rec.id,
        "ligand_id": lig.id,
        "box": {"center": list(rec.cocrystal_ligands[0].center), "size": [22.5, 22.5, 22.5]},
        "params": {"exhaustiveness": 2, "num_modes": 3, "cpu": 11, "seed": 42},
    }
    jid = client.post(f"/api/projects/{pid}/jobs", json=body).json()["job_id"]
    deadline = time.time() + 240
    while time.time() < deadline:
        j = client.get(f"/api/jobs/{jid}").json()
        if j["status"] == "completed":
            break
        time.sleep(2)
    assert j["status"] == "completed", j.get("error")

    resp = client.get(f"/api/jobs/{jid}/crystal_ligand")
    assert resp.status_code == 200
    assert "MK1" in resp.text or "CLS" in resp.text or "$$$$" in resp.text

    resp = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": 1})
    assert resp.status_code == 200
    data = resp.json()
    assert data["reference"] == "MK1:B"
    assert data["rmsd"] is not None
    assert 0.0 < data["rmsd"] < 5.0

    resp = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": 99})
    assert resp.status_code == 404


def test_crystal_rmsd_no_crystal(client, monkeypatch) -> None:
    """No co-crystal: rmsd is null with a null reference."""
    jid, _pid, rec, _lig = _seed(client, monkeypatch)
    rec_dir = Path(rec.pdbqt_path).parent
    meta = rec_dir / "cocrystals.json"
    meta.write_text("[]", encoding="utf-8")
    resp = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": 1})
    assert resp.status_code == 200
    data = resp.json()
    assert data["rmsd"] is None and data["reference"] is None


def test_cli_audit_e2e(client, monkeypatch, settings) -> None:
    jid, _pid, _rec, _lig = _seed(client, monkeypatch)
    import os

    from dockstudio.config import get_settings

    os.environ["DOCKSTUDIO_DATA_DIR"] = str(settings.data_dir)
    get_settings.cache_clear()
    import io
    from contextlib import redirect_stdout

    from dockstudio.cli import _cmd_audit

    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = _cmd_audit(jid)
    assert rc == 0
    out = buf.getvalue()
    assert "环境指纹" in out
    assert "SHA256" in out or "SHA256" in out.upper()
    assert "-8.50" in out

    # missing job -> exit 2
    with redirect_stdout(buf):
        rc2 = _cmd_audit("deadbeef")
    assert rc2 == 2
