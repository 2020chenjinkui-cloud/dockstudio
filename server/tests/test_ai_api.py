"""M4 API integration: capabilities -> submit -> poll -> cancel (mock engine + fake boltz)."""

from __future__ import annotations

import sys
import textwrap
import time
from pathlib import Path

import pytest

import dockstudio.core.queue as queue_mod
from dockstudio.engines.base import EngineManifest
from dockstudio.models import DockResult, Pose

_FAKE_PDBQT = """MODEL 1
REMARK VINA RESULT: -8.0 0.0 0.0
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
        time.sleep(0.1)
        (work_dir / "out.pdbqt").write_text(_FAKE_PDBQT, encoding="utf-8")
        return DockResult(
            engine="fake",
            engine_version="0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-8.0)],
            run_dir=str(work_dir),
        )


def _boltz_pkg(tmp_path: Path) -> Path:
    pkg = tmp_path / "boltzshim" / "boltz"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text('__version__ = "9.9.9-fake"\n', encoding="utf-8")
    (pkg / "__main__.py").write_text(
        textwrap.dedent(
            """\
            import json, pathlib, sys
            args = sys.argv[1:]
            out_dir = pathlib.Path(args[args.index("--out_dir") + 1])
            pred = out_dir / "predictions" / "x"
            pred.mkdir(parents=True, exist_ok=True)
            (pred / "affinity_x.json").write_text(json.dumps(
                {"affinity_pred_value": -6.5, "affinity_probability_binary": 0.93}))
            print("fake boltz ok")
            """
        ),
        encoding="utf-8",
    )
    return pkg.parent


@pytest.fixture()
def ai_client(client, monkeypatch, tmp_path):
    fake = _FakeEngine()
    original = queue_mod.get_engine
    monkeypatch.setattr(
        queue_mod, "get_engine", lambda name: fake if name == "fake" else original(name)
    )
    registry = client.app.state.registry

    class _Receptor:
        def __init__(self):
            self.id = "r1"
            self.original_path = str(Path(__file__).parent / "fixtures" / "1HSG.pdb")
            self.pdbqt_path = str(Path(__file__))
            self.cocrystal_ligands = []

    class _Ligand:
        def __init__(self):
            self.id = "l1"
            self.smiles = "CC(C)(C)Nc1ccccc1"
            self.pdbqt_path = str(Path(__file__))
            self.name = "mol"

    registry.receptors["r1"] = _Receptor()
    registry.ligands["l1"] = _Ligand()

    monkeypatch.setenv("DOCKSTUDIO_BOLTZ", sys.executable)
    monkeypatch.setenv("PYTHONPATH", str(_boltz_pkg(tmp_path)))
    monkeypatch.delenv("DOCKSTUDIO_BOLTZ_CMD", raising=False)
    return client


def _completed_job(client) -> str:
    pid = client.post("/api/projects", json={"name": "ai"}).json()["id"]
    body = {
        "engine": "fake",
        "receptor_id": "r1",
        "ligand_id": "l1",
        "box": {"center": [0, 0, 0], "size": [20, 20, 20]},
        "params": {"num_modes": 1},
    }
    jid = client.post(f"/api/projects/{pid}/jobs", json=body).json()["job_id"]
    deadline = time.time() + 15
    while time.time() < deadline:
        j = client.get(f"/api/jobs/{jid}").json()
        if j["status"] == "completed":
            return jid
        time.sleep(0.1)
    raise TimeoutError(jid)


def test_capabilities(ai_client) -> None:
    resp = ai_client.get("/api/capabilities")
    assert resp.status_code == 200
    ai = resp.json()["ai_rescore"]
    assert ai["available"] is True
    assert ai["mode"] == "venv"


def test_capabilities_unavailable(client, monkeypatch) -> None:
    monkeypatch.delenv("DOCKSTUDIO_BOLTZ", raising=False)
    monkeypatch.delenv("DOCKSTUDIO_BOLTZ_CMD", raising=False)
    ai = client.get("/api/capabilities").json()["ai_rescore"]
    assert ai["available"] is False
    assert "Boltz" in (ai["reason"] or "")
    jid_resp = client.post("/api/jobs/whatever1/rescore/ai", json={})
    assert jid_resp.status_code == 404


def test_rescore_lifecycle(ai_client) -> None:
    jid = _completed_job(ai_client)
    resp = ai_client.post(f"/api/jobs/{jid}/rescore/ai", json={"pose_indices": [1]})
    assert resp.status_code == 202, resp.text
    task_id = resp.json()["task_id"]

    deadline = time.time() + 30
    while time.time() < deadline:
        task = ai_client.get(f"/api/rescore/{task_id}").json()
        if task["status"] in ("completed", "failed", "canceled", "interrupted"):
            break
        time.sleep(0.2)
    assert task["status"] == "completed", task.get("error")
    assert task["job_id"] == jid
    assert len(task["results"]) == 1
    r = task["results"][0]
    assert r["log_ic50"] == pytest.approx(-6.5)
    assert r["binder_probability"] == pytest.approx(0.93)
    assert r["pose_index"] == 1
    assert task["msa_server_used"] is True

    tasks = ai_client.get(f"/api/jobs/{jid}/rescore").json()
    assert len(tasks) == 1 and tasks[0]["id"] == task_id


def test_rescore_cancel_queued(ai_client) -> None:
    jid = _completed_job(ai_client)
    resp = ai_client.post(f"/api/jobs/{jid}/rescore/ai", json={"pose_indices": [1]})
    task_id = resp.json()["task_id"]
    resp = ai_client.post(f"/api/rescore/{task_id}/cancel")
    assert resp.status_code == 200
    task = resp.json()
    assert task["id"] == task_id
    assert task["status"] in ("queued", "running", "canceled", "completed")


def test_rescore_404s(ai_client) -> None:
    assert ai_client.get("/api/rescore/nope0000").status_code == 404
    assert ai_client.post("/api/rescore/nope0000/cancel").status_code == 404
    assert ai_client.get("/api/jobs/nope0000/rescore").status_code == 404
