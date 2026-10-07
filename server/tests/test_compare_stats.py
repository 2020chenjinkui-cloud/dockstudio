"""GET /api/compare and /api/stats."""

from __future__ import annotations

import time
from pathlib import Path

import dockstudio.core.queue as queue_mod
from dockstudio.engines.base import EngineManifest
from dockstudio.models import DockResult, Pose

_FAKE_PDBQT = """MODEL 1
REMARK VINA RESULT: -8.1 0.0 0.0
ROOT
ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
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
            engine_version="0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-8.1)],
            run_dir=str(work_dir),
        )


def _wait(client, jid: str, timeout=15.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        j = client.get(f"/api/jobs/{jid}").json()
        if j["status"] in ("completed", "failed", "canceled", "interrupted"):
            return j
        time.sleep(0.1)
    raise TimeoutError(jid)


def _seed(client, monkeypatch) -> tuple[str, str]:
    fake = _FakeEngine()
    original = queue_mod.get_engine
    monkeypatch.setattr(
        queue_mod, "get_engine", lambda name: fake if name == "fake" else original(name)
    )
    registry = client.app.state.registry

    class _E:
        pdbqt_path = str(Path(__file__))
        name = "mol"
        smiles = "CCO"
        flex_pdbqt_path = None

    registry.receptors["r1"] = _E()
    registry.ligands["l1"] = _E()

    jids = []
    for _ in range(2):
        pid = client.post("/api/projects", json={"name": "cmp"}).json()["id"]
        body = {
            "engine": "fake",
            "receptor_id": "r1",
            "ligand_id": "l1",
            "box": {"center": [0, 0, 0], "size": [20, 20, 20]},
            "params": {"num_modes": 1},
        }
        jid = client.post(f"/api/projects/{pid}/jobs", json=body).json()["job_id"]
        jids.append(jid)
    for jid in jids:
        assert _wait(client, jid)["status"] == "completed"
    return jids


def test_compare_aggregation(client, monkeypatch) -> None:
    j1, j2 = _seed(client, monkeypatch)
    resp = client.get("/api/compare", params={"job_ids": f"{j1},{j2},missing1"})
    assert resp.status_code == 200
    rows = resp.json()["rows"]
    assert len(rows) == 2  # missing id skipped
    assert {r["job_id"] for r in rows} == {j1, j2}
    row = rows[0]
    assert row["top_affinity"] == -8.1
    assert row["engine"] == "fake"
    assert row["scoring"] == "vina"
    assert row["n_poses"] == 1
    assert row["status"] == "completed"


def test_stats(client, monkeypatch) -> None:
    _seed(client, monkeypatch)
    resp = client.get("/api/stats")
    assert resp.status_code == 200
    d = resp.json()
    assert d["projects"] >= 2
    assert d["jobs"] >= 2
    assert d["completed_jobs"] >= 2
    assert d["avg_dock_seconds"] is not None and d["avg_dock_seconds"] > 0
    assert d["by_engine"].get("fake", 0) >= 2
    assert isinstance(d["batches"], int)
