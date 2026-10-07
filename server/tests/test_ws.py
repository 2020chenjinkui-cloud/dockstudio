"""WebSocket integration test: job events are pushed to /ws/jobs/{jid}."""

from __future__ import annotations

import concurrent.futures
import time
from pathlib import Path

import pytest


def _recv_json_bounded(ws, timeout: float = 5.0):
    """receive_json with a hard timeout so event-delivery bugs fail, never hang."""
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        return pool.submit(ws.receive_json).result(timeout=timeout)
    finally:
        pool.shutdown(wait=False, cancel_futures=True)


import dockstudio.core.queue as queue_mod
from dockstudio.engines.base import EngineManifest
from dockstudio.models import DockResult, Pose

_FAKE_PDBQT = """MODEL 1
REMARK VINA RESULT: -7.0 0.0 0.0
ROOT
ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O
ENDROOT
ENDMDL
"""


class _WsFakeEngine:
    manifest = EngineManifest(name="fake", scoring_functions=["vina"])

    def __init__(self, delay=0.3):
        self.delay = delay

    def is_available(self):
        return True, None

    def cancel_current(self):
        pass

    def dock(self, receptor_pdbqt, ligand_pdbqt, box, params, work_dir, flex_pdbqt=None):
        work_dir.mkdir(parents=True, exist_ok=True)
        time.sleep(self.delay)
        (work_dir / "out.pdbqt").write_text(_FAKE_PDBQT, encoding="utf-8")
        return DockResult(
            engine="fake",
            engine_version="0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-7.0)],
            run_dir=str(work_dir),
        )


@pytest.fixture()
def ws_client(client, monkeypatch):
    """TestClient with a fake engine patched into the queue's engine lookup."""
    original = queue_mod.get_engine
    fake = _WsFakeEngine()
    monkeypatch.setattr(
        queue_mod, "get_engine", lambda name: fake if name == "fake" else original(name)
    )
    # register fake prepared entities directly in the app's registry
    registry = client.app.state.registry

    class _E:
        pdbqt_path = str(Path(__file__))

    registry.receptors["r1"] = _E()
    registry.ligands["l1"] = _E()
    return client


def _submit(client) -> str:
    pid = client.post("/api/projects", json={"name": "ws-test"}).json()["id"]
    body = {
        "engine": "fake",
        "receptor_id": "r1",
        "ligand_id": "l1",
        "box": {"center": [0, 0, 0], "size": [20, 20, 20]},
        "params": {"num_modes": 1},
    }
    return client.post(f"/api/projects/{pid}/jobs", json=body).json()["job_id"]


def test_ws_receives_events(ws_client) -> None:
    job_id = _submit(ws_client)
    events = []
    with ws_client.websocket_connect(f"/ws/jobs/{job_id}") as ws:
        deadline = time.time() + 15
        while time.time() < deadline:
            msg = _recv_json_bounded(ws)
            events.append(msg["type"])
            if msg["type"] == "finished":
                break
    assert "finished" in events, events
    assert "status" in events or "progress" in events
    # final event carries a completed job
    final_job = ws_client.get(f"/api/jobs/{job_id}").json()
    assert final_job["status"] == "completed"
    assert final_job["result"]["poses"][0]["affinity"] == -7.0


def test_poll_fallback_header(ws_client) -> None:
    job_id = _submit(ws_client)
    resp = ws_client.get(f"/api/jobs/{job_id}")
    assert resp.headers.get("x-poll-fallback") == "true"
