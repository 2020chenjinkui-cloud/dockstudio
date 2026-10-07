"""Batch screening: submit/aggregate/cancel/csv with mock engine + real e2e."""

from __future__ import annotations

import concurrent.futures
import csv
import io
import os
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
from dockstudio.core.orchestrator import PrepareRegistry
from dockstudio.engines.base import EngineManifest
from dockstudio.models import (
    BatchCreate,
    Box,
    DockParams,
    DockResult,
    Pose,
)

_FAKE_PDBQT = """MODEL 1
REMARK VINA RESULT: -7.0 0.0 0.0
ROOT
ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C
ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O
ENDROOT
ENDMDL
"""


def _engine_available() -> bool:
    from dockstudio.engines.vina import VinaEngine

    try:
        ok, _ = VinaEngine().is_available()
        return ok
    except Exception:  # noqa: BLE001
        return False


class _BatchFakeEngine:
    manifest = EngineManifest(name="fake", scoring_functions=["vina"])

    def __init__(self, delay: float = 0.1):
        self.delay = delay
        self.cancel_calls = 0

    def is_available(self):
        return True, None

    def cancel_current(self):
        self.cancel_calls += 1

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


def _batch_req(ligand_ids: list[str]) -> BatchCreate:
    return BatchCreate(
        name="测试批",
        engine="vina",
        receptor_id="r1",
        ligand_ids=ligand_ids,
        box=Box(center=(0, 0, 0), size=(20, 20, 20)),
        params=DockParams(num_modes=1),
    )


class _Lig:
    def __init__(self, lid: str, name: str, smiles: str):
        self.id = lid
        self.name = name
        self.smiles = smiles
        self.pdbqt_path = str(Path(__file__))


def _mock_client(client, monkeypatch, delay=0.1):
    fake = _BatchFakeEngine(delay=delay)
    original = queue_mod.get_engine
    monkeypatch.setattr(
        queue_mod, "get_engine", lambda name: fake if name == "fake" else original(name)
    )
    monkeypatch.setattr(
        "dockstudio.engines.get_engine", lambda name: fake if name == "fake" else original(name)
    )
    registry: PrepareRegistry = client.app.state.registry

    class _R:
        id = "r1"
        pdbqt_path = str(Path(__file__))
        flex_pdbqt_path = None

    registry.receptors["r1"] = _R()
    for lid in ["l1", "l2", "l3", "l4", "l5"]:
        registry.ligands[lid] = _Lig(lid, f"lig-{lid}", "CCO")
    return fake


def _submit_batch(client, ligands=("l1", "l2", "l3")) -> tuple[str, str]:
    pid = client.post("/api/projects", json={"name": "btest"}).json()["id"]
    body = {
        "name": "测试批",
        "engine": "fake",
        "receptor_id": "r1",
        "ligand_ids": list(ligands),
        "box": {"center": [0, 0, 0], "size": [20, 20, 20]},
        "params": {"num_modes": 1},
    }
    resp = client.post(f"/api/projects/{pid}/batches", json=body)
    assert resp.status_code == 201, resp.text
    return pid, resp.json()["batch_id"]


def _wait_batch(client, bid: str, timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        summary = client.get(f"/api/batches/{bid}").json()
        if summary["is_terminal"]:
            return summary
        time.sleep(0.2)
    raise TimeoutError(f"batch {bid} did not finish")


def test_batch_submit_summary_results_csv(client, monkeypatch) -> None:
    _mock_client(client, monkeypatch, delay=0.05)
    _, bid = _submit_batch(client)
    summary = _wait_batch(client, bid)
    assert summary["total"] == 3
    assert summary["done"] == 3
    assert summary["failed"] == 0

    resp = client.get(f"/api/batches/{bid}/results", params={"sort": "affinity"})
    data = resp.json()
    assert data["total"] == 3
    assert all(r["top_affinity"] == -7.0 for r in data["rows"])
    assert {r["ligand_id"] for r in data["rows"]} == {"l1", "l2", "l3"}

    resp = client.get(f"/api/batches/{bid}/results.csv")
    assert resp.status_code == 200
    assert resp.headers["content-type"].startswith("text/csv")
    rows = list(csv.reader(io.StringIO(resp.text)))
    assert rows[0] == ["rank", "ligand_name", "smiles", "top_affinity", "status", "job_id", "error"]
    assert len(rows) == 4
    assert rows[1][0] == "1"
    assert rows[1][3] == "-7.000"


def test_batch_cancel(client, monkeypatch) -> None:
    _mock_client(client, monkeypatch, delay=2.0)
    _, bid = _submit_batch(client, ligands=("l1", "l2", "l3", "l4", "l5"))
    time.sleep(0.4)  # first job is running, rest queued
    resp = client.post(f"/api/batches/{bid}/cancel")
    assert resp.status_code == 200
    summary = _wait_batch(client, bid)
    assert summary["is_terminal"]
    assert summary["canceled"] >= 3


def test_batch_404(client) -> None:
    assert client.get("/api/batches/nope0000").status_code == 404
    assert client.post("/api/batches/nope0000/cancel").status_code == 404
    assert client.get("/api/batches/nope0000/results").status_code == 404
    assert client.get("/api/batches/nope0000/results.csv").status_code == 404


def test_batch_ws_events(client, monkeypatch) -> None:
    _mock_client(client, monkeypatch, delay=0.1)
    _, bid = _submit_batch(client)
    types = []
    with client.websocket_connect(f"/ws/batches/{bid}") as ws:
        deadline = time.time() + 15
        while time.time() < deadline:
            msg = _recv_json_bounded(ws)
            types.append(msg["type"])
            data = msg["data"]
            assert "total" in data and "done" in data
            if msg["type"] == "finished":
                break
    assert "finished" in types, types


@pytest.mark.skipif(not _engine_available(), reason="vina binary not available")
def test_batch_real_e2e(client, hsg_pdb_bytes: bytes) -> None:
    """5 ligands through the real vina queue."""
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCN", "Oc1ccccc1"]
    resp = client.post("/api/prepare/receptor", files={"file": ("1HSG.pdb", hsg_pdb_bytes)})
    receptor = resp.json()
    mk1 = next(c for c in receptor["cocrystal_ligands"] if c["resname"] == "MK1")
    ligand_ids = []
    for s in smiles_list:
        r = client.post("/api/prepare/ligand", data={"smiles": s})
        assert r.status_code == 200, r.text
        ligand_ids.append(r.json()["id"])

    pid = client.post("/api/projects", json={"name": "batch-e2e"}).json()["id"]
    body = {
        "engine": "vina",
        "receptor_id": receptor["id"],
        "ligand_ids": ligand_ids,
        "box": {"center": list(mk1["center"]), "size": [22.5, 22.5, 22.5]},
        "params": {"exhaustiveness": 1, "num_modes": 2, "cpu": max(1, (os.cpu_count() or 2) - 1)},
    }
    resp = client.post(f"/api/projects/{pid}/batches", json=body)
    assert resp.status_code == 201, resp.text
    bid = resp.json()["batch_id"]
    summary = _wait_batch(client, bid, timeout=600)
    assert summary["done"] + summary["failed"] == 5
    assert summary["done"] >= 3, summary
    resp = client.get(f"/api/batches/{bid}/results")
    rows = resp.json()["rows"]
    good = [r for r in rows if r["top_affinity"] is not None]
    assert len(good) >= 3
    # ethanol/benzene scale ligands: expect weak-but-plausible affinities
    assert all(-10.0 <= r["top_affinity"] <= 0.0 for r in good)
    csv_text = client.get(f"/api/batches/{bid}/results.csv").text
    assert csv_text.splitlines()[0].startswith("rank,")
    assert len(csv_text.strip().splitlines()) == 6
