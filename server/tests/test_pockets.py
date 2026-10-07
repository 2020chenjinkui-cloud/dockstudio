"""Pocket detection: top pocket on 1HSG must sit near the MK1 centroid."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from dockstudio.analysis.pockets import detect_pockets, pocket_response
from dockstudio.prep.receptor import prepare_receptor


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, hsg_pdb_bytes: bytes):
    base = tmp_path_factory.mktemp("pockets")
    rec = prepare_receptor(hsg_pdb_bytes, "1HSG.pdb", base)
    protein = Path(rec.pdbqt_path).parent / "receptor_protein.pdb"
    return rec, protein


def test_top_pocket_near_mk1(prepared) -> None:
    rec, protein = prepared
    pockets = detect_pockets(protein, top_n=3)
    assert pockets, "no pockets detected on 1HSG"
    mk1 = rec.cocrystal_ligands[0]
    dist = math.dist(pockets[0].center, mk1.center)
    assert dist < 8.0, f"top pocket {dist:.1f} Å from MK1 centroid"
    p = pockets[0]
    assert all(14.9 <= s <= 40.1 for s in p.size)
    assert p.score == 1.0
    assert p.n_points >= 20


def test_pocket_cache_roundtrip(prepared, tmp_path) -> None:
    _, protein = prepared
    cache_dir = tmp_path / "cache"
    r1 = pocket_response("rec-test", protein, cache_dir, top_n=2)
    r2 = pocket_response("rec-test", protein, cache_dir, top_n=2)
    assert r1.model_dump() == r2.model_dump()
    assert r1.method == "builtin_grid"
    # wrong receptor_id must not reuse the cache
    r3 = pocket_response("other", protein, cache_dir, top_n=2)
    assert r3.receptor_id == "other"


def test_pockets_api(client, hsg_pdb_bytes: bytes) -> None:
    resp = client.post("/api/prepare/receptor", files={"file": ("1HSG.pdb", hsg_pdb_bytes)})
    rid = resp.json()["id"]
    resp = client.post("/api/pockets", json={"receptor_id": rid, "top_n": 3})
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["receptor_id"] == rid
    assert data["method"] == "builtin_grid"
    assert data["pockets"]
    assert data["pockets"][0]["score"] == 1.0
    # unknown receptor -> 400
    assert client.post("/api/pockets", json={"receptor_id": "nope"}).status_code == 400
