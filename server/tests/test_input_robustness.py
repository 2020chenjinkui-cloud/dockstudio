"""T14: upload/input robustness (size cap, empty/garbled inputs, unicode, symlink)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


def test_upload_over_50mb_rejected(client) -> None:
    big = b"ATOM  " * (55 * 1024 * 1024 // 6)
    resp = client.post("/api/prepare/receptor", files={"file": ("big.pdb", big)})
    assert resp.status_code == 413
    assert "detail" in resp.json()


def test_empty_file_rejected(client) -> None:
    resp = client.post("/api/prepare/receptor", files={"file": ("empty.pdb", b"")})
    assert resp.status_code == 400
    assert "detail" in resp.json()


def test_garbled_pdb_rejected(client) -> None:
    junk = b"\xff\xfe\x00\x01" + "垃圾数据".encode() * 100
    resp = client.post("/api/prepare/receptor", files={"file": ("junk.pdb", junk)})
    assert resp.status_code == 400


def test_wrong_extension_rejected(client, hsg_pdb_bytes: bytes) -> None:
    resp = client.post("/api/prepare/receptor", files={"file": ("rec.txt", hsg_pdb_bytes)})
    assert resp.status_code == 400
    assert "mmCIF" in resp.json()["detail"] or "unsupported" in resp.json()["detail"].lower()


def test_smiles_whitespace_truncated(client) -> None:
    resp = client.post("/api/prepare/ligand", data={"smiles": "  CCO ethanol  "})
    assert resp.status_code == 200
    assert resp.json()["smiles"] == "CCO"


def test_smiles_multi_component_hint(client) -> None:
    resp = client.post("/api/prepare/ligand", data={"smiles": "CCO.O"})
    assert resp.status_code == 400
    assert "拆分" in resp.json()["detail"] or "." in resp.json()["detail"]


def test_smiles_nonascii_rejected(client) -> None:
    resp = client.post("/api/prepare/ligand", data={"smiles": "ＣＣＯ"})  # full-width chars
    assert resp.status_code == 400


def test_unicode_filename_ligand(client) -> None:
    sdf = "乙醇\n  test\n\n  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"
    # a mol with zero atoms is invalid; build a tiny valid one instead
    sdf = (
        "吲哚\n"
        "  DockStudio\n\n"
        "  3  3  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.6000    1.0000    0.5000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "M  END\n$$$$\n"
    )
    resp = client.post(
        "/api/prepare/ligand",
        files={"file": ("乙醇配体.sdf", sdf.encode("utf-8"))},
    )
    assert resp.status_code in (200, 400)  # parse may fail on tiny mol, but never 500
    if resp.status_code != 200:
        assert "detail" in resp.json()


def test_unicode_project_name(client) -> None:
    resp = client.post("/api/projects", json={"name": "抑制剂筛选·测试🧪"})
    assert resp.status_code == 201
    pid = resp.json()["id"]
    assert client.get(f"/api/projects/{pid}").json()["name"] == "抑制剂筛选·测试🧪"


def test_file_endpoint_symlink_escape(client, settings, tmp_path: Path) -> None:
    """A symlink inside data_dir pointing outside must be rejected."""
    if not hasattr(os, "symlink"):
        pytest.skip("no symlink support")
    outside = tmp_path / "outside_secret.txt"
    outside.write_text("TOPSECRET", encoding="utf-8")
    link = settings.data_dir / "evil_link.txt"
    try:
        link.symlink_to(outside)
    except OSError as exc:
        pytest.skip(f"symlink creation not permitted on this system: {exc}")
    resp = client.get("/api/file", params={"path": str(link)})
    assert resp.status_code == 403


def test_file_endpoint_unicode_content(client, settings) -> None:
    target = settings.runs_dir / "中文内容.txt"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("中文内容 ✅", encoding="utf-8")
    resp = client.get("/api/file", params={"path": str(target)})
    assert resp.status_code == 200
    assert "中文内容 ✅" in resp.text
