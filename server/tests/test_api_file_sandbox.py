"""Unit tests for the /api/file path sandbox."""

from __future__ import annotations


def test_file_inside_data_dir(client, settings) -> None:
    target = settings.runs_dir / "hello.txt"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("hello world", encoding="utf-8")
    resp = client.get("/api/file", params={"path": str(target)})
    assert resp.status_code == 200
    assert resp.text == "hello world"


def test_file_outside_data_dir_forbidden(client) -> None:
    resp = client.get("/api/file", params={"path": "C:/Windows/win.ini"})
    assert resp.status_code == 403
    assert "detail" in resp.json()


def test_file_traversal_forbidden(client, settings) -> None:
    sneaky = str(settings.data_dir / ".." / "secret.txt")
    resp = client.get("/api/file", params={"path": sneaky})
    assert resp.status_code == 403
