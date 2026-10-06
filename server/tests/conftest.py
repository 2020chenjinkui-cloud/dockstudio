"""Shared fixtures: isolated data dir, API test client."""

from __future__ import annotations

import os
import uuid
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture()
def data_dir(tmp_path: Path) -> Path:
    return tmp_path / f"data_{uuid.uuid4().hex[:6]}"


@pytest.fixture()
def settings(data_dir: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("DOCKSTUDIO_DATA_DIR", str(data_dir))
    from dockstudio.config import Settings

    s = Settings(data_dir=data_dir)
    s.runs_dir.mkdir(parents=True, exist_ok=True)
    s.prepared_dir.mkdir(parents=True, exist_ok=True)
    return s


@pytest.fixture()
def client(settings):
    from fastapi.testclient import TestClient

    from dockstudio.api.app import create_app

    with TestClient(create_app(settings)) as c:
        yield c


@pytest.fixture(scope="session")
def indinavir_smiles() -> str:
    return (FIXTURES / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()


@pytest.fixture(scope="session")
def hsg_pdb_bytes() -> bytes:
    return (FIXTURES / "1HSG.pdb").read_bytes()


os.environ.setdefault("PYTHONIOENCODING", "utf-8")
