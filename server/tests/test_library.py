"""Ligand library ingestion: .sdf multi-mol / .smi / .csv + partial failures."""

from __future__ import annotations

import pytest

from dockstudio.core.orchestrator import PrepareRegistry
from dockstudio.models import PreparationError
from dockstudio.prep.library import prepare_library


def _sdf_block(name: str) -> str:
    return (
        f"{name}\n"
        "  DockStudio\n\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.4000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "M  END\n$$$$\n"
    )


def _registry(settings) -> PrepareRegistry:
    return PrepareRegistry(settings)


def test_smi_library(settings) -> None:
    content = b"CCO ethanol\nCC(=O)O acetic_acid\nnot_a_smiles((( bad\nc1ccccc1 benzene\n"
    lib = prepare_library(content, "test.smi", _registry(settings), "小分子库")
    assert lib.name == "小分子库"
    assert len(lib.ligand_ids) == 3
    assert len(lib.failed) == 1
    assert lib.failed[0]["index"] == 2
    assert "smiles" in lib.failed[0]


def test_csv_library(settings) -> None:
    content = "id,smiles,mw\na,CCO,46\nb,c1ccccc1,78\nbad,垃圾,,\n".encode()
    lib = prepare_library(content, "test.csv", _registry(settings))
    assert len(lib.ligand_ids) == 2
    assert len(lib.failed) == 1


def test_csv_missing_smiles_col(settings) -> None:
    with pytest.raises(PreparationError, match="smiles"):
        prepare_library(b"a,b,c\n1,2,3\n", "x.csv", _registry(settings))


def test_sdf_multi_mol(settings) -> None:
    content = (_sdf_block("molA") + _sdf_block("molB")).encode("utf-8")
    lib = prepare_library(content, "multi.sdf", _registry(settings))
    assert len(lib.ligand_ids) == 2
    assert not lib.failed


def test_bad_extension(settings) -> None:
    with pytest.raises(PreparationError, match=".sdf"):
        prepare_library(b"whatever", "x.pdb", _registry(settings))


def test_library_api(client) -> None:
    resp = client.post(
        "/api/prepare/library",
        files={"file": ("批量配体.smi", b"CCO ethanol\nc1ccccc1 benzene\n")},
    )
    assert resp.status_code == 201, resp.text
    lib = resp.json()
    assert len(lib["ligand_ids"]) == 2
    assert lib["failed"] == []
    # prepared ligands registered
    listing = client.get("/api/prepare/ligands").json()
    ids = {lig["id"] for lig in listing}
    assert set(lib["ligand_ids"]) <= ids
