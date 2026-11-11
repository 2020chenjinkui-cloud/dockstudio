"""Adversarial stress-testing suite for multi-component SMILES handling and desalting remediation.

Target commit: ef3b1b6b09a9af9eb9431d5b40b398fa8e44d8b6
Milestone: M4 Victory Remediation Challenger 1
"""

from __future__ import annotations

from pathlib import Path
import pytest
from fastapi.testclient import TestClient

from dockstudio.models import PreparationError
from dockstudio.prep.ligand import prepare_ligand_from_smiles


# ---------------------------------------------------------------------------
# Direct unit-level adversarial stress on prepare_ligand_from_smiles()
# ---------------------------------------------------------------------------


class TestMultiComponentAdversarialDirect:
    def test_mixture_two_distinct_drugs_aspirin_ibuprofen(self, tmp_path: Path) -> None:
        """Aspirin + Ibuprofen: CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O MUST raise PreparationError."""
        smiles = "CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "mix_asp_ibu", tmp_path)
        err_msg = str(excinfo.value)
        assert "多个独立组分" in err_msg or "拆分" in err_msg

    def test_mixture_two_distinct_drugs_paracetamol_caffeine(self, tmp_path: Path) -> None:
        """Paracetamol + Caffeine: CC(=O)Nc1ccc(O)cc1.Cn1c(=O)c2c(ncn2C)n(C)c1=O MUST raise PreparationError."""
        smiles = "CC(=O)Nc1ccc(O)cc1.Cn1c(=O)c2c(ncn2C)n(C)c1=O"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "mix_para_caff", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_mixture_three_distinct_molecules(self, tmp_path: Path) -> None:
        """Benzene + Ethanol + Ethanamine: c1ccccc1.CCO.CCN MUST raise PreparationError."""
        smiles = "c1ccccc1.CCO.CCN"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "mix_3mols", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_mixture_three_drugs(self, tmp_path: Path) -> None:
        """Aspirin + Ibuprofen + Paracetamol: 3 active drugs MUST raise PreparationError."""
        smiles = "CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O.CC(=O)Nc1ccc(O)cc1"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "mix_3drugs", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_mixture_solvated_water_benzene_ethanol(self, tmp_path: Path) -> None:
        """Water + Benzene + Ethanol: O.c1ccccc1.CCO (water < 3 heavy, but 2 distinct substantial) MUST raise PreparationError."""
        smiles = "O.c1ccccc1.CCO"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "water_mix", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_mixture_salt_counterion_two_drugs(self, tmp_path: Path) -> None:
        """Salt counterions + 2 distinct drugs: [Na+].[Cl-].c1ccccc1.CCO MUST strip salt and raise PreparationError for remaining 2."""
        smiles = "[Na+].[Cl-].c1ccccc1.CCO"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "salt_mix", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_mixture_polyhydrate_two_drugs(self, tmp_path: Path) -> None:
        """Polyhydrate + 2 drugs: O.O.O.CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O MUST raise PreparationError."""
        smiles = "O.O.O.CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "polyhydrate_mix", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_mixture_racemic_stereoisomers(self, tmp_path: Path) -> None:
        """Mixture of stereoisomers C[C@H](O)CC.C[C@@H](O)CC MUST raise PreparationError."""
        smiles = "C[C@H](O)CC.C[C@@H](O)CC"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "mix_stereoisomers", tmp_path)
        assert "多个独立组分" in str(excinfo.value)

    def test_single_molecule_benzene_hydrate(self, tmp_path: Path) -> None:
        """Benzene + Water: c1ccccc1.O MUST desalt cleanly and produce prepared ligand for benzene."""
        smiles = "c1ccccc1.O"
        res = prepare_ligand_from_smiles(smiles, "benzene_hydrate", tmp_path)
        assert res.smiles == "c1ccccc1"
        assert Path(res.pdbqt_path).exists()
        assert Path(res.pdbqt_path).stat().st_size > 0

    def test_single_molecule_aspirin_hydrate(self, tmp_path: Path) -> None:
        """Aspirin + Water: CC(=O)Oc1ccccc1C(=O)O.O MUST desalt cleanly and produce prepared ligand for aspirin."""
        smiles = "CC(=O)Oc1ccccc1C(=O)O.O"
        res = prepare_ligand_from_smiles(smiles, "aspirin_hydrate", tmp_path)
        assert "c1ccccc1" in res.smiles
        assert Path(res.pdbqt_path).exists()
        assert Path(res.pdbqt_path).stat().st_size > 0

    def test_single_molecule_aspirin_trihydrate(self, tmp_path: Path) -> None:
        """Aspirin + 3 Waters: CC(=O)Oc1ccccc1C(=O)O.O.O.O MUST desalt cleanly."""
        smiles = "CC(=O)Oc1ccccc1C(=O)O.O.O.O"
        res = prepare_ligand_from_smiles(smiles, "aspirin_trihydrate", tmp_path)
        assert "c1ccccc1" in res.smiles
        assert Path(res.pdbqt_path).exists()

    def test_single_molecule_salt_aniline_hydrochloride(self, tmp_path: Path) -> None:
        """Aniline hydrochloride: c1ccccc1[NH3+].[Cl-] MUST desalt cleanly."""
        smiles = "c1ccccc1[NH3+].[Cl-]"
        res = prepare_ligand_from_smiles(smiles, "aniline_hcl", tmp_path)
        assert "c1ccccc1" in res.smiles
        assert Path(res.pdbqt_path).exists()

    def test_single_molecule_methanol_solvate(self, tmp_path: Path) -> None:
        """Benzene + Methanol solvate: c1ccccc1.CO (methanol < 3 heavy) MUST desalt cleanly."""
        smiles = "c1ccccc1.CO"
        res = prepare_ligand_from_smiles(smiles, "benzene_methanol", tmp_path)
        assert res.smiles == "c1ccccc1"
        assert Path(res.pdbqt_path).exists()

    def test_stoichiometric_multication_berberine_sulfate(self, tmp_path: Path) -> None:
        """Berberine sulfate (2:1 cation/sulfate): 2 identical Berberine cations + sulfate.
        Deduplication MUST prevent false rejection."""
        berberine = "COc1ccc2c(c1OC)CC[N+]3=C2Cc4cc5c(cc43)OCO5"
        smiles = f"{berberine}.{berberine}.[O-]S(=O)(=O)[O-]"
        res = prepare_ligand_from_smiles(smiles, "berberine_sulfate_2_1", tmp_path)
        assert Path(res.pdbqt_path).exists()
        assert Path(res.pdbqt_path).stat().st_size > 0

    def test_homodimer_repeated_molecule(self, tmp_path: Path) -> None:
        """Identical repeated molecule c1ccccc1.c1ccccc1: Deduplication produces 1 unique species, desalted."""
        smiles = "c1ccccc1.c1ccccc1"
        res = prepare_ligand_from_smiles(smiles, "benzene_dimer", tmp_path)
        assert res.smiles == "c1ccccc1"
        assert Path(res.pdbqt_path).exists()

    def test_trivial_pure_solvate_rejected(self, tmp_path: Path) -> None:
        """Pure solvate mixture O.O (< min_heavy) MUST be rejected by desalt_and_filter."""
        smiles = "O.O"
        with pytest.raises(PreparationError) as excinfo:
            prepare_ligand_from_smiles(smiles, "pure_water", tmp_path)
        assert "脱盐" in str(excinfo.value) or "过滤失败" in str(excinfo.value) or "重原子" in str(excinfo.value)


# ---------------------------------------------------------------------------
# End-to-end FastAPI endpoint HTTP adversarial stress on /api/prepare/ligand
# ---------------------------------------------------------------------------


class TestMultiComponentAdversarialEndpoint:
    def test_http_mixture_aspirin_ibuprofen_returns_400(self, client) -> None:
        """Aspirin + Ibuprofen over HTTP MUST return HTTP 400 with detail message."""
        smiles = "CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        resp = client.post("/api/prepare/ligand", data={"smiles": smiles})
        assert resp.status_code == 400
        detail = resp.json().get("detail", "")
        assert "拆分" in detail or "." in detail

    def test_http_mixture_three_molecules_returns_400(self, client) -> None:
        """c1ccccc1.CCO.CCN over HTTP MUST return HTTP 400."""
        resp = client.post("/api/prepare/ligand", data={"smiles": "c1ccccc1.CCO.CCN"})
        assert resp.status_code == 400
        assert "拆分" in resp.json().get("detail", "")

    def test_http_mixture_water_benzene_ethanol_returns_400(self, client) -> None:
        """O.c1ccccc1.CCO over HTTP MUST return HTTP 400."""
        resp = client.post("/api/prepare/ligand", data={"smiles": "O.c1ccccc1.CCO"})
        assert resp.status_code == 400
        assert "拆分" in resp.json().get("detail", "")

    def test_http_solvate_benzene_hydrate_returns_200(self, client) -> None:
        """c1ccccc1.O over HTTP MUST return HTTP 200 with desalted benzene."""
        resp = client.post("/api/prepare/ligand", data={"smiles": "c1ccccc1.O"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["smiles"] == "c1ccccc1"
        assert data["source"] == "smiles"

    def test_http_solvate_aspirin_hydrate_returns_200(self, client) -> None:
        """CC(=O)Oc1ccccc1C(=O)O.O over HTTP MUST return HTTP 200 with desalted aspirin."""
        resp = client.post("/api/prepare/ligand", data={"smiles": "CC(=O)Oc1ccccc1C(=O)O.O"})
        assert resp.status_code == 200
        data = resp.json()
        assert "c1ccccc1" in data["smiles"]
        assert data["source"] == "smiles"
