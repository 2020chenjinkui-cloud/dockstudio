"""Ligand ↔ co-crystal identity: the red→green regression and the layer contract.

Background (measured, see docs/fixes/ligand-identity-2026-11-18.md):
a SMILES without ``@``/``@@`` does not mean "the molecule without stereo" — it
names a *set* of stereoisomers, and conformer generation picks one arbitrarily.
For 1HSG/MK1 that pick is a diastereomer of the crystal ligand, yet the platform
used to report a sub-2 Å redocking RMSD for it.  These tests pin the guard.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem

from dockstudio.models import LigandIdentity, PreparedLigand
from dockstudio.prep.identity import (
    IdentityCheck,
    check_ligand_against_cocrystal,
    check_ligand_identity,
    compare_inchikeys,
    has_defined_stereo,
    inchikey_layers,
    mol_inchikey,
)

FIXTURES = Path(__file__).parent / "fixtures"
GOLDEN = Path(__file__).parent / "golden" / "targets.json"

#: The stereo-complete MK1 SMILES (single source of truth in the repo).
MK1_STEREO = next(t for t in json.loads(GOLDEN.read_text())["targets"] if t["pdb_id"] == "1HSG")[
    "smiles"
]
#: The legacy stereo-undefined SMILES that used to ship in demo_data.
MK1_FLAT = (
    "CC(C)(C)NC(=O)C1CN(CCN1CC(CC(CC2=CC=CC=C2)C(=O)NC3C(CC4=CC=CC=C34)O)O)CC5=CN=CC=C5"
)

CRYSTAL_KEY = "CBVCZFGXHXORBI-PXQQMZJSSA-N"


# ---------------------------------------------------------------------------
# the shipped demo input must BE the crystal ligand
# ---------------------------------------------------------------------------
def test_demo_smiles_is_the_crystal_ligand() -> None:
    """demo_data must carry the same chemical individual as the MK1 co-crystal."""
    demo = (FIXTURES / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()
    assert demo == MK1_STEREO, "demo fixture drifted from golden/targets.json"
    verdict = check_ligand_identity(demo, MK1_STEREO)
    assert verdict.match is True
    assert verdict.level == "exact"
    assert verdict.ligand_inchikey == CRYSTAL_KEY


def test_demo_library_indinavir_is_the_crystal_ligand() -> None:
    """The library entry named indinavir must not be a different stereoisomer."""
    lib = Path(__file__).parents[1] / "dockstudio" / "demo_data" / "demo_library.smi"
    line = next(ln for ln in lib.read_text(encoding="utf-8").splitlines() if "indinavir" in ln)
    smiles = line.rsplit(" ", 1)[0].strip()
    assert smiles == MK1_STEREO
    assert check_ligand_identity(smiles, MK1_STEREO).match is True


def test_other_library_entries_are_not_treated_as_cocrystal() -> None:
    """The 9 non-co-crystal demo drugs stay as they are — they are not MK1."""
    lib = Path(__file__).parents[1] / "dockstudio" / "demo_data" / "demo_library.smi"
    others = [
        ln.rsplit(" ", 1)[0].strip()
        for ln in lib.read_text(encoding="utf-8").splitlines()
        if ln.strip() and "indinavir" not in ln
    ]
    assert len(others) == 9
    for smiles in others:
        # each is a valid molecule that is NOT the co-crystal ligand
        assert Chem.MolFromSmiles(smiles) is not None, smiles
        assert check_ligand_identity(smiles, MK1_STEREO).level == "mismatch", smiles


# ---------------------------------------------------------------------------
# RED: the legacy flat SMILES must be rejected
# ---------------------------------------------------------------------------
def test_flat_smiles_is_rejected() -> None:
    """The old demo SMILES has no stereo → must not pass as the crystal ligand."""
    verdict = check_ligand_identity(MK1_FLAT, MK1_STEREO)
    assert verdict.match is False
    assert verdict.level == "skeleton"  # same constitution, unproven stereo
    assert verdict.ligand_inchikey == "CBVCZFGXHXORBI-UHFFFAOYSA-N"
    assert verdict.ligand_stereo_defined is False
    assert verdict.warning and "立体" in verdict.warning


def test_flat_smiles_still_shares_the_skeleton() -> None:
    """It is the same constitution — which is exactly why the bug was invisible."""
    flat, stereo = mol_inchikey(Chem.MolFromSmiles(MK1_FLAT)), mol_inchikey(
        Chem.MolFromSmiles(MK1_STEREO)
    )
    assert inchikey_layers(flat)[0] == inchikey_layers(stereo)[0] == "CBVCZFGXHXORBI"
    assert inchikey_layers(flat)[1] != inchikey_layers(stereo)[1]
    assert has_defined_stereo(flat) is False
    assert has_defined_stereo(stereo) is True


def test_wrong_molecule_is_a_mismatch() -> None:
    """A genuinely different molecule reports the hard failure, not the soft one."""
    verdict = check_ligand_identity("CC(=O)Oc1ccccc1C(=O)O", MK1_STEREO)  # aspirin
    assert verdict.match is False
    assert verdict.level == "mismatch"
    assert verdict.warning and "不是同一分子" in verdict.warning


# ---------------------------------------------------------------------------
# the check rides on the prepared ligand (not just a log line)
# ---------------------------------------------------------------------------
def test_prepare_reports_identity_on_the_model(tmp_path: Path) -> None:
    """Preparation attaches the verdict to PreparedLigand.identity."""
    from dockstudio.prep.ligand import prepare_ligand_from_smiles

    lig = prepare_ligand_from_smiles(
        MK1_FLAT, "flat", tmp_path, cocrystal_smiles=MK1_STEREO
    )
    assert isinstance(lig, PreparedLigand)
    assert isinstance(lig.identity, LigandIdentity)
    assert lig.identity.match is False
    assert lig.identity.level == "skeleton"
    assert lig.identity.reference_inchikey == CRYSTAL_KEY


def test_prepare_accepts_the_golden_smiles(tmp_path: Path) -> None:
    from dockstudio.prep.ligand import prepare_ligand_from_smiles

    lig = prepare_ligand_from_smiles(
        MK1_STEREO, "stereo", tmp_path, cocrystal_smiles=MK1_STEREO
    )
    assert lig.identity is not None
    assert lig.identity.match is True
    assert lig.identity.level == "exact"


def test_prepare_without_reference_is_unchanged(tmp_path: Path) -> None:
    """Existing callers pass no reference: behaviour and shape stay as before."""
    from dockstudio.prep.ligand import prepare_ligand_from_smiles

    lig = prepare_ligand_from_smiles(MK1_STEREO, "plain", tmp_path)
    assert lig.identity is None
    assert Path(lig.pdbqt_path).is_file()


# ---------------------------------------------------------------------------
# comparing against the receptor's extracted (bond-order-free) SDF
# ---------------------------------------------------------------------------
def test_cocrystal_sdf_only_supports_skeleton_level(tmp_path: Path) -> None:
    """The extracted HETATM SDF has no bond orders, so it cannot prove stereo.

    This is the honest limit of the check: it must say "skeleton only" instead of
    pretending the (meaningless) InChIKey of the saturated skeleton is the answer.
    """
    from dockstudio.prep.receptor import prepare_receptor

    rec = prepare_receptor(
        (FIXTURES / "1HSG.pdb").read_bytes(), "1HSG.pdb", tmp_path / "rec"
    )
    sdf = Path(rec.pdbqt_path).parent / "cocrystal_MK1_B.sdf"
    assert sdf.is_file()

    verdict = check_ligand_against_cocrystal(MK1_STEREO, sdf)
    # Not "skeleton": InChI assigns bond orders, so the bond-order-free reference
    # hashes a *different* (saturated) molecule -- its skeleton layer differs too
    # (LCBGCXTW... vs CBVCZFGX...).  Claiming "skeleton" would assert evidence we
    # do not have; the honest verdict is "unverified".
    assert verdict.level == "unverified"
    assert verdict.match is False
    assert verdict.reference_stereo_known is False
    assert verdict.warning and "无法据此判定身份" in verdict.warning

    # with a curated SMILES the same pair upgrades to an exact match
    upgraded = check_ligand_against_cocrystal(MK1_STEREO, sdf, cocrystal_smiles=MK1_STEREO)
    assert upgraded.match is True and upgraded.level == "exact"


# ---------------------------------------------------------------------------
# layer helpers
# ---------------------------------------------------------------------------
def test_compare_inchikeys_handles_missing_keys() -> None:
    verdict = compare_inchikeys(None, CRYSTAL_KEY)
    assert verdict.match is False
    assert verdict.level == "mismatch"
    assert verdict.warning and "InChIKey" in verdict.warning


def test_identity_check_is_the_shared_model() -> None:
    """The verdict type is the pydantic model, so it serialises with the rest."""
    assert IdentityCheck is LigandIdentity
    verdict = check_ligand_identity(MK1_FLAT, MK1_STEREO)
    payload = verdict.model_dump()
    assert payload["match"] is False
    assert payload["ligand_skeleton"] == "CBVCZFGXHXORBI"
    assert json.dumps(payload)  # JSON-serialisable without a custom encoder