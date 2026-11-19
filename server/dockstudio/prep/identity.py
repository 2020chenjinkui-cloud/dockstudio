"""Ligand ↔ co-crystal identity: is the molecule we docked the molecule in the crystal?

Redocking only means something if the ligand and the co-crystal reference are the
**same chemical individual**.  They routinely are not:

* a SMILES without ``@``/``@@`` descriptors is *not* "the molecule without
  stereo" — ETKDG picks one arbitrary stereoisomer, so the docked species is
  whatever the embedder happened to build (for 1HSG/MK1 it is a diastereomer of
  the crystal ligand, InChIKey ``...-RXFGFDIZSA-N`` vs ``...-PXQQMZJSSA-N``);
* a co-crystal ligand extracted from a PDB ``HETATM`` block carries **no bond
  orders** (PDB has no bond-order column), so the reference the platform stores
  as ``cocrystal_*.sdf`` is a saturated skeleton whose own InChIKey
  (``LCBGCXTWZMACCR-...``) belongs to a different molecule than the crystal
  ligand really is.

So the check has to be layered, and each layer has to be honest about what it can
actually prove:

``exact``
    Full InChIKey equality — the **skeleton layer** *and* the **stereo layer**
    (plus protonation) agree.  Both molecules are the same chemical individual.
    Only available when the reference carries sound bond orders and stereo.
``skeleton``
    Only the InChIKey skeleton layer agrees: same constitution (same atoms,
    same connectivity), stereochemistry not proven equal — either the ligand has
    undefined stereo, or the reference cannot express it.
``mismatch``
    Even the skeleton differs: a different molecule.

Everything returns an :class:`IdentityCheck`; nothing is only logged.  Callers
get ``match``, the layer they matched at, both InChIKeys, and a human-readable
``warning`` (Chinese, matching the rest of the product's user-facing text).
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from rdkit import Chem

from ..logging import get_logger
from ..models import LigandIdentity

log = get_logger("prep.identity")

#: Marker used in InChIKey stereo blocks when the molecule has defined stereo.
_STEREO_UNDEFINED = "UHFFFAOYSA"

IdentityLevel = Literal["exact", "skeleton", "mismatch", "unverified"]


def inchikey_layers(inchikey: str | None) -> tuple[str, str]:
    """Split an InChIKey into (skeleton, stereo/protonation) layers.

    A standard InChIKey is ``AAAAAAAAAAAAAA-BBBBBBBBBB-C``.  The first block is
    the connectivity ("skeleton") hash; the second encodes stereochemistry,
    isotope and protonation.  ``UHFFFAOYSA`` in the second block means the
    structure has no defined stereocentres / no stereo layer at all.
    """
    if not inchikey:
        return "", ""
    parts = inchikey.split("-")
    skeleton = parts[0] if parts else ""
    stereo = parts[1] if len(parts) > 1 else ""
    return skeleton, stereo


def has_defined_stereo(inchikey: str | None) -> bool:
    """True when the InChIKey's stereo layer actually describes stereochemistry."""
    _, stereo = inchikey_layers(inchikey)
    return bool(stereo) and stereo != _STEREO_UNDEFINED


def mol_inchikey(mol: Chem.Mol | None) -> str | None:
    """InChIKey of a molecule, or ``None`` when it cannot be generated."""
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    try:
        key = Chem.MolToInchiKey(mol)
    except Exception as exc:  # 降级兜底
        log.warning("InChIKey generation failed: %s", exc)
        return None
    return key or None


#: The verdict type.  A pydantic model (not a local dataclass) so it can ride on
#: ``PreparedLigand.identity`` and be serialised by the normal model plumbing.
IdentityCheck = LigandIdentity


def _verdict(
    *,
    match: bool,
    level: IdentityLevel,
    ligand_inchikey: str | None,
    reference_inchikey: str | None,
    reference_stereo_known: bool,
    ligand_stereo_defined: bool,
    warning: str | None = None,
) -> LigandIdentity:
    """Build the verdict, filling the derived skeleton layers."""
    return LigandIdentity(
        match=match,
        level=level,
        ligand_inchikey=ligand_inchikey,
        reference_inchikey=reference_inchikey,
        ligand_skeleton=inchikey_layers(ligand_inchikey)[0],
        reference_skeleton=inchikey_layers(reference_inchikey)[0],
        ligand_stereo_defined=ligand_stereo_defined,
        reference_stereo_known=reference_stereo_known,
        warning=warning,
    )


def compare_inchikeys(
    ligand_inchikey: str | None,
    reference_inchikey: str | None,
    *,
    reference_stereo_known: bool = True,
) -> IdentityCheck:
    """Compare two InChIKeys layer by layer (see module docstring for the levels)."""
    lig_skel, lig_stereo = inchikey_layers(ligand_inchikey)
    ref_skel, ref_stereo = inchikey_layers(reference_inchikey)
    lig_defined = has_defined_stereo(ligand_inchikey)
    ref_defined = has_defined_stereo(reference_inchikey)

    if not lig_skel or not ref_skel:
        return _verdict(
            match=False,
            level="mismatch",
            ligand_inchikey=ligand_inchikey,
            reference_inchikey=reference_inchikey,
            reference_stereo_known=reference_stereo_known,
            ligand_stereo_defined=lig_defined,
            warning="无法生成 InChIKey，不能判定配体与共晶配体是否为同一化学个体。",
        )

    if not reference_stereo_known:
        # A bond-order-free reference has no trustworthy InChIKey: InChI assigns
        # bond orders, so its skeleton hash describes a different (saturated)
        # molecule.  Saying "mismatch" would be a false alarm and "skeleton"
        # would claim evidence we do not have, so refuse to judge.
        return _verdict(
            match=False,
            level="unverified",
            ligand_inchikey=ligand_inchikey,
            reference_inchikey=reference_inchikey,
            reference_stereo_known=False,
            ligand_stereo_defined=lig_defined,
            warning=(
                "共晶配体参考无键级（来自 PDB HETATM 抽取），其 InChIKey "
                f"{reference_inchikey} 描述的是饱和骨架而非真实配体，无法据此判定身份。"
                "请同时提供该共晶配体的立体 SMILES（cocrystal_smiles）作为参考。"
            ),
        )

    if lig_skel != ref_skel:
        return _verdict(
            match=False,
            level="mismatch",
            ligand_inchikey=ligand_inchikey,
            reference_inchikey=reference_inchikey,
            reference_stereo_known=reference_stereo_known,
            ligand_stereo_defined=lig_defined,
            warning=(
                f"配体与共晶配体不是同一分子：骨架层 {lig_skel} ≠ {ref_skel}。"
                "redocking 结论不成立（测的是另一个化合物）。"
            ),
        )

    # Skeleton agrees: same constitution.  Now the stereo layer decides.
    if reference_stereo_known and ref_defined:
        if ligand_inchikey == reference_inchikey:
            return _verdict(
                match=True,
                level="exact",
                ligand_inchikey=ligand_inchikey,
                reference_inchikey=reference_inchikey,
                reference_stereo_known=True,
                ligand_stereo_defined=True,
            )
        if not lig_defined:
            return _verdict(
                match=False,
                level="skeleton",
                ligand_inchikey=ligand_inchikey,
                reference_inchikey=reference_inchikey,
                reference_stereo_known=True,
                ligand_stereo_defined=False,
                warning=(
                    f"配体未指定立体化学（InChIKey {ligand_inchikey}，立体层 {lig_stereo}=未定义），"
                    f"而共晶配体有确定立体（{ref_stereo}）。三维构象生成会任选一个立体异构体，"
                    "实际对接的很可能是非对映体 —— 请在 SMILES 中补全 @/@@ 后再对接。"
                ),
            )
        return _verdict(
            match=False,
            level="skeleton",
            ligand_inchikey=ligand_inchikey,
            reference_inchikey=reference_inchikey,
            reference_stereo_known=True,
            ligand_stereo_defined=True,
            warning=(
                f"配体与共晶配体骨架相同但立体化学不同：{lig_stereo} ≠ {ref_stereo}，"
                "两者是立体异构体（非对映体/对映体），不是同一化学个体。"
            ),
        )


def check_ligand_identity(
    ligand: Chem.Mol | str | Path,
    reference: Chem.Mol | str | Path,
    *,
    reference_stereo_known: bool | None = None,
) -> IdentityCheck:
    """Compare a ligand (mol / SMILES / .sdf / .pdb) with a co-crystal reference.

    ``reference_stereo_known`` defaults to "believe the reference" — pass
    ``False`` when the reference came from a bond-order-free PDB ``HETATM``
    extraction (:func:`dockstudio.prep.receptor` writes such SDFs), otherwise a
    meaningless stereo layer would be compared as if it meant something.
    """
    ligand_mol = _coerce_mol(ligand)
    reference_mol = _coerce_mol(reference)
    if reference_stereo_known is None:
        reference_stereo_known = True
    return compare_inchikeys(
        mol_inchikey(ligand_mol),
        mol_inchikey(reference_mol),
        reference_stereo_known=reference_stereo_known,
    )


def check_ligand_against_cocrystal(
    ligand: Chem.Mol | str | Path,
    cocrystal_sdf: Path,
    *,
    cocrystal_smiles: str | None = None,
) -> IdentityCheck:
    """Compare a ligand against a receptor's extracted co-crystal SDF.

    The SDF written by :func:`dockstudio.prep.receptor.prepare_receptor` has no
    bond orders, so on its own it only supports a skeleton comparison.  Pass
    ``cocrystal_smiles`` (a curated, stereo-complete SMILES for the same CCD
    code) to upgrade the check to ``exact``.
    """
    if cocrystal_smiles:
        return check_ligand_identity(ligand, cocrystal_smiles, reference_stereo_known=True)
    reference_mol = _coerce_mol(cocrystal_sdf)
    return check_ligand_identity(ligand, reference_mol, reference_stereo_known=False)


def _coerce_mol(source: Chem.Mol | str | Path) -> Chem.Mol | None:
    """Turn a Mol / SMILES string / molecule file into a Mol (never raises)."""
    if isinstance(source, Chem.Mol):
        return source
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8", errors="replace")
        suffix = source.suffix.lower()
        if suffix in (".sdf", ".mol"):
            block = text.split("$$$$")[0]
            return Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
        if suffix in (".pdb", ".pdbqt", ".ent"):
            lines = [ln for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM"))]
            if lines:
                return Chem.MolFromPDBBlock(
                    "\n".join(lines) + "\nEND\n", sanitize=False, removeHs=False
                )
        # unknown extension: try molblock then SMILES
        mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)
        return mol if mol is not None else Chem.MolFromSmiles(text.strip())
    # str: SMILES
    mol = Chem.MolFromSmiles(str(source).strip())
    if mol is None:
        log.warning("could not parse ligand SMILES for identity check: %r", str(source)[:60])
    return mol
