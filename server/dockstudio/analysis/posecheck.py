"""Pose quality checks: PoseBusters when available, builtin geometry fallback."""

from __future__ import annotations

import math
from pathlib import Path

from rdkit import Chem

from ..logging import get_logger
from ..models import Box, PoseCheck, PoseCheckIssue

log = get_logger("analysis.posecheck")

_CLASH_DIST = 1.6  # Å: non-bonded heavy-atom pairs closer than this collide
_OUT_OF_BOX_MARGIN = 6.0  # Å beyond the box centroid counts as a warning

_RULE_LABELS = {
    "clash": "存在原子碰撞（非键原子对距离过近）",
    "bond_length": "存在异常键长",
    "out_of_box": "位姿质心超出对接盒过多",
}


def _load_pose(pose_sdf: Path) -> Chem.Mol:
    mol = Chem.MolFromMolBlock(pose_sdf.read_text(encoding="utf-8"), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not read pose {pose_sdf}")
    return mol


# ---------------------------------------------------------------------------
# PoseBusters route
# ---------------------------------------------------------------------------


def _posebusters_available() -> bool:
    try:
        import posebusters  # noqa: F401
    except Exception:  # noqa: BLE001  # ImportError, or OSError on non-ASCII rdkit paths
        return False
    return True


def _run_posebusters(pose_sdf: Path, pose_index: int, protein_pdb: Path | None = None) -> PoseCheck:
    from posebusters import PoseBusters

    # "redock" requires a crystal reference (mol_true); for a docked pose against
    # a receptor, "dock" is the correct PoseBusters mode.  Without a protein use
    # the molecule-only mode.
    buster = PoseBusters(config="dock" if protein_pdb is not None else "mol")
    # Pass in-memory Mol objects, not paths: RDKit's file suppliers cannot open
    # non-ASCII Windows paths (e.g. C:\Users\灰灰\...), which is common on the
    # target platform.  Parsing from text works regardless of the path encoding.
    pred = Chem.MolFromMolBlock(
        pose_sdf.read_text(encoding="utf-8", errors="replace"),
        sanitize=True,
        removeHs=False,
    )
    if pred is None:
        raise ValueError(f"RDKit could not read pose for PoseBusters: {pose_sdf.name}")
    cond = None
    if protein_pdb is not None and protein_pdb.is_file():
        cond = Chem.MolFromPDBBlock(
            protein_pdb.read_text(encoding="utf-8", errors="replace"),
            sanitize=False,
            removeHs=False,
        )
    # Passing the receptor enables the protein-ligand checks (clashes, volume
    # overlap); without it PoseBusters only validates ligand-internal geometry.
    if cond is not None:
        df = buster.bust(pred, mol_cond=cond)
    else:
        df = buster.bust(pred)
    issues: list[PoseCheckIssue] = []
    ok = True
    # bust() returns a DataFrame of rule -> True/False/NaN per molecule
    import math

    if df is not None and len(df.columns):
        row = df.iloc[0]
        for rule, value in row.items():
            is_nan = isinstance(value, float) and math.isnan(value)
            passed = True if is_nan else bool(value)  # NaN -> rule not applicable -> pass
            if not passed:
                ok = False
                issues.append(
                    PoseCheckIssue(
                        rule=str(rule),
                        severity="error",
                        message=f"PoseBusters 规则 {rule} 未通过：姿态几何异常，建议检查该姿态。",
                    )
                )
    return PoseCheck(pose_index=pose_index, ok=ok, source="posebusters", issues=issues)


# ---------------------------------------------------------------------------
# Builtin fallback
# ---------------------------------------------------------------------------


def _check_clashes(mol: Chem.Mol) -> list[PoseCheckIssue]:
    issues: list[PoseCheckIssue] = []
    conf = mol.GetConformer()
    bonded = {
        (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in mol.GetBonds()
    }
    # 1,3 (geminal) pairs share a neighbor and sit at ~2.2-2.5 Å legitimately,
    # so they get a relaxed threshold instead of being skipped entirely.
    neighbors: dict[int, set[int]] = {}
    for bond in mol.GetBonds():
        neighbors.setdefault(bond.GetBeginAtomIdx(), set()).add(bond.GetEndAtomIdx())
        neighbors.setdefault(bond.GetEndAtomIdx(), set()).add(bond.GetBeginAtomIdx())
    geminal = set()
    for nbs in neighbors.values():
        for a in nbs:
            for b in nbs:
                if a < b:
                    geminal.add((a, b))
    heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetSymbol() != "H"]
    pos = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    min_d = math.inf
    clash_pair: tuple[int, int] | None = None
    geminal_threshold = 1.9  # Å: geminal pairs are fine above this
    for a, i in enumerate(heavy):
        for j in heavy[a + 1 :]:
            if (i, j) in bonded:
                continue
            threshold = geminal_threshold if (i, j) in geminal else _CLASH_DIST
            d = math.dist((pos[i].x, pos[i].y, pos[i].z), (pos[j].x, pos[j].y, pos[j].z))
            if d < threshold and d < min_d:
                min_d = d
                clash_pair = (i, j)
    if clash_pair is not None:
        issues.append(
            PoseCheckIssue(
                rule="clash",
                severity="error",
                message=(
                    f"检测到原子碰撞：非键原子 {clash_pair[0]}-{clash_pair[1]} "
                    f"距离仅 {min_d:.2f} Å（<{_CLASH_DIST} Å）。"
                ),
            )
        )
    return issues


def _check_bond_lengths(mol: Chem.Mol) -> list[PoseCheckIssue]:
    issues: list[PoseCheckIssue] = []
    conf = mol.GetConformer()
    crazy: list[tuple[int, int, float]] = []
    # generous sanity band per bond type; only flags physically impossible geometry
    bounds = {1: (0.6, 1.9), 2: (1.0, 1.6), 3: (0.9, 1.4), 12: (1.2, 1.7)}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        a, b = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
        if a.GetSymbol() == "H" or b.GetSymbol() == "H":
            continue
        pi = conf.GetAtomPosition(i)
        pj = conf.GetAtomPosition(j)
        d = math.dist((pi.x, pi.y, pi.z), (pj.x, pj.y, pj.z))
        lo, hi = bounds.get(int(bond.GetBondType()), (0.6, 2.2))
        if not (lo <= d <= hi):
            crazy.append((i, j, d))
    if crazy:
        i, j, d = crazy[0]
        issues.append(
            PoseCheckIssue(
                rule="bond_length",
                severity="error",
                message=f"键长异常：原子 {i}-{j} 距离 {d:.2f} Å，超出合理范围。",
            )
        )
    return issues


def _check_out_of_box(mol: Chem.Mol, box: Box | None) -> list[PoseCheckIssue]:
    if box is None:
        return []
    conf = mol.GetConformer()
    heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetSymbol() != "H"]
    xs = [conf.GetAtomPosition(i) for i in heavy]
    cx = sum(p.x for p in xs) / len(xs)
    cy = sum(p.y for p in xs) / len(xs)
    cz = sum(p.z for p in xs) / len(xs)
    half = tuple(s / 2.0 + _OUT_OF_BOX_MARGIN for s in box.size)
    deltas = (
        abs(cx - box.center[0]),
        abs(cy - box.center[1]),
        abs(cz - box.center[2]),
    )
    if any(d > h for d, h in zip(deltas, half)):
        return [
            PoseCheckIssue(
                rule="out_of_box",
                severity="warning",
                message="位姿质心超出对接盒边界 6 Å 以上，可能采样不足或盒子偏移。",
            )
        ]
    return []


def check_pose(
    pose_sdf: Path,
    box: Box | None = None,
    use_posebusters: bool = True,
    pose_index: int | None = None,
    protein_pdb: Path | None = None,
) -> PoseCheck:
    """Run pose quality checks. Tries PoseBusters; falls back to builtin checks."""
    idx = pose_index
    if idx is None:
        stem = pose_sdf.stem
        idx = (
            int(stem.rsplit("_", 1)[-1]) if "_" in stem and stem.rsplit("_", 1)[-1].isdigit() else 1
        )
    if use_posebusters and _posebusters_available():
        try:
            return _run_posebusters(pose_sdf, idx, protein_pdb=protein_pdb)
        except Exception as exc:  # noqa: BLE001
            log.warning("PoseBusters failed, falling back to builtin: %s", exc)
    mol = _load_pose(pose_sdf)
    issues = _check_clashes(mol) + _check_bond_lengths(mol) + _check_out_of_box(mol, box)
    ok = not any(i.severity == "error" for i in issues)
    return PoseCheck(pose_index=idx, ok=ok, source="builtin", issues=issues)
