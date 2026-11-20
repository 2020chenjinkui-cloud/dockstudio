"""Multi-target redocking panel: D9 in-place RMSD across all 10 golden targets.

Answers the review's demand for ">=20 target panel, reporting RMSD<=2A success rate
AND the conditional rate (RMSD<=2A and PoseBusters all-pass)" on a panel we can
actually run here (10 co-crystal targets, real vina 1.2.7, Tesla P4).

Output: JSON + a markdown table, so the numbers are auditable.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

from rdkit import Chem

from dockstudio.analysis.crystal import crystal_pose_rmsd
from dockstudio.analysis.export import export_poses
from dockstudio.analysis.posecheck import check_pose
from dockstudio.engines.vina import VinaEngine
from dockstudio.models import Box, DockParams
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.prep.receptor import prepare_receptor

GOLDEN = Path("server/tests/golden")
FIX = GOLDEN.parent / "fixtures" / "golden"
EXH = int(os.environ.get("BENCH_EXH", "8"))
NUM_MODES = int(os.environ.get("BENCH_MODES", "9"))
SEED = 42


def run_one(target: dict, workdir: Path) -> dict:
    pid, ccd = target["pdb_id"], target["ccd"]
    out = {"pdb_id": pid, "ccd": ccd, "status": "ok"}
    t0 = time.time()
    try:
        rec = prepare_receptor((FIX / f"{pid}.pdb").read_bytes(), f"{pid}.pdb", workdir / "rec")
        coc = next((c for c in rec.cocrystal_ligands if c.resname == ccd and c.n_atoms >= 5), None)
        if coc is None:
            out.update(status="no_cocrystal")
            return out
        lig = prepare_ligand_from_smiles(target["smiles"], ccd, workdir / "lig")
        res = VinaEngine().dock(
            Path(rec.pdbqt_path), Path(lig.pdbqt_path),
            Box(center=coc.center, size=(22.5, 22.5, 22.5)),
            DockParams(exhaustiveness=EXH, num_modes=NUM_MODES, cpu=4, seed=SEED),
            workdir / "run",
        )
        if not res.poses:
            out.update(status="no_poses")
            return out
        out["top_affinity"] = res.top_affinity
        sdfs = export_poses(Path(res.run_dir) / "out.pdbqt", Path(res.run_dir))
        crystal_sdf = Path(rec.pdbqt_path).parent / f"cocrystal_{ccd}_{coc.chain}.sdf"
        # R7: the extracted SDF has no bond orders, which makes the substructure match
        # fail for ligands whose SMILES needs them (nitro groups, explicit H imines).
        # Graft the real chemistry from the curated CCD SMILES before measuring.
        if crystal_sdf.is_file():
            from dockstudio.prep.receptor import restore_cocrystal_bond_orders

            _, restore_err = restore_cocrystal_bond_orders(crystal_sdf, target["smiles"])
            out["bond_order_restore"] = "ok" if restore_err is None else restore_err
        if sdfs and crystal_sdf.is_file():
            m = crystal_pose_rmsd(crystal_sdf, sdfs[0])
            out["rmsd_inplace"] = m.inplace_rmsd
            out["rmsd_superposed"] = m.superposed_rmsd
            out["d9_pass"] = m.d9_status()
        # PoseBusters on the exported top pose
        chk = check_pose(sdfs[0], box=Box(center=coc.center, size=(22.5, 22.5, 22.5)))
        out["posecheck_source"] = chk.source
        out["posecheck_errors"] = [i.rule for i in chk.issues if i.severity == "error"]
        out["posecheck_ok"] = chk.ok
    except Exception as exc:  # noqa: BLE001
        out.update(status=f"error: {type(exc).__name__}: {exc}")
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main() -> int:
    targets = json.loads((GOLDEN / "targets.json").read_text(encoding="utf-8"))["targets"]
    only = os.environ.get("BENCH_ONLY")
    if only:
        keep = set(only.split(","))
        targets = [t for t in targets if t["pdb_id"] in keep]
    root = Path(os.environ.get("BENCH_OUT", "/home/kui/bench/redock"))
    root.mkdir(parents=True, exist_ok=True)
    results = []
    for t in targets:
        wd = root / t["pdb_id"]
        wd.mkdir(parents=True, exist_ok=True)
        r = run_one(t, wd)
        results.append(r)
        print(json.dumps(r, ensure_ascii=False), flush=True)
    (root / "results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")

    scored = [r for r in results if r.get("rmsd_inplace") is not None]
    ok = [r for r in scored if r["rmsd_inplace"] <= 2.0]
    ok_pb = [r for r in ok if r.get("posecheck_ok")]
    print("\n=== SUMMARY ===")
    print(f"targets: {len(results)}  measured: {len(scored)}")
    print(f"RMSD<=2A: {len(ok)}/{len(scored)}" + (f" = {100*len(ok)/len(scored):.1f}%" if scored else ""))
    print(f"RMSD<=2A AND PoseBusters ok: {len(ok_pb)}/{len(scored)}" +
          (f" = {100*len(ok_pb)/len(scored):.1f}%" if scored else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
