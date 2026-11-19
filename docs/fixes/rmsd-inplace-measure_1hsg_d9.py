#!/usr/bin/env python
"""Real 1HSG re-docking measured through the *fixed* D9 code path.

Runs on the P4 against real vina 1.2.7 and reports, per pose:

* in-place RMSD  (``rdMolAlign.CalcRMS``)   -> the D9 verdict
* superposed RMSD(``rdMolAlign.GetBestRMS``) -> auxiliary

Everything goes through the shipped modules, so the numbers are the numbers the
platform would report:

    prep.receptor.prepare_receptor
    prep.ligand.prepare_ligand_from_smiles
    engines.vina.VinaEngine.dock
    analysis.export.export_poses
    analysis.crystal.pose_rmsd_series / best_pose_rmsd

The ligand is the *stereochemically correct* MK1 SMILES from the repo's own
golden regression set (server/tests/golden/targets.json).  The demo SMILES is
deliberately NOT used here: it has no stereo descriptors and lands on a
diastereomer of MK1 — that is a separate, out-of-scope problem.

    flock /home/kui/.dockstudio-heavy.lock -c \
      "DOCKSTUDIO_VINA=... python measure_1hsg_d9.py <repo> <out>"
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        print("usage: measure_1hsg_d9.py <repo-root> <output-dir>")
        return 2
    repo = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(repo / "server"))

    import meeko
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    import dockstudio
    from dockstudio.analysis.crystal import (
        CRITERION_METRIC,
        D9_THRESHOLD_A,
        METRIC_LABELS_ZH,
        best_pose_rmsd,
        cocrystal_reference,
        crystal_pose_rmsd,
        pose_rmsd_series,
    )
    from dockstudio.analysis.export import export_poses
    from dockstudio.config import Settings, resolve_vina_binary
    from dockstudio.engines.vina import VinaEngine
    from dockstudio.models import Box, DockParams
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    settings = Settings(data_dir=out / "data", vina_binary=os.environ.get("DOCKSTUDIO_VINA") or None)
    for d in (settings.runs_dir, settings.prepared_dir):
        d.mkdir(parents=True, exist_ok=True)
    vina = resolve_vina_binary(settings)
    vina_sha = hashlib.sha256(Path(vina).read_bytes()).hexdigest() if vina else None

    report: dict = {
        "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "host": platform.node(),
        "python": sys.version.split()[0],
        "os": platform.platform(),
        "dockstudio": {"version": getattr(dockstudio, "__version__", None), "module": dockstudio.__file__},
        "rdkit": rdkit.__version__,
        "meeko": getattr(meeko, "__version__", None),
        "git_head": subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True
        ).stdout.strip(),
        "git_branch": subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "vina_binary": str(vina) if vina else None,
        "vina_sha256": vina_sha,
        "criterion_metric": CRITERION_METRIC,
        "threshold_a": D9_THRESHOLD_A,
        "metric_labels": METRIC_LABELS_ZH,
    }

    fixtures = repo / "server" / "tests" / "fixtures"
    targets = json.loads((repo / "server" / "tests" / "golden" / "targets.json").read_text())
    target = next(t for t in targets["targets"] if t["pdb_id"] == "1HSG")
    smiles = target["smiles"]

    # ---- receptor + co-crystal reference ---------------------------------
    rec = prepare_receptor(
        (fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", settings.prepared_dir / "receptors"
    )
    rec_dir = Path(rec.pdbqt_path).parent
    reference, crystal_sdf = cocrystal_reference(rec_dir)
    assert crystal_sdf is not None, "no co-crystal SDF for 1HSG"
    center = tuple(rec.cocrystal_ligands[0].center)

    report["receptor"] = {
        "id": rec.id,
        "pdbqt_atoms": sum(
            1 for line in Path(rec.pdbqt_path).read_text().splitlines() if line.startswith("ATOM")
        ),
        "warnings": rec.warnings,
        "cocrystals": [c.model_dump() for c in rec.cocrystal_ligands],
        "crystal_reference": reference,
        "crystal_sdf": str(crystal_sdf),
        "crystal_sdf_bytes": crystal_sdf.stat().st_size,
    }
    report["box"] = {
        "center_source": "PreparedReceptor.cocrystal_ligands[0].center (MK1:B heavy-atom centroid)",
        "center": [round(v, 3) for v in center],
        "size": [22.5, 22.5, 22.5],
    }

    # ---- ligand identity: what are we actually docking? -------------------
    lig = prepare_ligand_from_smiles(smiles, "indinavir_MK1", settings.prepared_dir / "ligands")

    def _ikey(smi: str) -> str | None:
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToInchiKey(m) if m is not None else None

    # identity of the crystal ligand: bond orders from the template, stereo from 3D
    from rdkit.Chem import AllChem

    pdb_text = (fixtures / "1HSG.pdb").read_text()
    het = "\n".join(
        l
        for l in pdb_text.splitlines()
        if l.startswith("HETATM") and l[17:20].strip() == "MK1" and l[21] == "B"
    )
    raw = Chem.MolFromPDBBlock(het + "\nEND\n", sanitize=False, removeHs=False)
    fixed = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles), Chem.RemoveHs(raw))
    Chem.AssignStereochemistryFrom3D(fixed)

    report["ligand"] = {
        "smiles_source": "server/tests/golden/targets.json (pdb_id=1HSG, ccd=MK1)",
        "input_smiles": smiles,
        "prepared_smiles": lig.smiles,
        "prepared_inchikey": _ikey(lig.smiles or ""),
        "crystal_inchikey_from_3d": Chem.MolToInchiKey(fixed),
        "identity_matches_crystal": _ikey(lig.smiles or "") == Chem.MolToInchiKey(fixed),
        "n_rotatable_bonds": lig.n_rotatable_bonds,
        "min_energy_kcal": round(lig.min_energy_kcal, 3) if lig.min_energy_kcal else None,
        "conformer_method": "ETKDGv3 numConfs=10 randomSeed=42 + MMFF94s(500 its), lowest-energy conf",
        "protonation": "meeko 0.8.0 PDBQTWriterLegacy (from ETKDG/MMFF hydrogens)",
    }

    # ---- dock -------------------------------------------------------------
    engine = VinaEngine(settings)
    ok, err = engine.is_available()
    params = DockParams(exhaustiveness=8, num_modes=9, energy_range=3.0, cpu=4, seed=42)
    run_dir = out / "run"
    work = run_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    box = Box(center=center, size=(22.5, 22.5, 22.5))

    t0 = time.time()
    result = engine.dock(Path(rec.pdbqt_path), Path(lig.pdbqt_path), box, params, work)
    wall = time.time() - t0
    pose_paths = export_poses(work / "out.pdbqt", work)

    report["engine"] = {
        "name": result.engine,
        "version": result.engine_version,
        "available": ok,
        "error": err,
        "params": params.model_dump(),
        "wall_s": round(wall, 1),
        "cmd": (
            f"vina --receptor receptor.pdbqt --ligand ligand.pdbqt "
            f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} "
            f"--size_x 22.5 --size_y 22.5 --size_z 22.5 "
            f"--exhaustiveness {params.exhaustiveness} --num_modes {params.num_modes} "
            f"--energy_range {params.energy_range} --cpu {params.cpu} "
            f"--scoring {params.scoring} --out out.pdbqt --seed {params.seed}"
        ),
        "affinities": [p.affinity for p in result.poses],
    }

    # ---- D9 measurements (via the shipped helpers) ------------------------
    series = pose_rmsd_series(crystal_sdf, pose_paths)
    best = best_pose_rmsd(series)
    rows = []
    for idx in sorted(series):
        r = series[idx]
        rows.append(
            {
                "pose_index": idx,
                "affinity": result.poses[idx - 1].affinity,
                "inplace_rmsd": r.inplace_rmsd,
                "superposed_rmsd": r.superposed_rmsd,
                "d9_status": r.d9_status(),
                "metric": r.metric,
                "auxiliary_metric": r.auxiliary_metric,
            }
        )

    best_row = min(
        (r for r in rows if r["inplace_rmsd"] is not None),
        key=lambda r: r["inplace_rmsd"],
        default=None,
    )
    report["d9"] = {
        "reference": reference,
        "poses": rows,
        "top1_inplace": rows[0]["inplace_rmsd"],
        "top1_superposed": rows[0]["superposed_rmsd"],
        "best_of_n_inplace": best_row["inplace_rmsd"] if best_row else None,
        "best_of_n_superposed": best_row["superposed_rmsd"] if best_row else None,
        "best_of_n_index": best_row["pose_index"] if best_row else None,
        "verdict": best_row["d9_status"] if best_row else "SKIP",
        "helper_agrees": best_row["inplace_rmsd"] == (best.inplace_rmsd if best else None),
    }

    # ---- translation-invariance probe (the regression, on real coordinates) --
    from rdkit import Chem as _C

    pose1 = pose_paths[0]
    mol = _C.MolFromMolBlock(pose1.read_text(), sanitize=False, removeHs=False)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, _C.rdGeometry.Point3D(p.x + 10.0, p.y, p.z))
    moved = out / "run" / "pose_1_translated_10A.sdf"
    moved.write_text(_C.MolToMolBlock(mol) + "\n$$$$\n", encoding="utf-8")
    probe = crystal_pose_rmsd(crystal_sdf, moved)
    report["translation_probe"] = {
        "note": "pose_1 translated +10 A along x; internal geometry unchanged",
        "inplace_rmsd": probe.inplace_rmsd,
        "superposed_rmsd": probe.superposed_rmsd,
        "d9_status": probe.d9_status(),
    }

    (out / "measure.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())