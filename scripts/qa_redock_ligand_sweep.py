#!/usr/bin/env python
"""Supplementary diagnostic: does the *ligand preparation* RNG (ETKDG seed)
explain the historical 4.226 A / 2.104 A readings?

The product's `prepare_ligand_from_smiles` fixes ETKDG `randomSeed=42`
(prep/tcm.py::generate_conformers_and_min_energy) and feeds a single
MMFF94s-minimised conformer into Vina. For a stereo-unspecified SMILES the
resulting 3D stereoisomer *and* its starting conformation therefore depend on
the embedder version/seed, while the Vina seed stays fixed.

This script varies ONLY the ETKDG seed (private pipeline steps reused verbatim
from prep/ligand.py, which the public function calls in the same order) and
re-docks with a fixed Vina seed, so the two RNGs are separated:

    _prepare_core_mol  ->  generate_conformers_and_min_energy(seed=X)  ->  _write_ligand

    flock /home/kui/.dockstudio-heavy.lock -c \
      "DOCKSTUDIO_VINA=... python qa_redock_ligand_sweep.py <repo> <out>"
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    sys.path.insert(0, str(repo / "server"))

    from rdkit import Chem

    from dockstudio.analysis.export import export_poses
    from dockstudio.config import Settings, resolve_vina_binary
    from dockstudio.engines.vina import VinaEngine
    from dockstudio.models import Box, DockParams
    from dockstudio.prep.ligand import _prepare_core_mol, _write_ligand
    from dockstudio.prep.receptor import prepare_receptor
    from dockstudio.prep.tcm import generate_conformers_and_min_energy

    from qa_redock_1hsg import (  # reuse the verified measurement code
        BOX_SIZE,
        GOLDEN_TARGETS,
        SMILES_FIXTURE,
        crystal_ligand_reference,
        rmsd_suite,
    )

    settings = Settings(data_dir=out / "data", vina_binary=os.environ.get("DOCKSTUDIO_VINA") or None)
    for d in (settings.runs_dir, settings.prepared_dir):
        d.mkdir(parents=True, exist_ok=True)
    vina = resolve_vina_binary(settings)

    fixtures = repo / "server" / "tests" / "fixtures"
    smiles = (repo / SMILES_FIXTURE).read_text().strip()
    pdb_text = (fixtures / "1HSG.pdb").read_text()

    rec = prepare_receptor(
        (fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", settings.prepared_dir / "receptors"
    )
    rec_dir = Path(rec.pdbqt_path).parent
    meta = json.loads((rec_dir / "cocrystals.json").read_text())
    entry = next(m for m in meta if m["resname"] == "MK1")
    crystal_sdf = Path(entry["sdf_path"])
    targets = json.loads((repo / GOLDEN_TARGETS).read_text())
    template_smiles = next(t["smiles"] for t in targets["targets"] if t["pdb_id"] == "1HSG")
    _xtal_raw, _xtal_fixed, _tmpl = crystal_ligand_reference(
        pdb_text, "MK1", entry["chain"], template_smiles
    )
    center = tuple(rec.cocrystal_ligands[0].center)

    engine = VinaEngine(settings)
    raw_mol = Chem.MolFromSmiles(smiles)
    core = _prepare_core_mol(raw_mol, 3, 120, "")

    report = {"vina": str(vina), "crystal_sdf": str(crystal_sdf), "runs": []}
    for etkdg_seed in (1, 7, 2026):
        mol_3d, min_e = generate_conformers_and_min_energy(core, seed=etkdg_seed)
        lig = _write_ligand(mol_3d, "smiles", f"etkdg{etkdg_seed}", settings.prepared_dir / "lig")
        run_dir = out / f"prepseed_{etkdg_seed}"
        run_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        result = engine.dock(
            Path(rec.pdbqt_path),
            Path(lig.pdbqt_path),
            Box(center=center, size=BOX_SIZE),
            DockParams(exhaustiveness=8, num_modes=9, energy_range=3.0, cpu=4, seed=42),
            run_dir / "work",
        )
        pose_paths = export_poses(run_dir / "work" / "out.pdbqt", run_dir / "work")
        rows = []
        for i, p in enumerate(pose_paths, start=1):
            m = rmsd_suite(crystal_sdf, p, _tmpl)
            m["index"] = i
            m["affinity"] = result.poses[i - 1].affinity
            rows.append(m)
        P = [r["platform_heavy_atom_rmsd"] for r in rows if r["platform_heavy_atom_rmsd"] is not None]
        C = [r["CalcRMS_sym_inplace"] for r in rows if r["CalcRMS_sym_inplace"] is not None]
        report["runs"].append(
            {
                "etkdg_seed": etkdg_seed,
                "prepared_smiles": lig.smiles,
                "prepared_3d_inchikey": rows[0]["pose_identity"].get("inchikey_from_3d"),
                "min_energy_kcal": round(min_e, 3),
                "wall_s": round(time.time() - t0, 1),
                "top1_affinity": result.poses[0].affinity,
                "top1_platform_rmsd": P[0],
                "best_platform_rmsd": min(P),
                "top1_inplace_rmsd": C[0],
                "best_inplace_rmsd": min(C),
                "poses": rows,
            }
        )
        print(
            "ETKDG seed %-5s minE=%8.3f top1_aff=%8.3f  platform top1=%.3f best=%.3f | "
            "inplace top1=%.3f  ikey=%s"
            % (
                etkdg_seed,
                min_e,
                result.poses[0].affinity,
                P[0],
                min(P),
                C[0],
                report["runs"][-1]["prepared_3d_inchikey"],
            )
        )

    (out / "ligand-prep-sweep.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
