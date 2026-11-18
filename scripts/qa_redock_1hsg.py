#!/usr/bin/env python
"""1HSG co-crystal redocking probe for DockStudio 2.0.5 (platform code path).

Measures the heavy-atom RMSD between docked poses and the co-crystal ligand
MK1 of PDB 1HSG, driving exactly the modules the product's HTTP API drives:

    dockstudio.prep.receptor.prepare_receptor        (POST /api/prepare/receptor)
    dockstudio.prep.ligand.prepare_ligand_from_smiles(POST /api/prepare/ligand)
    dockstudio.engines.vina.VinaEngine.dock          (core/queue.py job runner)
    dockstudio.analysis.export.export_poses          (core/queue.py pose export)
    dockstudio.analysis.crystal.heavy_atom_rmsd      (GET /api/jobs/{id}/crystal_rmsd)

Nothing here imports private helpers of the docking/analysis path: the only
non-public touch point is `_strip_bond_types`-equivalent re-implemented below
so that extra RMSD variants can be computed side by side with the shipped one.

Usage
-----
    flock /home/kui/.dockstudio-heavy.lock -c \
      "/home/kui/miniconda3/envs/ai-gpu/bin/python redock_1hsg.py \
         --repo /home/kui/wt/redock --out /home/kui/qa-runs/redock-1hsg"

The script is idempotent: it re-prepares from the PDB fixture every run and
overwrites the run directories under --out.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

# --------------------------------------------------------------------------
# constants: the exact inputs the repo itself uses
# --------------------------------------------------------------------------
PDB_FIXTURE = Path("server/tests/fixtures/1HSG.pdb")
SMILES_FIXTURE = Path("server/tests/fixtures/indinavir_smiles.txt")
# Curated identity of MK1 taken from the repo's own golden regression set
# (server/tests/golden/targets.json -> pdb_id 1HSG, ccd MK1).
GOLDEN_TARGETS = Path("server/tests/golden/targets.json")

# D9 acceptance criterion (team-wide): symmetry-corrected RMSD <= 2.0 A.
D9_CRITERION_A = 2.0

# Co-crystal autobox: identical to GET /api/demo/example and to
# server/tests/test_audit.py::test_crystal_ligand_and_rmsd
BOX_SIZE = (22.5, 22.5, 22.5)

RUNS = [
    # id                        ligand    exhaustiveness num_modes seed
    ("A_fixture_exh8",          "fixture", 8,            9,        42),
    ("A2_fixture_exh8_repeat",  "fixture", 8,            9,        42),
    ("B_golden_exh8",           "golden",  8,            9,        42),
    ("C_golden_exh32",          "golden",  32,           9,        42),
    # params of the repo's own automated check
    # (server/tests/test_audit.py::test_crystal_ligand_and_rmsd)
    ("D_fixture_exh2_seed42",   "fixture", 2,            3,        42),
    # product defaults: DockParams.seed == None -> vina picks a random seed
    ("E1_fixture_exh8_noseed",  "fixture", 8,            9,        None),
    ("E2_fixture_exh8_noseed",  "fixture", 8,            9,        None),
    ("F_fixture_exh8_seed7",    "fixture", 8,            9,        7),
    ("G_fixture_exh8_seed2026", "fixture", 8,            9,        2026),
    # redock the co-crystal geometry itself: ligand built from the receptor's
    # own extracted co-crystal SDF (no SMILES round trip)
    ("H_crystallsdf_exh8",      "crystal_sdf", 8,        9,        42),
]


# --------------------------------------------------------------------------
# RMSD helpers (all on heavy atoms)
# --------------------------------------------------------------------------
def _read_molblock(path: Path):
    from rdkit import Chem

    return Chem.MolFromMolBlock(
        path.read_text(encoding="utf-8"), sanitize=False, removeHs=False
    )


def strip_bond_types(mol):
    """Mirror of dockstudio.analysis.crystal._strip_bond_types (so the RMSD
    comparison is topology-only, independent of how bond orders were perceived
    by the PDB/SDF/OpenBabel reader)."""
    from rdkit import Chem

    mol = Chem.Mol(mol)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    return mol


def _heavy(mol):
    from rdkit import Chem

    return strip_bond_types(Chem.RemoveHs(mol))


def _centroid(mol):
    conf = mol.GetConformer()
    pts = [tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    n = len(pts)
    return tuple(sum(p[k] for p in pts) / n for k in range(3))


def _dist(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


def rmsd_suite(crystal_sdf: Path, pose_sdf: Path, tmpl_mol=None) -> dict:
    """All RMSD flavours for one pose, with the algorithm named for each."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    from dockstudio.analysis.crystal import heavy_atom_rmsd  # shipped code path

    out: dict = {}
    xtal_raw = _read_molblock(crystal_sdf)
    pose_raw = _read_molblock(pose_sdf)
    out["pose_atoms_all"] = pose_raw.GetNumAtoms()
    out["crystal_atoms_all"] = xtal_raw.GetNumAtoms()

    ref = _heavy(xtal_raw)
    prb = _heavy(pose_raw)
    out["crystal_heavy"] = ref.GetNumAtoms()
    out["pose_heavy"] = prb.GetNumAtoms()

    # 0) the shipped function -------------------------------------------------
    out["platform_heavy_atom_rmsd"] = heavy_atom_rmsd(crystal_sdf, pose_sdf)

    # 1) GetBestRMS  = symmetry-corrected + rigid-body superposition
    try:
        out["GetBestRMS_sym_superposed"] = round(
            float(rdMolAlign.GetBestRMS(Chem.Mol(prb), Chem.Mol(ref))), 3
        )
    except Exception as exc:  # noqa: BLE001
        out["GetBestRMS_sym_superposed"] = None
        out["GetBestRMS_error"] = str(exc)

    # 2) CalcRMS = symmetry-corrected, coordinates used as-is (in place)
    try:
        out["CalcRMS_sym_inplace"] = round(
            float(rdMolAlign.CalcRMS(Chem.Mol(prb), Chem.Mol(ref))), 3
        )
    except Exception as exc:  # noqa: BLE001
        out["CalcRMS_sym_inplace"] = None
        out["CalcRMS_error"] = str(exc)

    # 3) symmetry-free variants: first substructure match, then
    #    (a) in place  (b) after alignment
    match = prb.GetSubstructMatch(ref)
    if match:
        conf_p = prb.GetConformer()
        conf_r = ref.GetConformer()
        n = len(match)
        sq = 0.0
        for i, pi in enumerate(match):
            dp = conf_p.GetAtomPosition(pi)
            dr = conf_r.GetAtomPosition(i)
            sq += (dp.x - dr.x) ** 2 + (dp.y - dr.y) ** 2 + (dp.z - dr.z) ** 2
        out["raw_firstmatch_inplace"] = round((sq / n) ** 0.5, 3)
        amap = list(zip(range(len(match)), match))  # ref index -> probe index
        probe_aligned = Chem.Mol(prb)
        try:
            rdMolAlign.AlignMol(probe_aligned, ref, atomMap=amap)
            out["raw_firstmatch_superposed"] = round(
                float(rdMolAlign.AlignMol(Chem.Mol(prb), ref, atomMap=amap)), 3
            )
        except Exception as exc:  # noqa: BLE001
            out["raw_firstmatch_superposed"] = None
            out["raw_firstmatch_error"] = str(exc)
    else:
        out["raw_firstmatch_inplace"] = None
        out["raw_firstmatch_superposed"] = None
        out["raw_firstmatch_error"] = "no substructure match (graph mismatch)"
        match = ()

    # 4) graph identity evidence
    out["graph_isomorphic"] = bool(match) and bool(prb.GetSubstructMatch(ref)) and bool(
        ref.GetSubstructMatch(prb)
    )

    # 5) geometry of the pose relative to the crystal site
    pc = _centroid(prb)
    xc = _centroid(ref)
    out["pose_centroid"] = [round(v, 3) for v in pc]
    out["crystal_centroid"] = [round(v, 3) for v in xc]
    out["centroid_distance"] = round(_dist(pc, xc), 3)

    # 6) which stereoisomer is actually in this pose (3D-derived)
    if tmpl_mol is not None:
        out["pose_identity"] = pose_stereo_identity(pose_sdf, tmpl_mol)
    return out


def rmsd_from_pose_blocks(crystal_sdf: Path, pose_blocks: list[tuple[int, str]]) -> dict:
    """Same suite, but fed in-memory molblocks (used for raw MODEL chunks)."""
    import tempfile

    rows = {}
    with tempfile.TemporaryDirectory() as td:
        for idx, block in pose_blocks:
            p = Path(td) / f"pose_{idx}.sdf"
            p.write_text(block, encoding="utf-8")
            rows[idx] = rmsd_suite(crystal_sdf, p)
    return rows


# --------------------------------------------------------------------------
# identity helpers
# --------------------------------------------------------------------------
def _mol_to_inchikey(mol) -> str | None:
    from rdkit import Chem

    try:
        key = Chem.MolToInchiKey(mol)
    except Exception:  # noqa: BLE001
        return None
    return key or None


def crystal_ligand_reference(pdb_text: str, resname: str, chain: str, template_smiles: str):
    """Rebuild the co-crystal ligand with correct bond orders + stereo from 3D.

    The product itself stores the co-crystal ligand as a bond-order-free SDF
    (`Chem.MolFromPDBBlock(..., sanitize=False)` -> `MolToMolBlock`), which is
    enough for the bond-type-stripped RMSD but useless for an identity check.
    Here we restore bond orders from the curated template and read the
    stereochemistry back off the crystal coordinates.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    lines = [
        line
        for line in pdb_text.splitlines()
        if line.startswith("HETATM")
        and line[17:20].strip() == resname
        and line[21] == chain
    ]
    block = "\n".join(lines) + "\nEND\n"
    raw = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    heavy = Chem.RemoveHs(raw)
    tmpl = Chem.MolFromSmiles(template_smiles)
    fixed = AllChem.AssignBondOrdersFromTemplate(tmpl, heavy)
    Chem.AssignStereochemistryFrom3D(fixed)
    return raw, fixed, tmpl


def pose_stereo_identity(pose_sdf: Path, tmpl_mol) -> dict:
    """Which stereoisomer is actually in this pose, read off its 3D geometry.

    Bond orders are restored from the curated template (the co-crystal ligand),
    then CIP labels are assigned from the coordinates; the resulting InChIKey is
    comparable with the crystal one.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out: dict = {}
    pose_raw = _read_molblock(pose_sdf)
    if pose_raw is None:
        return {"error": "pose unreadable"}
    heavy = Chem.RemoveHs(pose_raw)
    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(Chem.Mol(tmpl_mol), heavy)
        Chem.AssignStereochemistryFrom3D(fixed)
        out["smiles_from_3d"] = Chem.MolToSmiles(fixed)
        out["inchikey_from_3d"] = _mol_to_inchikey(fixed)
        out["stereo_centers"] = Chem.FindMolChiralCenters(
            fixed, useLegacyImplementation=False, includeUnassigned=True
        )
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="worktree root (read-only use)")
    ap.add_argument("--out", required=True, help="artifact directory")
    ap.add_argument("--only", nargs="*", default=None, help="subset of run ids")
    ap.add_argument("--cpu", type=int, default=4)
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(repo / "server"))

    import rdkit
    import meeko
    from rdkit import Chem

    from dockstudio.analysis.crystal import heavy_atom_rmsd
    from dockstudio.analysis.export import export_poses
    from dockstudio.config import Settings, resolve_vina_binary
    from dockstudio.engines.vina import VinaEngine
    from dockstudio.models import Box, DockParams
    from dockstudio.prep.ligand import prepare_ligand_from_file, prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    import dockstudio

    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    git_head = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    git_branch = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
    ).stdout.strip()

    settings = Settings(
        data_dir=out / "data",
        vina_binary=os.environ.get("DOCKSTUDIO_VINA") or None,
    )
    for d in (settings.runs_dir, settings.prepared_dir):
        d.mkdir(parents=True, exist_ok=True)
    vina = resolve_vina_binary(settings)
    vina_sha = hashlib.sha256(Path(vina).read_bytes()).hexdigest() if vina else None

    summary: dict = {
        "started": started,
        "host": platform.node(),
        "python": sys.version.split()[0],
        "os": platform.platform(),
        "nproc": os.cpu_count(),
        "dockstudio_version": getattr(dockstudio, "__version__", None),
        "dockstudio_module": dockstudio.__file__,
        "rdkit": rdkit.__version__,
        "meeko": getattr(meeko, "__version__", None),
        "git_head": git_head,
        "git_branch": git_branch,
        "vina_binary": str(vina) if vina else None,
        "vina_sha256": vina_sha,
        "vina_version": None,
        "box_size": list(BOX_SIZE),
        "d9_criterion_A": D9_CRITERION_A,
    }

    # ---- 0. receptor -------------------------------------------------------
    rec = prepare_receptor(
        (repo / PDB_FIXTURE).read_bytes(),
        "1HSG.pdb",
        settings.prepared_dir / "receptors",
    )
    rec_dir = Path(rec.pdbqt_path).parent
    summary["receptor"] = {
        "id": rec.id,
        "pdbqt": rec.pdbqt_path,
        "chains": rec.chains,
        "warnings": rec.warnings,
        "cocrystal_ligands": [c.model_dump() for c in rec.cocrystal_ligands],
    }
    meta = json.loads((rec_dir / "cocrystals.json").read_text())
    summary["cocrystals_json"] = meta
    entry = next(m for m in meta if m["resname"] == "MK1")
    crystal_sdf = Path(entry["sdf_path"])
    summary["crystal_sdf"] = str(crystal_sdf)

    # ---- 1. co-crystal identity -------------------------------------------
    pdb_text = (repo / PDB_FIXTURE).read_text()
    targets = json.loads((repo / GOLDEN_TARGETS).read_text())
    golden_smiles = next(t["smiles"] for t in targets["targets"] if t["pdb_id"] == "1HSG")
    fixture_smiles = (repo / SMILES_FIXTURE).read_text().strip()

    xtal_raw, xtal_fixed, tmpl = crystal_ligand_reference(pdb_text, "MK1", entry["chain"], golden_smiles)
    xtal_fixed_smiles = Chem.MolToSmiles(xtal_fixed)
    summary["identity"] = {
        "golden_targets_smiles": golden_smiles,
        "fixture_smiles": fixture_smiles,
        "crystal_MK1_smiles_from_3d": xtal_fixed_smiles,
        "crystal_MK1_inchikey": _mol_to_inchikey(xtal_fixed),
        "golden_inchikey": _mol_to_inchikey(tmpl),
        "fixture_inchikey": _mol_to_inchikey(Chem.MolFromSmiles(fixture_smiles)),
        "crystal_and_golden_identical": (
            Chem.MolToSmiles(xtal_fixed) == Chem.MolToSmiles(tmpl)
        ),
        "fixture_matches_crystal": (
            _mol_to_inchikey(Chem.MolFromSmiles(fixture_smiles))
            == _mol_to_inchikey(tmpl)
        ),
        "crystal_MK1_inchikey_skeleton_only": None,
        "fixture_inchikey_skeleton_only": None,
    }
    def _skeleton(k):
        return k.split("-")[0] if k else None
    summary["identity"]["crystal_MK1_inchikey_skeleton_only"] = _skeleton(
        summary["identity"]["crystal_MK1_inchikey"]
    )
    summary["identity"]["fixture_inchikey_skeleton_only"] = _skeleton(
        summary["identity"]["fixture_inchikey"]
    )

    # ---- 2. autobox + coverage --------------------------------------------
    xtal_center = _centroid(xtal_raw)
    cocrystal_center = tuple(rec.cocrystal_ligands[0].center)
    summary["box"] = {
        "center_source": "PreparedReceptor.cocrystal_ligands[0].center (MK1:%s)" % entry["chain"],
        "center": [round(v, 3) for v in cocrystal_center],
        "crystal_centroid": [round(v, 3) for v in xtal_center],
        "center_vs_crystal_centroid_A": round(_dist(cocrystal_center, xtal_center), 3),
        "size": list(BOX_SIZE),
    }
    conf = xtal_raw.GetConformer()
    maxdev = {"x": 0.0, "y": 0.0, "z": 0.0}
    inside = 0
    for i in range(xtal_raw.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        ok = True
        for k, axis in enumerate("xyz"):
            dev = abs(getattr(p, axis) - cocrystal_center[k])
            maxdev[axis] = max(maxdev[axis], dev)
            if dev > BOX_SIZE[k] / 2:
                ok = False
        inside += int(ok)
    summary["box"]["crystal_atoms_inside_box"] = f"{inside}/{xtal_raw.GetNumAtoms()}"
    summary["box"]["max_axis_deviation_from_center"] = {k: round(v, 3) for k, v in maxdev.items()}

    # ---- 3. ligands --------------------------------------------------------
    lig_variants = {}
    for tag, smiles in (
        ("fixture", fixture_smiles),
        ("golden", golden_smiles),
        ("crystal_sdf", None),
    ):
        if tag == "crystal_sdf":
            lig = prepare_ligand_from_file(crystal_sdf, settings.prepared_dir / "ligands")
        else:
            lig = prepare_ligand_from_smiles(
                smiles, f"indinavir_{tag}", settings.prepared_dir / "ligands"
            )
        lig_variants[tag] = {
            "prepared_id": lig.id,
            "pdbqt": lig.pdbqt_path,
            "input_smiles": smiles,
            "prepared_smiles": lig.smiles,
            "prepared_inchikey": _mol_to_inchikey(Chem.MolFromSmiles(lig.smiles)),
            "n_rotatable_bonds": lig.n_rotatable_bonds,
            "min_energy_kcal": lig.min_energy_kcal,
            "derivative_type": lig.derivative_type,
        }
    summary["ligands"] = lig_variants

    # ---- 4. docking runs ---------------------------------------------------
    engine = VinaEngine(settings)
    ok, err = engine.is_available()
    summary["vina_available"] = {"ok": ok, "error": err, "version": engine.manifest.version}

    results: dict = {}
    for run_id, lig_tag, exh, nmodes, seed in RUNS:
        if args.only and run_id not in args.only and lig_tag not in args.only:
            continue
        run_dir = out / run_id
        work = run_dir / "work"
        work.mkdir(parents=True, exist_ok=True)
        box = Box(center=cocrystal_center, size=BOX_SIZE)
        params = DockParams(
            exhaustiveness=exh, num_modes=nmodes, energy_range=3.0, cpu=args.cpu, seed=seed
        )
        t0 = time.time()
        result = engine.dock(
            Path(rec.pdbqt_path), Path(lig_variants[lig_tag]["pdbqt"]), box, params, work
        )
        wall = time.time() - t0
        pose_paths = export_poses(work / "out.pdbqt", work)
        for pose, path in zip(result.poses, pose_paths):
            pose.file = str(path)

        # raw MODE->REMARK table, to prove export order == vina rank order
        remark = []
        for line in (work / "out.pdbqt").read_text(errors="replace").splitlines():
            if line.upper().lstrip().startswith("REMARK VINA RESULT"):
                remark.append(line.split(":", 1)[1].split())

        per_pose = []
        for i, path in enumerate(pose_paths, start=1):
            metrics = rmsd_suite(crystal_sdf, path, tmpl)
            metrics["index"] = i
            metrics["affinity"] = result.poses[i - 1].affinity
            metrics["file"] = str(path)
            per_pose.append(metrics)

        d9_key_platform = "platform_heavy_atom_rmsd"
        try:
            top1 = per_pose[0][d9_key_platform]
            best_platform = min(
                p[d9_key_platform] for p in per_pose if p[d9_key_platform] is not None
            )
            best_calc = min(
                p["CalcRMS_sym_inplace"] for p in per_pose if p["CalcRMS_sym_inplace"] is not None
            )
            best_best = min(
                p["GetBestRMS_sym_superposed"]
                for p in per_pose
                if p["GetBestRMS_sym_superposed"] is not None
            )
        except (ValueError, TypeError):
            top1 = best_platform = best_calc = best_best = None

        results[run_id] = {
            "ligand_variant": lig_tag,
            "ligand_prepared_smiles": lig_variants[lig_tag]["prepared_smiles"],
            "params": params.model_dump(),
            "box": box.model_dump(),
            "wall_s": round(wall, 1),
            "engine_version": result.engine_version,
            "n_poses": len(result.poses),
            "affinities": [p.affinity for p in result.poses],
            "vina_remark_rmsd_lb_ub": remark,
            "log_tail": "\n".join(result.log.strip().splitlines()[-6:]),
            "poses": per_pose,
            "summary": {
                "top1_platform_rmsd_A": top1,
                "best_of_N_platform_rmsd_A": best_platform,
                "best_of_N_CalcRMS_sym_inplace_A": best_calc,
                "best_of_N_GetBestRMS_sym_superposed_A": best_best,
                "D9_pass_platform": (best_platform is not None and best_platform <= D9_CRITERION_A),
                "D9_pass_CalcRMS": (best_calc is not None and best_calc <= D9_CRITERION_A),
            },
        }
        summary["runs"] = results
        (out / "results.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    (out / "results.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    # ---- 5. console report -------------------------------------------------
    print(json.dumps({k: v for k, v in summary.items() if k not in ("runs",)}, indent=2, ensure_ascii=False))
    for run_id, r in results.items():
        print(f"\n== {run_id} lig={r['ligand_variant']} exh={r['params']['exhaustiveness']} "
              f"seed={r['params']['seed']} n={r['n_poses']} wall={r['wall_s']}s")
        print("   affinities:", r["affinities"])
        for p in r["poses"]:
            print(
                f"   pose{p['index']}: aff={p['affinity']:>7.3f} "
                f"platform={p['platform_heavy_atom_rmsd']} "
                f"CalcRMS={p['CalcRMS_sym_inplace']} "
                f"GetBestRMS={p['GetBestRMS_sym_superposed']} "
                f"raw_inplace={p['raw_firstmatch_inplace']} "
                f"centroid_d={p['centroid_distance']}"
            )
        print("   summary:", json.dumps(r["summary"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
