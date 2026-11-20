"""Real enrichment benchmark on DUD-E (hivpr): dock a stratified sample and measure
EF1% / EF5% / AUROC / BEDROC with the platform's own compute_enrichment_metrics.

DUD-E ships ~1400 actives and ~36000 decoys per target; docking all of them is a
multi-day job on one P4, so we sample the DUD-E-recommended way -- keep every active
we can afford and a fixed multiple of decoys -- and report the sampling ratio, because
an enrichment factor is meaningless without knowing the class balance.
"""
from __future__ import annotations

import gzip
import json
import os
import random
import tempfile
from pathlib import Path

from rdkit import Chem

from dockstudio.analysis.scoring_rigor import compute_enrichment_metrics
from dockstudio.engines.vina import VinaEngine
from dockstudio.models import Box, DockParams
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.prep.receptor import prepare_receptor

WORK = Path("/home/kui/bench/enrich")
WORK.mkdir(parents=True, exist_ok=True)
N_ACTIVES = int(os.environ.get("ENRICH_ACTIVES", "60"))
N_DECOYS = int(os.environ.get("ENRICH_DECOYS", "180"))  # 1:3 -- disclosed, not hidden
SEED = 42


def load_smiles(path: Path, limit: int, rng: random.Random) -> list[tuple[str, str]]:
    """Read an SDF(.gz) and return (name, smiles), randomly subsampled."""
    opener = gzip.open if path.suffix in (".gz", ".tgz") else open
    out: list[tuple[str, str]] = []
    with opener(path, "rt", errors="replace") as fh:
        block: list[str] = []
        for line in fh:
            if line.startswith("$$$$"):
                if block:
                    txt = "".join(block)
                    mol = Chem.MolFromMolBlock(txt, sanitize=True)
                    if mol is not None:
                        smi = Chem.MolToSmiles(mol)
                        name = block[0].strip() or f"lig{len(out)}"
                        if smi:
                            out.append((name, smi))
                block = []
            else:
                block.append(line)
    rng.shuffle(out)
    return out[:limit]


def main() -> int:
    rng = random.Random(SEED)
    actives = load_smiles(Path("/tmp/dude_test.tgz"), N_ACTIVES, rng)
    decoys = load_smiles(Path("/tmp/dude_decoys.sdf.gz"), N_DECOYS, rng)
    print(f"sampled: {len(actives)} actives + {len(decoys)} decoys", flush=True)

    td = Path(tempfile.mkdtemp())
    rec = prepare_receptor(Path("/tmp/try_receptor.pdb").read_bytes(), "hivpr.pdb", td / "rec")

    # DUD-E ships a bare receptor with no HETATM ligand, so the binding site comes from
    # the separate crystal_ligand.mol2.  Compute its heavy-atom centroid rather than
    # inventing a box: a wrong site would make the whole enrichment number meaningless.
    site = None
    for c in rec.cocrystal_ligands:
        site = c
        break
    if site is None:
        from rdkit import Chem

        mol2 = Path("/tmp/dude_crystal_ligand.mol2")
        if not mol2.is_file():
            print("REFUSING to invent a box: no crystal ligand available.")
            return 2
        mol = Chem.MolFromMol2File(str(mol2), sanitize=False, removeHs=True)
        if mol is None or mol.GetNumAtoms() == 0:
            print("REFUSING to invent a box: crystal_ligand.mol2 unreadable.")
            return 2
        conf = mol.GetConformer()
        xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
        ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
        zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
        centre = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
        print(f"box centre from crystal_ligand.mol2: {centre}", flush=True)
        box = Box(center=centre, size=(22.5, 22.5, 22.5))
    else:
        box = Box(center=site.center, size=(22.5, 22.5, 22.5))

    eng = VinaEngine()
    rows: list[dict] = []
    pool = [(n, s, 1) for n, s in actives] + [(n, s, 0) for n, s in decoys]
    for idx, (name, smi, label) in enumerate(pool):
        try:
            lig = prepare_ligand_from_smiles(smi, f"lig{idx}", td / f"lig{idx}")
            res = eng.dock(
                Path(rec.pdbqt_path), Path(lig.pdbqt_path), box,
                DockParams(exhaustiveness=4, num_modes=3, cpu=4, seed=SEED),
                td / f"run{idx}",
            )
            rows.append({"name": name, "label": label, "affinity": res.top_affinity})
        except Exception as exc:  # noqa: BLE001
            rows.append({"name": name, "label": label, "affinity": None, "error": str(exc)[:120]})
        if (idx + 1) % 20 == 0:
            print(f"  docked {idx+1}/{len(pool)}", flush=True)

    scored = [r for r in rows if r.get("affinity") is not None]
    print(f"\nscored {len(scored)}/{len(rows)}", flush=True)
    if len(scored) < 10:
        print("too few scored to compute enrichment")
        return 1

    metrics = compute_enrichment_metrics(
        scores=[r["affinity"] for r in scored],
        is_active=[r["label"] for r in scored],
        higher_is_better=False,  # more negative affinity = better
    )
    n_act = sum(1 for r in scored if r["label"] == 1)
    print("\n=== ENRICHMENT (DUD-E hivpr) ===")
    print(f"actives in set: {n_act} / {len(scored)} (ratio 1:{len(scored)/max(1,n_act):.1f})")
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}")
    (WORK / "enrichment.json").write_text(
        json.dumps({"metrics": metrics, "n_scored": len(scored), "n_actives": n_act,
                    "rows": scored}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
