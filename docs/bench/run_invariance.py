"""Parameter-invariance benchmark (review item): does the ranking hold as you vary
exhaustiveness, and is a fixed seed actually reproducible?

Two experiments on the same target/ligand set:
  A) reproducibility: same seed twice -> identical affinities (T1 criterion)
  B) exhaustiveness 8 vs 16 -> how much does the top-N ranking move?
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from dockstudio.engines.vina import VinaEngine
from dockstudio.models import Box, DockParams
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.prep.receptor import prepare_receptor

FIX = Path("server/tests/fixtures")
GOLDEN = Path("server/tests/golden")
LIB = [
    ("indinavir", (FIX / "indinavir_smiles.txt").read_text().strip()),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("curcumin", "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"),
    ("baicalin", "OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccccc2)cc3=O)[C@H](O)[C@@H](O)[C@@H]1O"),
]


def dock_all(exh: int, seed: int | None, workdir: Path) -> dict[str, float]:
    rec = prepare_receptor((FIX / "1HSG.pdb").read_bytes(), "1HSG.pdb", workdir / "rec")
    mk1 = rec.cocrystal_ligands[0]
    box = Box(center=mk1.center, size=(22.5, 22.5, 22.5))
    eng = VinaEngine()
    out: dict[str, float] = {}
    for name, smi in LIB:
        lig = prepare_ligand_from_smiles(smi, name, workdir / f"lig_{name}")
        res = eng.dock(
            Path(rec.pdbqt_path), Path(lig.pdbqt_path), box,
            DockParams(exhaustiveness=exh, num_modes=5, cpu=4, seed=seed),
            workdir / f"run_{name}_e{exh}_s{seed}",
        )
        if res.top_affinity is not None:
            out[name] = round(float(res.top_affinity), 4)
    return out


def spearman(a: list[float], b: list[float]) -> float:
    n = len(a)
    if n < 2:
        return 1.0
    ra = {v: i for i, v in enumerate(sorted(a))}
    rb = {v: i for i, v in enumerate(sorted(b))}
    d2 = sum((ra[x] - rb[y]) ** 2 for x, y in zip(a, b))
    return 1 - (6 * d2) / (n * (n * n - 1))


def main() -> int:
    root = Path(os.environ.get("BENCH_OUT", "/home/kui/bench/invariance"))
    root.mkdir(parents=True, exist_ok=True)

    print("--- A) reproducibility: seed=42 twice, exh=8 ---", flush=True)
    a1 = dock_all(8, 42, Path(tempfile.mkdtemp()))
    a2 = dock_all(8, 42, Path(tempfile.mkdtemp()))
    identical = a1 == a2
    print("run1:", json.dumps(a1, ensure_ascii=False), flush=True)
    print("run2:", json.dumps(a2, ensure_ascii=False), flush=True)
    print("IDENTICAL:", identical, flush=True)

    print("\n--- B) exhaustiveness 8 vs 16 (seed=42) ---", flush=True)
    b8 = dock_all(8, 42, Path(tempfile.mkdtemp()))
    b16 = dock_all(16, 42, Path(tempfile.mkdtemp()))
    names = sorted(set(b8) & set(b16))
    rho = spearman([b8[n] for n in names], [b16[n] for n in names]) if names else 1.0
    print("exh=8 :", json.dumps(b8, ensure_ascii=False), flush=True)
    print("exh=16:", json.dumps(b16, ensure_ascii=False), flush=True)
    print(f"spearman(exh8, exh16) over {len(names)} ligands: {rho:.4f}", flush=True)

    print("\n--- C) random seed (None) twice, exh=8 ---", flush=True)
    c1 = dock_all(8, None, Path(tempfile.mkdtemp()))
    c2 = dock_all(8, None, Path(tempfile.mkdtemp()))
    print("random run1:", json.dumps(c1, ensure_ascii=False), flush=True)
    print("random run2:", json.dumps(c2, ensure_ascii=False), flush=True)
    print("IDENTICAL:", c1 == c2, flush=True)

    (root / "results.json").write_text(json.dumps(
        {"repro_seed42": [a1, a2], "identical_seed42": identical,
         "exh8": b8, "exh16": b16, "spearman_exh": rho,
         "random": [c1, c2], "identical_random": c1 == c2},
        indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
