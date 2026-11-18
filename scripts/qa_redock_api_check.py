#!/usr/bin/env python
"""Cross-check: drive the real HTTP API (FastAPI TestClient) exactly like
server/tests/test_audit.py::test_crystal_ligand_and_rmsd, but print the
per-pose crystal RMSD that GET /api/jobs/{jid}/crystal_rmsd returns.

Proves that the module-level driver (scripts/qa_redock_1hsg.py) reports the
same numbers as the product endpoint for identical inputs/params.

    flock /home/kui/.dockstudio-heavy.lock -c \
      "python redock_api_check.py /home/kui/wt/redock /home/kui/qa-runs/redock-1hsg"
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    sys.path.insert(0, str(repo / "server"))

    from fastapi.testclient import TestClient

    from dockstudio.api.app import create_app
    from dockstudio.config import Settings
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    data = out / "data_api"
    settings = Settings(data_dir=data, vina_binary=os.environ.get("DOCKSTUDIO_VINA") or None)
    settings.runs_dir.mkdir(parents=True, exist_ok=True)
    settings.prepared_dir.mkdir(parents=True, exist_ok=True)

    fixtures = repo / "server" / "tests" / "fixtures"
    smiles = (fixtures / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()

    report: dict = {"steps": []}
    with TestClient(create_app(settings)) as client:
        registry = client.app.state.registry
        rec = prepare_receptor((fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", registry.receptors_dir)
        registry.add_receptor(rec)
        report["steps"].append({"step": "prepare_receptor", "id": rec.id, "warnings": rec.warnings,
                                "cocrystals": [c.model_dump() for c in rec.cocrystal_ligands]})

        lig = prepare_ligand_from_smiles(smiles, "indinavir", registry.ligands_dir)
        registry.add_ligand(lig)
        report["steps"].append({"step": "prepare_ligand", "id": lig.id, "smiles": lig.smiles,
                                "rotB": lig.n_rotatable_bonds})

        pid = client.post("/api/projects", json={"name": "crystal"}).json()["id"]
        body = {
            "engine": "vina",
            "receptor_id": rec.id,
            "ligand_id": lig.id,
            "box": {"center": list(rec.cocrystal_ligands[0].center), "size": [22.5, 22.5, 22.5]},
            "params": {"exhaustiveness": 2, "num_modes": 3, "cpu": 4, "seed": 42},
        }
        t0 = time.time()
        jid = client.post(f"/api/projects/{pid}/jobs", json=body).json()["job_id"]
        deadline = time.time() + 600
        j = {}
        while time.time() < deadline:
            j = client.get(f"/api/jobs/{jid}").json()
            if j["status"] in ("completed", "failed"):
                break
            time.sleep(2)
        report["steps"].append({"step": "dock", "job": jid, "status": j["status"],
                                "wall_s": round(time.time() - t0, 1),
                                "error": j.get("error"),
                                "poses": [p["affinity"] for p in (j.get("result") or {}).get("poses", [])]})

        crystal = client.get(f"/api/jobs/{jid}/crystal_ligand")
        report["steps"].append({"step": "crystal_ligand", "status": crystal.status_code,
                                "bytes": len(crystal.text)})

        rows = []
        n = len((j.get("result") or {}).get("poses", [])) or 1
        for idx in range(1, n + 1):
            r = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": idx})
            rows.append({"pose_index": idx, "http": r.status_code, **r.json()})
        report["steps"].append({"step": "crystal_rmsd", "rows": rows})

        r404 = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": 99})
        report["steps"].append({"step": "crystal_rmsd_out_of_range", "http": r404.status_code})

    (out / "api-check.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
