"""DockStudio CLI: serve / engines / demo."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

from .config import Settings
from .engines import available_engines
from .models import Box, DockParams


def _cmd_engines() -> int:
    for engine in available_engines():
        ok, reason = engine.is_available()
        status = "可用" if ok else f"不可用: {reason}"
        print(f"{engine.manifest.name:8s} v{engine.manifest.version:12s} {status}")
        print(f"         homepage: {engine.manifest.homepage}")
        print(f"         scoring:  {', '.join(engine.manifest.scoring_functions)}")
    return 0


def _cmd_demo() -> int:
    """End-to-end demo on the bundled 1HSG fixture in a scratch data dir."""
    from .analysis.export import export_poses
    from .engines.vina import VinaEngine
    from .prep.ligand import prepare_ligand_from_smiles
    from .prep.receptor import prepare_receptor

    fixtures = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
    if not fixtures.is_dir():
        print(f"fixtures not found next to the package: {fixtures}", file=sys.stderr)
        return 1
    data_dir = Path(tempfile.mkdtemp(prefix="dockstudio_demo_"))
    print(f"demo data dir: {data_dir}")
    settings = Settings(data_dir=data_dir)

    receptor = prepare_receptor(
        (fixtures / "1HSG.pdb").read_bytes(), "1HSG.pdb", settings.prepared_dir / "receptors"
    )
    print(f"receptor prepared: {receptor.id} chains={receptor.chains}")
    smiles = (fixtures / "indinavir_smiles.txt").read_text(encoding="utf-8").strip()
    ligand = prepare_ligand_from_smiles(smiles, "indinavir", settings.prepared_dir / "ligands")
    print(f"ligand prepared: {ligand.id} rotatable={ligand.n_rotatable_bonds}")

    engine = VinaEngine(settings)
    ok, reason = engine.is_available()
    if not ok:
        print(f"vina unavailable: {reason}", file=sys.stderr)
        return 1
    mk1 = receptor.cocrystal_ligands[0]
    run_dir = data_dir / "runs" / "demo"
    result = engine.dock(
        Path(receptor.pdbqt_path),
        Path(ligand.pdbqt_path),
        Box(center=mk1.center, size=(22.5, 22.5, 22.5)),
        DockParams(exhaustiveness=4, num_modes=5, seed=42),
        run_dir,
    )
    export_poses(run_dir / "out.pdbqt", run_dir)
    print(f"dock finished: {len(result.poses)} poses in {run_dir}")
    for pose in result.poses:
        print(f"  pose {pose.index}: affinity = {pose.affinity:.2f} kcal/mol")
    print(f"top affinity: {result.poses[0].affinity:.2f} kcal/mol")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="dockstudio", description="DockStudio backend CLI")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("serve", help="start the uvicorn server")
    sub.add_parser("engines", help="list engines and availability")
    sub.add_parser("demo", help="run a full demo docking on the bundled 1HSG fixture")
    args = parser.parse_args()

    if args.command == "serve":
        from .api.app import main as serve

        serve()
        return 0
    if args.command == "engines":
        return _cmd_engines()
    if args.command == "demo":
        return _cmd_demo()
    return 2


if __name__ == "__main__":
    sys.exit(main())
