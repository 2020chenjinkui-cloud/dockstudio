"""DockStudio Unified Command-Line Interface and Terminal Agent Package.

Provides:
- Subcommand parsing and dispatch (serve, engines, demo, prep, dock, screen, verify, audit, view, agent)
- Backward compatibility exports (_cmd_audit, _cmd_engines, _cmd_demo, main)
- Interactive Agent REPL launch on bare invocation in TTY
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from ..config import Settings, get_settings
from ..engines import available_engines
from ..models import Box, DockParams
from .agent import run_agent_repl


def _cmd_engines() -> int:
    """Print engines and availability to stdout (backward compatible)."""
    for engine in available_engines():
        ok, reason = engine.is_available()
        status = "可用" if ok else f"不可用: {reason}"
        print(f"{engine.manifest.name:8s} v{engine.manifest.version:12s} {status}")
        print(f"         homepage: {engine.manifest.homepage}")
        print(f"         scoring:  {', '.join(engine.manifest.scoring_functions)}")
    return 0


def _cmd_demo() -> int:
    """End-to-end demo on the bundled 1HSG fixture in a scratch data dir (backward compatible)."""
    from ..analysis.export import export_poses
    from ..engines.vina import VinaEngine
    from ..prep.ligand import prepare_ligand_from_smiles
    from ..prep.receptor import prepare_receptor

    fixtures = Path(__file__).resolve().parents[2] / "tests" / "fixtures"
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


def _cmd_audit(job_id: str) -> int:
    """Print the audit Markdown chain for one job (scans all projects). Backward compatible."""
    import json

    from ..analysis.audit import build_audit_markdown
    from ..core.events import EventLog
    from ..core.project import ProjectStore

    settings = get_settings()
    store = ProjectStore(settings)
    job = store.find_job(job_id)
    if job is None:
        print(
            f"未找到任务 {job_id}：请确认 job_id 来自 /api/jobs 或项目目录，"
            f"当前 data_dir={settings.data_dir}",
            file=sys.stderr,
        )
        return 2
    events = [
        e
        for e in EventLog(store.root / job.project_id).read()
        if (e.get("subject") or {}).get("id") in (job.id, job.project_id)
    ]
    sidecar = store.job_dir(job.project_id, job.id) / "inputs.sha256.json"
    input_hashes: dict = {}
    if sidecar.is_file():
        try:
            input_hashes = json.loads(sidecar.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            input_hashes = {}
    project = store.get(job.project_id)
    print(
        build_audit_markdown(
            job, events, input_hashes, project_name=(project.name if project else job.project_id)
        )
    )
    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Construct the unified argparse argument parser."""
    parser = argparse.ArgumentParser(
        prog="dockstudio",
        description="DockStudio: 面向无头 Linux 服务器与 SSH 终端的分子对接智能体 (CLI Terminal Agent)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 2.0.5")

    sub = parser.add_subparsers(dest="command", help="子命令清单")

    # 1. serve
    serve_p = sub.add_parser("serve", help="启动 FastAPI + Uvicorn 后端服务")
    serve_p.add_argument("--host", default="127.0.0.1", help="监听地址 (默认 127.0.0.1)")
    serve_p.add_argument("--port", type=int, default=8000, help="监听端口 (默认 8000)")
    serve_p.add_argument("--reload", action="store_true", help="启用代码热重载 (开发模式)")

    # 2. engines
    sub.add_parser("engines", help="探查硬件算力与计算引擎可用性 (CPU/GPU, Vina, Uni-Dock, GNINA 等)")

    # 3. demo
    sub.add_parser("demo", help="在内置 1HSG 蛋白酶上执行全流程端到端演示对接")

    # 4. prep
    prep_p = sub.add_parser("prep", help="受体 / 配体预处理向导")
    prep_p.add_argument("target", nargs="?", choices=["receptor", "ligand"], help="准备目标 (receptor 或 ligand)")
    prep_p.add_argument("-i", "--input", help="输入文件路径 (PDB 或 SDF/SMILES)")

    # 5. dock
    dock_p = sub.add_parser("dock", help="启动单分子高精度对接五步向导")
    dock_p.add_argument("--demo", action="store_true", help="使用内置 1HSG 示例运行向导")

    # 6. screen
    screen_p = sub.add_parser("screen", help="启动高通量虚拟筛选向导 (GPU Uni-Dock / 多进程)")
    screen_p.add_argument("--demo", action="store_true", help="使用内置 1HSG + 示例药物分子库运行向导")

    # 7. verify
    verify_p = sub.add_parser("verify", help="运行 T1~T5 自动化交付物数值验收套件")
    verify_p.add_argument("path", help="待核算交付包路径 (目录或 .zip 压缩包)")

    # 8. audit
    audit_p = sub.add_parser("audit", help="查看不可篡改事件链与 SHA256 完整性指纹")
    audit_p.add_argument("job_id_or_path", help="任务 ID (8位十六进制) 或项目/工作目录路径")
    audit_p.add_argument("--markdown", action="store_true", help="以原始 Markdown 格式输出 (向后兼容模式)")

    # 9. view
    view_p = sub.add_parser("view", help="启动即用即焚内网 Web 3D 快速预览 / 终端 ASCII 查看")
    view_p.add_argument("path", help="待预览的结构文件路径 (PDB/PDBQT/SDF)")
    view_p.add_argument("--port", type=int, default=0, help="HTTP 监听端口 (默认 0 自动分配)")
    view_p.add_argument("--ascii", action="store_true", help="仅在终端以 ASCII / 相互作用树展示")

    # 10. agent
    sub.add_parser("agent", help="启动交互式 REPL 智能体")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Unified entrypoint for dockstudio command-line tool."""
    parser = _build_parser()

    args_list = list(argv) if argv is not None else sys.argv[1:]
    if not args_list:
        if sys.stdin.isatty():
            return run_agent_repl()
        parser.print_help()
        return 0

    args = parser.parse_args(args_list)

    if args.command == "agent":
        return run_agent_repl()

    if args.command == "serve":
        from ..api.app import main as serve
        serve()
        return 0

    if args.command == "engines":
        from .commands.engines import run_engines_command
        return run_engines_command()

    if args.command == "demo":
        return _cmd_demo()

    if args.command == "prep":
        from .commands.prep import run_prep_ligand_interactive, run_prep_receptor_interactive
        target = (args.target or "").lower()
        if target == "receptor":
            in_p = Path(args.input) if args.input else None
            run_prep_receptor_interactive(in_p)
            return 0
        if target == "ligand":
            run_prep_ligand_interactive()
            return 0
        from .ui import select_one
        choice = select_one("选择准备目标:", [
            ("受体蛋白准备 (Receptor Prep)", "清洗、加氢、指派电荷、转 PDBQT"),
            ("小分子配体准备 (Ligand Prep)", "SMILES/SDF 质子化、低能 3D 构象生成"),
        ])
        if choice == 0:
            run_prep_receptor_interactive()
        else:
            run_prep_ligand_interactive()
        return 0

    if args.command == "dock":
        from .commands.dock import run_dock_interactive
        return run_dock_interactive(demo=args.demo)

    if args.command == "screen":
        from .commands.screen import run_screen_interactive
        return run_screen_interactive(demo=args.demo)

    if args.command == "verify":
        from .commands.verify import run_verify
        target_path = Path(args.path).expanduser().resolve()
        return run_verify(target_path)

    if args.command == "audit":
        if args.markdown:
            return _cmd_audit(args.job_id_or_path)
        from .commands.audit import run_audit
        return run_audit(args.job_id_or_path)

    if args.command == "view":
        from .commands.view import run_view
        target_path = Path(args.path).expanduser().resolve()
        return run_view(target_path, port=args.port, ascii_mode=args.ascii)

    parser.print_help()
    return 2


__all__ = [
    "_cmd_audit",
    "_cmd_demo",
    "_cmd_engines",
    "main",
    "run_agent_repl",
]
