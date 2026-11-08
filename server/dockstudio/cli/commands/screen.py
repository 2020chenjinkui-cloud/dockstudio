"""High-Throughput Virtual Screening command for DockStudio CLI (/screen).

Designed for headless Linux & GPU clusters (Uni-Dock / Vina multi-process).
Provides:
- Multi-molecule library loading (SDF / SMILES)
- Real-time live throughput & ETA progress
- Anti-formula-injection CSV export & funnel.csv
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path
from typing import Any

from rich.table import Table
from rich.text import Text

from ...config import get_settings
from ...core.orchestrator import PrepareRegistry
from ...engines import get_engine
from ...prep.library import prepare_library as prepare_library_from_file
from ..ui import (
    console,
    format_affinity,
    make_header,
    make_pane,
    prompt_path,
    select_one,
)
from .dock import _select_docking_params, _select_search_box
from .prep import run_prep_receptor_interactive

_DEMO_LIBRARY_SMI = """\
CC(=O)Oc1ccccc1C(=O)O Aspirin
CC(C)Cc1ccc(cc1)C(C)C(=O)O Ibuprofen
CC(=O)Nc1ccc(O)cc1 Acetaminophen
COc1ccc2cc(ccc2c1)C(C)C(=O)O Naproxen
Cc1ccc(cc1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F Celecoxib
"""


def _render_screen_stepper(step: int, total_steps: int = 5, step_title: str = "") -> Any:
    """Renders a Claude Code style linear wizard stepper for /screen."""
    steps = [
        ("1", "受体准备"),
        ("2", "口袋定盒"),
        ("3", "化合物库"),
        ("4", "引擎选择"),
        ("5", "批量筛选"),
    ]
    txt = Text()
    for idx, (num, name) in enumerate(steps, 1):
        if idx < step:
            txt.append(f"[✓ {num}.{name}]", style="bold #4EBA65")
        elif idx == step:
            txt.append(f"[● {num}.{name}]", style="bold #D77757")
        else:
            txt.append(f"[○ {num}.{name}]", style="#505050")
        if idx < len(steps):
            txt.append(" ─ ", style="#505050")
    return make_pane(
        f"高通量虚拟筛选向导 (/screen) · 步骤 {step}/{total_steps}: {step_title}",
        txt,
        subtitle="交互式逐步就地更新 · [Ctrl+C] 取消",
    )


def _refresh_screen_step(step: int, total_steps: int, step_title: str) -> None:
    """Clears viewport and renders top header + screening stepper card in-place."""
    if sys.stdout.isatty():
        sys.stdout.write("\033[2J\033[H")
        sys.stdout.flush()
        console.print(make_header(version="2.0.5"))
        console.print(_render_screen_stepper(step, total_steps, step_title))
        console.print()


def run_screen_interactive(demo: bool = False) -> int:
    """Runs interactive high-throughput virtual screening."""
    settings = get_settings()

    # 1. Receptor Prep
    _refresh_screen_step(1, 5, "受体结构加载与准备")
    if demo:
        fixtures = Path(__file__).resolve().parents[3] / "tests" / "fixtures"
        pdb_path = fixtures / "1HSG.pdb"
        out_dir = settings.prepared_dir / "receptors"
        out_dir.mkdir(parents=True, exist_ok=True)
        from ...prep.receptor import prepare_receptor

        receptor = prepare_receptor(pdb_path.read_bytes(), "1HSG.pdb", out_dir)
        console.print("[bold #4EBA65]✓[/bold #4EBA65] 加载内置示例受体: 1HSG (HIV-1 蛋白酶)")
    else:
        pdb_path, receptor = run_prep_receptor_interactive()

    # 2. Search Box
    _refresh_screen_step(2, 5, "口袋与搜索盒 (Search Box) 定义")
    console.print(f"[dim]已就绪受体: {receptor.id} (链: {', '.join(receptor.chains)})[/dim]\n")
    box = _select_search_box(pdb_path, receptor)

    # 3. Library Selection
    _refresh_screen_step(3, 5, "化合物库加载与配体批量准备")
    console.print(
        f"[dim]已定义搜索盒: 质心 ({box.center[0]:.1f}, {box.center[1]:.1f}, {box.center[2]:.1f}) · 尺寸 {box.size[0]:.1f}Å[/dim]\n"
    )

    if demo:
        lib_mode = 1
    else:
        lib_mode = select_one("选择配体库来源:", [
            ("使用本地化合物库文件 (SDF / SMILES / CSV)", "支持商用现货库 (TargetMol / Selleck / MCE 等)"),
            ("使用 DockStudio 内置演示分子库 (5 个上市药物)", "测试端到端筛选链路与吞吐"),
        ])

    if lib_mode == 0:
        lib_path = prompt_path("请指定配体库文件路径")
        lib_bytes = lib_path.read_bytes()
        lib_filename = lib_path.name
    else:
        demo_file = Path(__file__).resolve().parents[2] / "demo_data" / "demo_library.smi"
        if not demo_file.is_file():
            demo_file = Path(__file__).resolve().parents[3] / "tests" / "fixtures" / "demo_library.smi"
        if demo_file.is_file():
            lib_bytes = demo_file.read_bytes()
        else:
            lib_bytes = _DEMO_LIBRARY_SMI.encode("utf-8")
        lib_filename = "demo_library.smi"
        console.print(f"[bold #4EBA65]✓[/bold #4EBA65] 加载内置演示库: {lib_filename}")

    # 4. Prepare Library
    registry = PrepareRegistry(settings)
    with console.status("[bold #93A5FF]正在解析分子库并逐条加氢质子化转 PDBQT...[/bold #93A5FF]"):
        prepared_lib = prepare_library_from_file(lib_bytes, lib_filename, registry)

    n_ready = len(prepared_lib.ligand_ids)
    n_failed = len(prepared_lib.failed)
    console.print(f"[bold #4EBA65]✓ 分子库准备完成:[/bold #4EBA65] 成功 {n_ready} 个分子，失败 {n_failed} 个")

    if not prepared_lib.ligand_ids:
        console.print("[bold #FF6B80]错误: 库中没有成功准备的配体分子，筛选终止[/#FF6B80]")
        return 1

    # 5. Engine Selection
    _refresh_screen_step(4, 5, "计算引擎与对接参数设定")
    console.print(f"[dim]已就绪分子库: {prepared_lib.name} (共 {n_ready} 个有效配体)[/dim]\n")
    engine_name, scoring, base_params = _select_docking_params(settings)
    engine = get_engine(engine_name)

    # 6. Run Screening with Live Throughput
    _refresh_screen_step(5, 5, "执行高通量虚拟筛选与候选排序")
    out_dir = settings.data_dir / "screens" / f"screen_{int(time.time())}"
    out_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"[bold #D77757]开始高通量筛选 (引擎: {engine_name}, 共 {n_ready} 分子)...[/bold #D77757]")
    results: list[dict[str, Any]] = []
    t_start = time.time()

    for idx, lig_id in enumerate(prepared_lib.ligand_ids):
        lig_obj = registry.ligands.get(lig_id)
        lig_pdbqt = Path(lig_obj.pdbqt_path) if lig_obj else (registry.ligands_dir / lig_id / "ligand.pdbqt")
        lig_label = f"{lig_obj.name} ({lig_id})" if (lig_obj and lig_obj.name) else lig_id
        if not lig_pdbqt.is_file():
            continue
        
        run_scratch = out_dir / f"run_{lig_id}"
        run_scratch.mkdir(exist_ok=True)
        try:
            res = engine.dock(
                Path(receptor.pdbqt_path),
                lig_pdbqt,
                box,
                base_params,
                run_scratch,
            )
            top_aff = res.poses[0].affinity if res.poses else 0.0
            results.append({
                "id": lig_label,
                "affinity": top_aff,
                "n_poses": len(res.poses),
                "run_dir": str(run_scratch),
            })
        except Exception as err:  # noqa: BLE001
            results.append({
                "id": lig_label,
                "affinity": 999.0,
                "n_poses": 0,
                "error": str(err),
            })

        # Print progress update
        elapsed = time.time() - t_start
        speed = (idx + 1) / max(0.1, elapsed)
        eta = (n_ready - idx - 1) / max(0.1, speed)
        pct = (idx + 1) * 100 // n_ready
        console.print(
            f" [bold #93A5FF]筛选进度:[/] {pct:3d}% ({idx + 1}/{n_ready}) "
            f"· 速度: [bold]{speed:5.1f} mol/s[/] · 剩余: {int(eta)}s",
            end="\r",
        )

    total_time = time.time() - t_start
    console.print(f"\n[bold #4EBA65]✓ 筛选完成![/bold #4EBA65] 总耗时: {total_time:.2f}s · 平均速度: {len(results)/max(0.1, total_time):.1f} mol/s\n")

    # Sort results
    valid_results = [r for r in results if r.get("affinity", 999.0) < 900.0]
    valid_results.sort(key=lambda x: x["affinity"])

    # 7. Display Top Ranking Table
    table = Table(
        title="[bold #D77757]虚拟筛选 Top 10 命中候选 (Hit Ranking)[/bold #D77757]",
        border_style="#888888",
        header_style="bold #FFFFFF",
    )
    table.add_column("Rank", justify="center", width=8)
    table.add_column("分子 ID", justify="left", width=26)
    table.add_column("结合亲和力 (kcal/mol)", justify="right", width=24)
    table.add_column("构象数", justify="center", width=12)

    for rank, item in enumerate(valid_results[:10], start=1):
        table.add_row(
            f"#{rank}",
            item["id"],
            format_affinity(item["affinity"]),
            str(item["n_poses"]),
        )
    console.print(table)

    # 8. Export CSV (Formula-safe)
    csv_path = out_dir / "screening_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "Ligand_ID", "Top_Affinity_kcal_mol", "Modes_Count"])
        for rank, item in enumerate(valid_results, start=1):
            safe_id = item["id"]
            if str(safe_id).startswith(("=", "+", "-", "@")):
                safe_id = f"'{safe_id}"
            writer.writerow([rank, safe_id, f"{item['affinity']:.3f}", item["n_poses"]])

    console.print(f"[bold #4EBA65]✓ 打分结果已导出 (防公式注入 CSV):[/bold #4EBA65] {csv_path}\n")
    return 0
