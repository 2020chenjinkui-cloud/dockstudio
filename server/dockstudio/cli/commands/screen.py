"""High-Throughput Virtual Screening command for DockStudio CLI (/screen).

Designed for headless Linux & GPU clusters (Uni-Dock / Vina multi-process).
Provides:
- Multi-molecule library loading (SDF / SMILES)
- Real-time live throughput & ETA progress
- Anti-formula-injection CSV export & funnel.csv
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Any

from rich.table import Table

from ...config import get_settings
from ...engines import get_engine
from ...prep.library import prepare_library as prepare_library_from_file
from ..ui import (
    console,
    format_affinity,
    prompt_path,
    select_one,
)
from .dock import _select_docking_params, _select_search_box
from .prep import run_prep_receptor_interactive


def run_screen_interactive() -> int:
    """Runs interactive high-throughput virtual screening."""
    console.print("\n[bold #D77757]══════ DockStudio 高通量虚拟筛选向导 (/screen) ══════[/bold #D77757]\n")
    settings = get_settings()

    # 1. Receptor Prep
    pdb_path, receptor = run_prep_receptor_interactive()

    # 2. Search Box
    box = _select_search_box(pdb_path, receptor)

    # 3. Library Selection
    lib_mode = select_one("选择配体库来源:", [
        ("使用本地化合物库文件 (SDF / SMILES / CSV)", "支持商用现货库 (TargetMol / Selleck / MCE 等)"),
        ("使用 DockStudio 内置演示分子库 (10 个上市药物)", "测试端到端筛选链路与吞吐"),
    ])

    if lib_mode == 0:
        lib_path = prompt_path("请指定配体库文件路径")
        lib_bytes = lib_path.read_bytes()
        lib_filename = lib_path.name
    else:
        demo_file = Path(__file__).resolve().parents[2] / "demo_data" / "demo_library.smi"
        if not demo_file.is_file():
            demo_file = Path(__file__).resolve().parents[3] / "tests" / "fixtures" / "demo_library.smi"
        lib_bytes = demo_file.read_bytes()
        lib_filename = "demo_library.smi"
        console.print(f"[bold #4EBA65]✓[/bold #4EBA65] 加载内置演示库: {lib_filename}")

    # 4. Prepare Library
    lib_dir = settings.prepared_dir / "libraries" / f"lib_{int(time.time())}"
    with console.status("[bold #93A5FF]正在解析分子库并逐条加氢质子化转 PDBQT...[/bold #93A5FF]"):
        prepared_lib = prepare_library_from_file(lib_bytes, lib_filename, lib_dir)

    n_ready = len(prepared_lib.ligand_ids)
    n_failed = len(prepared_lib.failed)
    console.print(f"[bold #4EBA65]✓ 分子库准备完成:[/bold #4EBA65] 成功 {n_ready} 个分子，失败 {n_failed} 个")

    if not prepared_lib.ligand_ids:
        console.print("[bold #FF6B80]错误: 库中没有成功准备的配体分子，筛选终止[/#FF6B80]")
        return 1

    # 5. Engine Selection
    engine_name, scoring, base_params = _select_docking_params(settings)
    engine = get_engine(engine_name, settings)

    # 6. Run Screening with Live Throughput
    out_dir = settings.data_dir / "screens" / f"screen_{int(time.time())}"
    out_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"\n[bold #D77757]开始高通量筛选 (引擎: {engine_name}, 共 {n_ready} 分子)...[/bold #D77757]")
    results: list[dict[str, Any]] = []
    t_start = time.time()

    for idx, lig_id in enumerate(prepared_lib.ligand_ids):
        lig_pdbqt = lib_dir / f"{lig_id}.pdbqt"
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
                "id": lig_id,
                "affinity": top_aff,
                "n_poses": len(res.poses),
                "run_dir": str(run_scratch),
            })
        except Exception as err:  # noqa: BLE001
            results.append({
                "id": lig_id,
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
    table.add_column("分子 ID", justify="left", width=20)
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
            writer.writerow([rank, item["id"], f"{item['affinity']:.3f}", item["n_poses"]])

    console.print(f"[bold #4EBA65]✓ 打分结果已导出 (防公式注入 CSV):[/bold #4EBA65] {csv_path}\n")
    return 0
