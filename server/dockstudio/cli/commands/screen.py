"""High-Throughput Virtual Screening command for DockStudio CLI (/screen).

Provides:
- 5-Step interactive wizard with Claude Code Open-ClaudeCode design
- Dual-mode execution (Interactive Terminal Wizard & Non-Interactive / Headless CLI)
- Multi-tier live Rich progress display (ScreenProgressManager)
- Multi-dimensional consensus leaderboard (render_screening_leaderboard)
- TCM dual-track deglycosylation and side-by-side aglycone comparison
- Key residue pharmacophore anchoring
- Full delivery package assembly for /verify T1~T5 compliance
"""

from __future__ import annotations

import csv
import io
import math
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import click
from prompt_toolkit import prompt
from prompt_toolkit.formatted_text import HTML
from rich.box import ROUNDED
from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TaskProgressColumn,
    TextColumn,
    TimeRemainingColumn,
)
from rich.table import Table
from rich.text import Text

from ...analysis.pockets import detect_pockets
from ...config import get_settings
from ...core.orchestrator import PrepareRegistry
from ...core.screen_runner import ScreeningProgressEvent, ScreeningRunner
from ...engines import available_engines, get_engine
from ...models import (
    Box,
    DockParams,
    PreparedLibrary,
    ScreeningCandidateRow,
    ScreeningConfig,
    ScreeningFunnelReport,
)
from ...prep.library import prepare_library
from ...prep.receptor import prepare_receptor
from ..ui import (
    console,
    format_affinity,
    make_header,
    make_pane,
)
from ..ui import confirm as _ui_confirm
from ..ui import prompt_path as _ui_prompt_path
from ..ui import select_one as _ui_select_one

_DEMO_LIBRARY_SMI = """\
O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O Baicalin
O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12 Quercetin
O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1 Curcumin
CC(C)Cc1ccc(cc1)C(C)C(=O)O Ibuprofen
CC(=O)Oc1ccccc1C(=O)O Aspirin
"""


# ---------------------------------------------------------------------------
# Cross-Platform Safe Input Handlers (supporting both TTY and CliRunner)
# ---------------------------------------------------------------------------

def select_one(
    title: str,
    choices: Any,
    default: int = 0,
) -> int:
    """Selects one option; delegates to interactive UI when TTY or uses stream input in tests."""
    if sys.stdin.isatty():
        try:
            return _ui_select_one(title, choices, default=default)
        except Exception:
            pass

    options: list[str] = [(c[0] if isinstance(c, tuple) else str(c)) for c in choices]
    content = Text()
    for idx, label in enumerate(options):
        is_default = idx == default
        marker = "❯" if is_default else " "
        idx_str = f"[{idx + 1}]"
        if is_default:
            content.append(f" {marker} {idx_str:4s} {label}\n", style="bold #D77757")
        else:
            content.append(f"   {idx_str:4s} {label}\n", style="#CCCCCC")
    console.print(make_pane(title, content, subtitle=f"默认: [{default + 1}]"))

    while True:
        try:
            choice_str = input(f"请选择 [1-{len(options)}] (回车默认={default + 1}): ").strip()
        except (EOFError, OSError):
            return default
        except KeyboardInterrupt:
            raise
        if not choice_str:
            return default
        try:
            num = int(choice_str)
            if 1 <= num <= len(options):
                return num - 1
            console.print(f"[#FF6B80]输入超出范围，请输入 1 到 {len(options)}[/]")
        except ValueError:
            console.print(f"[#FF6B80]无效输入 '{choice_str}'，请输入数字[/]")


def prompt_path(
    message: str = "请输入文件路径",
    default: str = "",
    must_exist: bool = True,
) -> Path:
    """Prompts user for path; uses interactive UI when TTY or input() in tests."""
    if sys.stdin.isatty():
        try:
            return _ui_prompt_path(message=message, default=default, must_exist=must_exist)
        except Exception:
            pass

    while True:
        try:
            val = input(f"{message} [{default}]: ").strip()
        except (EOFError, OSError):
            if default:
                val = default
            else:
                raise
        except KeyboardInterrupt:
            raise
        if not val:
            if default:
                val = default
            else:
                continue
        p = Path(val).expanduser().resolve()
        if must_exist and not p.exists():
            console.print(f"[#FF6B80]文件不存在: {p}，请重新输入[/]")
            continue
        return p


def confirm(message: str, default: bool = True) -> bool:
    """Prompts for Y/n confirmation; safe across both TTY and test streams."""
    if sys.stdin.isatty():
        try:
            return _ui_confirm(message, default=default)
        except Exception:
            pass

    suffix = "[Y/n]" if default else "[y/N]"
    try:
        ans = input(f"{message} {suffix}: ").strip().lower()
        if not ans:
            return default
        return ans in ("y", "yes", "true", "1")
    except (EOFError, OSError):
        return default
    except KeyboardInterrupt:
        raise


def _prompt_text(message: str, default: str = "") -> str:
    """Prompts for a string input with default value."""
    if sys.stdin.isatty():
        try:
            return prompt(HTML(f"<ansicyan>{message}:</ansicyan> "), default=default).strip()
        except (KeyboardInterrupt, EOFError):
            raise
        except Exception:
            pass
    try:
        val = input(f"{message} [{default}]: ").strip()
        return val if val else default
    except (EOFError, OSError):
        return default
    except KeyboardInterrupt:
        raise


def _select_docking_params(settings: Any) -> tuple[str, str, DockParams]:
    """Select engine and parameters (for backward compatibility and test monkeypatching)."""
    return "vina", "vina", DockParams(exhaustiveness=8, num_modes=9, cpu=2, seed=42, scoring="vina")


# ---------------------------------------------------------------------------
# Stepper & Header Visuals
# ---------------------------------------------------------------------------

def _render_screen_stepper(step: int, total_steps: int = 5, step_title: str = "") -> Any:
    """Renders a Claude Code style linear wizard stepper for /screen."""
    steps = [
        ("1", "受体定盒"),
        ("2", "中药双态"),
        ("3", "药效锚定"),
        ("4", "三级漏斗"),
        ("5", "引擎并发"),
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


def _resolve_demo_files() -> tuple[Path, Path]:
    """Locates built-in demo receptor and ligand library fixtures."""
    base_dir = Path(__file__).resolve().parents[2]  # dockstudio/
    fixtures_dir = Path(__file__).resolve().parents[3] / "tests" / "fixtures"

    demo_rec = base_dir / "demo_data" / "1HSG.pdb"
    if not demo_rec.is_file():
        demo_rec = fixtures_dir / "1HSG.pdb"

    demo_lib = base_dir / "demo_data" / "demo_library.smi"
    if not demo_lib.is_file():
        demo_lib = fixtures_dir / "demo_library.smi"

    return demo_rec, demo_lib


# ---------------------------------------------------------------------------
# Live Progress Manager
# ---------------------------------------------------------------------------

class ScreenProgressManager:
    """Manages 3-tier cascade progress display for HTVS via rich.progress.Progress."""

    def __init__(self, console_obj: Console | None = None, is_interactive: bool = True):
        self.console = console_obj or console
        self.is_interactive = is_interactive and sys.stdout.isatty()
        self.progress = Progress(
            SpinnerColumn(spinner_name="dots", style="bold #93A5FF"),
            TextColumn("[bold #FFFFFF]{task.description}"),
            BarColumn(bar_width=24, complete_style="#D77757", finished_style="#4EBA65"),
            MofNCompleteColumn(),
            TaskProgressColumn(),
            TimeRemainingColumn(),
            TextColumn("[dim]({task.fields[speed]:.1f} mol/s)[/dim]"),
            TextColumn("[#93A5FF]{task.fields[current_mol]}[/#93A5FF]"),
            console=self.console,
            transient=False,
            disable=not self.is_interactive,
        )
        self.tasks: dict[int, TaskID] = {}

    def __enter__(self) -> ScreenProgressManager:
        if self.is_interactive:
            self.progress.start()
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if self.is_interactive:
            self.progress.stop()

    def callback(self, event: ScreeningProgressEvent) -> None:
        """Callback invoked by ScreeningRunner._emit_progress."""
        if not self.is_interactive:
            if event.current == 1 or event.current == event.total or event.current % max(1, event.total // 10) == 0:
                pct = int(event.current * 100 // max(1, event.total))
                self.console.print(
                    f"  [Tier {event.tier}] {event.stage_name}: {event.current}/{event.total} "
                    f"({pct}%) · {event.speed_mol_per_sec:.1f} mol/s"
                )
            return

        tier = event.tier
        if tier not in self.tasks:
            for prev_tier, prev_tid in self.tasks.items():
                if prev_tier < tier:
                    t_task = self.progress.tasks[prev_tid]
                    self.progress.update(prev_tid, completed=t_task.total)

            tid = self.progress.add_task(
                f"Tier {tier}: {event.stage_name}",
                total=event.total,
                speed=event.speed_mol_per_sec,
                current_mol=event.current_ligand_id,
            )
            self.tasks[tier] = tid
        else:
            tid = self.tasks[tier]
            self.progress.update(
                tid,
                completed=event.current,
                total=event.total,
                speed=event.speed_mol_per_sec,
                current_mol=event.current_ligand_id,
            )


# ---------------------------------------------------------------------------
# Post-Screening Rich Summary Leaderboard
# ---------------------------------------------------------------------------

def render_screening_leaderboard(
    candidates: list[ScreeningCandidateRow],
    n_raw: int,
    n_ready: int,
    top_k: int = 20,
    out_dir: Path | None = None,
) -> Panel:
    """Renders the comprehensive multi-dimensional consensus ranking table."""
    table = Table(
        box=ROUNDED,
        header_style="bold #FFFFFF",
        border_style="#888888",
        expand=True,
    )

    table.add_column("Rank", justify="center", width=7, style="bold #D77757")
    table.add_column("配体名称 (ID)", justify="left", width=22)
    table.add_column("形态 (Derivative)", justify="center", width=14)
    table.add_column("ΔG ± σ (kcal/mol)", justify="right", width=18)
    table.add_column("SILE", justify="right", width=8)
    table.add_column("ΔE_strain", justify="right", width=10)
    table.add_column("ECR共识 (ECR)", justify="right", width=12)
    table.add_column("骨架簇 (Cluster)", justify="center", width=12)
    table.add_column("警示徽章 (Badges)", justify="left", width=22)

    display_candidates = candidates[:top_k]
    for rank_idx, c in enumerate(display_candidates, start=1):
        if rank_idx == 1:
            rank_str = "🥇 #1"
        elif rank_idx == 2:
            rank_str = "🥈 #2"
        elif rank_idx == 3:
            rank_str = "🥉 #3"
        else:
            rank_str = f"#{rank_idx}"

        if c.ligand_name and c.ligand_name != c.ligand_id:
            name_id = f"{c.ligand_name[:14]} [dim]({c.ligand_id[:10]})[/dim]"
        else:
            name_id = f"{c.ligand_id[:20]}"

        if c.derivative_type == "aglycone":
            type_badge = "[#93A5FF]衍生苷元[/#93A5FF]"
        elif c.derivative_type == "stereoisomer":
            type_badge = "[#FFC107]异构体[/#FFC107]"
        else:
            type_badge = "[#4EBA65]原形单体[/#4EBA65]"

        if c.affinity_std > 0.0:
            aff_text = f"{c.affinity_mean:6.2f} ± {c.affinity_std:4.2f}"
        else:
            aff_text = f"{c.affinity_mean:6.2f}"

        if c.affinity_mean <= -9.0:
            aff_styled = f"[bold #4EBA65]{aff_text}[/bold #4EBA65]"
        elif c.affinity_mean <= -7.0:
            aff_styled = f"[bold #FFC107]{aff_text}[/bold #FFC107]"
        else:
            aff_styled = f"[#FF6B80]{aff_text}[/#FF6B80]"

        if c.sile >= 2.0:
            sile_styled = f"[bold #4EBA65]{c.sile:.2f}[/bold #4EBA65]"
        elif c.sile >= 1.5:
            sile_styled = f"[#FFFFFF]{c.sile:.2f}[/#FFFFFF]"
        else:
            sile_styled = f"[dim]{c.sile:.2f}[/dim]"

        if c.strain_energy > 6.0:
            strain_styled = f"[bold #FF6B80]⚠️ {c.strain_energy:.1f}[/bold #FF6B80]"
        elif c.strain_energy <= 3.0:
            strain_styled = f"[#4EBA65]{c.strain_energy:.1f}[/#4EBA65]"
        else:
            strain_styled = f"[#FFFFFF]{c.strain_energy:.1f}[/#FFFFFF]"

        ecr_styled = f"[bold]{c.ecr_score:.3f}[/]"

        if c.is_cluster_rep:
            cluster_str = f"#{c.cluster_id} [bold #93A5FF]★[/]"
        else:
            cluster_str = f"#{c.cluster_id}"

        if not c.alerts:
            badge_str = "[bold #4EBA65]● Clean[/bold #4EBA65]"
        else:
            badge_items = []
            for alert in c.alerts:
                al_lower = alert.lower()
                if "catechol" in al_lower:
                    badge_items.append("[bold #FFC107]Catechol[/]")
                elif "quinone" in al_lower:
                    badge_items.append("[bold #FF6B80]Quinone[/]")
                elif "michael" in al_lower:
                    badge_items.append("[bold #FF6B80]Michael[/]")
                elif "tannin" in al_lower:
                    badge_items.append("[bold #FFC107]Tannin[/]")
                elif "pains" in al_lower or "brenk" in al_lower:
                    badge_items.append("[bold #FF6B80]PAINS[/]")
                elif "posecheck" in al_lower:
                    # Real pose-geometry QC diagnostics (P1-5) instead of a hard-coded
                    # "always pass" verdict; without this branch the 8-char fallback
                    # would render the useless "POSECHEC".
                    badge_items.append("[bold #FF6B80]PoseQC[/]")
                elif "high_strain" in al_lower or "unconverged" in al_lower:
                    badge_items.append("[bold #FF6B80]QC-降级[/]")
                else:
                    badge_items.append(f"[#FFC107]{alert[:8]}[/]")
            badge_str = " ".join(badge_items)

        table.add_row(
            rank_str,
            name_id,
            type_badge,
            aff_styled,
            sile_styled,
            strain_styled,
            ecr_styled,
            cluster_str,
            badge_str,
        )

    subtitle = (
        f"展示 Top {min(top_k, len(candidates))} / 最终入围 {len(candidates)} "
        f"(筛选就绪 {n_ready}，原始输入 {n_raw})"
    )
    pane = make_pane(
        "DockStudio HTVS 综合共识排行榜 (Consensus Leaderboard)",
        table,
        subtitle=subtitle,
    )
    return pane


# ---------------------------------------------------------------------------
# Coordinate & Option Helpers
# ---------------------------------------------------------------------------

def parse_box_coords(box_center_str: str | None, box_size_str: str) -> Box:
    """Parses center and size coordinates from CLI string flags."""
    if box_center_str:
        clean_center = box_center_str.replace(";", ",").replace("，", ",")
        parts = [float(v.strip()) for v in clean_center.split(",") if v.strip()]
        if len(parts) != 3:
            parts = [float(v.strip()) for v in box_center_str.split() if v.strip()]
        if len(parts) != 3:
            raise ValueError(
                f"box-center 必须为 3 个由逗号分隔的浮点数，如 '15.2,24.8,-5.1'，当前为: '{box_center_str}'"
            )
        cx, cy, cz = parts
    else:
        cx, cy, cz = (16.0, 25.0, 4.0)

    clean_size = box_size_str.replace(";", ",").replace("，", ",")
    size_parts = [float(v.strip()) for v in clean_size.split(",") if v.strip()]
    if len(size_parts) == 1:
        sx = sy = sz = size_parts[0]
    elif len(size_parts) == 3:
        sx, sy, sz = size_parts
    else:
        size_parts = [float(v.strip()) for v in box_size_str.split() if v.strip()]
        if len(size_parts) == 1:
            sx = sy = sz = size_parts[0]
        elif len(size_parts) == 3:
            sx, sy, sz = size_parts
        else:
            sx, sy, sz = (22.5, 22.5, 22.5)

    return Box(center=(cx, cy, cz), size=(sx, sy, sz))


def _resolve_receptor_files(receptor_path: Path, settings: Any) -> tuple[Path, Path | None, str, Any | None]:
    """Resolves receptor PDBQT and PDB files, calling prepare_receptor if necessary."""
    rec_out = settings.prepared_dir / "receptors"
    rec_out.mkdir(parents=True, exist_ok=True)

    if receptor_path.suffix.lower() == ".pdbqt":
        rec_pdbqt = receptor_path
        rec_dir = receptor_path.parent
        if (rec_dir / "receptor_protein.pdb").is_file():
            rec_pdb = rec_dir / "receptor_protein.pdb"
        elif (rec_dir / "receptor_clean.pdb").is_file():
            rec_pdb = rec_dir / "receptor_clean.pdb"
        elif receptor_path.with_suffix(".pdb").is_file():
            rec_pdb = receptor_path.with_suffix(".pdb")
        else:
            rec_pdb = None
        rec_id = receptor_path.stem
        rec_obj = None
    elif receptor_path.with_suffix(".pdbqt").is_file():
        rec_pdbqt = receptor_path.with_suffix(".pdbqt")
        rec_dir = rec_pdbqt.parent
        if (rec_dir / "receptor_protein.pdb").is_file():
            rec_pdb = rec_dir / "receptor_protein.pdb"
        elif (rec_dir / "receptor_clean.pdb").is_file():
            rec_pdb = rec_dir / "receptor_clean.pdb"
        else:
            rec_pdb = receptor_path
        rec_id = receptor_path.stem
        rec_obj = None
    else:
        rec_obj = prepare_receptor(receptor_path.read_bytes(), receptor_path.name, rec_out)
        rec_pdbqt = Path(rec_obj.pdbqt_path)
        rec_dir = rec_pdbqt.parent
        if (rec_dir / "receptor_protein.pdb").is_file():
            rec_pdb = rec_dir / "receptor_protein.pdb"
        elif (rec_dir / "receptor_clean.pdb").is_file():
            rec_pdb = rec_dir / "receptor_clean.pdb"
        else:
            rec_pdb = receptor_path
        rec_id = rec_obj.id

    return rec_pdbqt, rec_pdb, rec_id, rec_obj


# ---------------------------------------------------------------------------
# Interactive 5-Step Wizard Execution
# ---------------------------------------------------------------------------

def run_screen_interactive(demo: bool = False, **overrides: Any) -> int:
    """Runs the 5-step interactive HTVS wizard with Claude Code style panels."""
    settings = get_settings()
    registry = PrepareRegistry(settings)

    try:
        # Step 1: Receptor Selection & Pocket Box Definition
        _refresh_screen_step(1, 5, "受体结构加载与结合口袋定盒")

        if demo:
            demo_rec, demo_lib = _resolve_demo_files()
            rec_path = demo_rec
            rec_pdbqt, rec_pdb, rec_id, rec_obj = _resolve_receptor_files(rec_path, settings)
            console.print(f"[bold #4EBA65]✓[/bold #4EBA65] 加载内置演示受体: 1HSG (HIV-1 蛋白酶催化二聚体)")
            dock_box = Box(center=(16.0, 25.0, 4.0), size=(22.5, 22.5, 22.5))
            console.print(
                f"[bold #4EBA65]✓ 搜索盒已设定:[/bold #4EBA65] 中心 ({dock_box.center[0]:.1f}, {dock_box.center[1]:.1f}, {dock_box.center[2]:.1f}) · 尺寸 {dock_box.size[0]:.1f}Å"
            )
        else:
            rec_mode = select_one("选择受体来源:", [
                ("指定受体文件路径 (PDB / mmCIF / PDBQT)", "从本地文件系统加载"),
                ("使用内置演示受体: 1HSG", "HIV-1 蛋白酶催化二聚体"),
            ], default=0)

            if rec_mode == 0:
                rec_path = prompt_path("请指定受体结构文件路径")
            else:
                demo_rec, _ = _resolve_demo_files()
                rec_path = demo_rec

            rec_pdbqt, rec_pdb, rec_id, rec_obj = _resolve_receptor_files(rec_path, settings)
            console.print(f"[bold #4EBA65]✓ 受体就绪:[/bold #4EBA65] {rec_id} ({rec_pdbqt.name})")

            # Pocket Box Selection
            box_choices = [
                ("围绕受体催化位点自动外包盒", "自动计算结合位点质心构建 22.5Å 搜索盒"),
                ("运行腔体预测算法 (fpocket / 几何网格)", "在受体全表面自动寻找潜在活性口袋"),
                ("手动指定中心坐标与尺寸", "手动输入 Center X/Y/Z 与 Size"),
            ]
            box_mode = select_one("选择结合位点与搜索盒 (Search Box):", box_choices, default=0)

            if box_mode == 0:
                if rec_obj and getattr(rec_obj, "cocrystal_ligands", None):
                    c = rec_obj.cocrystal_ligands[0]
                    dock_box = Box(center=c.center, size=(22.5, 22.5, 22.5))
                else:
                    dock_box = Box(center=(16.0, 25.0, 4.0), size=(22.5, 22.5, 22.5))
                console.print(
                    f"[bold #4EBA65]✓ 搜索盒已设定:[/bold #4EBA65] 中心 ({dock_box.center[0]:.1f}, {dock_box.center[1]:.1f}, {dock_box.center[2]:.1f}) · 尺寸 {dock_box.size[0]:.1f}Å"
                )
            elif box_mode == 1:
                with console.status("[bold #93A5FF]正在分析受体表面几何腔体...[/bold #93A5FF]"):
                    pockets = detect_pockets(rec_pdb or rec_pdbqt, top_n=3)
                if pockets:
                    p_choices = [
                        (f"口袋 #{p.id} (评分: {p.score:.1f})", f"中心: ({p.center[0]:.1f}, {p.center[1]:.1f}, {p.center[2]:.1f}) · 尺寸: ({p.size[0]:.1f}, {p.size[1]:.1f}, {p.size[2]:.1f})")
                        for p in pockets
                    ]
                    p_idx = select_one("选择预测口袋:", p_choices, default=0)
                    dock_box = Box(center=pockets[p_idx].center, size=pockets[p_idx].size)
                else:
                    dock_box = Box(center=(16.0, 25.0, 4.0), size=(22.5, 22.5, 22.5))
            else:
                console.print("[dim]请输入盒子中心与尺寸 (以回车确认):[/dim]")
                try:
                    cx_str = _prompt_text("Center X", default="10.0")
                    cx = float(cx_str) if cx_str else 10.0
                except (ValueError, Exception):
                    cx = 10.0
                try:
                    cy_str = _prompt_text("Center Y", default="10.0")
                    cy = float(cy_str) if cy_str else 10.0
                except (ValueError, Exception):
                    cy = 10.0
                try:
                    cz_str = _prompt_text("Center Z", default="10.0")
                    cz = float(cz_str) if cz_str else 10.0
                except (ValueError, Exception):
                    cz = 10.0
                try:
                    size_str = _prompt_text("Box Size (Å)", default="20.0")
                    size_val = float(size_str) if size_str else 20.0
                except (ValueError, Exception):
                    size_val = 20.0
                dock_box = Box(center=(cx, cy, cz), size=(size_val, size_val, size_val))

        # Step 2: Library Selection & TCM Dual-Track Mode
        _refresh_screen_step(2, 5, "化合物库加载与中药单体双轨衍生")

        if demo:
            lib_bytes = _DEMO_LIBRARY_SMI.encode("utf-8")
            lib_name = "demo_library.smi"
            tcm_mode = True
            console.print("[bold #4EBA65]✓[/bold #4EBA65] 加载内置中药单体及药物库 (黄芩苷、槲皮素、姜黄素、布洛芬、阿司匹林)")
        else:
            lib_choices = [
                ("指定化合物库文件路径 (SMILES / SDF / CSV)", "从本地文件系统加载化合物库"),
                ("使用内置演示化合物库 (含糖苷单体与已上市药物)", "测试端到端中药筛选双轨比对链路"),
            ]
            lib_mode = select_one("选择配体库来源:", lib_choices, default=0)
            if lib_mode == 0:
                lib_path = prompt_path("请指定化合物库文件路径 (SMILES / SDF / CSV)")
                lib_bytes = lib_path.read_bytes()
                lib_name = lib_path.name
            else:
                lib_bytes = _DEMO_LIBRARY_SMI.encode("utf-8")
                lib_name = "demo_library.smi"

            tcm_info = Text()
            tcm_info.append("💡 中药单体双轨比对机制 (TCM Dual-Track Screening):\n", style="bold #93A5FF")
            tcm_info.append("• 糖苷假阳性消除: 自动水解可切断糖苷键衍生苷元，解决原形糖链非特异氢键假阳性。\n", style="#CCCCCC")
            tcm_info.append("• 双轨同屏比对: 原形糖苷与衍生苷元同时进入筛选池，排行榜同屏对比活性。\n", style="#CCCCCC")
            tcm_info.append("• 碳苷精密保护: 葛根素、牡荆素等 C-C 键糖苷结构得到精确保护不断裂。", style="#CCCCCC")
            console.print(make_pane("中药单体预处理特性", tcm_info))

            tcm_mode = confirm("是否启用中药单体双轨模式 (自动水解衍生苷元同屏比对)?", default=True)

        with console.status("[bold #93A5FF]正在解析分子库、执行脱盐水解与 PDBQT 准备...[/bold #93A5FF]"):
            prepared_lib = prepare_library(
                content=lib_bytes,
                filename=lib_name,
                registry=registry,
                tcm_mode=tcm_mode,
            )

        n_ready = len(prepared_lib.ligand_ids)
        console.print(
            f"[bold #4EBA65]✓ 配体库就绪:[/bold #4EBA65] 原始分子 {prepared_lib.n_raw} -> 筛选就绪 {n_ready} (原形+苷元双轨) "
            f"· 损耗(脱盐{prepared_lib.n_desalt_fail}, 去重{prepared_lib.n_dedup}, 异常{prepared_lib.n_prep_fail + prepared_lib.n_chiral_fail})"
        )

        # Step 3: Key Residue Pharmacophore Anchoring
        _refresh_screen_step(3, 5, "关键残基相互作用药效团锚定")

        if demo:
            parsed_anchors: list[str] = []
            console.print("[dim]内置演示模式：未设定锚定残基，执行全口袋自由对接[/dim]")
        else:
            anchor_choices = [
                ("不设锚定残基 (标准全口袋对接)", "允许配体在整个搜索盒内部自由结合"),
                ("手动输入关键残基代码", "输入活性中心关键催化残基 (例如: ASP25:A,THR80:A)"),
            ]
            anchor_mode = select_one("选择锚定残基设定方式:", anchor_choices, default=0)
            parsed_anchors = []
            if anchor_mode == 1:
                raw_anchors = _prompt_text("请输入锚定残基代码 (逗号分隔，如 ASP25:A,THR80:A)", default="ASP25:A,THR80:A")
                parsed_anchors = [a.strip().upper() for a in raw_anchors.replace(";", ",").split(",") if a.strip()]
                console.print(f"[bold #4EBA65]✓ 药效团锚定残基:[/bold #4EBA65] {', '.join(parsed_anchors)}")
            else:
                console.print("[dim]未设定锚定残基，执行全口袋自由对接[/dim]")

        # Step 4: 3-Tier Hierarchical Funnel Configuration
        _refresh_screen_step(4, 5, "三级级联漏斗参数与严谨性打分配置")

        custom_engine = "vina"
        custom_scoring = "vina"
        custom_params = None
        try:
            custom_engine, custom_scoring, custom_params = _select_docking_params(settings)
        except Exception:
            pass

        if demo:
            tier1_ratio = overrides.get("tier1_ratio", 0.5)
            tier2_ratio = overrides.get("tier2_ratio", 0.5)
            parsed_seeds = [custom_params.seed] if (custom_params and custom_params.seed is not None) else [42]
            strain_threshold = overrides.get("strain_threshold", 8.0)
            enable_gnina = False
            console.print("[bold #4EBA65]✓ 演示模式漏斗:[/bold #4EBA65] Tier1 50% -> Tier2 单种子 -> 快速完成")
        else:
            preset_choices = [
                ("标准三级漏斗 (推荐)", "Tier1 15% SILE -> Tier2 3-Seed 收敛性/张力能 (2%) -> Tier3 GNINA/ECR 聚类"),
                ("极速高通量模式 (Fast)", "Tier1 10% -> Tier2 2-Seed (1%) -> Tier3 ECR (跳过 GNINA 重打分)"),
                ("严谨药物发现模式 (Strict)", "Tier1 20% -> Tier2 3-Seed (5%, 张力能 <= 6.0) -> Tier3 GNINA/ECR"),
                ("自定义高级参数 (Custom)", "手动指定各级漏斗晋级比例、随机种子数与张力能阈值"),
            ]
            preset_idx = select_one("选择漏斗配置方案:", preset_choices, default=0)

            if preset_idx == 0:
                tier1_ratio = 0.15
                tier2_ratio = 0.02
                parsed_seeds = [42, 1024, 2026]
                strain_threshold = 8.0
                enable_gnina = True
            elif preset_idx == 1:
                tier1_ratio = 0.10
                tier2_ratio = 0.01
                parsed_seeds = [42, 1024]
                strain_threshold = 10.0
                enable_gnina = False
            elif preset_idx == 2:
                tier1_ratio = 0.20
                tier2_ratio = 0.05
                parsed_seeds = [42, 1024, 2026]
                strain_threshold = 6.0
                enable_gnina = True
            else:
                t1_str = _prompt_text("Tier 1 粗筛晋级比例 (0.01~1.0)", default="0.15")
                tier1_ratio = float(t1_str) if t1_str else 0.15
                t2_str = _prompt_text("Tier 2 精筛晋级比例 (0.001~0.5)", default="0.02")
                tier2_ratio = float(t2_str) if t2_str else 0.02
                seeds_str = _prompt_text("Tier 2 随机种子列表 (逗号分隔)", default="42,1024,2026")
                parsed_seeds = [int(s.strip()) for s in seeds_str.split(",") if s.strip()] or [42, 1024, 2026]
                strain_str = _prompt_text("配体张力能上限 (kcal/mol)", default="8.0")
                strain_threshold = float(strain_str) if strain_str else 8.0
                enable_gnina = confirm("Tier 3 是否启用 GNINA CNN 正交重打分?", default=True)

        # Step 5: Engine, Concurrency & Checkpoint Resumption
        _refresh_screen_step(5, 5, "计算引擎、多核并发分配与筛选启动")

        if demo:
            engine_name = custom_engine if custom_engine else "vina"
            resolved_engine = _resolve_engine(engine_name)
            out_dir = overrides.get("out") or (settings.data_dir / "screens" / f"screen_demo_{int(time.time())}")
            out_dir = Path(out_dir).expanduser().resolve()
            out_dir.mkdir(parents=True, exist_ok=True)
            resume = False
            total_cpus = overrides.get("cpu") or (custom_params.cpu if custom_params and custom_params.cpu else 2)
            cpu_per_job = overrides.get("cpu_per_job") or (custom_params.cpu if custom_params and custom_params.cpu else 2)
            console.print(f"[bold #4EBA65]✓ 演示模式引擎:[/bold #4EBA65] {resolved_engine.manifest.name} (CPU {total_cpus}) -> 立即启动")
        else:
            eng_choices = [
                ("自动推荐 (Auto: GPU 可用首选 Uni-Dock，否则 Vina)", "根据系统硬件自适应选择最高吞吐后端"),
                ("Vina (CPU 多进程)", "跨平台通用标准对接引擎"),
                ("Uni-Dock (GPU 批处理加速)", "NVIDIA CUDA 批量高通量加速"),
                ("GNINA (深度学习打分)", "含 CNN 打分模型的高精度引擎"),
            ]
            eng_choice_idx = select_one("选择对接引擎:", eng_choices, default=0)

            if eng_choice_idx == 1:
                engine_name = "vina"
            elif eng_choice_idx == 2:
                engine_name = "unidock"
            elif eng_choice_idx == 3:
                engine_name = "gnina"
            else:
                engine_name = "auto"

            resolved_engine = _resolve_engine(engine_name)

            out_dir = overrides.get("out") or (settings.data_dir / "screens" / f"screen_{int(time.time())}")
            out_dir = Path(out_dir).expanduser().resolve()
            out_dir.mkdir(parents=True, exist_ok=True)

            ckpt_file = out_dir / "checkpoint.jsonl"
            resume = False
            if ckpt_file.is_file() and ckpt_file.stat().st_size > 0:
                resume = confirm("检测到已有断点文件 checkpoint.jsonl，是否从断点继续筛选 (跳过已完成分子)?", default=True)

            total_cpus = overrides.get("cpu") or max(2, min(os.cpu_count() or 4, 16))
            cpu_per_job = overrides.get("cpu_per_job") or 2

            preflight = Text()
            preflight.append(f"• 受体结构: {rec_id} · 搜索盒中心 ({dock_box.center[0]:.1f}, {dock_box.center[1]:.1f}, {dock_box.center[2]:.1f}) · 尺寸 {dock_box.size[0]:.1f}Å\n", style="#CCCCCC")
            preflight.append(f"• 配体分子: {prepared_lib.name} · 就绪分子 {n_ready} 个 (中药双态: {tcm_mode})\n", style="#CCCCCC")
            preflight.append(f"• 药效锚定: {', '.join(parsed_anchors) if parsed_anchors else '全口袋自由结合'}\n", style="#CCCCCC")
            preflight.append(f"• 漏斗配置: Tier1 留存 {int(tier1_ratio*100)}% -> Tier2 留存 {int(tier2_ratio*100)}% (Seeds: {parsed_seeds}, 张力能 <= {strain_threshold} kcal/mol)\n", style="#CCCCCC")
            preflight.append(f"• 计算后端: 引擎 {resolved_engine.manifest.name} · 总 CPU {total_cpus} (每任务 {cpu_per_job} 线程)\n", style="#CCCCCC")
            preflight.append(f"• 交付路径: {out_dir} (断点续筛: {resume})", style="#CCCCCC")
            console.print(make_pane("筛选启动预检确认 (Pre-flight Summary)", preflight))

            start_confirm = confirm("确认以上配置，立即启动高通量筛选?", default=True)
            if not start_confirm:
                console.print("[dim]用户取消启动[/dim]")
                return 0

        config = ScreeningConfig(
            receptor_id=rec_id,
            library_id=prepared_lib.name or "library",
            box=dock_box,
            tcm_mode=tcm_mode,
            tcm_deglycosylation=tcm_mode,
            enable_tiered_funnel=True,
            tier1_keep_ratio=tier1_ratio,
            tier1_max_candidates=2 if demo else 1000,
            tier2_seeds=parsed_seeds,
            tier2_keep_ratio=tier2_ratio,
            tier2_max_candidates=2 if demo else 100,
            max_strain_kcal=strain_threshold,
            anchor_residues=parsed_anchors,
            enable_gnina_rescore=enable_gnina,
            engine=resolved_engine.manifest.name,
            cpu=total_cpus,
        )

        lib_ligands = {
            lid: registry.ligands[lid]
            for lid in prepared_lib.ligand_ids
            if lid in registry.ligands
        }

        with ScreenProgressManager(console, is_interactive=True) as progress_mgr:
            runner = ScreeningRunner(
                config=config,
                receptor_pdbqt=rec_pdbqt,
                protein_pdb=rec_pdb,
                ligands=lib_ligands,
                prepared_library=prepared_lib,
                out_dir=out_dir,
                engine=resolved_engine,
                resume=resume,
                cpu_per_job=cpu_per_job,
                callback=progress_mgr.callback,
            )
            report = runner.run()

        console.print()
        leaderboard_panel = render_screening_leaderboard(
            report.candidates,
            report.n_raw,
            report.n_ready,
            top_k=20,
            out_dir=out_dir,
        )
        console.print(leaderboard_panel)

        console.print(f"\n[bold #4EBA65]✓ 高通量筛选全部完成![/bold #4EBA65] 产出交付包目录: {out_dir}")
        console.print("[dim]交付文件: funnel.csv · screening_results.csv · conf.txt · receptor.pdbqt · poses/[/dim]\n")
        return 0

    except (KeyboardInterrupt, EOFError):
        console.print("\n[bold #FFC107]操作已取消 (User Cancelled)[/bold #FFC107]")
        return 130


# ---------------------------------------------------------------------------
# Headless / Non-Interactive Execution
# ---------------------------------------------------------------------------

def _resolve_engine(engine_name: str) -> Any:
    """Resolves engine instance, with auto-fallback to unidock GPU or vina."""
    eng_str = (engine_name or "auto").lower()
    if eng_str == "auto":
        unidock = get_engine("unidock")
        if unidock is not None and unidock.is_available()[0]:
            return unidock
        vina = get_engine("vina")
        if vina is not None:
            return vina
        engs = available_engines()
        return engs[0] if engs else get_engine("vina")
    eng = get_engine(eng_str)
    if eng is not None:
        return eng
    return get_engine("vina")


def run_screen_headless(
    receptor_path: Path | None,
    library_path: Path | None,
    box_center_str: str | None,
    box_size_str: str,
    tcm_mode: bool,
    anchor_residues_str: str,
    tier1_ratio: float,
    tier2_ratio: float,
    seeds_str: str,
    strain_threshold: float,
    engine_name: str,
    cpu: int | None,
    cpu_per_job: int,
    resume: bool | None,
    out_dir: Path | None,
    enable_gnina_rescore: bool,
    demo: bool,
) -> int:
    """Executes high-throughput screening in non-interactive headless mode."""
    settings = get_settings()
    registry = PrepareRegistry(settings)

    if demo:
        demo_rec, demo_lib = _resolve_demo_files()
        if receptor_path is None:
            receptor_path = demo_rec
        if library_path is None:
            library_path = demo_lib

    if receptor_path is None or not receptor_path.is_file():
        console.print(f"[bold #FF6B80]错误: 受体文件不存在或未指定: {receptor_path}[/]")
        return 2

    if library_path is None or not library_path.is_file():
        console.print(f"[bold #FF6B80]错误: 化合物库文件不存在或未指定: {library_path}[/]")
        return 2

    if out_dir is None:
        out_dir = settings.data_dir / "screens" / f"screen_{int(time.time())}"
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ckpt_file = out_dir / "checkpoint.jsonl"
    if resume is None:
        resume_flag = ckpt_file.is_file() and ckpt_file.stat().st_size > 0
    else:
        resume_flag = bool(resume)

    dock_box = parse_box_coords(box_center_str, box_size_str)
    rec_pdbqt, rec_pdb, rec_id, _ = _resolve_receptor_files(receptor_path, settings)

    lib_bytes = library_path.read_bytes()
    prepared_lib = prepare_library(
        content=lib_bytes,
        filename=library_path.name,
        registry=registry,
        tcm_mode=tcm_mode,
    )

    n_ready = len(prepared_lib.ligand_ids)
    if n_ready == 0:
        console.print("[bold #FF6B80]错误: 化合物库中未提取到有效小分子结构[/]")
        return 1

    parsed_seeds = [int(s.strip()) for s in seeds_str.split(",") if s.strip()] or [42, 1024, 2026]
    parsed_anchors = [a.strip().upper() for a in anchor_residues_str.replace(";", ",").split(",") if a.strip()]

    resolved_engine = _resolve_engine(engine_name)
    total_cpu = cpu or max(2, min(os.cpu_count() or 4, 16))

    config = ScreeningConfig(
        receptor_id=rec_id,
        library_id=prepared_lib.name or library_path.stem,
        box=dock_box,
        tcm_mode=tcm_mode,
        tcm_deglycosylation=tcm_mode,
        enable_tiered_funnel=True,
        tier1_keep_ratio=tier1_ratio,
        tier2_seeds=parsed_seeds,
        tier2_keep_ratio=tier2_ratio,
        max_strain_kcal=strain_threshold,
        anchor_residues=parsed_anchors,
        enable_gnina_rescore=enable_gnina_rescore,
        engine=resolved_engine.manifest.name,
        cpu=total_cpu,
    )

    console.print(
        f"[bold #D77757]启动 HTVS 三级漏斗筛选:[/] 受体={rec_id}, 库规模={n_ready} (中药双态={tcm_mode}), "
        f"引擎={resolved_engine.manifest.name}, 锚点={parsed_anchors or '无'}"
    )

    lib_ligands = {
        lid: registry.ligands[lid]
        for lid in prepared_lib.ligand_ids
        if lid in registry.ligands
    }

    with ScreenProgressManager(console, is_interactive=False) as progress_mgr:
        runner = ScreeningRunner(
            config=config,
            receptor_pdbqt=rec_pdbqt,
            protein_pdb=rec_pdb,
            ligands=lib_ligands,
            prepared_library=prepared_lib,
            out_dir=out_dir,
            engine=resolved_engine,
            resume=resume_flag,
            cpu_per_job=cpu_per_job,
            callback=progress_mgr.callback,
        )
        report = runner.run()

    leaderboard_panel = render_screening_leaderboard(
        report.candidates,
        report.n_raw,
        report.n_ready,
        top_k=20,
        out_dir=out_dir,
    )
    console.print(leaderboard_panel)
    console.print(f"[bold #4EBA65]✓ 筛选完成![/bold #4EBA65] 产出目录: {out_dir}")
    return 0


# ---------------------------------------------------------------------------
# Click CLI Entrypoint
# ---------------------------------------------------------------------------

@click.command(name="screen", help="高通量虚拟筛选 (HTVS) 三级漏斗调度引擎")
@click.option(
    "--receptor", "-r",
    type=click.Path(path_type=Path),
    default=None,
    help="受体蛋白结构文件路径 (PDB / mmCIF / PDBQT)",
)
@click.option(
    "--library", "-l",
    type=click.Path(path_type=Path),
    default=None,
    help="待筛选小分子化合物库文件 (SMILES / SDF / CSV)",
)
@click.option(
    "--box-center",
    type=str,
    default=None,
    help="对接搜索盒中心坐标，格式: 'x,y,z' (例如 '15.2,24.8,-5.1')",
)
@click.option(
    "--box-size",
    type=str,
    default="22.5,22.5,22.5",
    help="对接搜索盒三维边长，格式: 'sx,sy,sz' 或单一数值 '22.5'",
)
@click.option(
    "--tcm-mode/--no-tcm-mode",
    default=True,
    help="是否启用中药单体双轨比对模式 (自动水解衍生苷元同屏比对，默认启用)",
)
@click.option(
    "--anchor-residues",
    type=str,
    default="",
    help="口袋关键药效团锚定残基列表，逗号分隔 (例如 'ASP25:A,THR80:A')",
)
@click.option(
    "--tier1-ratio",
    type=click.FloatRange(min=0.001, max=1.0),
    default=0.15,
    show_default=True,
    help="Tier 1 快速粗筛晋级比例 (按 SILE 排序)",
)
@click.option(
    "--tier2-ratio",
    type=click.FloatRange(min=0.001, max=1.0),
    default=0.02,
    show_default=True,
    help="Tier 2 平衡精筛晋级比例 (按收敛性与张力能过滤)",
)
@click.option(
    "--seeds",
    type=str,
    default="42,1024,2026",
    show_default=True,
    help="Tier 2 多 seed 随机种子列表，逗号分隔",
)
@click.option(
    "--strain-threshold",
    type=float,
    default=8.0,
    show_default=True,
    help="配体内部构象张力能惩罚上限 (kcal/mol)",
)
@click.option(
    "--engine", "-e",
    type=str,
    default="auto",
    show_default=True,
    help="对接计算引擎 (auto: GPU 可用时自动首选 unidock，否则 vina)",
)
@click.option(
    "--cpu", "-c",
    type=click.IntRange(min=1),
    default=None,
    help="分配的总 CPU 逻辑核心数 (默认自动探测系统最大核心)",
)
@click.option(
    "--cpu-per-job",
    type=click.IntRange(min=1),
    default=2,
    show_default=True,
    help="每个 Worker 进程分配的 CPU 线程数",
)
@click.option(
    "--resume/--no-resume",
    default=None,
    help="是否从断点文件 checkpoint.jsonl 恢复筛选 (默认自动检测是否存在)",
)
@click.option(
    "--out", "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="筛选产出与交付包目录 (默认在 data_dir/screens/screen_<timestamp>)",
)
@click.option(
    "--enable-gnina-rescore/--no-gnina-rescore",
    default=True,
    show_default=True,
    help="Tier 3 中是否启用 GNINA CNN 正交重打分 (需 GNINA 可用)",
)
@click.option(
    "--demo", "-d",
    is_flag=True,
    default=False,
    help="使用内置 1HSG 蛋白酶与上市药物库运行端到端演示",
)
@click.option(
    "--non-interactive",
    is_flag=True,
    default=False,
    help="非交互/无头模式：跳过所有终端提问，参数不全直接报错退出",
)
def screen_cli(
    receptor: Path | None,
    library: Path | None,
    box_center: str | None,
    box_size: str,
    tcm_mode: bool,
    anchor_residues: str,
    tier1_ratio: float,
    tier2_ratio: float,
    seeds: str,
    strain_threshold: float,
    engine: str,
    cpu: int | None,
    cpu_per_job: int,
    resume: bool | None,
    out: Path | None,
    enable_gnina_rescore: bool,
    demo: bool,
    non_interactive: bool,
) -> int:
    """Click entrypoint for `/screen` command."""
    try:
        if non_interactive and not demo:
            if not receptor:
                console.print("[bold #FF6B80]错误: 非交互模式必须提供 --receptor 参数[/]")
                sys.exit(2)
            if not library:
                console.print("[bold #FF6B80]错误: 非交互模式必须提供 --library 参数[/]")
                sys.exit(2)

        if non_interactive or (receptor and library) or (demo and not sys.stdin.isatty()):
            rc = run_screen_headless(
                receptor_path=receptor,
                library_path=library,
                box_center_str=box_center,
                box_size_str=box_size,
                tcm_mode=tcm_mode,
                anchor_residues_str=anchor_residues,
                tier1_ratio=tier1_ratio,
                tier2_ratio=tier2_ratio,
                seeds_str=seeds,
                strain_threshold=strain_threshold,
                engine_name=engine,
                cpu=cpu,
                cpu_per_job=cpu_per_job,
                resume=resume,
                out_dir=out,
                enable_gnina_rescore=enable_gnina_rescore,
                demo=demo,
            )
            if rc != 0:
                sys.exit(rc)
            return rc

        rc = run_screen_interactive(
            demo=demo,
            out=out,
            cpu=cpu,
            cpu_per_job=cpu_per_job,
        )
        if rc != 0:
            sys.exit(rc)
        return rc

    except (KeyboardInterrupt, EOFError):
        console.print("\n[bold #FFC107]操作已取消 (User Cancelled)[/bold #FFC107]")
        sys.exit(130)
