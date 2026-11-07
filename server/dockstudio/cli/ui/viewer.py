"""Terminal molecular viewer & interaction tree visualizer.

Provides:
- render_interaction_tree: Rich Tree rendering of PLIP protein-ligand contacts
- render_posecheck_table: Table of PoseBusters 18 geometry/chemical checks
- render_ligand_summary: Molecular summary & 2D depiction
"""

from __future__ import annotations

from typing import Any

from rich.table import Table
from rich.tree import Tree

from .components import format_affinity
from .theme import (
    COLOR_BORDER,
    COLOR_CLAUDE,
    ICON_ERROR,
    ICON_SUCCESS,
    ICON_WARNING,
)


def render_interaction_tree(
    receptor_name: str,
    ligand_name: str,
    interactions: list[dict[str, Any]],
    pose_index: int = 1,
    affinity: float = 0.0,
) -> Tree:
    """Renders a rich textual tree of PLIP contacts for a docked pose."""
    header = (
        f"[bold {COLOR_CLAUDE}]构象 #{pose_index}[/bold {COLOR_CLAUDE}] "
        f"({format_affinity(affinity)} kcal/mol) "
        f"• 靶点: [bold]{receptor_name}[/bold] ── 配体: [bold]{ligand_name}[/bold]"
    )
    root = Tree(header)

    if not interactions:
        root.add("[#999999]未检测到显著强相互作用力 (可能主要由非特异性空间互补主导)[/#999999]")
        return root

    # Group by interaction type
    groups: dict[str, list[dict[str, Any]]] = {
        "hbond": [],
        "salt_bridge": [],
        "pi_stacking": [],
        "hydrophobic": [],
        "halogen": [],
        "water_bridge": [],
    }
    type_labels = {
        "hbond": ("氢键相互作用 (Hydrogen Bonds)", "#4EBA65"),
        "salt_bridge": ("盐桥作用 (Salt Bridges)", "#FFC107"),
        "pi_stacking": ("π-π 堆叠 (Pi-Stacking)", "#B1B9F9"),
        "hydrophobic": ("疏水接触 (Hydrophobic Contacts)", "#93A5FF"),
        "halogen": ("卤键 (Halogen Bonds)", "#48968C"),
        "water_bridge": ("水桥作用 (Water Bridges)", "#999999"),
    }

    for item in interactions:
        t = item.get("type", "hydrophobic")
        if t in groups:
            groups[t].append(item)
        else:
            groups.setdefault("other", []).append(item)

    for itype, (label, color) in type_labels.items():
        items = groups.get(itype, [])
        if not items:
            continue
        branch = root.add(f"[bold {color}]{label}[/bold {color}] ([#FFFFFF]{len(items)} 条[/#FFFFFF])")
        for contact in items:
            res = contact.get("residue", "Unknown")
            dist = contact.get("distance")
            dist_str = f"{dist:.2f} Å" if dist is not None else ""
            detail = contact.get("detail", "")
            
            line = f"[{color}]•[/{color}] [bold #FFFFFF]{res:10s}[/bold #FFFFFF]"
            if dist_str:
                line += f" 距离: [{color}]{dist_str}[/{color}]"
            if detail:
                line += f"  [dim]({detail})[/dim]"
            branch.add(line)

    return root


def render_posecheck_table(pose_check: dict[str, Any]) -> Table:
    """Renders PoseBusters 18 checks into a compact terminal table."""
    table = Table(
        title="[bold #D77757]PoseBusters 物理化学构象质检 (18项)[/bold #D77757]",
        border_style=COLOR_BORDER,
        header_style="bold #FFFFFF",
    )
    table.add_column("检查项 (Check Name)", justify="left")
    table.add_column("状态", justify="center")
    table.add_column("详情 / 阈值", justify="left")

    results = pose_check.get("results", {})
    if not results:
        table.add_row("质检状态", ICON_WARNING, "PoseBusters 未生成逐项报告 (使用内置几何筛查)")
        return table

    for name, detail in results.items():
        passed = detail.get("passed", True)
        icon = ICON_SUCCESS if passed else ICON_ERROR
        desc = detail.get("message", "正常") if not passed else "通过"
        table.add_row(name, icon, desc)

    return table
