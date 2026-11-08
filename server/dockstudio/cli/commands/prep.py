"""Receptor and Ligand Preparation commands for DockStudio CLI."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from rich.text import Text

from ...config import get_settings
from ...prep.ligand import prepare_ligand_from_file, prepare_ligand_from_smiles
from ...prep.receptor import prepare_receptor
from ..ui import console, make_pane, prompt_path, select_one


def run_prep_receptor_interactive(pdb_path: Path | None = None) -> tuple[Path, Any]:
    """Interactive receptor preparation wizard."""
    if not pdb_path:
        pdb_path = prompt_path("请指定受体 PDB/mmCIF 文件路径")

    settings = get_settings()
    out_dir = settings.prepared_dir / "receptors"
    out_dir.mkdir(parents=True, exist_ok=True)

    with console.status(f"[bold #93A5FF]正在准备受体: {pdb_path.name} (加氢、去杂原子)...[/bold #93A5FF]"):
        receptor = prepare_receptor(pdb_path.read_bytes(), pdb_path.name, out_dir)

    content = Text()
    content.append(f"✓ 受体文件: {pdb_path.name}\n", style="bold #4EBA65")
    content.append(f"• 生成 PDBQT: {Path(receptor.pdbqt_path).name}\n", style="#FFFFFF")
    content.append(f"• 包含多肽链: {', '.join(receptor.chains) if receptor.chains else '单链'}\n", style="#999999")

    if receptor.cocrystal_ligands:
        content.append(f"\n检测到 {len(receptor.cocrystal_ligands)} 个共晶小分子配体:\n", style="bold #D77757")
        for idx, c in enumerate(receptor.cocrystal_ligands):
            cx, cy, cz = c.center
            content.append(
                f"  [{idx + 1}] {c.resname}:{c.chain} (原子数: {c.n_atoms}, 质心: {cx:.2f}, {cy:.2f}, {cz:.2f})\n",
                style="#B1B9F9",
            )
    else:
        content.append("\n• 未检测到共晶配体 (后续需使用腔体预测或手动定盒)\n", style="#999999")

    if receptor.warnings:
        content.append("\n注意事项 / 警告:\n", style="bold #FFC107")
        for w in receptor.warnings:
            content.append(f"  ▲ {w}\n", style="#FFC107")

    console.print(make_pane("受体准备完成 (Receptor Prepared)", content))
    return pdb_path, receptor


def run_prep_ligand_interactive() -> tuple[str, Any]:
    """Interactive ligand preparation wizard."""
    mode = select_one("选择配体输入方式:", [
        ("输入 SMILES 字符串", "粘贴或输入单个分子的标准 SMILES"),
        ("指定配体文件路径 (SDF / MOL2)", "从本地已有 2D/3D 配体文件准备"),
    ])

    settings = get_settings()
    out_dir = settings.prepared_dir / "ligands"
    out_dir.mkdir(parents=True, exist_ok=True)

    if mode == 0:
        from prompt_toolkit import prompt
        from prompt_toolkit.formatted_text import HTML
        
        smiles = prompt(HTML("<ansicyan>请输入配体 SMILES:</ansicyan> ")).strip()
        name = prompt(HTML("<ansibrightblack>配体名称 (可选，回车默认 ligand):</ansibrightblack> ")).strip() or "ligand"
        
        with console.status(f"[bold #93A5FF]正在由 SMILES 构建 3D 构象并质子化 (pH 7.4): {name}...[/bold #93A5FF]"):
            ligand = prepare_ligand_from_smiles(smiles, name, out_dir)
    else:
        lig_path = prompt_path("请指定配体文件路径 (SDF / MOL2)")
        with console.status(f"[bold #93A5FF]正在解析配体文件并加氢转 PDBQT: {lig_path.name}...[/bold #93A5FF]"):
            ligand = prepare_ligand_from_file(lig_path, out_dir)

    content = Text()
    content.append(f"✓ 配体准备成功: {ligand.name or 'ligand'}\n", style="bold #4EBA65")
    content.append(f"• 生成 PDBQT: {Path(ligand.pdbqt_path).name}\n", style="#FFFFFF")
    if ligand.n_rotatable_bonds is not None:
        content.append(f"• 可旋转键数目: {ligand.n_rotatable_bonds}\n", style="#B1B9F9")
    if ligand.smiles:
        content.append(f"• SMILES: {ligand.smiles[:60]}...\n", style="#999999")

    console.print(make_pane("配体准备完成 (Ligand Prepared)", content))
    return ligand.name or "ligand", ligand
