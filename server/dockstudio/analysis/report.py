"""Self-contained HTML report for a finished docking job (3Dmol.js via CDN).

M3 ("report v2"): optional section filtering (`sections`) and custom title."""

from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path

from ..models import DockJobRequest, DockResult, Interaction, PoseAnalysis, Project

ALL_SECTIONS = ("snapshot", "scores", "interactions", "diagram", "viewer3d", "disclaimer")


def _esc(text: object) -> str:
    return html.escape(str(text))


def _fmt(value: object, ndigits: int = 2) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.{ndigits}f}"
    return _esc(value)


_TYPE_LABELS = {
    "hbond": "氢键",
    "hydrophobic": "疏水接触",
    "pi_stacking": "π 堆积",
    "salt_bridge": "盐桥",
    "halogen": "卤键",
    "water_bridge": "水桥",
}


def _interactions_rows(interactions: list[Interaction]) -> str:
    if not interactions:
        return '<tr><td colspan="4" class="muted">未检测到相互作用</td></tr>'
    rows = []
    for i in interactions:
        rows.append(
            "<tr>"
            f"<td>{_esc(_TYPE_LABELS.get(i.type, i.type))}</td>"
            f"<td>{_esc(i.residue)}</td>"
            f"<td>{_fmt(i.distance)} Å</td>"
            f'<td class="muted">{_esc(i.detail)}</td>'
            "</tr>"
        )
    return "\n".join(rows)


def _poses_rows(result: DockResult) -> str:
    rows = []
    for pose in result.poses:
        rows.append(
            "<tr>"
            f"<td>{pose.index}</td>"
            f"<td><b>{pose.affinity:.2f}</b></td>"
            f"<td>{_fmt(pose.rmsd_lb)}</td>"
            f"<td>{_fmt(pose.rmsd_ub)}</td>"
            "</tr>"
        )
    return "\n".join(rows)


def _pose_check_html(analysis: PoseAnalysis | None) -> str:
    if analysis is None or analysis.pose_check is None:
        return ""
    pc = analysis.pose_check
    badge = "✅ 通过" if pc.ok else "❌ 未通过"
    rows = (
        "".join(
            f"<tr><td>{_esc(i.rule)}</td><td>{_esc(i.severity)}</td><td>{_esc(i.message)}</td></tr>"
            for i in pc.issues
        )
        or '<tr><td colspan="3" class="muted">未发现问题</td></tr>'
    )
    return (
        f"<h3>姿态质检（{_esc(pc.source)}）：{badge}</h3>"
        f"<table><tr><th>规则</th><th>级别</th><th>说明</th></tr>{rows}</table>"
    )


def _inline_file(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def generate_report(
    project: Project,
    request: DockJobRequest,
    result: DockResult,
    pose_analyses: list[PoseAnalysis],
    receptor_pdb_path: Path,
    pose_sdf_paths: list[Path],
    receptor_name: str,
    ligand_name: str,
    sections: list[str] | None = None,
    title: str | None = None,
    diagram_svg: str = "",
) -> str:
    """Render a self-contained HTML report (structures inlined as text).

    `sections` selects a subset of ALL_SECTIONS (default: all)."""
    want = set(sections) if sections else set(ALL_SECTIONS)
    top_analysis = pose_analyses[0] if pose_analyses else None
    receptor_pdb = _inline_file(receptor_pdb_path)
    top_pose_sdf = _inline_file(pose_sdf_paths[0]) if pose_sdf_paths else ""
    report_title = (title or "").strip()[:80] or f"{project.name} 分子对接报告"

    box = result.box
    params = result.params
    engine_line = (
        f"{_esc(result.engine)} {_esc(result.engine_version)} · 打分函数 {_esc(params.scoring)}"
    )
    snapshot_rows = [
        ("项目", project.name),
        ("生成时间", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ("引擎", engine_line),
        ("受体", receptor_name),
        ("配体", ligand_name),
        ("受体 ID", request.receptor_id),
        ("配体 ID", request.ligand_id),
        ("盒子中心", f"({box.center[0]:.2f}, {box.center[1]:.2f}, {box.center[2]:.2f}) Å"),
        ("盒子尺寸", f"{box.size[0]:.1f} × {box.size[1]:.1f} × {box.size[2]:.1f} Å"),
        (
            "采样参数",
            f"exhaustiveness={params.exhaustiveness}, num_modes={params.num_modes}, "  # noqa: ISC004
            f"energy_range={params.energy_range}, cpu={params.cpu}, seed={params.seed}",
        ),
    ]
    snapshot_html = "\n".join(
        f"<tr><th>{_esc(k)}</th><td>{_esc(v)}</td></tr>" for k, v in snapshot_rows
    )
    NOW_TAG = datetime.now().strftime("%Y-%m-%d %H:%M")

    snapshot_block = ""
    if "snapshot" in want:
        snapshot_block = f"<h2>任务快照</h2>\n<table>\n{snapshot_html}\n</table>\n"

    scores_block = ""
    if "scores" in want:
        scores_block = (
            "<h2>姿态打分</h2>\n<table>\n"
            '<tr><th style="width:auto">#</th><th style="width:auto">打分 (kcal/mol)</th>\n'
            '    <th style="width:auto">RMSD 下界 (Å)</th><th style="width:auto">RMSD 上界 (Å)</th></tr>\n'
            f"{_poses_rows(result)}\n</table>\n"
            f'<p class="muted">共 {len(result.poses)} 个姿态，按打分升序。</p>\n'
        )

    interactions_block = ""
    if "interactions" in want:
        warnings_html = "".join(
            f'<p class="muted">警告: {_esc(w)}</p>'
            for w in (top_analysis.warnings if top_analysis else [])
        )
        interactions_block = (
            f"<h2>Top 姿态相互作用（{_esc(top_analysis.source if top_analysis else 'n/a')} 分析）</h2>\n"
            "<table>\n"
            '<tr><th style="width:auto">类型</th><th style="width:auto">残基</th>\n'
            '    <th style="width:auto">距离</th><th style="width:auto">备注</th></tr>\n'
            f"{_interactions_rows(top_analysis.interactions if top_analysis else [])}\n</table>\n"
            f"{warnings_html}\n{_pose_check_html(top_analysis)}\n"
        )

    diagram_block = ""
    if "diagram" in want and diagram_svg:
        diagram_block = f"<h2>2D 相互作用示意图</h2>\n{diagram_svg}\n"

    viewer_block = ""
    if "viewer3d" in want:
        receptor_json = json.dumps(receptor_pdb)
        pose_json = json.dumps(top_pose_sdf)
        viewer_block = (
            "<h2>3D 结合模式</h2>\n"
            '<div id="viewer"></div>\n'
            "<script>\n(function() {\n"
            f"  const receptorPDB = {receptor_json};\n"
            f"  const poseSDF = {pose_json};\n"
            '  const viewer = $3Dmol.createViewer("viewer", {backgroundColor: "white"});\n'
            "  if (receptorPDB) {\n"
            '    viewer.addModel(receptorPDB, "pdb");\n'
            '    viewer.setStyle({model: 0}, {cartoon: {color: "spectrum"}});\n'
            "  }\n"
            "  if (poseSDF) {\n"
            '    viewer.addModel(poseSDF, "sdf");\n'
            '    viewer.setStyle({model: 1}, {stick: {colorscheme: "greenCarbon"}});\n'
            "  }\n"
            "  viewer.zoomTo();\n"
            "  viewer.render();\n"
            "})();\n</script>\n"
        )

    disclaimer_block = ""
    if "disclaimer" in want:
        disclaimer_block = (
            '<div class="disclaimer">\n<b>免责声明：</b>对接打分（kcal/mol）是经验打分函数的相对排序指标，\n'
            "<b>不是</b>实验结合亲和力实测值，仅可用于同一靶点同一配体系列内的相对比较。\n</div>\n"
        )

    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>{_esc(report_title)} — DockStudio</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  :root {{ color-scheme: light; }}
  body {{ font-family: "Segoe UI", "Microsoft YaHei", system-ui, sans-serif;
         max-width: 960px; margin: 0 auto; padding: 24px; color: #1f2328; }}
  h1 {{ font-size: 1.6rem; border-bottom: 2px solid #2563eb; padding-bottom: 8px; }}
  h2 {{ font-size: 1.15rem; margin-top: 28px; color: #1e40af; }}
  h3 {{ font-size: 1.0rem; margin-top: 18px; color: #334; }}
  table {{ border-collapse: collapse; width: 100%; margin-top: 8px; }}
  th, td {{ border: 1px solid #d0d7de; padding: 6px 10px; text-align: left;
           font-size: 0.92rem; }}
  th {{ background: #f6f8fa; width: 160px; }}
  .muted {{ color: #6e7781; }}
  #viewer {{ width: 100%; height: 420px; position: relative; border: 1px solid #d0d7de;
             border-radius: 6px; }}
  .disclaimer {{ background: #fff8e6; border: 1px solid #f0c36d; border-radius: 6px;
                 padding: 10px 14px; margin-top: 24px; font-size: 0.9rem; }}
  .refs {{ font-size: 0.85rem; color: #444; }}
  .refs li {{ margin: 4px 0; }}
  @media print {{ #viewer {{ page-break-inside: avoid; }} }}
</style>
</head>
<body>
<h1>{_esc(report_title)}</h1>

{snapshot_block}
{scores_block}
{interactions_block}
{diagram_block}
{viewer_block}
{disclaimer_block}

<h2>引用</h2>
<ul class="refs">
  <li>Eberhardt, J. et al. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field,
      and Python Bindings. <i>J. Chem. Inf. Model.</i> 2021.</li>
  <li>Trott, O.; Olson, A. J. AutoDock Vina: Improving the Speed and Accuracy of Docking.
      <i>J. Comput. Chem.</i> 2010.</li>
  <li>Meeko — Forli Lab, Scripps Research. https://github.com/forlilab/meeko</li>
  <li>RDKit: Open-source cheminformatics. https://www.rdkit.org</li>
  <li>Salentin, S. et al. PLIP: fully automated protein-ligand interaction profiler.
      <i>Nucleic Acids Res.</i> 2015.</li>
  <li>Rego, N.; Koes, D. 3Dmol.js: molecular visualization with WebGL.
      <i>Bioinformatics</i> 2015.</li>
</ul>

<p class="muted">DockStudio v0.4 生成 · {NOW_TAG}</p>
</body>
</html>
"""
