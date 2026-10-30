"""Audit bundle: snapshot + events + engine log + report.md for agent review."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

from ..core.events import sha256_file
from ..models import DockResult, Job

_DISCLAIMER = (
    "对接打分（kcal/mol）是经验打分函数的相对排序指标，不是实验结合亲和力实测值；"
    "AI 复核结果为独立预测模型输出，仅供交叉参考。"
)


def _hline(text: str) -> str:
    return f"\n## {text}\n"


def build_audit_markdown(
    job: Job,
    events: list[dict],
    input_hashes: dict | None = None,
    project_name: str = "",
) -> str:
    """事理链 Markdown：何时何地 -> 输入 -> 环境指纹 -> 参数 -> 结果 -> 事件序列."""
    parts: list[str] = []
    parts.append(f"# DockStudio 对接任务审计报告 — job {job.id}\n")

    parts.append(_hline("起始终起"))
    parts.append(
        f"- 项目: {project_name or job.project_id}（{job.project_id}）\n"
        f"- 任务: {job.id}\n"
        f"- 创建: {job.created_at.isoformat()}\n"
        f"- 开始: {job.started_at.isoformat() if job.started_at else '—'}\n"
        f"- 结束: {job.finished_at.isoformat() if job.finished_at else '—'}\n"
        f"- 状态: {job.status.value}\n" + (f"- 错误: {job.error}\n" if job.error else "")
    )

    req = job.request
    parts.append(_hline("输入"))
    parts.append(f"- receptor_id: `{req.receptor_id}`\n- ligand_id: `{req.ligand_id}`\n")
    if input_hashes:
        for key, digest in input_hashes.items():
            parts.append(f"- {key} SHA256: `{digest}`\n")
            parts.append(f"  - 前缀: `{digest[:16]}`\n")

    parts.append(_hline("环境指纹"))
    env = job.env
    if env is None:
        parts.append("- 无（该任务记录于指纹功能上线之前）\n")
    else:
        parts.append(
            f"- python: {env.python}\n"
            f"- rdkit: {env.rdkit}\n"
            f"- meeko: {env.meeko}\n"
            f"- openbabel: {env.openbabel}\n"
            f"- platform: {env.platform}\n"
        )

    parts.append(_hline("盒子与参数"))
    box = req.box
    parts.append(
        f"- center: ({box.center[0]:.3f}, {box.center[1]:.3f}, {box.center[2]:.3f}) Å\n"
        f"- size: {box.size[0]:.1f} × {box.size[1]:.1f} × {box.size[2]:.1f} Å\n"
    )
    params = req.params
    parts.append(
        f"- engine: {req.engine}\n"
        f"- scoring: {params.scoring}, exhaustiveness: {params.exhaustiveness}, "
        f"num_modes: {params.num_modes}, energy_range: {params.energy_range}, "
        f"cpu: {params.cpu}, seed: {params.seed}\n"
    )

    parts.append(_hline("姿态打分"))
    result: DockResult | None = job.result
    if result is None or not result.poses:
        parts.append("- 无结果（任务未完成或失败）\n")
    else:
        parts.append("| # | affinity (kcal/mol) | rmsd_lb | rmsd_ub |\n")
        parts.append("|---|---|---|---|\n")
        for p in result.poses:
            lb = f"{p.rmsd_lb:.2f}" if p.rmsd_lb is not None else "—"
            ub = f"{p.rmsd_ub:.2f}" if p.rmsd_ub is not None else "—"
            parts.append(f"| {p.index} | {p.affinity:.2f} | {lb} | {ub} |\n")

    parts.append(_hline("事件序列"))
    if not events:
        parts.append("- （该项目尚无 job 相关事件）\n")
    else:
        parts.append("| ts | action | actor | payload 摘要 |\n|---|---|---|---|\n")
        for e in events:
            payload = json.dumps(e.get("payload", {}), ensure_ascii=False)[:120]
            parts.append(
                f"| {e.get('ts', '')} | {e.get('action', '')} | {e.get('actor', '')} | {payload} |\n"
            )

    parts.append(_hline("免责声明"))
    parts.append(_DISCLAIMER + "\n")
    return "".join(parts)


def input_hashes_for_job(job: Job, receptor_pdbqt: Path | None, ligand_pdbqt: Path | None) -> dict:
    out: dict[str, str] = {}
    for name, path in (("receptor", receptor_pdbqt), ("ligand", ligand_pdbqt)):
        if path is not None and Path(path).is_file():
            out[name] = sha256_file(Path(path))
    return out


def build_audit_zip(
    job: Job,
    project_dir: Path,
    events: list[dict],
    out_path: Path,
    include_files: dict[str, Path] | None = None,
    project_name: str = "",
    input_hashes: dict | None = None,
    prepared: dict | None = None,
) -> Path:
    """Pack snapshot.json + events.jsonl + engine.log + artifacts + report.md."""
    snapshot = {
        "job": json.loads(job.model_dump_json()),
        "inputs": input_hashes or {},
    }
    if prepared:
        snapshot["prepared"] = prepared
    report_md = build_audit_markdown(job, events, input_hashes, project_name)
    job_dir = Path(job.result.run_dir) if job.result else None

    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("snapshot.json", json.dumps(snapshot, ensure_ascii=False, indent=2))
        zf.writestr(
            "events.jsonl",
            "\n".join(json.dumps(e, ensure_ascii=False) for e in events) + "\n",
        )
        if job.result and job.result.log:
            zf.writestr("engine.log", job.result.log)
        zf.writestr("report.md", report_md)
        seen = set()
        for name, path in (include_files or {}).items():
            p = Path(path)
            if p.is_file() and name not in seen:
                zf.write(p, arcname=name)
                seen.add(name)
        if job_dir is not None and job_dir.is_dir():
            for name in ("receptor.pdbqt", "ligand.pdbqt", "out.pdbqt"):
                candidate = job_dir / name
                if candidate.is_file() and name not in seen:
                    zf.write(candidate, arcname=name)
                    seen.add(name)
    return out_path
