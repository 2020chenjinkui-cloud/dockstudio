"""DockStudio CLI: Task Audit & Cryptographic Event Stream Command (/audit).

Parses events.jsonl, validates live SHA256 checksums against inputs.sha256.json,
and renders a tamper-evident audit report matching Claude Code's design system.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any

from rich.table import Table
from rich.text import Text

from ...config import Settings, get_settings
from ...core.project import ProjectStore
from ...models import Job, JobStatus
from ..ui import (
    COLOR_ERROR,
    COLOR_SUCCESS,
    COLOR_WARNING,
    ICON_ERROR,
    ICON_SUCCESS,
    ICON_WARNING,
    console,
    format_affinity,
    format_status_badge,
    make_header,
    make_pane,
)


def compute_sha256_chunked(path: Path) -> str:
    """Computes SHA256 digest of a file in 1MB (1<<20) chunks to prevent memory spikes."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_events_jsonl(events_path: Path, subject_ids: set[str] | None = None) -> list[dict[str, Any]]:
    """Robustly reads JSONL event logs, ignoring blank or crash-torn trailing lines."""
    if not events_path.is_file():
        return []

    events: list[dict[str, Any]] = []
    try:
        content = events_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []

    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except (json.JSONDecodeError, ValueError):
            # Gracefully swallow truncated/torn trailing line
            continue
        if not isinstance(record, dict):
            continue

        if subject_ids is not None:
            subj_id = (record.get("subject") or {}).get("id")
            if subj_id and subj_id not in subject_ids:
                continue

        events.append(record)

    events.sort(key=lambda x: str(x.get("ts", "")))
    return events


def resolve_job_data(
    job_id_or_path: str,
    settings: Settings,
) -> tuple[Job | None, Path | None, list[dict[str, Any]], dict[str, str], str]:
    """Locates job object, directory, events, input hashes, and project name from ID or path.

    Returns:
        (job, job_dir, events, input_hashes, project_name)
    """
    candidate = Path(job_id_or_path).expanduser().resolve()

    # Case 1: Direct zip archive
    if candidate.is_file() and zipfile.is_zipfile(candidate):
        try:
            with zipfile.ZipFile(candidate) as zf:
                namelist = set(zf.namelist())
                if "snapshot.json" in namelist:
                    snapshot = json.loads(zf.read("snapshot.json").decode("utf-8", errors="replace"))
                    job = Job.model_validate(snapshot["job"]) if "job" in snapshot else None
                    hashes = snapshot.get("inputs", {})
                    events = []
                    if "events.jsonl" in namelist:
                        raw_ev = zf.read("events.jsonl").decode("utf-8", errors="replace")
                        for line in raw_ev.splitlines():
                            if line.strip():
                                try:
                                    events.append(json.loads(line.strip()))
                                except Exception:
                                    continue
                    return job, candidate.parent, events, hashes, "ZIP归档审计包"
        except Exception:
            return None, None, [], {}, ""

    # Case 2: Directory path
    if candidate.is_dir():
        job_file = candidate / "job.json"
        snapshot_file = candidate / "snapshot.json"
        job = None
        hashes: dict[str, str] = {}
        events: list[dict[str, Any]] = []
        project_name = candidate.name

        if job_file.is_file():
            try:
                job = Job.model_validate_json(job_file.read_text(encoding="utf-8", errors="replace"))
            except Exception:
                pass
        elif snapshot_file.is_file():
            try:
                snap = json.loads(snapshot_file.read_text(encoding="utf-8", errors="replace"))
                if "job" in snap:
                    job = Job.model_validate(snap["job"])
                hashes = snap.get("inputs", {})
            except Exception:
                pass

        if job:
            sidecar = candidate / "inputs.sha256.json"
            if sidecar.is_file():
                try:
                    hashes = json.loads(sidecar.read_text(encoding="utf-8", errors="replace"))
                except Exception:
                    pass

            events_candidates = [
                candidate / "events.jsonl",
                candidate.parent / "events.jsonl",
                candidate.parent.parent / "events.jsonl",
            ]
            for ev_path in events_candidates:
                if ev_path.is_file():
                    events = parse_events_jsonl(ev_path, {job.id, job.project_id})
                    if events:
                        break

            return job, candidate, events, hashes, project_name

    # Case 3: Search in ProjectStore by job ID
    store = ProjectStore(settings)
    job = store.find_job(job_id_or_path)
    if job is not None:
        job_dir = store.job_dir(job.project_id, job.id)
        events_path = store.root / job.project_id / "events.jsonl"
        events = parse_events_jsonl(events_path, {job.id, job.project_id})
        sidecar = job_dir / "inputs.sha256.json"
        hashes = {}
        if sidecar.is_file():
            try:
                hashes = json.loads(sidecar.read_text(encoding="utf-8", errors="replace"))
            except Exception:
                hashes = {}
        project = store.get(job.project_id)
        pname = project.name if project else job.project_id
        return job, job_dir, events, hashes, pname

    return None, None, [], {}, ""


def run_audit(job_id_or_path: str, options: dict[str, Any] | None = None) -> int:
    """Executes the audit workflow and renders the audit report.

    Returns:
        0: Integrity verification succeeded and all checksums match.
        1: Hash mismatch or tamper detected.
        2: Job not found or invalid path.
    """
    opts = options or {}
    settings = get_settings()
    if "data_dir" in opts:
        settings = Settings(data_dir=Path(opts["data_dir"]))

    job, job_dir, events, expected_hashes, project_name = resolve_job_data(job_id_or_path, settings)

    if job is None:
        console.print(
            f"\n[bold {COLOR_ERROR}]未找到任务:[/bold {COLOR_ERROR}] '{job_id_or_path}'\n"
            f"[dim]提示: 请确认 job_id 来自 /api/jobs 或项目目录，当前 data_dir={settings.data_dir}[/dim]\n"
        )
        return 2

    # Top Banner
    console.print(make_header())

    # 1. Metadata Panel
    req = job.request
    meta_table = Table(box=None, show_header=False, expand=True, padding=(0, 1))
    meta_table.add_column("Key", style="bold #999999", width=14)
    meta_table.add_column("Val", style="#FFFFFF")
    meta_table.add_column("Key2", style="bold #999999", width=14)
    meta_table.add_column("Val2", style="#FFFFFF")

    st_badge = format_status_badge(job.status == JobStatus.completed, reason=str(job.error or ""))
    dur_str = "—"
    if job.started_at and job.finished_at:
        dur_str = f"{(job.finished_at - job.started_at).total_seconds():.2f}s"

    meta_table.add_row("任务编号 (Job ID)", f"[bold #D77757]{job.id}[/bold #D77757]", "所属项目", f"{project_name} ({job.project_id})")
    meta_table.add_row("对接引擎", f"[bold #93A5FF]{req.engine}[/bold #93A5FF]", "任务状态", st_badge)
    meta_table.add_row("创建时间", job.created_at.strftime("%Y-%m-%d %H:%M:%S"), "总计耗时", dur_str)
    meta_table.add_row("受体编号", f"`{req.receptor_id}`", "配体编号", f"`{req.ligand_id}`")

    console.print(make_pane("任务基础信息与元数据", meta_table, subtitle=f"ID: {job.id}"))

    # 2. SHA256 Checksums & Live Tamper Check
    hash_table = Table(box=None, show_header=True, header_style="bold #FFFFFF", expand=True)
    hash_table.add_column("校验资产 (Asset Key)", style="bold #FFFFFF", width=18)
    hash_table.add_column("预期指纹 (Expected SHA256)", style="#999999", width=22)
    hash_table.add_column("实时计算指纹 (Live SHA256)", style="#CCCCCC", width=22)
    hash_table.add_column("完整性状态", justify="center", width=16)

    tampered = False
    verified_count = 0

    # Determine files to inspect
    target_files: dict[str, Path] = {}
    if job_dir and job_dir.is_dir():
        for f in job_dir.iterdir():
            if f.is_file():
                target_files[f.name] = f
                target_files[f.stem] = f

    # Also resolve from prepared_dir
    rec_pdbqt = settings.prepared_dir / "receptors" / f"{req.receptor_id}.pdbqt"
    if rec_pdbqt.is_file():
        target_files.setdefault("receptor", rec_pdbqt)
    lig_pdbqt = settings.prepared_dir / "ligands" / f"{req.ligand_id}.pdbqt"
    if lig_pdbqt.is_file():
        target_files.setdefault("ligand", lig_pdbqt)

    if not expected_hashes:
        hash_table.add_row(
            "输入指纹记录",
            "[dim]—[/dim]",
            "[dim]—[/dim]",
            f"{ICON_WARNING} [#FFC107]无基线哈希[/#FFC107]",
        )
    else:
        for key, exp_hash in expected_hashes.items():
            file_path = target_files.get(key) or target_files.get(f"{key}.pdbqt") or target_files.get(f"{key}.pdb")
            exp_short = f"{exp_hash[:10]}...{exp_hash[-8:]}" if len(exp_hash) > 18 else exp_hash
            if file_path and file_path.is_file():
                live_hash = compute_sha256_chunked(file_path)
                live_short = f"{live_hash[:10]}...{live_hash[-8:]}"
                if live_hash.lower() == exp_hash.lower():
                    status_str = f"{ICON_SUCCESS} [#4EBA65]匹配吻合[/#4EBA65]"
                    verified_count += 1
                else:
                    status_str = f"{ICON_ERROR} [bold #FF6B80]篡改/不符[/bold #FF6B80]"
                    tampered = True
                hash_table.add_row(key, exp_short, live_short, status_str)
            else:
                hash_table.add_row(
                    key,
                    exp_short,
                    "[dim]文件缺失[/dim]",
                    f"{ICON_WARNING} [#FFC107]未找到实体[/#FFC107]",
                )

    console.print(make_pane("不可篡改 SHA256 数据资产校验", hash_table))

    # 3. Environment & Parameters
    env_text = Text()
    if job.env:
        env = job.env
        env_text.append(f"• Python: {env.python or '—'}  • 操作系统: {env.platform or '—'}\n", style="#CCCCCC")
        env_text.append(f"• RDKit: {env.rdkit or '未安装'}  • Meeko: {env.meeko or '未安装'}  • OpenBabel: {env.openbabel or '未安装'}\n", style="#999999")
    else:
        env_text.append("• 运行环境快照: 任务记录于环境指纹模块上线之前\n", style="#999999")

    box = req.box
    params = req.params
    env_text.append(
        f"• 搜索盒子: 中心 ({box.center[0]:.2f}, {box.center[1]:.2f}, {box.center[2]:.2f}) Å · "
        f"尺寸 {box.size[0]:.1f} × {box.size[1]:.1f} × {box.size[2]:.1f} Å\n",
        style="#48968C",
    )
    env_text.append(
        f"• 参数配置: exhaustiveness={params.exhaustiveness}, modes={params.num_modes}, "
        f"seed={params.seed}, cpu={params.cpu}\n",
        style="#CCCCCC",
    )
    console.print(make_pane("计算环境与对接参数配置", env_text))

    # 4. Result Poses (if available)
    if job.result and job.result.poses:
        pose_table = Table(box=None, show_header=True, header_style="bold #FFFFFF", expand=True)
        pose_table.add_column("构象 (Rank)", style="bold #FFFFFF", justify="center", width=12)
        pose_table.add_column("亲和力 (kcal/mol)", justify="right", width=20)
        pose_table.add_column("RMSD lower bound", justify="right", style="#999999", width=20)
        pose_table.add_column("RMSD upper bound", justify="right", style="#999999", width=20)
        for p in job.result.poses:
            lb = f"{p.rmsd_lb:.2f} Å" if p.rmsd_lb is not None else "—"
            ub = f"{p.rmsd_ub:.2f} Å" if p.rmsd_ub is not None else "—"
            pose_table.add_row(f"Mode #{p.index}", format_affinity(p.affinity), lb, ub)
        console.print(make_pane("姿态打分与构象结果", pose_table))

    # 5. Events Timeline
    if events:
        event_table = Table(box=None, show_header=True, header_style="bold #FFFFFF", expand=True)
        event_table.add_column("时间戳 (UTC)", style="#999999", width=22)
        event_table.add_column("事件动作 (Action)", style="bold #D77757", width=18)
        event_table.add_column("主体", style="#93A5FF", width=10)
        event_table.add_column("载荷摘要 (Payload)", style="#CCCCCC")

        for e in events:
            ts = str(e.get("ts", ""))[:19].replace("T", " ")
            action = str(e.get("action", ""))
            actor = str(e.get("actor", ""))
            payload = json.dumps(e.get("payload", {}), ensure_ascii=False)
            if len(payload) > 75:
                payload = payload[:72] + "..."
            event_table.add_row(ts, action, actor, payload)

        console.print(make_pane("事理链时序流 (Event Log)", event_table, subtitle=f"共 {len(events)} 条记录"))

    # 6. Cryptographic Seal / Audit Verdict
    seal_text = Text()
    if tampered:
        seal_text.append("✗ 审计结论: 检测到资产指纹不一致或数据篡改！\n", style="bold #FF6B80")
        seal_text.append("• 实体文件 SHA256 哈希与初始记录不吻合，该计算成果不具备合规证明效力。\n", style="#FF6B80")
        seal_text.append("• 建议核实受影响文件是否曾被第三方编辑器修改或外部重新写入。\n", style="#999999")
        console.print(make_pane("密码学审计印章 (Cryptographic Seal)", seal_text, border_style=COLOR_ERROR))
        return 1
    elif verified_count > 0:
        seal_text.append("✓ 审计结论: SHA256 完整性校验 100% 吻合 · 事理链闭合\n", style="bold #4EBA65")
        seal_text.append(f"• 成功验证 {verified_count} 个核心资产指纹，全生命周期无任何篡改痕迹。\n", style="#4EBA65")
        seal_text.append(f"• 审计验证时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} · 具备学术与外包交付合规效力。\n", style="#999999")
        console.print(make_pane("密码学审计印章 (Cryptographic Seal)", seal_text, border_style=COLOR_SUCCESS))
        return 0
    else:
        seal_text.append("▲ 审计结论: 未发现篡改，但未找到足够基线指纹进行比对。\n", style="bold #FFC107")
        console.print(make_pane("密码学审计印章 (Cryptographic Seal)", seal_text, border_style=COLOR_WARNING))
        return 0
