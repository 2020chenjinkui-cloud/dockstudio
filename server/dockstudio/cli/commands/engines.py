"""Engine & Hardware Probe command for DockStudio CLI."""

from __future__ import annotations

import os
import shutil
import subprocess

from rich.table import Table

from ...config import get_settings
from ...engines import available_engines
from ..ui import console, format_status_badge, make_header, make_pane


def detect_gpu_info() -> str | None:
    """Best-effort detection of NVIDIA GPU model and VRAM."""
    # 1. Try nvidia-smi command
    if shutil.which("nvidia-smi"):
        try:
            out = subprocess.check_output(
                ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
                text=True,
                timeout=2,
            ).strip()
            if out:
                first_gpu = out.splitlines()[0]
                parts = [p.strip() for p in first_gpu.split(",")]
                if len(parts) >= 2:
                    return f"{parts[0]} ({int(parts[1]) // 1024}GB VRAM)"
                return parts[0]
        except Exception:  # noqa: BLE001
            pass

    # 2. Try torch if available
    try:
        import torch
        if torch.cuda.is_available():
            name = torch.cuda.get_device_name(0)
            mem = torch.cuda.get_device_properties(0).total_memory // (1024**3)
            return f"{name} ({mem}GB VRAM)"
    except Exception:  # noqa: BLE001
        pass

    return None


def build_engines_pane() -> Any:
    """Construct the rich Panel containing the engines probe table."""
    settings = get_settings()
    cpu_count = os.cpu_count() or 1

    table = Table(
        box=None,
        show_header=True,
        header_style="bold #FFFFFF",
        expand=True,
    )
    table.add_column("引擎 (Engine)", justify="left", style="bold #FFFFFF", width=12)
    table.add_column("版本", justify="left", style="#999999", width=10)
    table.add_column("状态", justify="left", width=22)
    table.add_column("打分函数", justify="left", style="#48968C", width=18)
    table.add_column("特性 / 算力类型", justify="left", style="#999999")

    # Get available engines
    engines = available_engines()
    for eng in engines:
        ok, reason = eng.is_available()
        is_gpu = getattr(eng.manifest, "requires_gpu", False) or "unidock" in eng.manifest.name.lower()
        status_str = format_status_badge(ok, is_gpu=is_gpu, reason=reason)
        scoring = ", ".join(eng.manifest.scoring_functions)

        # Tags / note
        notes = []
        if is_gpu:
            notes.append("[bold #93A5FF]GPU 加速[/bold #93A5FF]")
        if "vina" in eng.manifest.name:
            notes.append("内置官方二进制 · 黄金标准")
        elif "unidock" in eng.manifest.name:
            notes.append("大规模虚拟筛选推荐 (1000x)")
        elif "gnina" in eng.manifest.name:
            notes.append("CNN 深度学习打分")
        elif "smina" in eng.manifest.name:
            notes.append("Vinardo 打分 / 自定义力场")

        note_str = " · ".join(notes) if notes else getattr(eng.manifest, "homepage", "")[:30]
        table.add_row(
            eng.manifest.name,
            f"v{eng.manifest.version}",
            status_str,
            scoring,
            note_str,
        )

    # Check boltz AI rescore capability
    try:
        from ...analysis.ai_rescore import get_availability as get_boltz_availability
        boltz_info = get_boltz_availability()
        boltz_ok = boltz_info.available
        boltz_reason = boltz_info.reason or ""
    except Exception as exc:
        boltz_ok = False
        boltz_reason = str(exc)

    boltz_status = format_status_badge(boltz_ok, is_gpu=True, reason=boltz_reason)
    table.add_row(
        "boltz-2",
        "v2.3.0",
        boltz_status,
        "log IC50 / binder%",
        "[#48968C]AI 亲和力复核 (全复合物结构与 log IC50)[/#48968C]",
    )

    return make_pane(
        "对接引擎与算力探测",
        table,
        subtitle=f"CPU 核心数: {cpu_count} · 数据目录: {settings.data_dir.name}",
    )


def run_engines_command(include_header: bool = True) -> int:
    """Inspect and display engines, hardware, and GPU availability."""
    gpu = detect_gpu_info()
    if include_header:
        console.print(make_header(gpu_info=gpu))
    console.print(build_engines_pane())
    return 0
