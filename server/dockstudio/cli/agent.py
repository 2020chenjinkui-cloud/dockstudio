"""DockStudio CLI Interactive Terminal Agent (REPL).

Provides:
- Claude Code Terracotta prompt (dockstudio ❯)
- Slash command auto-completion (/dock, /screen, /verify, /audit, etc.)
- Dynamic bottom status toolbar
- History management with InMemory fallback
- Robust signal handling & subprocess cleanup
- Non-interactive pipe/stream detection
"""

from __future__ import annotations

import shlex
import subprocess
import sys
import time
from collections.abc import Iterable
from pathlib import Path
from typing import Any, ClassVar

from prompt_toolkit import PromptSession
from prompt_toolkit.auto_suggest import AutoSuggestFromHistory
from prompt_toolkit.completion import CompleteEvent, Completer, Completion, PathCompleter
from prompt_toolkit.document import Document
from prompt_toolkit.formatted_text import HTML
from prompt_toolkit.history import FileHistory, InMemoryHistory
from rich.table import Table
from rich.text import Text

from ..config import get_settings
from .commands.engines import detect_gpu_info
from .ui import console, make_header, make_pane

# Global registry for active background subtasks/subprocesses
ACTIVE_PROCESSES: set[subprocess.Popen] = set()


def register_process(proc: subprocess.Popen) -> None:
    """Register a running child process for signal protection."""
    ACTIVE_PROCESSES.add(proc)


def unregister_process(proc: subprocess.Popen) -> None:
    """Unregister a completed child process."""
    ACTIVE_PROCESSES.discard(proc)


def cleanup_active_processes() -> None:
    """Safely terminate and reap all registered active child processes."""
    for proc in list(ACTIVE_PROCESSES):
        try:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    proc.kill()
        except Exception:  # noqa: BLE001
            pass
    ACTIVE_PROCESSES.clear()


class SlashCompleter(Completer):
    """Completer for DockStudio slash commands with parameter path delegation."""

    COMMANDS: ClassVar[list[tuple[str, str]]] = [
        ("/dock", "单分子高精度对接五步向导 (受体 → 口袋 → 配体 → 对接 → 质检)"),
        ("/screen", "大通量虚拟筛选向导 (GPU Uni-Dock / 多进程 Vina)"),
        ("/engines", "硬件算力与对接引擎可用性探针"),
        ("/prep", "受体 / 配体准备向导 (加氢、转 PDBQT、质子化)"),
        ("/verify", "运行 T1-T5 自动化交付物验收套件 (复算外包交付包)"),
        ("/audit", "查看不可篡改事件链与 SHA256 指纹"),
        ("/view", "即用即焚内网 Web 3D 快速预览 / 终端 ASCII 查看"),
        ("/help", "查看命令指南与热键提示"),
        ("/clear", "清空终端屏幕"),
        ("/exit", "安全退出交互式 Agent"),
        ("/quit", "安全退出交互式 Agent (同 /exit)"),
    ]

    def __init__(self) -> None:
        self.path_completer = PathCompleter(expanduser=True)

    def get_completions(
        self, document: Document, complete_event: CompleteEvent
    ) -> Iterable[Completion]:
        text = document.text_before_cursor.lstrip()

        # When typing arguments for commands that accept file paths
        if " " in text:
            parts = text.split(maxsplit=1)
            cmd = parts[0].lower()
            if cmd in ("/verify", "/view", "verify", "view"):
                arg_text = parts[1] if len(parts) > 1 else ""

                # Handle optional flags for /view (e.g. /view --ascii <path> or /view -a <path>)
                if cmd in ("/view", "view"):
                    if arg_text.startswith("-") and " " not in arg_text:
                        for opt, desc in [
                            ("--ascii", "终端 ASCII 预览模式"),
                            ("-a", "终端 ASCII 预览模式 (简写)"),
                        ]:
                            if opt.startswith(arg_text):
                                yield Completion(
                                    opt,
                                    start_position=-len(arg_text),
                                    display=opt,
                                    display_meta=desc,
                                )
                        return

                    for flag in ("--ascii", "-a"):
                        if arg_text.startswith(flag + " "):
                            arg_text = arg_text[len(flag) :].lstrip()
                            break

                sub_doc = Document(arg_text, len(arg_text))
                yield from self.path_completer.get_completions(sub_doc, complete_event)
            return

        # Complete slash command
        prefix = text.lower()
        for cmd, desc in self.COMMANDS:
            bare = cmd.lstrip("/")
            if cmd.startswith(prefix) or (not prefix.startswith("/") and bare.startswith(prefix)):
                yield Completion(
                    cmd,
                    start_position=-len(text),
                    display=cmd,
                    display_meta=desc,
                )


def _print_help() -> None:
    """Print the formatted command guide card matching Claude Code style."""
    table = Table(box=None, padding=(0, 2), show_header=False)
    table.add_column("Command", style="bold #D77757", width=12)
    table.add_column("Description", style="#FFFFFF")

    for cmd, desc in SlashCompleter.COMMANDS:
        if cmd == "/quit":
            continue
        table.add_row(cmd, desc)

    hints = Text()
    hints.append("• 快捷热键: ", style="#999999")
    hints.append("[Tab]", style="bold #D77757")
    hints.append(" 自动补全  ", style="#999999")
    hints.append("[↑/↓]", style="bold #D77757")
    hints.append(" 历史命令  ", style="#999999")
    hints.append("[Ctrl+C]", style="bold #D77757")
    hints.append(" 取消当前输入 / 双击退出", style="#999999")

    console.print(make_pane("DockStudio 终端命令指南", table, subtitle="输入 / 即可呼出命令菜单"))
    console.print(hints)
    console.print()


def _dispatch_command(line: str, store: Any = None) -> int:
    """Parse and dispatch a single command line. Returns exit code."""
    line = line.strip()
    if not line or line.startswith("#"):
        return 0

    try:
        parts = shlex.split(line, posix=(sys.platform != "win32"))
        if sys.platform == "win32":
            parts = [p.strip("\"'") for p in parts]
    except ValueError:
        parts = line.split()

    raw_cmd = parts[0]
    cmd = raw_cmd.lower()
    args = parts[1:]

    # Normalize slash vs bare commands
    if not cmd.startswith("/"):
        cmd = f"/{cmd}"

    if cmd in ("/exit", "/quit"):
        return 999  # Signal to terminate REPL loop

    if cmd == "/help":
        _print_help()
        return 0

    if cmd == "/clear":
        console.clear()
        return 0

    if cmd == "/engines":
        from .commands.engines import run_engines_command

        return run_engines_command()

    if cmd == "/dock":
        from .commands.dock import run_dock_interactive

        demo_flag = "--demo" in args or "-d" in args
        return run_dock_interactive(demo=demo_flag)

    if cmd == "/screen":
        from .commands.screen import run_screen_interactive

        return run_screen_interactive()

    if cmd == "/prep":
        from .commands.prep import run_prep_ligand_interactive, run_prep_receptor_interactive

        sub = args[0].lower() if args else ""
        if sub in ("receptor", "rec"):
            run_prep_receptor_interactive()
            return 0
        if sub in ("ligand", "lig"):
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

    if cmd == "/verify":
        from .commands.verify import run_verify
        from .ui import prompt_path

        path_str = args[0] if args else None
        if not path_str:
            target_path = prompt_path("请输入待核算交付包路径 (目录或 .zip)")
        else:
            target_path = Path(path_str).expanduser().resolve()
        return run_verify(target_path)

    if cmd == "/audit":
        from prompt_toolkit import prompt as pt_prompt

        from .commands.audit import run_audit

        target = args[0] if args else None
        if not target:
            try:
                target = pt_prompt("请输入任务 ID 或项目目录: ").strip()
            except (KeyboardInterrupt, EOFError):
                return 0
        if not target:
            console.print("[#FF6B80]任务标识不能为空[/#FF6B80]")
            return 2
        return run_audit(target)

    if cmd == "/view":
        from .commands.view import run_view
        from .ui import prompt_path

        ascii_mode = "--ascii" in args or "-a" in args
        port = 0
        target_path: Path | None = None

        clean_args = [a for a in args if a not in ("--ascii", "-a")]
        if clean_args:
            target_path = Path(clean_args[0]).expanduser().resolve()
        else:
            target_path = prompt_path("请输入结构文件路径 (PDB/PDBQT/SDF)")

        return run_view(target_path, port=port, ascii_mode=ascii_mode)

    console.print(f"[#FF6B80]未知命令: {raw_cmd}[/#FF6B80] (输入 [bold #D77757]/help[/bold #D77757] 查看帮助)")
    return 1


def run_agent_repl(store: Any = None) -> int:
    """Runs the main interactive REPL agent loop."""
    settings = get_settings()
    gpu_info = detect_gpu_info()

    # Welcome banner
    console.print()
    console.print(make_header(version="2.0.5", gpu_info=gpu_info))
    console.print(
        "[dim]输入 [bold #D77757]/help[/bold #D77757] 查看命令列表，"
        "输入 [bold #D77757]/exit[/bold #D77757] 退出交互环境[/dim]\n"
    )

    # 1. Non-interactive pipe/stream execution (e.g. cat commands.txt | dockstudio)
    if not sys.stdin.isatty():
        try:
            for raw_line in sys.stdin:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                rc = _dispatch_command(line, store=store)
                if rc == 999:
                    break
                if rc != 0:
                    cleanup_active_processes()
                    return rc
        except Exception as e:
            console.print(f"[#FF6B80]管道执行发生异常: {e}[/#FF6B80]")
            cleanup_active_processes()
            return 1
        cleanup_active_processes()
        return 0

    # 2. History setup
    history_file = settings.data_dir / ".agent_history"
    try:
        history_file.parent.mkdir(parents=True, exist_ok=True)
        # Verify writable before instantiating FileHistory
        with open(history_file, "a", encoding="utf-8"):
            pass
        history = FileHistory(str(history_file))
    except Exception:  # noqa: BLE001
        history = InMemoryHistory()

    # 3. Dynamic toolbar renderer
    mode_str = "GPU加速" if gpu_info else "CPU多核"

    def _bottom_toolbar() -> HTML:
        return HTML(
            f'<style bg="#2C323E" fg="#FFFFFF"><b> DockStudio </b></style> '
            f'<style fg="#999999">[?] /help 帮助  [/] 命令补全  [Tab] 补全  [Ctrl+C] 取消/退出</style> '
            f'<style fg="#48968C"> | SSH无头模式 ({mode_str})</style>'
        )

    # 4. PromptSession initialization with Claude Code Terracotta prompt
    prompt_str = HTML(
        '<style fg="#D77757"><b>dockstudio</b></style> <style fg="#D77757"><b>❯</b></style> '
    )
    session: PromptSession = PromptSession(
        message=prompt_str,
        completer=SlashCompleter(),
        history=history,
        bottom_toolbar=_bottom_toolbar,
        auto_suggest=AutoSuggestFromHistory(),
        enable_history_search=True,
    )

    last_ctrl_c: float = 0.0

    while True:
        try:
            user_input = session.prompt()
            last_ctrl_c = 0.0  # Reset Ctrl+C counter on successful input
        except KeyboardInterrupt:
            now = time.time()
            if now - last_ctrl_c < 1.5:
                console.print("\n[dim]检测到连续 Ctrl+C，正在退出 DockStudio Agent...[/dim]")
                cleanup_active_processes()
                return 0
            last_ctrl_c = now
            console.print("\n[dim]^C (再次按 Ctrl+C 或输入 /exit 退出)[/dim]")
            continue
        except EOFError:
            console.print("\n[dim]检测到 EOF 输入，正在退出 DockStudio Agent...[/dim]")
            cleanup_active_processes()
            return 0

        if not user_input.strip():
            continue

        try:
            rc = _dispatch_command(user_input, store=store)
            if rc == 999:
                console.print("[dim]感谢使用 DockStudio，再见！[/dim]")
                cleanup_active_processes()
                return 0
        except KeyboardInterrupt:
            console.print("\n[dim]当前操作已被中断[/dim]")
            cleanup_active_processes()
        except Exception as err:  # noqa: BLE001
            console.print(f"[bold #FF6B80]执行出错: {err}[/bold #FF6B80]")

    cleanup_active_processes()
    return 0
