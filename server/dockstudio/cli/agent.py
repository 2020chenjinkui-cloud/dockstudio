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
from rich.console import Group, RenderableType
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


def _build_help_pane() -> Any:
    """Construct the formatted command guide card matching Claude Code style."""
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

    return make_pane(
        "DockStudio 终端命令指南",
        Group(table, Text(), hints),
        subtitle="输入 / 即可呼出命令菜单 · 按 [Tab] 自动补全",
    )


def _print_help() -> None:
    """Print the formatted command guide card to stdout (backward compatible)."""
    console.print(_build_help_pane())
    console.print()


def _build_welcome_pane(gpu_info: str | None = None) -> Any:
    """Construct the initial welcome dashboard view matching Open-ClaudeCode."""
    body = Text()
    body.append("欢迎使用 ", style="#999999")
    body.append("DockStudio Terminal Agent", style="bold #D77757")
    body.append("。专为无头 Linux 服务器与 SSH 终端设计的交互式分子对接工作台。\n\n", style="#999999")

    table = Table(box=None, padding=(0, 2), show_header=False)
    table.add_column("Command", style="bold #D77757", width=12)
    table.add_column("Action", style="#FFFFFF", width=26)
    table.add_column("Description", style="#999999")

    table.add_row("/dock", "单分子高精度对接向导", "受体 → 口袋定盒 → 配体 → Vina/GNINA 对接 → 质检")
    table.add_row("/screen", "高通量虚拟筛选向导", "GPU Uni-Dock (1000x) / 多核并行海量小分子库打分")
    table.add_row("/engines", "硬件算力与引擎探针", "探查 CPU、NVIDIA Tesla/GeForce GPU 与对接引擎就绪状态")
    table.add_row("/prep", "受体 / 配体预处理向导", "PDB/SDF 质子化、加氢、指派 Gasteiger 电荷与转 PDBQT")
    table.add_row("/verify", "外包交付物合规验收", "运行 T1~T5 自动化数值验收套件，复算重现性与数据链自洽")
    table.add_row("/audit", "事件审计与防篡改验证", "检验 inputs.sha256 指纹与 events.jsonl 不可篡改哈希链")
    table.add_row("/view", "3D 结构即用即焚查看", "启动轻量内网 Web 3D 快速预览或纯终端 ASCII 构象轮廓")
    table.add_row("/help", "命令手册与热键指南", "查看全部斜杠指令、参数格式与快捷操作热键")
    table.add_row("/exit", "退出 Agent 工作台", "安全终止后台子任务并返回宿主系统 Shell")

    hints = Text()
    hints.append("• 快速开始: ", style="#999999")
    hints.append("输入 ", style="#999999")
    hints.append("/dock --demo", style="bold #D77757")
    hints.append(" 可直接在 1HSG 蛋白酶上运行示例对接；输入 ", style="#999999")
    hints.append("/engines", style="bold #D77757")
    hints.append(" 可检测 GPU 加速状态。", style="#999999")

    return make_pane(
        "交互式分子对接智能体工作台",
        Group(body, table, Text(), hints),
        subtitle="输入 / 即可呼出命令菜单 · 按 [Tab] 自动补全 · 界面实时就地刷新",
    )


class AgentWorkspace:
    """Manages the interactive terminal workspace canvas for DockStudio Agent.

    Provides:
    - Alternate Screen Buffer (DEC 1049: \\033[?1049h / \\033[?1049l) matching Open-ClaudeCode
    - In-place viewport refresh: clears canvas, resets cursor, prints fixed Top Header,
      active View Card, optional flash notices, and bottom prompt area
    - Eliminates scrolling clutter, duplicate stacked headers, and stale terminal logs
    """

    def __init__(self, gpu_info: str | None = None, use_alt_screen: bool = True) -> None:
        self.gpu_info = gpu_info
        self.use_alt_screen = use_alt_screen and sys.stdout.isatty()
        self.in_alt_screen = False
        self.current_view: Any = None
        self.flash_notice: tuple[str, str] | None = None  # (message, style)

    def enter(self) -> None:
        """Enter alternate screen buffer and initialize canvas."""
        if self.use_alt_screen and not self.in_alt_screen:
            sys.stdout.write("\033[?1049h\033[2J\033[H")
            sys.stdout.flush()
            self.in_alt_screen = True

    def exit(self) -> None:
        """Exit alternate screen buffer and restore host terminal shell."""
        if self.in_alt_screen:
            sys.stdout.write("\033[?1049l")
            sys.stdout.flush()
            self.in_alt_screen = False

    def clear(self) -> None:
        """Clear the terminal screen and reset cursor."""
        if sys.stdout.isatty():
            sys.stdout.write("\033[2J\033[H")
            sys.stdout.flush()

    def set_view(self, view: Any) -> None:
        """Update active view and clear any previous flash notice."""
        self.current_view = view
        self.flash_notice = None

    def prepare_command_canvas(self, show_header: bool = True) -> None:
        """Clear canvas for an interactive/streaming command and optionally render top header."""
        self.current_view = None
        self.flash_notice = None
        if sys.stdout.isatty():
            self.clear()
            if show_header:
                console.print(make_header(version="2.0.5", gpu_info=self.gpu_info))

    def refresh(self) -> None:
        """Redraws the entire workspace canvas in place."""
        if not sys.stdout.isatty():
            if self.current_view is not None:
                console.print(self.current_view)
            return

        if self.current_view is not None:
            self.clear()
            console.print(make_header(version="2.0.5", gpu_info=self.gpu_info))
            console.print(self.current_view)

        if self.flash_notice:
            msg, style = self.flash_notice
            console.print(f"[{style}]{msg}[/{style}]")
            self.flash_notice = None


def _dispatch_command(
    line: str,
    store: Any = None,
    workspace: AgentWorkspace | None = None,
) -> int:
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
        if workspace:
            workspace.set_view(_build_help_pane())
            workspace.refresh()
            return 0
        _print_help()
        return 0

    if cmd == "/clear":
        if workspace:
            workspace.set_view(_build_welcome_pane(gpu_info=workspace.gpu_info))
            workspace.refresh()
            return 0
        console.clear()
        return 0

    if cmd == "/engines":
        from .commands.engines import build_engines_pane, run_engines_command

        if workspace:
            workspace.set_view(build_engines_pane())
            workspace.refresh()
            return 0
        return run_engines_command()

    if cmd == "/dock":
        from .commands.dock import run_dock_interactive

        demo_flag = "--demo" in args or "-d" in args
        if workspace:
            workspace.prepare_command_canvas(show_header=False)
        return run_dock_interactive(demo=demo_flag)

    if cmd == "/screen":
        from .commands.screen import run_screen_interactive

        demo_flag = "--demo" in args or "-d" in args
        if workspace:
            workspace.prepare_command_canvas(show_header=False)
        return run_screen_interactive(demo=demo_flag)

    if cmd == "/prep":
        from .commands.prep import run_prep_ligand_interactive, run_prep_receptor_interactive

        sub = args[0].lower() if args else ""
        if sub in ("receptor", "rec"):
            if workspace:
                workspace.prepare_command_canvas(show_header=True)
            run_prep_receptor_interactive()
            return 0
        if sub in ("ligand", "lig"):
            if workspace:
                workspace.prepare_command_canvas(show_header=True)
            run_prep_ligand_interactive()
            return 0
        from .ui import select_one

        if workspace:
            workspace.prepare_command_canvas(show_header=True)
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
        if workspace:
            workspace.prepare_command_canvas(show_header=True)
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
            if workspace:
                workspace.flash_notice = ("任务标识不能为空", "bold #FF6B80")
                workspace.refresh()
            else:
                console.print("[#FF6B80]任务标识不能为空[/#FF6B80]")
            return 2
        if workspace:
            workspace.prepare_command_canvas(show_header=False)
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

        if workspace:
            workspace.prepare_command_canvas(show_header=False)
        return run_view(target_path, port=port, ascii_mode=ascii_mode)

    if workspace:
        workspace.flash_notice = (f"未知命令: {raw_cmd} (输入 /help 查看帮助)", "bold #FF6B80")
        workspace.refresh()
    else:
        console.print(f"[#FF6B80]未知命令: {raw_cmd}[/#FF6B80] (输入 [bold #D77757]/help[/bold #D77757] 查看帮助)")
    return 1


def run_agent_repl(store: Any = None) -> int:
    """Runs the main interactive REPL agent loop."""
    settings = get_settings()
    gpu_info = detect_gpu_info()

    # 1. Non-interactive pipe/stream execution (e.g. cat commands.txt | dockstudio)
    if not sys.stdin.isatty():
        console.print()
        console.print(make_header(version="2.0.5", gpu_info=gpu_info))
        console.print(
            "[dim]输入 [bold #D77757]/help[/bold #D77757] 查看命令列表，"
            "输入 [bold #D77757]/exit[/bold #D77757] 退出交互环境[/dim]\n"
        )
        try:
            for raw_line in sys.stdin:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                rc = _dispatch_command(line, store=store, workspace=None)
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

    # 2. Interactive Workspace & History setup
    workspace = AgentWorkspace(gpu_info=gpu_info)
    workspace.enter()
    workspace.set_view(_build_welcome_pane(gpu_info=gpu_info))

    history_file = settings.data_dir / ".agent_history"
    try:
        history_file.parent.mkdir(parents=True, exist_ok=True)
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

    try:
        workspace.refresh()

        while True:
            try:
                user_input = session.prompt()
                last_ctrl_c = 0.0
            except KeyboardInterrupt:
                now = time.time()
                if now - last_ctrl_c < 1.5:
                    workspace.exit()
                    console.print("\n[dim]检测到连续 Ctrl+C，正在退出 DockStudio Agent...[/dim]")
                    cleanup_active_processes()
                    return 0
                last_ctrl_c = now
                workspace.flash_notice = ("^C (再次按 Ctrl+C 或输入 /exit 退出)", "dim")
                workspace.refresh()
                continue
            except EOFError:
                workspace.exit()
                console.print("\n[dim]检测到 EOF 输入，正在退出 DockStudio Agent...[/dim]")
                cleanup_active_processes()
                return 0

            if not user_input.strip():
                if sys.stdout.isatty():
                    sys.stdout.write("\033[1A\033[2K\r")
                    sys.stdout.flush()
                continue

            try:
                rc = _dispatch_command(user_input, store=store, workspace=workspace)
                if rc == 999:
                    workspace.exit()
                    console.print("[dim]感谢使用 DockStudio，再见！[/dim]")
                    cleanup_active_processes()
                    return 0
            except KeyboardInterrupt:
                workspace.flash_notice = ("当前操作已被中断", "dim")
                workspace.refresh()
                cleanup_active_processes()
            except Exception as err:  # noqa: BLE001
                workspace.flash_notice = (f"执行出错: {err}", "bold #FF6B80")
                workspace.refresh()
    finally:
        workspace.exit()
        cleanup_active_processes()

    return 0
