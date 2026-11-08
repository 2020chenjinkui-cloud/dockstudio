"""Interactive terminal picker and input dialogs.

Provides:
- select_one: Arrow-key / number picker matching Claude Code FuzzyPicker
- prompt_path: Path input with filesystem auto-completion (PathCompleter)
- confirm: Y/N confirmation prompt
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from prompt_toolkit import prompt
from prompt_toolkit.completion import PathCompleter
from prompt_toolkit.formatted_text import HTML
from rich.text import Text

from .components import console, make_pane


def _read_single_key() -> str:
    """Read a single key or arrow sequence portably across Windows & POSIX."""
    if sys.platform == "win32":
        import msvcrt

        ch = msvcrt.getch()
        if ch in (b"\x00", b"\xe0"):
            ch2 = msvcrt.getch()
            if ch2 == b"H":
                return "up"
            if ch2 == b"P":
                return "down"
            return ""
        if ch in (b"\r", b"\n"):
            return "enter"
        if ch == b"\x03":
            raise KeyboardInterrupt
        if ch == b"\x1b":
            return "escape"
        try:
            return ch.decode("utf-8", errors="ignore")
        except Exception:
            return ""
    else:
        import select
        import termios
        import tty

        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setcbreak(fd)
            ch = sys.stdin.read(1)
            if ch == "\x1b":
                if select.select([sys.stdin], [], [], 0.05)[0]:
                    ch2 = sys.stdin.read(1)
                    if ch2 in ("[", "O") and select.select([sys.stdin], [], [], 0.05)[0]:
                        ch3 = sys.stdin.read(1)
                        seq = ch + ch2 + ch3
                        if seq in ("\x1b[A", "\x1bOA"):
                            return "up"
                        if seq in ("\x1b[B", "\x1bOB"):
                            return "down"
                return "escape"
            if ch in ("\r", "\n"):
                return "enter"
            if ch == "\x03":
                raise KeyboardInterrupt
            return ch
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)


def prompt_path(
    message: str = "请输入文件路径",
    default: str = "",
    must_exist: bool = True,
) -> Path:
    """Prompt user for a file or directory path with Tab auto-completion."""
    completer = PathCompleter(expanduser=True)
    while True:
        try:
            val = prompt(
                HTML(f"<ansibrightblack>[</ansibrightblack><ansiyellow>📁</ansiyellow><ansibrightblack>]</ansibrightblack> <ansicyan>{message}:</ansicyan> "),
                default=default,
                completer=completer,
            ).strip()
        except (KeyboardInterrupt, EOFError):
            console.print("\n[dim]操作已取消[/dim]")
            raise

        if not val:
            if default:
                val = default
            else:
                console.print("[#FF6B80]路径不能为空，请重新输入[/#FF6B80]")
                continue

        p = Path(val).expanduser().resolve()
        if must_exist and not p.exists():
            console.print(f"[#FF6B80]文件不存在: {p}，请重新输入[/#FF6B80]")
            continue
        return p


def confirm(message: str, default: bool = True) -> bool:
    """Prompt user for Y/n confirmation."""
    suffix = "[Y/n]" if default else "[y/N]"
    try:
        ans = prompt(
            HTML(f"<ansiyellow>?</ansiyellow> <b>{message}</b> <ansibrightblack>{suffix}:</ansibrightblack> ")
        ).strip().lower()
    except (KeyboardInterrupt, EOFError):
        console.print("\n[dim]操作已取消[/dim]")
        raise

    if not ans:
        return default
    return ans in ("y", "yes", "true", "1")


def select_one(
    title: str,
    choices: Sequence[str | tuple[str, str]],
    default: int = 0,
) -> int:
    """Select one option from a list with arrow-key in-place selection matching Claude Code.

    Supports:
    - [↑/↓] or [k/j] to navigate in place without terminal spam
    - [1-9] direct numeric jump
    - [Enter] or [Space] to confirm selection
    - Automatic in-place erasure of the menu upon confirmation, leaving a single concise summary:
      ✓ {title}: {chosen_label}
    - Graceful fallback to prompt() if stdin/stdout is non-interactive or piped.
    """
    options: list[tuple[str, str]] = []
    for item in choices:
        if isinstance(item, tuple):
            options.append((str(item[0]), str(item[1])))
        else:
            options.append((str(item), ""))

    if not options:
        return 0

    current_idx = max(0, min(default, len(options) - 1))

    # Check if we can run interactive live arrow-key picker
    is_interactive = sys.stdin.isatty() and sys.stdout.isatty()
    if is_interactive:
        try:
            from rich.live import Live

            def _build_pane(sel: int) -> Any:
                content = Text()
                for idx, (label, sub) in enumerate(options):
                    is_sel = idx == sel
                    marker = "❯" if is_sel else " "
                    idx_str = f"[{idx + 1}]"
                    if is_sel:
                        content.append(f" {marker} ", style="bold #D77757")
                        content.append(f"{idx_str:4s} ", style="bold #D77757")
                        content.append(f"{label}\n", style="bold #FFFFFF")
                        if sub:
                            content.append(f"       {sub}\n", style="#999999")
                    else:
                        content.append(f" {marker} ", style="#505050")
                        content.append(f"{idx_str:4s} ", style="#999999")
                        content.append(f"{label}\n", style="#CCCCCC")
                        if sub:
                            content.append(f"       {sub}\n", style="#505050")
                return make_pane(
                    title,
                    content,
                    subtitle="[↑/↓] 移动光标 · [Enter] 确认 · [1-9] 直达",
                )

            with Live(
                _build_pane(current_idx),
                console=console,
                auto_refresh=False,
                transient=True,
            ) as live:
                while True:
                    live.update(_build_pane(current_idx))
                    live.refresh()
                    key = _read_single_key()
                    if key in ("up", "k"):
                        current_idx = (current_idx - 1) % len(options)
                    elif key in ("down", "j"):
                        current_idx = (current_idx + 1) % len(options)
                    elif key in ("enter", " "):
                        break
                    elif key.isdigit():
                        num = int(key)
                        if 1 <= num <= len(options):
                            current_idx = num - 1
                            break
                    elif key in ("escape", "q"):
                        raise KeyboardInterrupt

            chosen_label = options[current_idx][0]
            console.print(f"[bold #4EBA65]✓[/bold #4EBA65] {title} [bold #FFFFFF]{chosen_label}[/bold #FFFFFF]")
            return current_idx
        except (KeyboardInterrupt, EOFError):
            console.print("\n[dim]操作已取消[/dim]")
            raise
        except Exception:
            # Fallback to standard prompt below if raw keyboard mode is unavailable
            pass

    # Non-interactive / pipe / fallback path
    content = Text()
    for idx, (label, sub) in enumerate(options):
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
            choice_str = prompt(
                HTML(f"<ansibrightblack>请选择 [1-{len(options)}] (回车默认={default + 1}):</ansibrightblack> ")
            ).strip()
        except (KeyboardInterrupt, EOFError):
            console.print("\n[dim]操作已取消[/dim]")
            raise

        if not choice_str:
            return default

        try:
            num = int(choice_str)
            if 1 <= num <= len(options):
                return num - 1
            console.print(f"[#FF6B80]输入超出范围，请输入 1 到 {len(options)}[/#FF6B80]")
        except ValueError:
            console.print(f"[#FF6B80]无效输入 '{choice_str}'，请输入数字[/#FF6B80]")
