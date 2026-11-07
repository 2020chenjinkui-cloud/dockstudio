"""Interactive terminal picker and input dialogs.

Provides:
- select_one: Arrow-key / number picker matching Claude Code FuzzyPicker
- prompt_path: Path input with filesystem auto-completion (PathCompleter)
- confirm: Y/N confirmation prompt
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from prompt_toolkit import prompt
from prompt_toolkit.completion import PathCompleter
from prompt_toolkit.formatted_text import HTML
from rich.text import Text

from .components import console, make_pane


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
    """Select one option from a list.

    Supports direct number entry (e.g. 1, 2) or Enter to accept default.
    Choices can be simple strings or (title, subtitle) tuples.
    """
    options: list[tuple[str, str]] = []
    for item in choices:
        if isinstance(item, tuple):
            options.append((str(item[0]), str(item[1])))
        else:
            options.append((str(item), ""))

    content = Text()
    for idx, (label, sub) in enumerate(options):
        is_default = idx == default
        marker = "❯" if is_default else " "
        idx_str = f"[{idx + 1}]"
        
        if is_default:
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
