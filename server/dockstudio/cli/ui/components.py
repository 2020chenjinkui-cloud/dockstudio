"""Terminal UI components modeled after Claude Code's Design System.

Provides Pane, Dialog, Header, KeyHint, and scientific formatting helpers.
"""

from __future__ import annotations

import platform
import shutil

from rich.box import ROUNDED
from rich.console import Console, RenderableType
from rich.panel import Panel
from rich.text import Text

from .theme import (
    COLOR_PERMISSION,
    DOCKSTUDIO_THEME,
    ICON_ERROR,
    ICON_RUNNING,
    ICON_SUCCESS,
    ICON_UNAVAILABLE,
    ICON_WARNING,
)

# Global rich console initialized with Claude Code theme
console = Console(theme=DOCKSTUDIO_THEME, highlight=False)


def get_terminal_width(default: int = 90) -> int:
    try:
        width = shutil.get_terminal_size((default, 24)).columns
        return max(60, min(width, 120))
    except Exception:  # noqa: BLE001
        return default


def make_header(version: str = "2.0.5", gpu_info: str | None = None) -> Panel:
    """Claude Code style top title banner."""
    title_text = Text()
    title_text.append(" DockStudio ", style="bold #D77757")
    title_text.append("Terminal Agent ", style="bold #FFFFFF")
    title_text.append(f"v{version}", style="#999999")

    content = Text()
    content.append("• 平台: ", style="#999999")
    content.append(f"{platform.system()} {platform.machine()}  ", style="#FFFFFF")

    if gpu_info:
        content.append("• GPU: ", style="#999999")
        content.append(f"{gpu_info}  ", style="bold #93A5FF")
    else:
        content.append("• 计算: ", style="#999999")
        content.append("CPU 多核  ", style="#FFFFFF")

    content.append("• 模式: ", style="#999999")
    content.append("无头服务器/SSH终端", style="#48968C")

    return Panel(
        content,
        title=title_text,
        title_align="left",
        border_style="#888888",
        box=ROUNDED,
        padding=(0, 1),
    )


def make_pane(
    title: str,
    content: RenderableType,
    subtitle: str | None = None,
    border_style: str = "#888888",
    title_style: str = "bold #D77757",
) -> Panel:
    """Claude Code <Pane> component: Bounded region with top border title and padding."""
    formatted_title = f"[{title_style}] {title} [/{title_style}]"
    formatted_subtitle = f"[#505050]{subtitle}[/#505050]" if subtitle else None

    return Panel(
        content,
        title=formatted_title,
        title_align="left",
        subtitle=formatted_subtitle,
        subtitle_align="right",
        border_style=border_style,
        box=ROUNDED,
        padding=(0, 2),
    )


def make_dialog(
    title: str,
    content: RenderableType,
    action_hint: str = "[Enter 确认]  [Esc 取消]",
    color: str = COLOR_PERMISSION,
) -> Panel:
    """Claude Code <Dialog> component: Modal card with keyboard action hint."""
    return Panel(
        content,
        title=f"[bold {color}] {title} [/bold {color}]",
        subtitle=f"[dim]{action_hint}[/dim]",
        subtitle_align="right",
        border_style=color,
        box=ROUNDED,
        padding=(1, 2),
    )


def format_status_badge(available: bool, is_gpu: bool = False, reason: str = "") -> str:
    """Renders a status badge matching Claude Code style."""
    if available:
        if is_gpu:
            return f"{ICON_RUNNING} [bold #93A5FF]可用 (GPU)[/bold #93A5FF]"
        return f"{ICON_SUCCESS} [#4EBA65]可用 (CPU)[/#4EBA65]"
    msg = f" ({reason})" if reason else ""
    return f"{ICON_UNAVAILABLE} [#999999]不可用{msg}[/#999999]"


def format_affinity(affinity: float) -> str:
    """Color-code docking affinity (kcal/mol)."""
    if affinity <= -9.0:
        return f"[bold #4EBA65]{affinity:6.2f}[/bold #4EBA65]"
    if affinity <= -7.0:
        return f"[bold #FFC107]{affinity:6.2f}[/bold #FFC107]"
    return f"[#FF6B80]{affinity:6.2f}[/#FF6B80]"


def format_posecheck(passed_count: int, total_count: int = 18, details: str = "") -> str:
    """Color-code PoseBusters quality score."""
    if passed_count >= total_count:
        return f"{ICON_SUCCESS} [#4EBA65]{passed_count}/{total_count} 全部合规[/#4EBA65]"
    if passed_count >= total_count - 2:
        extra = f" ({details})" if details else ""
        return f"{ICON_WARNING} [#FFC107]{passed_count}/{total_count} 基本合格{extra}[/#FFC107]"
    extra = f" ({details})" if details else ""
    return f"{ICON_ERROR} [#FF6B80]{passed_count}/{total_count} 异常{extra}[/#FF6B80]"


def make_key_hint(hints: list[tuple[str, str]]) -> Text:
    """Renders keyboard shortcut hints like: [Enter] 确认  [Esc] 返回."""
    line = Text()
    for key, desc in hints:
        line.append(f"[{key}]", style="bold #D77757")
        line.append(f" {desc}  ", style="#999999")
    return line
