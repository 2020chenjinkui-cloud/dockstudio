"""DockStudio CLI UI Theme & Color Palette.

Directly adapted from Anthropic Claude Code (Open-ClaudeCode darkTheme)
using explicit 24-bit TrueColor RGB values to prevent terminal palette corruption.
"""

from __future__ import annotations

from rich.style import Style
from rich.theme import Theme

# Explicit RGB hex colors matching Open-ClaudeCode darkTheme
COLOR_CLAUDE = "#D77757"  # Claude Terracotta / Primary brand accent
COLOR_SPINNER = "#93A5FF"  # Lavender blue for waiting state and spinners
COLOR_SPINNER_SHIMMER = "#B1C3FF"  # Shimmer effect
COLOR_PERMISSION = "#B1B9F9"  # Interactive prompt / permission highlight
COLOR_PLAN_MODE = "#48968C"  # Sage green for tasks / plan mode
COLOR_SUCCESS = "#4EBA65"  # Soft bright green for passed checks / good scores
COLOR_WARNING = "#FFC107"  # Bright amber for warnings / fallback
COLOR_ERROR = "#FF6B80"  # Coral red for errors / clashes / failures
COLOR_TEXT = "#FFFFFF"  # Pure white for main text
COLOR_INACTIVE = "#999999"  # Medium gray for secondary text
COLOR_SUBTLE = "#505050"  # Dark gray for subtle hints and dividers
COLOR_BORDER = "#888888"  # Medium gray for panel / card borders
COLOR_CARD_BG = "#2C323E"  # Cool dark blue-gray for card selection / action background
COLOR_SELECTION_BG = "#264F78"  # VS Code dark selection blue
COLOR_USER_MSG_BG = "#373737"  # Neutral gray background

# Semantic symbols
ICON_SUCCESS = "[bold #4EBA65]✓[/bold #4EBA65]"
ICON_WARNING = "[bold #FFC107]▲[/bold #FFC107]"
ICON_ERROR = "[bold #FF6B80]✗[/bold #FF6B80]"
ICON_RUNNING = "[bold #93A5FF]●[/bold #93A5FF]"
ICON_UNAVAILABLE = "[#505050]○[/#505050]"
ICON_BULLET = "[#D77757]❯[/#D77757]"

DOCKSTUDIO_THEME = Theme({
    # Brand & Panels
    "brand": Style(color=COLOR_CLAUDE, bold=True),
    "panel.border": Style(color=COLOR_BORDER),
    "panel.title": Style(color=COLOR_CLAUDE, bold=True),
    "panel.subtitle": Style(color=COLOR_INACTIVE),
    "spinner": Style(color=COLOR_SPINNER, bold=True),
    "prompt.border": Style(color=COLOR_BORDER),
    # State Semantics
    "success": Style(color=COLOR_SUCCESS, bold=True),
    "warning": Style(color=COLOR_WARNING, bold=True),
    "error": Style(color=COLOR_ERROR, bold=True),
    "muted": Style(color=COLOR_INACTIVE),
    "dim": Style(color=COLOR_SUBTLE),
    "highlight": Style(color=COLOR_PERMISSION, bold=True),
    "plan": Style(color=COLOR_PLAN_MODE, bold=True),
    # Scientific & Affinity scoring
    "energy.excellent": Style(color=COLOR_SUCCESS, bold=True),  # < -9.0 kcal/mol
    "energy.moderate": Style(color=COLOR_WARNING, bold=True),  # -7.0 ~ -9.0 kcal/mol
    "energy.poor": Style(color=COLOR_ERROR),  # > -7.0 kcal/mol
    # Hardware badges
    "badge.gpu": Style(color="#FFFFFF", bgcolor="#264F78", bold=True),
    "badge.cpu": Style(color="#FFFFFF", bgcolor="#373737", bold=True),
})
