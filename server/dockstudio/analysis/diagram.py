"""2D interaction diagram as SVG: RDKit 2D ligand drawing + residue labels
colored by interaction type (self-drawn; no extra plotting dependency)."""

from __future__ import annotations

import html
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw

from ..logging import get_logger
from ..models import Interaction

log = get_logger("analysis.diagram")

_TYPE_COLORS = {
    "hbond": "#2563eb",
    "salt_bridge": "#dc2626",
    "pi_stacking": "#9333ea",
    "hydrophobic": "#059669",
    "halogen": "#d97706",
    "water_bridge": "#0891b2",
}
_TYPE_LABELS_ZH = {
    "hbond": "氢键",
    "salt_bridge": "盐桥",
    "pi_stacking": "π 堆积",
    "hydrophobic": "疏水",
    "halogen": "卤键",
    "water_bridge": "水桥",
}

_LABEL_W = 78
_LABEL_H = 26
_LIGAND_W = 420
_LIGAND_H = 320


def _ligand_svg(pose_sdf: Path, width: int, height: int) -> str:
    """Render the pose's ligand as a 2D depiction (SVG inner content)."""
    mol = Chem.MolFromMolBlock(pose_sdf.read_text(encoding="utf-8"), sanitize=False, removeHs=True)
    if mol is None:
        raise ValueError(f"RDKit could not read pose {pose_sdf}")
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _place_labels(interactions: list[Interaction], cx: float, cy: float) -> list[str]:
    """Lay residue pills out on two arcs around the ligand, colored by type."""
    by_residue: dict[str, list[Interaction]] = {}
    for it in interactions:
        if it.residue:
            by_residue.setdefault(it.residue, []).append(it)
    groups = sorted(by_residue.items(), key=lambda kv: len(kv[1]), reverse=True)[:14]
    import math

    n = len(groups)
    parts: list[str] = []
    for i, (residue, its) in enumerate(groups):
        angle = -math.pi / 2 + (2 * math.pi * i / max(n, 1))
        rx, ry = 240.0, 165.0
        x = cx + rx * math.cos(angle) - _LABEL_W / 2
        y = cy + ry * math.sin(angle) - _LABEL_H / 2
        main_type = its[0].type
        color = _TYPE_COLORS.get(main_type, "#6b7280")
        line_x1, line_y1 = cx, cy
        line_x2, line_y2 = x + _LABEL_W / 2, y + _LABEL_H / 2
        parts.append(
            f'<line x1="{line_x1:.1f}" y1="{line_y1:.1f}" x2="{line_x2:.1f}" y2="{line_y2:.1f}"'
            f' stroke="{color}" stroke-width="1.2" stroke-dasharray="4 3" opacity="0.65"/>'
        )
        parts.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{_LABEL_W}" height="{_LABEL_H}" rx="13"'
            f' fill="{color}" fill-opacity="0.12" stroke="{color}" stroke-width="1.2"/>'
        )
        dist = its[0].distance
        detail = f"{dist:.1f}Å" if dist is not None else ""
        parts.append(
            f'<text x="{x + _LABEL_W / 2:.1f}" y="{y + 11:.1f}" text-anchor="middle"'
            f' font-size="9.5" font-family="monospace" fill="{color}" font-weight="bold">'
            f"{html.escape(residue)}</text>"
        )
        parts.append(
            f'<text x="{x + _LABEL_W / 2:.1f}" y="{y + 21:.1f}" text-anchor="middle"'
            f' font-size="8" font-family="monospace" fill="{color}" opacity="0.8">'
            f"{html.escape(_TYPE_LABELS_ZH.get(main_type, main_type))} {html.escape(detail)}</text>"
        )
    return parts


def generate_diagram_svg(
    pose_sdf: Path,
    interactions: list[Interaction],
    title: str = "",
    width: int = 920,
    height: int = 640,
) -> str:
    """Compose the full SVG: ligand drawing center-left, residue pills around it."""
    ligand_inner = _ligand_svg(pose_sdf, _LIGAND_W, _LIGAND_H)
    cx = 300.0
    cy = height / 2.0 + 20
    labels = _place_labels(interactions, cx, cy)
    used_types = sorted({i.type for i in interactions if i.residue})
    legend = []
    lx = 16
    for t in used_types:
        color = _TYPE_COLORS.get(t, "#6b7280")
        legend.append(
            f'<rect x="{lx}" y="{height - 44}" width="12" height="12" rx="3" fill="{color}"/>'
            f'<text x="{lx + 18}" y="{height - 34}" font-size="11px" '
            f'font-family="sans-serif" fill="#333">{html.escape(_TYPE_LABELS_ZH.get(t, t))}</text>'
        )
        lx += 90
    # RDKit's MolDraw2DSVG returns a full <svg>...</svg> document; strip its
    # header and keep only the drawing body, wrapped in a translated group.
    body = re.sub(r"<\?xml[^>]*\?>", "", ligand_inner)
    body = re.sub(r"<!DOCTYPE[^>]*>", "", body)
    body = re.sub(r"</?svg[^>]*>", "", body).strip()
    title_text = (
        f'<text x="{width / 2:.0f}" y="30" text-anchor="middle" font-size="16"'
        f' font-weight="bold" font-family="sans-serif" fill="#1f2328">'
        f"{html.escape(title)}</text>"
        if title
        else ""
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"'
        f' viewBox="0 0 {width} {height}" font-family="sans-serif">\n'
        f'<rect width="{width}" height="{height}" fill="white"/>\n'
        f"{title_text}\n"
        f'<g transform="translate({cx - _LIGAND_W / 2:.0f},{cy - _LIGAND_H / 2:.0f})">\n'
        f"{body}\n</g>\n" + "\n".join(labels) + "\n" + "\n".join(legend) + "\n</svg>"
    )
