"""Render a laid-out architecture diagram to standalone SVG."""
from __future__ import annotations

import re
from xml.sax.saxutils import escape

from . import icons as ic
from .layout import FONT, HEADER, ICON, NOTE_FONT, NOTE_LINE, card_height, layout, text_width
from .spec import Block, Diagram, SpecError

STYLE = {  # kind -> (fill, stroke, label colour)
    "module": ("#F3F6FA", "#41566B", ic.INK),
    "data": ("#FFFFFF", "#9AA8B6", ic.INK),
}
NESTED_FILLS = ("#F3F6FA", "#E7EDF5", "#DCE5F0")
RADIUS = 12.0
CORNER = 14.0
_REF = re.compile(r"^\s*([^.\s]+)\.(left|right|top|bottom|cx|cy)\s*(?:([+-])\s*([0-9.]+))?\s*$")


def _f(value: float) -> str:
    return f"{value:.1f}".rstrip("0").rstrip(".")


def _anchor(block: Block, side: str) -> tuple[float, float]:
    return {"left": (block.x, block.y + block.h / 2), "right": (block.x + block.w, block.y + block.h / 2),
            "top": (block.x + block.w / 2, block.y), "bottom": (block.x + block.w / 2, block.y + block.h)}[side]


def _resolve(value, blocks: dict[str, Block]) -> float:
    """A waypoint coordinate: a number, or 'block.anchor[+-offset]'."""
    if isinstance(value, (int, float)):
        return float(value)
    match = _REF.match(str(value))
    if not match or match.group(1) not in blocks:
        raise SpecError(f"cannot resolve waypoint coordinate {value!r}")
    b = blocks[match.group(1)]
    base = {"left": b.x, "right": b.x + b.w, "top": b.y, "bottom": b.y + b.h,
            "cx": b.x + b.w / 2, "cy": b.y + b.h / 2}[match.group(2)]
    offset = float(match.group(4) or 0.0)
    return base - offset if match.group(3) == "-" else base + offset


def _route(edge, blocks) -> list[tuple[float, float]]:
    start = _anchor(blocks[edge.source], edge.source_side)
    end = _anchor(blocks[edge.target], edge.target_side)
    points = [start]
    cursor = start
    for waypoint in edge.route:  # each waypoint may fix x, y or both; the other keeps the cursor's
        x = _resolve(waypoint["x"], blocks) if "x" in waypoint else cursor[0]
        y = _resolve(waypoint["y"], blocks) if "y" in waypoint else cursor[1]
        cursor = (x, y)
        points.append(cursor)
    horizontal_end = edge.target_side in ("left", "right")
    if not edge.route:
        if edge.source_side in ("left", "right") and horizontal_end:
            mid = (start[0] + end[0]) / 2
            points += [(mid, start[1]), (mid, end[1])]
        elif not horizontal_end and edge.source_side in ("top", "bottom"):
            mid = (start[1] + end[1]) / 2
            points += [(start[0], mid), (end[0], mid)]
        else:
            points.append((end[0], start[1]) if horizontal_end is False else (start[0], end[1]))
    else:  # close the route orthogonally into the target side
        points.append((cursor[0], end[1]) if horizontal_end else (end[0], cursor[1]))
    points.append(end)
    clean = [points[0]]
    for p in points[1:]:
        if abs(p[0] - clean[-1][0]) > 1e-6 or abs(p[1] - clean[-1][1]) > 1e-6:
            clean.append(p)
    return clean


def _rounded_path(points: list[tuple[float, float]]) -> str:
    """Orthogonal polyline with rounded corners."""
    d = [f"M{_f(points[0][0])} {_f(points[0][1])}"]
    for i in range(1, len(points) - 1):
        (x0, y0), (x1, y1), (x2, y2) = points[i - 1], points[i], points[i + 1]
        len_in, len_out = abs(x1 - x0) + abs(y1 - y0), abs(x2 - x1) + abs(y2 - y1)
        r = min(CORNER, len_in / 2, len_out / 2)
        ux, uy = ((x1 - x0) / len_in, (y1 - y0) / len_in) if len_in else (0, 0)
        vx, vy = ((x2 - x1) / len_out, (y2 - y1) / len_out) if len_out else (0, 0)
        d.append(f"L{_f(x1 - ux * r)} {_f(y1 - uy * r)}")
        d.append(f"Q{_f(x1)} {_f(y1)} {_f(x1 + vx * r)} {_f(y1 + vy * r)}")
    d.append(f"L{_f(points[-1][0])} {_f(points[-1][1])}")
    return " ".join(d)


def _text(x, y, text, size, colour=ic.INK, anchor="middle", weight="400", mono=False, math=True) -> str:
    family = "Consolas, 'Cascadia Mono', monospace" if mono else "'Segoe UI', Helvetica, Arial, sans-serif"
    return (f'<text x="{_f(x)}" y="{_f(y)}" font-family="{family}" font-size="{_f(size)}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{colour}">{_math(text) if math else escape(text)}</text>')


def _math(text: str) -> str:
    """Escape text; 'B_ψ', 'ψ_j' or 'x_{ab}' become subscripts, identifiers like 'krom_extract' stay literal."""
    out, i = [], 0
    for match in re.finditer(r"_\{([^}]*)\}|(?<=[A-Za-zͰ-Ͽ])_([A-Za-z0-9Ͱ-Ͽ])(?![A-Za-z0-9_])", text):
        out.append(escape(text[i:match.start()]))
        sub = match.group(1) if match.group(1) is not None else match.group(2)
        out.append(f'<tspan baseline-shift="sub" font-size="75%">{escape(sub)}</tspan>')
        i = match.end()
    out.append(escape(text[i:]))
    return "".join(out)


def _card(block: Block, top: float, lines: list[str], fill: str, stroke: str, ink: str) -> str:
    if not lines:
        return ""
    out = [f'<rect x="{_f(block.x + 8)}" y="{_f(top)}" width="{_f(block.w - 16)}" height="{_f(card_height(lines) - 4)}" '
           f'rx="5" fill="{fill}" stroke="{stroke}" stroke-width="1" stroke-dasharray="3 2"/>']
    for i, line in enumerate(lines):
        out.append(_text(block.x + 16, top + 15 + NOTE_LINE * i, line, NOTE_FONT, ink, "start", mono=True, math=False))
    return "".join(out)


def _notes(block: Block, top: float) -> str:
    """Settings (gold card), then audit flags (red card)."""
    return (_card(block, top, block.notes, "#FFF8E6", ic.GOLD, "#5A4A1A")
            + _card(block, top + card_height(block.notes), block.flags, "#FDECEC", ic.RED, "#8A1C1F"))


def _block(block: Block, depth: int) -> str:
    out = []
    if block.kind == "group":
        for child in block.children:
            out.append(_block(child, depth))
        return "".join(out)
    fill, stroke, colour = STYLE[block.kind]
    if block.container and block.kind == "module":
        fill = NESTED_FILLS[min(depth, len(NESTED_FILLS) - 1)]
    dash = ' stroke-dasharray="6 4"' if block.kind == "data" and block.container else ""
    out.append(f'<rect x="{_f(block.x)}" y="{_f(block.y)}" width="{_f(block.w)}" height="{_f(block.h)}" '
               f'rx="{RADIUS}" fill="{fill}" stroke="{stroke}" stroke-width="{2.2 if block.kind == "module" else 1.6}"'
               f'{dash} filter="url(#shadow)"/>')
    if block.container:
        x = block.x + 16
        if block.icon:
            out.append(ic.icon_markup(block.icon, x, block.y + 7, 30))
            x += 38
        out.append(_text(x, block.y + 28, block.label, FONT + 3, colour, "start", "700"))
        out.append(_notes(block, block.y + HEADER))
        for child in block.children:
            out.append(_block(child, depth + 1))
        return "".join(out)
    y = block.y + 14
    if block.icon:
        out.append(ic.icon_markup(block.icon, block.x + (block.w - ICON) / 2, y, ICON))
        y += ICON
    weight = "700" if block.kind == "module" else "600"
    out.append(_text(block.x + block.w / 2, y + 22, block.label, FONT, colour, weight=weight))
    out.append(_notes(block, y + 34))
    return "".join(out)


def _edge(edge, blocks) -> str:
    points = _route(edge, blocks)
    dash = ' stroke-dasharray="7 5"' if edge.style in ("dashed", "flagged") else ""
    colour, marker = (ic.RED, "arrow-flagged") if edge.style == "flagged" else (ic.INK, "arrow")
    out = [f'<path d="{_rounded_path(points)}" fill="none" stroke="{colour}" stroke-width="2.4"{dash} '
           f'marker-end="url(#{marker})"/>']
    if edge.label:
        colour = ic.RED if edge.style == "flagged" else ic.MUTED
        segments = list(zip(points, points[1:]))
        horizontal = [s for s in segments if abs(s[0][1] - s[1][1]) < 1e-6]
        need = text_width(edge.label, NOTE_FONT + 1) + 12
        roomy = [s for s in horizontal if abs(s[1][0] - s[0][0]) >= need]
        if roomy or horizontal:  # centre on the longest horizontal stretch
            (x0, y0), (x1, _) = max(roomy or horizontal, key=lambda s: abs(s[1][0] - s[0][0]))
            out.append(_text((x0 + x1) / 2, y0 - 8, edge.label, NOTE_FONT + 1, colour))
        else:  # beside the longest vertical stretch
            (x0, y0), (_, y1) = max(segments, key=lambda s: abs(s[1][1] - s[0][1]))
            out.append(_text(x0 + 8, (y0 + y1) / 2 + 4, edge.label, NOTE_FONT + 1, colour, "start"))
    return "".join(out)


def render(diagram: Diagram) -> str:
    title_band = 46.0 if diagram.title else 0.0
    width, height = layout(diagram.root, title_band=title_band)
    defs = ('<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{ic.INK}"/></marker>'
            '<marker id="arrow-flagged" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{ic.RED}"/></marker>{ic.ICON_DEFS}'
            '<filter id="shadow" x="-5%" y="-5%" width="110%" height="115%"><feDropShadow dx="0" dy="1.5" '
            'stdDeviation="1.8" flood-color="#1B2733" flood-opacity="0.16"/></filter></defs>')
    for edge in diagram.edges:  # grow the canvas for edges routed outside the blocks
        for x, y in _route(edge, diagram.blocks):
            width, height = max(width, x + 40.0), max(height, y + 40.0)
    body = [_block(diagram.root, 0)] + [_edge(edge, diagram.blocks) for edge in diagram.edges]
    title = _text(40, 46, diagram.title, 24, ic.INK, "start", "700") if diagram.title else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{_f(width)}" height="{_f(height)}" '
            f'viewBox="0 0 {_f(width)} {_f(height)}">{defs}<rect width="100%" height="100%" fill="#FFFFFF"/>'
            f'{title}{"".join(body)}</svg>\n')
