"""Nested block layout: measure every block bottom-up, then place it top-down."""
from __future__ import annotations

from .spec import Block

FONT = 15.0           # block label size
NOTE_FONT = 11.5      # settings/notes size
ICON = 104.0          # icon edge length for leaf blocks
LEAF_MIN_WIDTH = 168.0
PAD = 22.0            # container padding
HEADER = 44.0         # container title band
GAP = 64.0            # default spacing between siblings (room for arrows)
NOTE_LINE = 15.0


def text_width(text: str, size: float, mono: bool = False) -> float:
    """Approximate rendered width; good enough for box sizing without a font engine."""
    return len(text) * size * (0.61 if mono else 0.56)


def card_height(lines: list[str]) -> float:
    return 10.0 + NOTE_LINE * len(lines) if lines else 0.0


def _notes_height(block: Block) -> float:
    """Settings card plus audit-flag card."""
    return card_height(block.notes) + card_height(block.flags)


def _notes_width(block: Block) -> float:
    return max((text_width(line, NOTE_FONT, True) for line in block.notes + block.flags), default=0.0) + 20.0


def measure(block: Block) -> None:
    if not block.container:
        label = text_width(block.label, FONT) + 24.0
        block.w = max(LEAF_MIN_WIDTH, label, _notes_width(block))
        block.h = 14.0 + (ICON if block.icon else 0.0) + 34.0 + _notes_height(block)
        return
    for child in block.children:
        measure(child)
    gap = GAP if block.gap is None else float(block.gap)
    spans = [(c.w, c.h) for c in block.children]
    if block.layout == "row":
        inner_w = sum(w for w, _ in spans) + gap * (len(spans) - 1)
        inner_h = max(h for _, h in spans)
    else:
        inner_w = max(w for w, _ in spans)
        inner_h = sum(h for _, h in spans) + gap * (len(spans) - 1)
    framed = block.kind != "group"
    pad = PAD if framed else 0.0
    header = (HEADER if block.label or block.icon else PAD) if framed else 0.0
    block.w = max(inner_w + 2 * pad, text_width(block.label, FONT + 3) + 90.0 if framed else 0.0, _notes_width(block))
    block.h = inner_h + header + pad + _notes_height(block)
    block._inner = (inner_w, inner_h, gap, pad, header)


def place(block: Block, x: float, y: float) -> None:
    block.x, block.y = x, y
    if not block.container:
        return
    inner_w, inner_h, gap, pad, header = block._inner
    left = x + (block.w - inner_w) / 2.0
    top = y + header + _notes_height(block)
    cursor = left if block.layout == "row" else top
    for child in block.children:
        if block.layout == "row":
            offset = {"start": 0.0, "center": (inner_h - child.h) / 2.0, "end": inner_h - child.h}[block.align]
            place(child, cursor, top + offset)
            cursor += child.w + gap
        else:
            offset = {"start": 0.0, "center": (inner_w - child.w) / 2.0, "end": inner_w - child.w}[block.align]
            place(child, left + offset, cursor)
            cursor += child.h + gap


def layout(root: Block, margin: float = 40.0, title_band: float = 0.0) -> tuple[float, float]:
    measure(root)
    place(root, margin, margin + title_band)
    return root.w + 2 * margin, root.h + 2 * margin + title_band
