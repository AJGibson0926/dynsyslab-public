"""DynSysLab architecture diagrams: a JSON spec of nested blocks rendered to SVG.

    from architecture_diagrams import render_file
    render_file("docs/diagrams/pipeline-overview.diagram.json", "out.svg")
"""
from __future__ import annotations

from pathlib import Path

from .icons import ICONS
from .render import render
from .spec import SCHEMA, Diagram, SpecError, load, parse

__all__ = ["ICONS", "SCHEMA", "Diagram", "SpecError", "load", "parse", "render", "render_file"]


def render_file(spec_path: str | Path, svg_path: str | Path) -> Path:
    diagram = load(spec_path, set(ICONS))
    out = Path(svg_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(diagram), encoding="utf-8", newline="\n")
    return out
