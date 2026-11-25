"""Versioned architecture-diagram spec: dynsyslab.architecture-diagram/v1.

A diagram is a tree of blocks plus a list of edges. Blocks nest: a module may
contain other modules and data blocks, laid out in a row or a column.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA = "dynsyslab.architecture-diagram/v1"
KINDS = ("module", "data", "group")
LAYOUTS = ("row", "column")
SIDES = ("left", "right", "top", "bottom")
EDGE_STYLES = ("solid", "dashed", "flagged")


class SpecError(ValueError):
    """The diagram spec violates the v1 contract."""


@dataclass
class Block:
    id: str
    kind: str = "module"
    label: str = ""
    icon: str | None = None
    notes: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)  # audit findings, drawn in red
    layout: str = "row"
    gap: float | None = None
    align: str = "center"
    children: list["Block"] = field(default_factory=list)
    # Filled in by layout.
    x: float = 0.0
    y: float = 0.0
    w: float = 0.0
    h: float = 0.0

    @property
    def container(self) -> bool:
        return bool(self.children)


@dataclass
class Edge:
    source: str
    target: str
    source_side: str = "right"
    target_side: str = "left"
    label: str = ""
    style: str = "solid"
    route: list[dict] = field(default_factory=list)


@dataclass
class Diagram:
    title: str
    root: Block
    edges: list[Edge]
    blocks: dict[str, Block]


def _notes(value, where: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, dict):
        return [f"{key}: {json.dumps(item) if not isinstance(item, str) else item}" for key, item in value.items()]
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return list(value)
    raise SpecError(f"{where}: notes must be a list of strings or an object")


def _block(raw: dict, blocks: dict[str, Block], icons: set[str], where: str) -> Block:
    if not isinstance(raw, dict):
        raise SpecError(f"{where}: a block must be an object")
    unknown = set(raw) - {"id", "kind", "label", "icon", "notes", "flags", "layout", "gap", "align", "children"}
    if unknown:
        raise SpecError(f"{where}: unknown block fields {sorted(unknown)}")
    ident = raw.get("id")
    if not isinstance(ident, str) or not ident or "." in ident:
        raise SpecError(f"{where}: id must be a nonempty string without '.'")
    if ident in blocks:
        raise SpecError(f"{where}: duplicate block id {ident!r}")
    kind = raw.get("kind", "module")
    if kind not in KINDS:
        raise SpecError(f"{ident}: kind must be one of {KINDS}")
    layout = raw.get("layout", "row")
    if layout not in LAYOUTS:
        raise SpecError(f"{ident}: layout must be one of {LAYOUTS}")
    align = raw.get("align", "center")
    if align not in ("start", "center", "end"):
        raise SpecError(f"{ident}: align must be start, center or end")
    icon = raw.get("icon")
    if icon is not None and icon not in icons:
        raise SpecError(f"{ident}: unknown icon {icon!r}; known: {sorted(icons)}")
    block = Block(ident, kind, str(raw.get("label", "")), icon, notes=_notes(raw.get("notes"), ident),
                  layout=layout, gap=raw.get("gap"), align=align)
    block.flags = _notes(raw.get("flags"), ident)
    blocks[ident] = block
    children = raw.get("children", [])
    if not isinstance(children, list):
        raise SpecError(f"{ident}: children must be a list")
    block.children = [_block(child, blocks, icons, f"{ident}.children[{i}]") for i, child in enumerate(children)]
    if kind == "group" and not block.children:
        raise SpecError(f"{ident}: a group needs children")
    return block


def _edge(raw: dict, blocks: dict[str, Block], where: str) -> Edge:
    if not isinstance(raw, dict):
        raise SpecError(f"{where}: an edge must be an object")
    unknown = set(raw) - {"from", "to", "from_side", "to_side", "label", "style", "route"}
    if unknown:
        raise SpecError(f"{where}: unknown edge fields {sorted(unknown)}")
    edge = Edge(raw.get("from"), raw.get("to"), raw.get("from_side", "right"), raw.get("to_side", "left"),
                str(raw.get("label", "")), raw.get("style", "solid"), raw.get("route", []))
    for ident in (edge.source, edge.target):
        if ident not in blocks:
            raise SpecError(f"{where}: unknown block {ident!r}")
    if edge.source_side not in SIDES or edge.target_side not in SIDES:
        raise SpecError(f"{where}: sides must be one of {SIDES}")
    if edge.style not in EDGE_STYLES:
        raise SpecError(f"{where}: style must be one of {EDGE_STYLES}")
    if not isinstance(edge.route, list) or not all(isinstance(p, dict) and set(p) <= {"x", "y"} for p in edge.route):
        raise SpecError(f"{where}: route must be a list of {{x, y}} waypoints")
    return edge


def parse(raw: dict, icons: set[str]) -> Diagram:
    if not isinstance(raw, dict) or raw.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA}")
    unknown = set(raw) - {"schema", "title", "description", "root", "edges"}
    if unknown:
        raise SpecError(f"unknown diagram fields {sorted(unknown)}")
    blocks: dict[str, Block] = {}
    root = _block(raw.get("root"), blocks, icons, "root")
    edges = [_edge(item, blocks, f"edges[{i}]") for i, item in enumerate(raw.get("edges", []))]
    return Diagram(str(raw.get("title", "")), root, edges, blocks)


def load(path: str | Path, icons: set[str]) -> Diagram:
    with open(path, encoding="utf-8") as handle:
        return parse(json.load(handle), icons)
