"""Starting specs for new diagrams. Each renders as-is; edit labels, notes and steps from there."""
from __future__ import annotations

from .spec import SCHEMA


def _module(title: str) -> dict:
    return {
        "schema": SCHEMA, "title": title,
        "description": "Cite the source files and the job whose settings the notes repeat.",
        "root": {"id": "page", "kind": "group", "layout": "row", "gap": 70, "children": [
            {"id": "input", "kind": "data", "label": "Input archive", "icon": "archive", "notes": ["arrays it reads"]},
            {"id": "module", "kind": "module", "label": "Module", "icon": "monitor", "layout": "row", "gap": 60,
             "flags": ["audit findings go here"],
             "children": [
                 {"id": "step-1", "kind": "module", "label": "Step 1", "icon": "time_series", "notes": {"setting": "value"}},
                 {"id": "step-2", "kind": "module", "label": "Step 2", "icon": "matrix", "notes": {"setting": "value"}},
                 {"id": "step-3", "kind": "data", "label": "Result", "icon": "heatmap"}]},
            {"id": "output", "kind": "data", "label": "Output archive", "icon": "archive", "notes": ["arrays it writes"]},
            {"id": "next", "kind": "module", "label": "Next module", "icon": "monitor"}]},
        "edges": [{"from": "input", "to": "module"}, {"from": "step-1", "to": "step-2"},
                  {"from": "step-2", "to": "step-3"}, {"from": "module", "to": "output"},
                  {"from": "output", "to": "next"}]}


def _overview(title: str) -> dict:
    return {
        "schema": SCHEMA, "title": title,
        "root": {"id": "pipeline", "kind": "group", "layout": "row", "gap": 70, "children": [
            {"id": "solver", "kind": "module", "label": "Solver", "icon": "monitor"},
            {"id": "solution", "kind": "data", "label": "Solution", "icon": "time_series"},
            {"id": "postprocessing", "kind": "module", "label": "Post-Processing", "icon": "monitor",
             "layout": "column", "children": [
                 {"id": "method", "kind": "module", "label": "Method", "icon": "matrix"}]},
            {"id": "krom", "kind": "data", "label": "K-ROM", "icon": "krom"}]},
        "edges": [{"from": "solver", "to": "solution"}, {"from": "solution", "to": "postprocessing"},
                  {"from": "postprocessing", "to": "krom"}]}


TEMPLATES = {"module": _module, "overview": _overview}


def scaffold(template: str, title: str) -> dict:
    return TEMPLATES[template](title)
