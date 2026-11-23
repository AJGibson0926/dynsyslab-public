"""Small Python CLI that delegates every operation to installed native tools."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Sequence

from .workflow import installation_root, library_candidates


def _binary(name: str) -> Path:
    suffix = ".exe" if os.name == "nt" and not name.endswith(".exe") else ""
    return installation_root() / "bin" / f"{name}{suffix}"


def _require(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(
            f"Compiled DynSysLab component is missing: {path}. Install the native suite "
            "or set DYNSYSLAB_INSTALL_DIR.")
    return path


def parser() -> argparse.ArgumentParser:
    """Build the non-executing optional-client parser."""
    result = argparse.ArgumentParser(
        prog="dynsyslab", description="Optional Python client for compiled DynSysLab")
    commands = result.add_subparsers(dest="command")
    commands.add_parser("doctor", help="report compiled installation status")
    commands.add_parser("ui", help="launch Solver/Post-processor")
    commands.add_parser("viewer", help="launch Viewer")
    commands.add_parser("docs", help="print the installed start document")
    for name in ("systems", "validate", "plan", "run", "archive", "queue", "koopman-spectral", "lie-algebra", "krom-audit", "mathematics"):
        command = commands.add_parser(name, help=f"delegate {name} to the compiled CLI")
        command.add_argument("arguments", nargs=argparse.REMAINDER)
    return result


def diagnostics() -> dict[str, object]:
    """Return dependency-free diagnostics without requiring a working install."""
    root = installation_root()
    paths = {"cli": _binary("dynsyslab-cli"), "ui": _binary("DynSysLab"),
             "viewer": _binary("DynSysLabViewer"),
             "workflow": library_candidates()[0]}
    return {"schema": "dynsyslab.python-client-diagnostics/v1",
            "python_required_by_native": False, "install_root": str(root),
            "ok": all(path.is_file() for path in paths.values()),
            "paths": {name: {"path": str(path), "exists": path.is_file()}
                      for name, path in paths.items()}}


def main(argv: Sequence[str] | None = None) -> int:
    """Launch an app or pass an operation to the canonical compiled CLI."""
    selected = parser().parse_args(argv)
    if selected.command is None:
        parser().print_help(); return 0
    try:
        if selected.command == "doctor":
            print(json.dumps(diagnostics(), indent=2, sort_keys=True)); return 0
        if selected.command == "docs":
            path = installation_root() / "resources" / "docs" / "START_HERE.md"
            print(_require(path)); return 0
        if selected.command in {"ui", "viewer"}:
            name = "DynSysLab" if selected.command == "ui" else "DynSysLabViewer"
            subprocess.Popen([str(_require(_binary(name)))], cwd=installation_root() / "bin")
            return 0
        command = [str(_require(_binary("dynsyslab-cli"))), selected.command,
                   *selected.arguments]
        return subprocess.call(command)
    except (FileNotFoundError, OSError) as error:
        print(f"dynsyslab: {error}", file=sys.stderr); return 2


def compute_main() -> int:
    """Compatibility executable name for the primary native workbench."""
    return main(["ui"])


def viewer_main() -> int:
    """Compatibility executable name for native Viewer."""
    return main(["viewer"])
