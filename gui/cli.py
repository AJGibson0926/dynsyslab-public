"""Command-line boundary for the compatibility Tk desktop interface."""
from __future__ import annotations

import argparse
from importlib.metadata import PackageNotFoundError, version
from typing import Sequence


def _version() -> str:
    try:
        return version("dynsyslab")
    except PackageNotFoundError:
        return "source-tree"


def parser() -> argparse.ArgumentParser:
    """Build the non-launching CLI parser for the legacy Tk interface."""
    result = argparse.ArgumentParser(
        prog="dynsyslab-gui",
        description=("Launch the compatibility Tk interface. The native "
                     "Solver/Post-processor is the primary DynSysLab UI."),
    )
    result.add_argument("--version", action="version", version=f"%(prog)s {_version()}")
    commands = result.add_subparsers(dest="command")
    commands.add_parser("launch", help="launch the compatibility Tk application")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Parse CLI arguments and launch Tk only after an explicit command."""
    selected = parser(); arguments = selected.parse_args(argv)
    if arguments.command != "launch":
        selected.print_help()
        return 0
    from .desktop import launch
    launch()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
