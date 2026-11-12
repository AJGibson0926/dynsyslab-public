"""CLI and launcher for the native C++ DynSysLab Viewer workbench."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
from typing import Sequence

from .ui import write_catalog, write_payload


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _cache_root() -> Path:
    base = Path(os.environ.get("LOCALAPPDATA", str(PROJECT_ROOT / "build")))
    return base / "DynSysLab" / "Viewer"


def _build_dir(root: Path) -> Path:
    return root / "build" / "viewer-ui"


def _executable(root: Path) -> Path:
    directory = _build_dir(root)
    installed = Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "DynSysLab" / "bin"
    candidates = (
        installed / "DynSysLabViewer.exe",
        root / "build" / "native-suite" / "bin" / "Release" / "DynSysLabViewer.exe",
        directory / "bin" / "Release" / "DynSysLabViewer.exe",
        directory / "Release" / "DynSysLabViewer.exe",
        directory / "Release" / "dynsyslab_viewer.exe",
    )
    return next((path for path in candidates if path.is_file()), candidates[0])


def build_native(root: Path, *, configuration: str = "Release") -> Path:
    source = root / "viewer" / "native"
    build = _build_dir(root)
    subprocess.run(["cmake", "-S", str(source), "-B", str(build)], check=True,
                   cwd=root)
    subprocess.run(["cmake", "--build", str(build), "--config", configuration,
                    "--parallel"], check=True, cwd=root)
    executable = _executable(root)
    if not executable.is_file():
        raise FileNotFoundError(f"Native Viewer build did not produce {executable}")
    return executable


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(prog="dynsyslab-viewer",
        description="Catalog, inspect, build, and launch the native DynSysLab Viewer UI.")
    result.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    commands = result.add_subparsers(dest="command", required=True)
    catalog = commands.add_parser("catalog", help="write a compact metadata catalog")
    catalog.add_argument("--output", type=Path, required=True)
    payload = commands.add_parser("payload", help="write one lazy native plot payload")
    payload.add_argument("--source", type=Path, required=True)
    payload.add_argument("--view", required=True)
    payload.add_argument("--observable", default="state_0")
    payload.add_argument("--output", type=Path, required=True)
    native = commands.add_parser("build", help="build the C++ Viewer application")
    native.add_argument("--configuration", default="Release")
    launch = commands.add_parser("launch", help="launch the C++ Viewer application")
    launch.add_argument("--rebuild", action="store_true")
    launch.add_argument("--wait", action="store_true")
    launch.add_argument("--source", type=Path)
    return result


def _launch(root: Path, args: argparse.Namespace) -> Path:
    executable = _executable(root)
    if args.rebuild or not executable.is_file():
        executable = build_native(root)
    cache = _cache_root()
    cache.mkdir(parents=True, exist_ok=True)
    from archiver.gateway import VIEWER_ARCHIVE
    command = [str(executable), "--project-root", str(root),
               "--archive-root", str(VIEWER_ARCHIVE.archive_root),
               "--cache", str(cache)]
    if args.source:
        command.extend(("--source", str(args.source.expanduser().resolve())))
    process = subprocess.Popen(command, cwd=root)
    if args.wait:
        return_code = process.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, command)
    return executable


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    root = args.project_root.expanduser().resolve()
    if args.command == "catalog":
        target = write_catalog(root, args.output)
    elif args.command == "payload":
        source = args.source if args.source.is_absolute() else root / args.source
        target = write_payload(source, args.view, args.output,
                               observable=args.observable)
    elif args.command == "build":
        target = build_native(root, configuration=args.configuration)
    else:
        target = _launch(root, args)
    print(json.dumps({"ok": True, "command": args.command, "path": str(target)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
