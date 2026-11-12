"""One discoverable command-line front door for humans and automation."""
from __future__ import annotations

import argparse
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
try:
    import tomllib
except ModuleNotFoundError:  # Historical reference lane supports Python 3.10.
    import tomli as tomllib
from typing import Sequence

from project_paths import PROJECT_ROOT


_PASSTHROUGH = {
    "archive": ("archiver.cli", "main"),
    "krom-training": ("scripts.design_krom_training_campaign", "main"),
}


def _distribution_version(name: str, source_project: Path | None = None) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        if source_project and source_project.is_file():
            return tomllib.loads(source_project.read_text(
                encoding="utf-8"))["project"]["version"] + "+source"
        return "not-installed"


def _version() -> str:
    return _distribution_version("dynsyslab", PROJECT_ROOT / "pyproject.toml")


def _native_executable(name: str) -> Path:
    """Resolve the stable installation first, then supported developer builds."""
    install = Path(os.environ.get(
        "DYNSYSLAB_INSTALL_DIR",
        str(Path(os.environ.get("LOCALAPPDATA", PROJECT_ROOT))
            / "Programs" / "DynSysLab"),
    )).expanduser()
    candidates = (
        install / "bin" / name,
        PROJECT_ROOT / "build" / "native-suite-fast" / "bin" / name,
        PROJECT_ROOT / "build" / "native-suite" / "bin" / "Release" / name,
    )
    return next((path.resolve() for path in candidates if path.is_file()),
                candidates[0].resolve())


def parser() -> argparse.ArgumentParser:
    """Return the non-launching DynSysLab front-door parser."""
    result = argparse.ArgumentParser(
        prog="dynsyslab",
        description=("DynSysLab front door: inspect the workspace, launch the "
                     "native workbenches, or discover supported systems."),
    )
    result.add_argument("--version", action="version",
                        version=f"%(prog)s {_version()}")
    commands = result.add_subparsers(dest="command")
    commands.add_parser("doctor", help="report paths, versions, and capabilities")
    commands.add_parser("systems", help="list registered dynamical systems")
    ui = commands.add_parser("ui", help="launch the primary Solver/Post-processor")
    ui.add_argument("--rebuild", action="store_true")
    ui.add_argument("--refresh-data", action="store_true")
    viewer = commands.add_parser("viewer", help="launch the native Viewer")
    viewer.add_argument("--rebuild", action="store_true")
    legacy = commands.add_parser("legacy-gui", help="launch the compatibility Tk UI")
    legacy.add_argument("--confirm", action="store_true",
                        help="acknowledge that the native UI is preferred")
    commands.add_parser("docs", help="print the human/agent start-document path")
    audit = commands.add_parser("audit", help="run the architecture-policy audit")
    audit.add_argument("--format", choices=("json", "markdown"), default="markdown")
    audit.add_argument("--strict", action="store_true")
    commands.add_parser("solve", help="plan or run a versioned Solver workbench job")
    commands.add_parser("postprocess", help="plan or run archive-safe post-processing")
    commands.add_parser("archive", help="audit or conservatively deduplicate archives")
    commands.add_parser("jobs", help="inspect or run Campaigner job manifests")
    commands.add_parser("estimate", help="run observability and KLQE workflows")
    commands.add_parser("krom-training", help="plan generic KROM domain-expansion training")
    return result


def workspace_diagnostics() -> dict:
    """Return a machine-readable, dependency-light workspace diagnosis."""
    archive = Path(os.environ.get("DYNSYSLAB_ARCHIVE_DIR",
                                  str(PROJECT_ROOT / "archives"))).expanduser().resolve()
    paths = {
        "archive": archive,
        "compute_executable": _native_executable("DynSysLab.exe"),
        "viewer_executable": _native_executable("DynSysLabViewer.exe"),
        "native_cli": _native_executable("dynsyslab-cli.exe"),
        "start_document": PROJECT_ROOT / "docs" / "START_HERE.md",
    }
    return {
        "schema": "dynsyslab.workspace-diagnostics/v1",
        "ok": _version() != "not-installed",
        "project_root": str(PROJECT_ROOT),
        "python": {"executable": sys.executable, "version": platform.python_version()},
        "versions": {
            "python_suite": _version(),
            "native_kernels": _distribution_version(
                "dynsyslab-native", PROJECT_ROOT / "native" / "pyproject.toml"),
        },
        "paths": {name: {"path": str(path), "exists": path.exists()}
                  for name, path in paths.items()},
        "next_commands": {
            "primary_ui": str(paths["compute_executable"]),
            "viewer": str(paths["viewer_executable"]),
            "native_doctor": f'{paths["native_cli"]} doctor',
            "architecture_audit": "python -m quality.architecture --format markdown",
        },
    }


def system_registry() -> dict:
    """Return the canonical catalog reported by the compiled registry."""
    executable = _native_executable("dynsyslab-cli.exe")
    completed = subprocess.run(
        [str(executable), "systems"], cwd=PROJECT_ROOT, capture_output=True,
        text=True, check=True,
    )
    return json.loads(completed.stdout)


def _launch_module(module: str, arguments: list[str]) -> int:
    return subprocess.call([sys.executable, "-m", module, *arguments], cwd=PROJECT_ROOT)


def _native_cli(arguments: list[str]) -> int:
    """Run an accepted compiled command from the optional Python front door."""
    executable = _native_executable("dynsyslab-cli.exe")
    if not executable.is_file():
        raise FileNotFoundError(
            "the compiled DynSysLab CLI is unavailable; build or install the native suite"
        )
    return subprocess.call([str(executable), *arguments], cwd=PROJECT_ROOT)


def _require_native_source(application: str, source: Path, executable: Path,
                           selected: argparse.ArgumentParser) -> None:
    if not executable.is_file() and not source.is_file():
        selected.error(
            f"the native {application} is not bundled in this wheel; "
            "launch it from a DynSysLab source checkout"
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Execute one explicitly selected front-door operation."""
    raw = list(sys.argv[1:] if argv is None else argv)
    if raw == ["systems"]:
        return _native_cli(raw)
    if raw and raw[0] in {"solve", "postprocess", "estimate"}:
        from .workflow_front_door import main as workflow_main
        return workflow_main(raw[0], raw[1:])
    if raw and raw[0] == "jobs":
        from .job_front_door import main as jobs_main
        return jobs_main(raw[1:])
    if (len(raw) >= 2 and raw[0] == "archive"
            and raw[1] in {"audit", "deduplicate"}
            and "--help" not in raw):
        return _native_cli(raw)
    if raw and raw[0] in _PASSTHROUGH:
        from importlib import import_module
        module, function = _PASSTHROUGH[raw[0]]
        return int(getattr(import_module(module), function)(raw[1:]))
    selected = parser()
    args = selected.parse_args(raw)
    if args.command is None:
        selected.print_help()
        return 0
    if args.command == "doctor":
        print(json.dumps(workspace_diagnostics(), indent=2, sort_keys=True))
    elif args.command == "systems":
        return _native_cli(["systems"])
    elif args.command == "docs":
        path = PROJECT_ROOT / "docs" / "START_HERE.md"
        print(path if path.is_file() else
              "START_HERE.md is not bundled in this installation; use the package README.")
    elif args.command == "ui":
        executable = _native_executable("DynSysLab.exe")
        _require_native_source(
            "Solver/Post-processor",
            PROJECT_ROOT / "compute_workbench" / "native" / "CMakeLists.txt",
            executable,
            selected,
        )
        if not args.rebuild and not args.refresh_data and executable.is_file():
            subprocess.Popen([str(executable)], cwd=executable.parent)
            return 0
        options = ["launch"]
        if args.rebuild:
            options.append("--rebuild")
        if args.refresh_data:
            options.append("--refresh-data")
        return _launch_module("compute_workbench.cli", options)
    elif args.command == "viewer":
        executable = _native_executable("DynSysLabViewer.exe")
        _require_native_source(
            "Viewer", PROJECT_ROOT / "viewer" / "native" / "CMakeLists.txt",
            executable,
            selected,
        )
        if not args.rebuild and executable.is_file():
            subprocess.Popen([str(executable)], cwd=executable.parent)
            return 0
        options = ["launch"] + (["--rebuild"] if args.rebuild else [])
        return _launch_module("viewer.ui_cli", options)
    elif args.command == "legacy-gui":
        if not args.confirm:
            selected.error("legacy-gui requires --confirm; use `dynsyslab ui` by default")
        return _launch_module("gui", ["launch"])
    else:
        from quality.architecture import main as audit_main
        options = ["--format", args.format]
        if args.strict:
            options.append("--strict")
        return audit_main(options)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
