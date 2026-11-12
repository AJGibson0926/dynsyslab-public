"""CLI and native launcher for the Solver and Post-processor workbench."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
from typing import Sequence

from .contracts import load_job, postprocess_template, solve_template, validate_job

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _cache_root() -> Path:
    base = Path(os.environ.get("LOCALAPPDATA", str(PROJECT_ROOT / "build")))
    return base / "DynSysLab" / "Compute"


def _build_dir(root: Path) -> Path:
    return root / "build" / "compute-ui"


def _executable(root: Path) -> Path:
    base = _build_dir(root)
    installed = Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "DynSysLab" / "bin"
    candidates = (
        installed / "DynSysLab.exe",
        root / "build" / "native-suite" / "bin" / "Release" / "DynSysLab.exe",
        base / "bin" / "Release" / "DynSysLab.exe",
        base / "Release" / "DynSysLab.exe",
        base / "Release" / "dynsyslab_compute.exe",
    )
    return next((path for path in candidates if path.is_file()), candidates[0])


def build_native(root: Path, configuration: str = "Release") -> Path:
    """Configure and compile the native Solver/Post-processor executable."""
    subprocess.run(["cmake", "-S", str(root / "compute_workbench" / "native"),
                    "-B", str(_build_dir(root))], cwd=root, check=True)
    subprocess.run(["cmake", "--build", str(_build_dir(root)), "--config",
                    configuration, "--parallel"], cwd=root, check=True)
    executable = _executable(root)
    if not executable.is_file():
        raise FileNotFoundError(f"Native compute build did not produce {executable}")
    return executable


def parser() -> argparse.ArgumentParser:
    """Build the non-executing public command-line contract."""
    result = argparse.ArgumentParser(prog="dynsyslab-compute",
        description="Plan, execute, build, and launch the DynSysLab Solver and Post-processor.")
    result.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    commands = result.add_subparsers(dest="command", required=True)
    registry = commands.add_parser("registry"); registry.add_argument("--output", type=Path)
    template = commands.add_parser("template"); template.add_argument("kind", choices=("solve", "postprocess")); template.add_argument("--output", type=Path)
    validate = commands.add_parser("validate"); validate.add_argument("--job", type=Path, required=True); validate.add_argument("--output", type=Path); validate.add_argument("--backend", choices=("native", "python-reference"), default="native"); validate.add_argument("--native-library", type=Path)
    plan = commands.add_parser("plan"); plan.add_argument("--job", type=Path, required=True); plan.add_argument("--output", type=Path); plan.add_argument("--backend", choices=("native", "python-reference"), default="native"); plan.add_argument("--native-library", type=Path)
    run = commands.add_parser("run"); run.add_argument("--job", type=Path, required=True); run.add_argument("--execute", action="store_true"); run.add_argument("--catalog-output", type=Path); run.add_argument("--progress-output", type=Path); run.add_argument("--output", type=Path)
    run.add_argument("--backend", choices=("native", "python-reference"), default="native")
    run.add_argument("--native-library", type=Path)
    run.add_argument("--max-threads", type=int, default=0)
    run.add_argument("--max-memory-mib", type=int, default=0)
    build = commands.add_parser("build"); build.add_argument("--configuration", default="Release")
    launch = commands.add_parser("launch"); launch.add_argument("--rebuild", action="store_true"); launch.add_argument("--refresh-data", action="store_true"); launch.add_argument("--wait", action="store_true")
    return result


def _emit(payload: dict, output: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if output:
        output = output.expanduser().resolve(); output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
        print(json.dumps({"ok": True, "path": str(output)}))
    else:
        print(text)


def _launch(root: Path, args: argparse.Namespace) -> Path:
    executable = _executable(root)
    if args.rebuild or not executable.is_file():
        executable = build_native(root)
    cache = _cache_root(); cache.mkdir(parents=True, exist_ok=True)
    from archiver.gateway import VIEWER_ARCHIVE
    if args.refresh_data:
        from archiver.gateway import ArchiveGateway
        ArchiveGateway("manage", VIEWER_ARCHIVE.archive_root).refresh_index()
    launch_config = {
        "schema": "dynsyslab.compute-launch/v1",
        "project_root": str(root),
        "archive_root": str(VIEWER_ARCHIVE.archive_root),
        "cache": str(cache),
    }
    (cache / "launch.json").write_text(
        json.dumps(launch_config, indent=2, sort_keys=True), encoding="utf-8")
    command = [str(executable), "--project-root", str(root),
               "--archive-root", str(VIEWER_ARCHIVE.archive_root),
               "--cache", str(cache)]
    process = subprocess.Popen(command, cwd=root)
    if args.wait:
        code = process.wait()
        if code: raise subprocess.CalledProcessError(code, command)
    return executable


def main(argv: Sequence[str] | None = None) -> int:
    """Dispatch one explicitly selected compute, build, or launch operation."""
    args = parser().parse_args(argv); root = args.project_root.expanduser().resolve()
    try:
        if args.command == "registry":
            from .registry import build_registry
            _emit(build_registry(root), args.output)
        elif args.command == "template": _emit(solve_template() if args.kind == "solve" else postprocess_template(), args.output)
        elif args.command == "validate":
            job = load_job(args.job)
            if args.backend == "native":
                from .native_workflow import plan_job as native_plan
                plan = native_plan(job, project_root=root, library=args.native_library)
                _emit({"ok": True, "job": job, "native_plan": plan}, args.output)
            else:
                _emit({"ok": True, "job": validate_job(job)}, args.output)
        elif args.command == "plan":
            job = load_job(args.job)
            if args.backend == "native":
                from .native_workflow import plan_job
                payload = plan_job(job, project_root=root, library=args.native_library)
            else:
                from .execution import plan_reference_job
                payload = plan_reference_job(validate_job(job), root)
            _emit(payload, args.output)
        elif args.command == "run":
            if not args.execute: raise ValueError("run requires --execute; plan the job first.")
            job = load_job(args.job)
            if args.backend == "python-reference":
                from .execution import execute_reference_job
                payload = execute_reference_job(
                    validate_job(job), root, catalog_output=args.catalog_output,
                    progress_output=args.progress_output)
            else:
                from .native_workflow import run_job
                archive = Path(os.environ.get(
                    "DYNSYSLAB_ARCHIVE_DIR", str(root / "archives")))
                payload = run_job(
                    job, archive, project_root=root, library=args.native_library,
                    progress_path=args.progress_output,
                    maximum_threads=args.max_threads,
                    maximum_memory_bytes=args.max_memory_mib * 1024 * 1024,
                )
                if args.catalog_output:
                    from viewer.ui.catalog import write_catalog
                    payload["catalog"] = str(write_catalog(root, args.catalog_output))
            _emit(payload, args.output)
        elif args.command == "build": _emit({"ok": True, "path": str(build_native(root, args.configuration))}, None)
        else: _emit({"ok": True, "path": str(_launch(root, args))}, None)
    except (FileNotFoundError, KeyError, TypeError, ValueError, RuntimeError) as error:
        parser().exit(2, f"dynsyslab-compute: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
