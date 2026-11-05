"""Thin Python CLI over the canonical compiled DynSysLab workflow ABI."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Sequence

from project_paths import PROJECT_ROOT


_JOB_KINDS = {"solve": "solve", "postprocess": "postprocess",
              "estimate": "postprocess"}
_ESTIMATION_METHODS = {"observability", "klqe"}


def _parser(scope: str) -> argparse.ArgumentParser:
    subject = {
        "solve": "Solver and Post-processor workflows",
        "postprocess": "post-processing workflows",
        "estimate": "observability and KLQE workflows",
    }[scope]
    parser = argparse.ArgumentParser(
        prog=f"dynsyslab {scope}",
        description=f"Optional Python front door for compiled DynSysLab {subject}.",
    )
    parser.add_argument("--native-library", type=Path)
    commands = parser.add_subparsers(dest="operation", required=True)
    template = commands.add_parser("template", help="emit a versioned job template")
    template.add_argument("--output", type=Path)
    for name in ("validate", "plan"):
        command = commands.add_parser(name, help=f"{name} a job without execution")
        command.add_argument("job", type=Path)
        command.add_argument("--output", type=Path)
    run = commands.add_parser("run", help="execute through the compiled workflow")
    run.add_argument("job", type=Path)
    run.add_argument("--execute", action="store_true",
                     help="required acknowledgement for scientific execution")
    run.add_argument("--archive-root", type=Path)
    run.add_argument("--progress", type=Path)
    run.add_argument("--max-threads", type=int, default=0)
    run.add_argument("--max-memory-mib", type=int, default=0)
    run.add_argument("--output", type=Path)
    return parser


def _emit(payload: dict[str, Any], output: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if output is None:
        print(text)
        return
    destination = output.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text + "\n", encoding="utf-8")
    print(json.dumps({"ok": True, "path": str(destination)}))


def _template(scope: str) -> dict[str, Any]:
    from compute_workbench.contracts import postprocess_template, solve_template
    if scope == "solve":
        return solve_template()
    if scope == "postprocess":
        return postprocess_template()
    return {
        "schema": "dynsyslab.compute-job/v1",
        "kind": "postprocess",
        "label": "interactive-estimation",
        "sources": [],
        "methods": [{
            "name": "observability", "enabled": True,
            "settings": {
                "krom": "", "measurement": "echo",
                "state_reconstruction": "state",
                "observability_method": "continuous-gramian",
                "observability_horizon": 10.0,
            },
        }],
        "archive": {"enabled": True},
    }


def _job(path: Path, scope: str) -> dict[str, Any]:
    from compute_workbench.contracts import load_job
    result = load_job(path)
    kind = _JOB_KINDS[scope]
    if result.get("kind") != kind:
        raise ValueError(f"{scope} front door requires a {kind!r} job document")
    if scope == "estimate":
        methods = {str(item.get("name", "")).lower()
                   for item in result.get("methods", ()) if item.get("enabled", True)}
        if not methods or not methods.issubset(_ESTIMATION_METHODS):
            raise ValueError(
                "estimate jobs may contain only observability and klqe stages")
    return result


def main(scope: str, argv: Sequence[str] | None = None) -> int:
    """Run a kind-scoped adapter without implementing numerical work."""
    if scope not in _JOB_KINDS:
        raise ValueError(f"unknown native workflow scope: {scope}")
    parser = _parser(scope)
    args = parser.parse_args(argv)
    try:
        if args.operation == "template":
            _emit(_template(scope), args.output)
            return 0
        job = _job(args.job, scope)
        from compute_workbench.native_workflow import NativeWorkflow
        workflow = NativeWorkflow(args.native_library, PROJECT_ROOT)
        if args.operation in {"validate", "plan"}:
            plan = workflow.plan(job)
            payload = ({"schema": "dynsyslab.job-validation/v1", "ok": True,
                        "scope": scope, "kind": job["kind"], "plan": plan}
                       if args.operation == "validate" else plan)
            _emit(payload, args.output)
            return 0
        if not args.execute:
            parser.error("run requires --execute; inspect `plan` first")
        if args.max_threads < 0 or args.max_memory_mib < 0:
            raise ValueError("thread and memory limits must be nonnegative")
        archive = args.archive_root or Path(os.environ.get(
            "DYNSYSLAB_ARCHIVE_DIR", str(PROJECT_ROOT / "archives")))
        result = workflow.run(
            job, archive, progress_path=args.progress,
            maximum_threads=args.max_threads,
            maximum_memory_bytes=args.max_memory_mib * 1024 * 1024,
        )
        _emit(result, args.output)
    except (FileNotFoundError, KeyError, TypeError, ValueError, RuntimeError) as error:
        parser.exit(2, f"dynsyslab {scope}: {error}\n")
    return 0


def solve_main(argv: Sequence[str] | None = None) -> int:
    """Console-script entry point for native solve jobs."""
    return main("solve", argv)


def postprocess_main(argv: Sequence[str] | None = None) -> int:
    """Console-script entry point for native post-processing jobs."""
    return main("postprocess", argv)


def estimate_main(argv: Sequence[str] | None = None) -> int:
    """Console-script entry point for native observability and KLQE jobs."""
    return main("estimate", argv)


__all__ = ["estimate_main", "main", "postprocess_main", "solve_main"]
