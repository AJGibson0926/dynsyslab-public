"""Thin optional Python front door for the compiled durable job queue."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence


def parser() -> argparse.ArgumentParser:
    """Return the non-executing native queue client parser."""
    result = argparse.ArgumentParser(
        prog="dynsyslab jobs",
        description="Manage native campaigns and workflow jobs without a Python scheduler.",
    )
    commands = result.add_subparsers(dest="operation", required=True)
    submit = commands.add_parser("submit", help="validate and enqueue a workflow job")
    submit.add_argument("job", type=Path)
    commands.add_parser("list", help="list durable queue records")
    status = commands.add_parser("status", help="inspect one queue record")
    status.add_argument("job_id")
    cancel = commands.add_parser("cancel", help="cancel one queued or running record")
    cancel.add_argument("job_id")
    commands.add_parser("run-next", help="claim and execute the next queued job")
    recover = commands.add_parser("recover", help="recover stale running records")
    recover.add_argument("age_seconds", type=float)
    for command in commands.choices.values():
        command.add_argument("--archive-root", type=Path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Translate queue arguments without owning scheduling or execution."""
    selected = parser().parse_args(argv)
    arguments = ["queue", selected.operation]
    if selected.operation == "submit":
        arguments.append(str(selected.job.expanduser().resolve()))
    elif selected.operation in {"status", "cancel"}:
        arguments.append(selected.job_id)
    elif selected.operation == "recover":
        if selected.age_seconds < 0:
            parser().error("recover age must be nonnegative")
        arguments.append(str(selected.age_seconds))
    if selected.archive_root is not None:
        arguments.extend(["--archive-root", str(selected.archive_root.expanduser().resolve())])
    from .front_door import _native_cli
    return _native_cli(arguments)


__all__ = ["main", "parser"]
