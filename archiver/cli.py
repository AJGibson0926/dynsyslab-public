"""Command-line entry point for archive audits and conservative cleanup."""
from __future__ import annotations
import argparse
import json
from typing import Sequence
from .core import ArchiveScope, Archiver

def main(argv: Sequence[str] | None = None) -> int:
    """Run one explicitly selected archive audit or cleanup operation."""
    parser=argparse.ArgumentParser(prog="dynsyslab-archive")
    parser.add_argument("command", choices=("audit", "deduplicate", "index"))
    parser.add_argument("--days", type=float, default=2.0)
    parser.add_argument("--section", action="append", default=["solutions"])
    parser.add_argument("--system", action="append", default=[])
    parser.add_argument("--apply", action="store_true", help="Perform conservative deletion; default is dry-run.")
    args=parser.parse_args(argv); archiver=Archiver(); scope=ArchiveScope(days=args.days,sections=tuple(args.section),systems=tuple(args.system))
    if args.command=="audit": result=archiver.audit(scope)
    elif args.command=="index": result=archiver.reconcile_index()
    else: result=archiver.deduplicate(scope,dry_run=not args.apply)
    payload=result if isinstance(result,dict) else (result.payload() if hasattr(result,"payload") else result.__dict__)
    print(json.dumps(payload,indent=2,default=str)); return 0

if __name__=="__main__": raise SystemExit(main())
