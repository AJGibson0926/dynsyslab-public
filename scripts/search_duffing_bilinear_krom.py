"""Bounded search CLI for the explicit Duffing bilinear-control KROM candidate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from operations.duffing_bilinear import (
    DEFAULT_CONTROL_SOURCE, DEFAULT_MANIFEST, DEFAULT_OUTPUT, _jsonable, run,
)


def main(argv: list[str] | None = None) -> int:
    """Parse bounded-campaign options and delegate to the candidate API."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--controlled-data", type=Path, default=DEFAULT_CONTROL_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-attempts", type=int, default=25)
    parser.add_argument("--maximum-mode-count", type=int, default=6)
    parser.add_argument("--output-dt", type=float, default=0.001)
    parser.add_argument("--grid-size", type=int, default=81)
    parser.add_argument("--case-timeout-seconds", type=float, default=300.0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    if args.maximum_mode_count < 1:
        parser.error("--maximum-mode-count must be positive")
    if args.output_dt <= 0.0 or args.grid_size < 2 or args.case_timeout_seconds <= 0.0:
        parser.error("--output-dt must be positive and --grid-size must be at least 2")
    report = run(
        args.manifest, args.controlled_data, args.output,
        max_attempts=args.max_attempts,
        maximum_mode_count=args.maximum_mode_count,
        output_dt=args.output_dt,
        grid_size=args.grid_size,
        case_timeout_seconds=args.case_timeout_seconds,
        dry_run=args.dry_run,
        resume=args.resume,
    )
    print(json.dumps(_jsonable({
        "dry_run": bool(report.get("dry_run", False)),
        "attempt_count": report.get("attempt_count", 0),
        "stop_reason": report.get("stop_reason"),
        "successful_attempt_number": report.get("successful_attempt_number"),
        "campaign_report": str(args.output.resolve() / "campaign-report.json"),
        "campaign_ledger": str(args.output.resolve() / "campaign-ledger.json"),
        "candidate_screening": str(args.output.resolve() / "candidate-screening.json"),
    }), indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
