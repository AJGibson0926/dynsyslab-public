"""Generic command-line operations for any dynsyslab campaign manifest."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from .manager import JobManager


def _campaign(path: str | Path) -> tuple[Path, dict[str, Any]]:
    manifest = Path(path).expanduser().resolve()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    if "campaign_id" not in payload or "jobs" not in payload:
        raise ValueError("Expected a dynsyslab Campaigner campaign manifest.")
    return manifest, payload


def _runnable_campaign_jobs(manager: JobManager, campaign_id: str) -> list[str]:
    return [job_id for job_id, row in manager.status().items()
            if row.get("tags", {}).get("campaign_id") == campaign_id
            and (row.get("status") in {"pending", "blocked"}
                 or (row.get("status")=="failed" and row.get("retry_eligible",False)))]


def run_campaign(manifest_path: str | Path, *, workers: int = 4,
                 adaptive_policy: str = "adaptive_local") -> dict[str, Any]:
    """Safely run all currently runnable Jobs belonging to one Campaign."""
    manifest, campaign = _campaign(manifest_path)
    manager = JobManager()
    parents = list(campaign["jobs"])
    transitioned = [job_id for job_id in parents if manager.status(job_id).get("parallel_transition_manifest")]
    if transitioned:
        review = manager.review_parallel_transition(parents, campaign_path=manifest)
        if not review["safe_to_execute"]:
            raise RuntimeError(f"Post-transition review failed; execution held. Review: {review['report_path']}")
    jobs = _runnable_campaign_jobs(manager, campaign["campaign_id"])
    result = manager.run(jobs, mode="parallel", backend="local_process", max_workers=workers,
                         adaptive_policy=adaptive_policy)
    return {"campaign_id": campaign["campaign_id"], "submitted_job_ids": jobs, "result": result}


def main(argv: Sequence[str] | None = None) -> int:
    """Run one explicitly selected campaign scheduling operation."""
    parser = argparse.ArgumentParser(prog="python -m job_manager", description="Run and review dynsyslab campaigns.")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run-campaign", "status-campaign", "review-transition", "repair-transition"):
        sub = commands.add_parser(name)
        sub.add_argument("campaign_manifest", type=Path)
    run = commands.choices["run-campaign"]
    run.add_argument("--workers", type=int, default=4)
    run.add_argument("--adaptive-policy", choices=("fixed", "adaptive_local", "adaptive_portable"), default="adaptive_local")
    tray_list=commands.add_parser("tray-list",help="List reserved, unscheduled plans.")
    tray_list.add_argument("--include-removed",action="store_true")
    tray_add=commands.add_parser("tray-add",help="Reserve a campaign, collection, or job manifest without submitting it.")
    tray_add.add_argument("manifest",type=Path); tray_add.add_argument("--note",default="")
    tray_remove=commands.add_parser("tray-remove",help="Remove a reserved plan without touching scheduled work.")
    tray_remove.add_argument("entry_id"); tray_remove.add_argument("--keep-snapshot",action="store_true")
    args = parser.parse_args(argv)
    manager = JobManager()
    if args.command=="tray-list":
        print(json.dumps(manager.side_tray(include_removed=args.include_removed),indent=2,default=str)); return 0
    if args.command=="tray-add":
        print(json.dumps(manager.reserve(args.manifest,note=args.note),indent=2,default=str)); return 0
    if args.command=="tray-remove":
        print(json.dumps(manager.remove_reservation(args.entry_id,delete_snapshot=not args.keep_snapshot),indent=2,default=str)); return 0
    manifest, campaign = _campaign(args.campaign_manifest)
    if args.command == "status-campaign":
        print(json.dumps(manager.campaign_status(manifest), indent=2, default=str)); return 0
    if args.command == "review-transition":
        print(json.dumps(manager.review_parallel_transition(campaign["jobs"], campaign_path=manifest), indent=2, default=str)); return 0
    if args.command == "repair-transition":
        print(json.dumps(manager.repair_transition_interruptions(campaign["jobs"], campaign_path=manifest), indent=2, default=str)); return 0
    print(json.dumps(run_campaign(manifest, workers=args.workers, adaptive_policy=args.adaptive_policy), indent=2, default=str)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
