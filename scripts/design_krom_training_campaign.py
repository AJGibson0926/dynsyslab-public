"""CLI for generic, screened KROM training-domain expansion campaigns."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from campaigner import design_krom_training_plan
from operations.krom_training_screening import (KROMTrainingScreeningPolicy,
    publish_screened_krom_training_campaign, save_krom_training_screening_report,
    screen_krom_training_plan)
from rom import AxisAlignedStateDomain, load_krom_archive


def parser() -> argparse.ArgumentParser:
    """Return the non-executing-by-default campaign planner parser."""
    result = argparse.ArgumentParser(description=(
        "Design a system-agnostic KROM training-domain expansion campaign"))
    result.add_argument("krom", type=Path, help="source KROM archive")
    result.add_argument("--required-min", type=float, nargs="+", required=True)
    result.add_argument("--required-max", type=float, nargs="+", required=True)
    result.add_argument("--points-per-axis", type=int, nargs="+", required=True)
    result.add_argument("--anchor", type=float, nargs="+", action="append", default=[])
    result.add_argument("--preset")
    result.add_argument("--system-parameters-json", type=Path)
    result.add_argument("--margin-fraction", type=float, default=.1)
    result.add_argument("--model-time-end", type=float, default=25.)
    result.add_argument("--model-time-step", type=float, default=.001)
    result.add_argument("--label")
    result.add_argument("--screen", action="store_true",
        help="run bounded in-memory safety screening")
    result.add_argument("--screen-time-end", type=float, default=.25)
    result.add_argument("--screen-maximum-step", type=float, default=.002)
    result.add_argument("--state-minimum", type=float, nargs="+")
    result.add_argument("--state-maximum", type=float, nargs="+")
    result.add_argument("--maximum-absolute-state", type=float, nargs="+")
    result.add_argument("--screen-report", type=Path)
    result.add_argument("--save-screen-report", action="store_true")
    result.add_argument("--publish", action="store_true",
        help="publish immutable job plans after a successful --screen")
    return result


def main(argv: list[str] | None = None) -> int:
    """Plan, optionally screen, and explicitly publish one campaign."""
    args = parser().parse_args(argv)
    if args.publish and not args.screen:
        raise ValueError("--publish requires --screen.")
    if args.save_screen_report and not args.screen:
        raise ValueError("--save-screen-report requires --screen.")
    krom = load_krom_archive(args.krom.resolve())
    coordinates = AxisAlignedStateDomain.from_krom(krom).coordinate_system
    required = AxisAlignedStateDomain(tuple(args.required_min), tuple(args.required_max), coordinates)
    parameters = (None if args.system_parameters_json is None else
        json.loads(args.system_parameters_json.read_text(encoding="utf-8")))
    plan = design_krom_training_plan(source_krom=args.krom,
        required_domain=required, points_per_axis=args.points_per_axis,
        anchors=args.anchor, margin_fraction=args.margin_fraction,
        model_time_end=args.model_time_end, model_time_step=args.model_time_step,
        preset=args.preset, system_parameters=parameters)
    payload = {**plan.summary(), "screening": None, "published": False,
        "campaign_manifest": None, "job_plan_paths": []}
    report = None
    if args.screen:
        policy = KROMTrainingScreeningPolicy(model_time_end=args.screen_time_end,
            maximum_model_time_step=args.screen_maximum_step,
            state_minimum=None if args.state_minimum is None else tuple(args.state_minimum),
            state_maximum=None if args.state_maximum is None else tuple(args.state_maximum),
            maximum_absolute_state=(None if args.maximum_absolute_state is None
                else tuple(args.maximum_absolute_state)))
        report = screen_krom_training_plan(plan, policy=policy)
        payload["screening"] = report.summary(include_results=False)
        if args.screen_report is not None:
            destination = args.screen_report.resolve()
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(json.dumps(report.summary(), indent=2, sort_keys=True),
                encoding="utf-8")
            payload["screening_report_path"] = str(destination)
        if args.save_screen_report:
            payload["immutable_screening_report_path"] = str(
                save_krom_training_screening_report(plan, report))
    if args.publish:
        campaign, collection = publish_screened_krom_training_campaign(
            plan, report, label=args.label)
        payload.update({"published": True,
            "campaign_manifest": collection.tags["campaign_manifest"],
            "job_plan_paths": campaign.metadata.get("job_plan_paths", [])})
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
