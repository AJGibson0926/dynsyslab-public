"""Lightweight command line interface for reduced-order-model operations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .selection import ObservableModeSelectionPolicy, select_observable_modes
from .observation import ObservationMapQualityPolicy, fit_state_observation_map
from .storage import load_krom_archive, reserve_krom_archive, reserve_observation_map_archive


def _parser() -> argparse.ArgumentParser:
    parser=argparse.ArgumentParser(prog="python -m rom")
    commands=parser.add_subparsers(dest="command",required=True)
    select=commands.add_parser("select-observable",help="select polynomial KROM modes by grouped state observability")
    select.add_argument("archive",type=Path)
    select.add_argument("--r2-threshold",type=float,default=.99)
    select.add_argument("--max-modes",type=int,default=30)
    select.add_argument("--minimum-relative-improvement",type=float,default=1e-10)
    select.add_argument("--control-authority-weight",type=float,default=.05)
    select.add_argument("--label",default="observable-selected-krom")
    select.add_argument("--report-only",action="store_true")
    observe=commands.add_parser("fit-observation",help="fit and validate an affine state decoder for a KROM")
    observe.add_argument("archive",type=Path)
    observe.add_argument("--validation-solution",type=Path)
    observe.add_argument("--validation-fraction",type=float,default=.2)
    observe.add_argument("--split",choices=("interleaved","tail"),default="interleaved")
    observe.add_argument("--ridge",type=float,default=0.)
    observe.add_argument("--minimum-r2",type=float,default=.99)
    observe.add_argument("--maximum-normalized-rmse",type=float,default=.10)
    observe.add_argument("--maximum-condition-number",type=float,default=1e12)
    observe.add_argument("--allow-quality-failure",action="store_true",
        help="publish a clearly unvalidated diagnostic map when thresholds fail")
    observe.add_argument("--label",default="observed-krom")
    observe.add_argument("--report-only",action="store_true")
    return parser


def main(argv: list[str]|None=None) -> int:
    args=_parser().parse_args(argv)
    source=args.archive.resolve()
    krom=load_krom_archive(source)
    with np.load(source,allow_pickle=False) as data:
        if "training_states" not in data.files:
            raise KeyError("KROM archive has no training_states array.")
        states=np.asarray(data["training_states"],dtype=float)
    if args.command=="fit-observation":
        validation=None
        if args.validation_solution is not None:
            from analysis_data import load_solution
            validation=np.asarray(load_solution(args.validation_solution.resolve()).values,dtype=float)
        policy=ObservationMapQualityPolicy(args.minimum_r2,args.maximum_normalized_rmse,
            args.maximum_condition_number)
        observation=fit_state_observation_map(krom,states,validation_states=validation,
            validation_fraction=args.validation_fraction,split=args.split,ridge=args.ridge,
            quality_policy=policy,require_quality=not args.allow_quality_failure)
        payload={"source_archive":str(source),"observation_map":observation.summary(),
            "validation_solution":None if args.validation_solution is None else str(args.validation_solution.resolve()),
            "published_archive":None}
        if not args.report_only:
            payload["published_archive"]=str(reserve_observation_map_archive(source,observation,label=args.label))
        print(json.dumps(payload,indent=2,sort_keys=True))
        return 0
    policy=ObservableModeSelectionPolicy(args.r2_threshold,args.max_modes,
        args.minimum_relative_improvement,args.control_authority_weight)
    report=select_observable_modes(krom,states,policy)
    payload={"source_archive":str(source),"selection":report.summary(),"published_archive":None}
    if not args.report_only:
        published=reserve_krom_archive(source,report.selected_retained_indices,label=args.label,
            selection_metadata=report.summary())
        payload["published_archive"]=str(published)
    print(json.dumps(payload,indent=2,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
