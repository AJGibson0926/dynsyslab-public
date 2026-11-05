"""Lightweight CLI for persistent KROM/eigenfunction validation records."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from postprocessing.core import integration_coordinates, load_solution
from postprocessing.koopman import FreeTrajectoryDataset
from rom import (load_controlled_dataset, load_krom_archive,
                 load_krom_error_evaluations, save_krom_error_evaluations)


def _parser() -> argparse.ArgumentParser:
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest="command",required=True)
    free=commands.add_parser("free",help="evaluate eigenfunctions on free trajectories")
    free.add_argument("--krom",required=True,type=Path)
    free.add_argument("--trajectory",required=True,action="append",type=Path)
    free.add_argument("--role",required=True,choices=("training","test"))
    free.add_argument("--label",required=True)
    free.add_argument("--horizon",type=float,action="append")
    free.add_argument("--output",required=True,type=Path)
    control=commands.add_parser("control",help="evaluate one stored B_psi approximation")
    control.add_argument("--krom",required=True,type=Path)
    control.add_argument("--controlled-data",required=True,type=Path)
    control.add_argument("--representation",required=True,choices=("constant","bilinear"))
    control.add_argument("--label",required=True)
    control.add_argument("--bounds",required=True,nargs=4,type=float,
        metavar=("XMIN","XMAX","YMIN","YMAX"))
    control.add_argument("--grid-shape",nargs=2,type=int,default=(81,81),metavar=("NY","NX"))
    control.add_argument("--output",required=True,type=Path)
    compare=commands.add_parser("compare",help="compare evaluations from one or more bundles")
    compare.add_argument("--krom",required=True,type=Path)
    compare.add_argument("--bundle",required=True,action="append",type=Path)
    compare.add_argument("--metric",default="relative.mean")
    return parser


def _free(arguments) -> dict:
    model=load_krom_archive(arguments.krom); trajectories=[]; sources=[]
    for path in arguments.trajectory:
        source=load_solution(path)
        if model.coordinate_system=="dimensionless": source=integration_coordinates(source)
        trajectories.append(source.time_series); sources.append(str(path.resolve()))
    dataset=FreeTrajectoryDataset(trajectories,label=arguments.label,role=arguments.role,
        coordinate_system="model",provenance={"source_archives":sources})
    evaluations=model.evaluate_eigenfunction_errors(dataset,horizons=arguments.horizon)
    destination=save_krom_error_evaluations(arguments.output,model,
        metadata={"command":"free","krom_archive":str(arguments.krom.resolve())})
    return {"bundle":str(destination),"dataset":dataset.summary(),
        "eigenfunction_evaluations":{str(index):item.summary() for index,item in evaluations.items()}}


def _control(arguments) -> dict:
    model=load_krom_archive(arguments.krom); dataset=load_controlled_dataset(arguments.controlled_data)
    if arguments.representation not in (model.control_representations or {}):
        raise ValueError(f"KROM archive has no stored {arguments.representation!r} representation.")
    xmin,xmax,ymin,ymax=arguments.bounds
    evaluation=model.evaluate_control_approximation(
        model.control_representations[arguments.representation],dataset,
        ((xmin,ymin),(xmax,ymax)),label=arguments.label,
        grid_shape=tuple(arguments.grid_shape),
        metadata={"controlled_dataset_archive":str(arguments.controlled_data.resolve())})
    destination=save_krom_error_evaluations(arguments.output,model,
        metadata={"command":"control","krom_archive":str(arguments.krom.resolve())})
    return {"bundle":str(destination),"evaluation":evaluation.summary()}


def _compare(arguments) -> dict:
    model=load_krom_archive(arguments.krom)
    for path in arguments.bundle:
        load_krom_error_evaluations(path,model,attach=True)
    return model.compare_control_error_evaluations(metric=arguments.metric)


def main(argv=None) -> int:
    """Evaluate and archive free/test trajectory errors for one KROM."""
    arguments=_parser().parse_args(argv)
    result=_free(arguments) if arguments.command=="free" else (
        _control(arguments) if arguments.command=="control" else _compare(arguments))
    print(json.dumps(result,indent=2,sort_keys=True,default=str)); return 0


if __name__=="__main__": raise SystemExit(main())
