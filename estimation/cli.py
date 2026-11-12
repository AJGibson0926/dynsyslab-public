"""Lightweight CLI for observability and offline KLQE experiments."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def build_parser() -> argparse.ArgumentParser:
    parser=argparse.ArgumentParser(prog="dynsyslab-estimation")
    commands=parser.add_subparsers(dest="command",required=True)
    observe=commands.add_parser("observability",help="analyze the observable C attached to a KROM archive")
    observe.add_argument("krom",type=Path)
    observe.add_argument("--measurement",help="named KROM observable projection; defaults to the archived default")
    observe.add_argument("--method",choices=("matrix","discrete-gramian","continuous-gramian",
                         "continuous-infinite-gramian"),default="continuous-gramian")
    observe.add_argument("--horizon",type=float,default=10.)
    observe.add_argument("--sample-time",type=float)
    observe.add_argument("--output",type=Path)
    observe.add_argument("--viewer",type=Path,help="optional native Viewer HTML report")
    design=commands.add_parser("design",help="design and audit KLQE from named KROM observable projections")
    design.add_argument("krom",type=Path)
    design.add_argument("--measurement",required=True)
    design.add_argument("--state-output",default="state")
    design.add_argument("--process-variance",type=float,default=1.)
    design.add_argument("--measurement-variance",type=float,default=1.)
    design.add_argument("--output",type=Path,help="explicit immutable design artifact path")
    benchmark=commands.add_parser("benchmark",help="run KLQE against an archived free trajectory")
    benchmark.add_argument("krom",type=Path)
    benchmark.add_argument("solution",type=Path)
    benchmark.add_argument("--measurement",help="named measurement projection; defaults to the archived default")
    benchmark.add_argument("--state-output",default="state",help="named full-state projection or validated state observation map")
    benchmark.add_argument("--process-variance",type=float,default=1.)
    benchmark.add_argument("--measurement-variance",type=float,default=1.)
    benchmark.add_argument("--initial-estimate",help="comma-separated model-coordinate initial estimate")
    benchmark.add_argument("--output",type=Path)
    benchmark.add_argument("--viewer",type=Path,help="optional native Viewer HTML report")
    return parser


def _initial_estimate(argument: str | None, krom) -> np.ndarray:
    if argument is not None:
        result=np.fromstring(argument,sep=",",dtype=float)
        if result.size!=np.asarray(krom.system.initial_state).size:
            raise ValueError("--initial-estimate does not match the system state dimension.")
        return result
    physical=np.asarray(krom.system.initial_state,dtype=float)
    if krom.coordinate_system=="dimensionless" and hasattr(krom.system,"to_dimensionless_state"):
        return np.asarray(krom.system.to_dimensionless_state(physical),dtype=float)
    return physical


def main(argv: list[str]|None=None) -> int:
    args=build_parser().parse_args(argv)
    from rom import load_krom_archive
    krom=load_krom_archive(args.krom.resolve())
    if args.command=="observability":
        from rom.observability import analyze_krom_observability, save_observability_analysis
        horizon=None if args.method=="continuous-infinite-gramian" else args.horizon
        if args.method in {"matrix","discrete-gramian"} and horizon is not None:
            if not float(horizon).is_integer():
                raise ValueError("Matrix and discrete-Gramian horizons must be integer step counts.")
            horizon=int(horizon)
        reconstruction=None if args.measurement is None else args.measurement
        analysis=analyze_krom_observability(krom,reconstruction,method=args.method,horizon=horizon,
                                            sample_time=args.sample_time)
        output=None
        if args.output is not None:
            output=save_observability_analysis(args.output.resolve(),analysis)
        viewer_output=None
        if args.viewer is not None:
            from importlib import import_module
            presentation=import_module("viewer")
            ObservabilitySource,PlotRequest,Viewer=(presentation.ObservabilitySource,
                presentation.PlotRequest,presentation.Viewer)
            report=Viewer(renderer="plotly").add(PlotRequest("observability",ObservabilitySource(analysis))).render()
            viewer_output=Viewer(renderer="plotly").write_html(report,args.viewer.resolve())
        print(json.dumps({**analysis.summary(),"artifact":None if output is None else str(output),
                          "viewer":None if viewer_output is None else str(viewer_output)},indent=2,sort_keys=True))
        return 0
    if args.command=="design":
        from estimation.koopman import design_klqe_from_krom,save_klqe_design
        if args.process_variance<0. or args.measurement_variance<=0.:
            raise ValueError("Process variance must be nonnegative and measurement variance positive.")
        measurement=krom.observable_map(args.measurement)
        split_dimension=krom.A(output="split").shape[0]
        design=design_klqe_from_krom(krom,args.measurement,args.state_output,
            args.process_variance*np.eye(split_dimension),
            args.measurement_variance*np.eye(measurement.observable_dimension))
        output=None
        if args.output is not None:
            output=save_klqe_design(args.output.resolve(),design,krom_archive=args.krom.resolve(),
                                    measurement_name=args.measurement,state_output_name=args.state_output)
        print(json.dumps({"design":design.summary(),"artifact":None if output is None else str(output)},
                         indent=2,sort_keys=True))
        return 0
    measurement_name=(args.measurement or krom.default_observable_reconstruction)
    if measurement_name is None:
        raise ValueError("KLQE benchmark requires an attached measurement ObservableReconstruction.")
    measurement_reconstruction=krom.observable_map(measurement_name)
    if args.state_output not in krom.observable_reconstructions and krom.observation_map is None:
        raise ValueError("KLQE benchmark requires a named full-state reconstruction or attached state observation map.")
    from estimation.koopman import design_klqe_from_krom, run_klqe, save_klqe_result
    from postprocessing.core import integration_coordinates, load_solution
    trajectory=integration_coordinates(load_solution(args.solution.resolve()))
    states=np.asarray(trajectory.values,dtype=float); time=np.asarray(trajectory.time,dtype=float)
    measurement=measurement_reconstruction.observable.evaluate(states)
    split_dimension=krom.A(output="split").shape[0]
    measurement_dimension=measurement_reconstruction.observable_dimension
    if args.process_variance<0. or args.measurement_variance<=0.:
        raise ValueError("Process variance must be nonnegative and measurement variance positive.")
    design=design_klqe_from_krom(krom,measurement_name,args.state_output,
        args.process_variance*np.eye(split_dimension),args.measurement_variance*np.eye(measurement_dimension))
    result=run_klqe(design.estimator(_initial_estimate(args.initial_estimate,krom)),time,measurement,
                    true_states=states)
    output=None
    if args.output is not None:
        output=save_klqe_result(args.output.resolve(),result)
    viewer_output=None
    if args.viewer is not None:
        from importlib import import_module
        presentation=import_module("viewer")
        KLQESource,PlotRequest,TilePlan,Viewer=(presentation.KLQESource,
            presentation.PlotRequest,presentation.TilePlan,presentation.Viewer)
        viewer=Viewer(tile_plan=TilePlan("klqe",columns=2),renderer="plotly")
        report=viewer.add(PlotRequest("klqe",KLQESource(result))).render()
        viewer_output=viewer.write_html(report,args.viewer.resolve())
    print(json.dumps({"design":design.summary(),"result":result.summary(),
                      "artifact":None if output is None else str(output),
                      "viewer":None if viewer_output is None else str(viewer_output)},
                     indent=2,sort_keys=True))
    return 0


__all__=["build_parser","main"]
