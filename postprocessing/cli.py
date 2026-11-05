"""Lightweight command line for inspectable DynSysLab post-processing workflows."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from .planning import plan_trajectory_hdmd_sindy, summarize_processing_runs
from .processor import PostprocessingPipeline
from .linalg_backends import available_linear_algebra_backends
from dataclasses import asdict


def _observable_argument(value: str):
    """Accept either a state index or a named system-observable key."""
    try:
        return int(value)
    except ValueError:
        return value


def _workflow_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--matrix-size",type=int)
    parser.add_argument("--sampling-window",type=int)
    parser.add_argument("--observable",type=_observable_argument,
                        help="Scalar state index or named system observable (for example, far_field_acoustic_echo).")
    parser.add_argument("--rank",type=int,help="Explicit positive truncation rank; omit to retain the full tolerance-filtered spectrum.")
    parser.add_argument("--svd-tolerance",type=float)
    parser.add_argument("--svd-driver",choices=("gesvd","gesdd"))
    parser.add_argument("--svd-threads",type=int,help="Explicit BLAS thread limit during the H-DMD SVD; omit to inherit the caller/Job policy.")
    parser.add_argument("--eig-threads",dest="eigendecomposition_threads",type=int,help="Explicit BLAS thread limit during the reduced eigendecomposition.")
    parser.add_argument("--linalg-backend",dest="linear_algebra_backend",choices=("auto","scipy","torch-cuda"),
                        help="Dense H-DMD backend; auto retains the admitted SciPy/LAPACK reference.")
    parser.add_argument("--sample-stride",type=int)
    parser.add_argument("--max-memory-mb",type=float)
    parser.add_argument("--null-space-backend",choices=("gram-eigh","svd"))
    parser.add_argument("--no-conjugate-reuse",action="store_true")
    parser.add_argument("--sindy-storage",dest="archive_storage",choices=("coefficient-first","full"),
                        help="Derived eigenfunction storage; coefficient-first reconstructs values lazily.")
    parser.add_argument("--sindy-coordinates",nargs="+",
                        help="Named observable-coordinate preset or ordered observable names used to represent SINDy eigenfunctions.")
    parser.add_argument("--coordinate-system",choices=("dimensionless","dimensional"),
                        help="State/time coordinate system of the SINDy training data.")


def _workflow(args, sources=None, *, archive=True):
    names=("matrix_size","sampling_window","observable","rank","svd_tolerance","svd_driver","svd_threads","eigendecomposition_threads","linear_algebra_backend","sample_stride","max_memory_mb")
    hdmd={name:getattr(args,name) for name in names if getattr(args,name,None) is not None}
    if hdmd.get("rank",1)<=0: raise ValueError("--rank must be positive when supplied.")
    if hdmd.get("svd_threads",1)<=0: raise ValueError("--svd-threads must be positive when supplied.")
    if hdmd.get("eigendecomposition_threads",1)<=0: raise ValueError("--eig-threads must be positive when supplied.")
    sindy={}
    if getattr(args,"null_space_backend",None) is not None: sindy["null_space_backend"]=args.null_space_backend
    if getattr(args,"no_conjugate_reuse",False): sindy["reuse_conjugate_fits"]=False
    if getattr(args,"max_memory_mb",None) is not None: sindy["max_memory_mb"]=args.max_memory_mb
    if getattr(args,"archive_storage",None) is not None: sindy["archive_storage"]=args.archive_storage
    selected=getattr(args,"sindy_coordinates",None)
    if selected is not None: sindy["coordinate_observable"]=selected[0] if len(selected)==1 else tuple(selected)
    if getattr(args,"coordinate_system",None) is not None: sindy["coordinate_system"]=args.coordinate_system
    return PostprocessingPipeline.preset("trajectory-hdmd-sindy",sources=sources,archive=archive,sindy_settings=sindy,**hdmd)


def build_parser() -> argparse.ArgumentParser:
    parser=argparse.ArgumentParser(prog="dynsyslab-postprocess",description="Plan or run archive-safe DynSysLab post-processing workflows.")
    commands=parser.add_subparsers(dest="command",required=True)
    describe=commands.add_parser("describe",help="Print the standard workflow without loading trajectory data."); _workflow_options(describe)
    plan=commands.add_parser("plan",help="Validate one source and report matrix/resource geometry without computing."); plan.add_argument("source",type=Path); _workflow_options(plan)
    run=commands.add_parser("run",help="Execute the standard trajectory H-DMD to SINDy KROM workflow.")
    run.add_argument("sources",type=Path,nargs="+"); run.add_argument("--execute",action="store_true",help="Required confirmation for the expensive numerical run.")
    run.add_argument("--no-archive",action="store_true",help="Return in-memory results without publishing derived archives."); _workflow_options(run)
    observables=commands.add_parser("observables",help="List the archived system's observable catalog and H-DMD eligibility without computing.")
    observables.add_argument("source",type=Path)
    compare=commands.add_parser("compare-observables",help="Run identical H-DMD settings for several scalar observables and compare their spectra.")
    compare.add_argument("source",type=Path)
    compare.add_argument("--observables",nargs="+",required=True,type=_observable_argument,
                         help="Two or more named scalar observables or state indices.")
    compare.add_argument("--execute",action="store_true",help="Required confirmation for the numerical decompositions.")
    compare.add_argument("--no-archive",action="store_true",help="Do not publish the comparison archive.")
    compare.add_argument("--viewer-output",type=Path,help="Optional interactive Viewer HTML output.")
    compare.add_argument("--matrix-size",type=int)
    compare.add_argument("--sampling-window",type=int)
    compare.add_argument("--rank",type=int)
    compare.add_argument("--svd-tolerance",type=float)
    compare.add_argument("--svd-driver",choices=("gesvd","gesdd"))
    compare.add_argument("--svd-threads",type=int)
    compare.add_argument("--eig-threads",dest="eigendecomposition_threads",type=int)
    compare.add_argument("--linalg-backend",dest="linear_algebra_backend",choices=("auto","scipy","torch-cuda"))
    compare.add_argument("--sample-stride",type=int)
    compare.add_argument("--max-memory-mb",type=float)
    commands.add_parser("backends",help="Report available dense linear-algebra runtimes without computing.")
    return parser


def main(argv: list[str]|None=None) -> int:
    parser=build_parser(); args=parser.parse_args(argv)
    try:
        if args.command=="backends":
            payload={"format":"dynsyslab.postprocessing-backends/v1","linear_algebra":{
                name:asdict(backend) for name,backend in available_linear_algebra_backends().items()}}
        elif args.command=="describe": payload=_workflow(args,sources=[],archive=True).describe()
        elif args.command=="plan":
            workflow=_workflow(args,sources=[args.source],archive=True)
            payload=plan_trajectory_hdmd_sindy(args.source,workflow.methods[0].settings,workflow.methods[1].settings)
        elif args.command=="observables":
            from dynsyslib import resolve_archived_system
            from .core import integration_coordinates, load_solution
            source=integration_coordinates(load_solution(args.source))
            system=resolve_archived_system(source.metadata.get("system"),
                                           source.metadata.get("system_parameters"))
            coordinates=(source.metadata.get("postprocessing_preprocessing") or {}).get("coordinate_system")
            catalog=system.observable_catalog(coordinates)
            summary=catalog.summary()
            payload={"format":"dynsyslab.hdmd-observables/v1","system":system.name,
                     "coordinate_system":coordinates,"aliases":summary["aliases"],
                     "observables":[]}
            for name,item in summary["observables"].items():
                payload["observables"].append({"key":name,**item,
                    "hdmd_scalar_eligible":item["output_dimension"]==1 and item["value_type"]=="real"})
        elif args.command=="run":
            if not args.execute: parser.error("run requires --execute; use plan first to inspect geometry and memory policy.")
            workflow=_workflow(args,sources=args.sources,archive=not args.no_archive)
            payload=summarize_processing_runs(workflow.run())
        else:
            if not args.execute:
                parser.error("compare-observables requires --execute because it runs one H-DMD decomposition per observable.")
            from .koopman import (compare_hankel_observables_archive,
                                  observable_spectrum_comparison_report)
            names=("matrix_size","sampling_window","rank","svd_tolerance","svd_driver",
                   "svd_threads","eigendecomposition_threads","linear_algebra_backend",
                   "sample_stride","max_memory_mb")
            settings={name:getattr(args,name) for name in names if getattr(args,name,None) is not None}
            comparison=compare_hankel_observables_archive(
                args.source,args.observables,archive=not args.no_archive,**settings)
            payload=comparison.summary()
            if args.viewer_output is not None:
                from importlib import import_module
                Viewer=import_module("viewer").Viewer
                output=Viewer(renderer="plotly").write_html(
                    observable_spectrum_comparison_report(comparison),args.viewer_output)
                payload["viewer_output"]=str(output)
    except (FileNotFoundError,KeyError,TypeError,ValueError,MemoryError) as error:
        parser.exit(2,f"dynsyslab-postprocess: {error}\n")
    print(json.dumps(payload,indent=2,sort_keys=True)); return 0


if __name__=="__main__": raise SystemExit(main())
