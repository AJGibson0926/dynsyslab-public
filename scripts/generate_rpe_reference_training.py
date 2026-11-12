"""Generate the retained dimensionless RPE training trajectory for parity work."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from dynsyslib import load_system
from solver import InterpolationOptions, SolveRequest, SolverOptions, solve


RPE_REFERENCE_INITIAL_STATE_DIMENSIONLESS=np.asarray((2.3333865986017348,0.),dtype=float)
RPE_REFERENCE_HORIZON_TAU=25.
RPE_REFERENCE_DTAU=.001
RPE_REFERENCE_RTOL=1e-9
RPE_REFERENCE_ATOL=1e-9


def build_request(*, archive_root: str|None=None) -> SolveRequest:
    """Build the exact, reusable first-stage parity request.

    This autonomous training request carries its exact seeded non-equilibrium
    IC explicitly. Free solves otherwise retain their system/request defaults;
    the solver supplies equilibrium only for controlled RPE without an IC.
    """
    system=load_system("rpe",preset="air_water_20c")
    scale=float(system.dimensionless_parameters.time_scale_s)
    initial_state=np.asarray(system.to_dimensional_state(RPE_REFERENCE_INITIAL_STATE_DIMENSIONLESS),dtype=float)
    return SolveRequest(
        system=system,
        kind="ode",
        time_span=(0.,RPE_REFERENCE_HORIZON_TAU*scale),
        initial_state=initial_state,
        options=SolverOptions(
            method="BDF",
            rtol=RPE_REFERENCE_RTOL,
            atol=RPE_REFERENCE_ATOL,
            use_dimensionless_form=True,
            execution_backend="python",
        ),
        interpolation=InterpolationOptions(
            enabled=True,
            time_step=RPE_REFERENCE_DTAU*scale,
            method="pchip",
        ),
        archive_root=archive_root,
        label="rpe-reference-training-bdf-tau25-dtau0p001",
    )


def main(argv: list[str]|None=None) -> int:
    """Generate the retained RPE reference free-training ensemble."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root",type=Path)
    args=parser.parse_args(argv)
    request=build_request(archive_root=None if args.archive_root is None else str(args.archive_root))
    result=solve(request)
    payload={"status":result.status,"message":result.message,"archive_path":result.archive_path}
    if result.archive_path:
        path=Path(result.archive_path)
        payload["sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
        payload["output_samples"]=result.time_series.sample_count
        payload["raw_samples"]=None if result.raw_time_series is None else result.raw_time_series.sample_count
    print(json.dumps(payload,indent=2,sort_keys=True))
    return 0 if result.status=="complete" else 1


if __name__=="__main__":
    raise SystemExit(main())
