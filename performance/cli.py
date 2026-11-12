"""Command-line performance baseline for representative numerical kernels."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from dynsyslib import load_system
from dynsyslib.lattice import IsingModel
from postprocessing.koopman.dmd import dmd, hankel_dmd
from postprocessing.koopman.koopman_sindy import (
    _generator_gram,
    _generator_gram_terms,
    _sparse_null_vector,
    _sparse_null_vector_from_gram,
)
from postprocessing.koopman.virtual_arrays import polynomial_eigenfunction_values
from solver.models import serializable_system
from solver.models import HistoryOptions
from solver.native_ode import integrate_fixed_native
from solver.backends import sundials_backend
from solver.sundials import integrate_adaptive_sundials
from solver.ode import _fixed_step
from postprocessing.linalg_backends import torch_cuda_backend
from .core import BenchmarkCase, run_benchmarks
from .gates import evaluate_relative_speedups, evaluate_speedups


def _cases(iterations: int, warmup: int) -> tuple[BenchmarkCase, ...]:
    duffing = load_system("duffing", preset="double_well_reference")
    rpe = load_system("rpe", preset="air_water_20c")
    state = np.asarray(duffing.initial_state, dtype=float)
    rpe_state = np.asarray((2.3334, 0.0), dtype=float)
    time = np.linspace(0.0, 20.0, 4001)
    reference = solve_ivp(duffing.rhs, (float(time[0]), float(time[-1])), state, t_eval=time)
    if not reference.success:
        raise RuntimeError(f"Could not construct the benchmark trajectory: {reference.message}")
    signal = reference.y.T
    generator = np.random.default_rng(20260819)
    krom_eigenvalues = -generator.random(64) + 1j * generator.normal(size=64)
    krom_modes = generator.normal(size=(64, 64)) + 1j * generator.normal(size=(64, 64))
    krom_initial = generator.normal(size=64) + 1j * generator.normal(size=64)
    krom_time = np.linspace(0.0, 20.0, 4001)
    ising_reference = IsingModel(size=32, seed=20260819)
    ising_initial = np.asarray(ising_reference.spins, dtype=np.int64)
    sindy_theta = generator.normal(size=(5000, 15))
    sindy_gamma = generator.normal(size=(5000, 15))
    sindy_eigenvalues = np.linspace(-1.0, 1.0, 20) + 1j * np.linspace(0.2, 2.0, 20)
    sindy_states = generator.normal(size=(5000, 2))
    sindy_coefficients = generator.normal(size=(256, 15)) + 1j * generator.normal(size=(256, 15))
    sindy_dictionary = {"encoding":"polynomial-dictionary/v1","class":"PolynomialDictionary",
                        "degree":4,"include_constant":True}
    sindy_preprocessing = {"center":[0.,0.],"scale":[1.,1.],"multiplier":[1.,1.]}

    def numpy_krom_rollout():
        lifted = np.exp(np.multiply.outer(krom_time - krom_time[0], krom_eigenvalues)) * krom_initial
        return lifted, lifted @ krom_modes.T

    def sindy_null_svd():
        return [_sparse_null_vector(sindy_gamma-value*sindy_theta,.05,20) for value in sindy_eigenvalues]

    def sindy_null_gram_eigh():
        terms=_generator_gram_terms(sindy_theta,sindy_gamma); diagnostics={}
        return [_sparse_null_vector_from_gram(_generator_gram(terms,value),.05,20,
            matrix_factory=lambda value=value: sindy_gamma-value*sindy_theta,diagnostics=diagnostics)
            for value in sindy_eigenvalues]
    metadata = {
        "system": duffing.name,
        "system_parameters": serializable_system(duffing),
        "uniform_time_steps": True,
    }

    class _Session:
        request = type("Request", (), {"history": HistoryOptions()})()

    def rhs_calls():
        value = None
        for index in range(100_000):
            value = duffing.rhs(index * 1e-4, state)
        return value

    cases = [
        BenchmarkCase(
            "duffing-rhs-100k",
            rhs_calls,
            "rhs",
            iterations,
            warmup,
        ),
        BenchmarkCase("numpy-krom-rollout-4001x64x64", numpy_krom_rollout, "krom", iterations, warmup),
        BenchmarkCase(
            "python-duffing-rk4-20k-steps",
            lambda: _fixed_step(duffing.rhs, 0.0, 20.0, state, 0.001, "RK4", _Session()),
            "ode",
            iterations,
            warmup,
        ),
        BenchmarkCase(
            "python-rpe-dimensionless-rk4-25k-steps",
            lambda: _fixed_step(rpe.dimensionless_rhs, 0.0, 25.0, rpe_state, 0.001, "RK4", _Session()),
            "ode",
            iterations,
            warmup,
        ),
        BenchmarkCase(
            "scipy-duffing-rk45-20s",
            lambda: solve_ivp(duffing.rhs, (0.0, 20.0), state, method="RK45", rtol=1e-8, atol=1e-10),
            "ode-adaptive",
            iterations,
            warmup,
        ),
        BenchmarkCase(
            "python-ising2d-32x32-20-sweeps",
            lambda: IsingModel(size=32, seed=20260819).run(20),
            "lattice",
            iterations,
            warmup,
        ),
        BenchmarkCase("dmd-4001x2", lambda: dmd(signal, time=time, rank=2, metadata=metadata), "koopman", iterations, warmup),
        BenchmarkCase(
            "hankel-dmd-4001x2-delay32",
            lambda: hankel_dmd(signal,time=time,matrix_size=32,rank=8,metadata=metadata,
                               linear_algebra_backend="scipy"),
            "koopman",
            iterations,
            warmup,
        ),
        BenchmarkCase(
            "scipy-hdmd-full-3000x256",
            lambda: hankel_dmd(signal,time=time,matrix_size=256,sampling_window=3000,
                               observable=0,matrix_layout="window-snapshots",rank=None,
                               svd_tolerance=0.,svd_driver="gesvd",metadata=metadata,
                               linear_algebra_backend="scipy"),
            "koopman-full-spectrum",iterations,warmup),
        BenchmarkCase("sindy-null-svd-5000x15x20",sindy_null_svd,"koopman-sindy",iterations,warmup),
        BenchmarkCase("sindy-null-gram-eigh-5000x15x20",sindy_null_gram_eigh,"koopman-sindy",iterations,warmup),
        BenchmarkCase("sindy-lazy-polynomial-values-5000x256",
            lambda: polynomial_eigenfunction_values(sindy_states,sindy_coefficients,
                dictionary=sindy_dictionary,state_preprocessing=sindy_preprocessing),
            "koopman-sindy-view",iterations,warmup),
    ]
    if torch_cuda_backend().available:
        cases.append(BenchmarkCase(
            "torch-cuda-hdmd-full-3000x256",
            lambda: hankel_dmd(signal,time=time,matrix_size=256,sampling_window=3000,
                               observable=0,matrix_layout="window-snapshots",rank=None,
                               svd_tolerance=0.,svd_driver="gesvd",metadata=metadata,
                               linear_algebra_backend="torch-cuda"),
            "koopman-full-spectrum",iterations,warmup))
    if importlib.util.find_spec("_dynsyslab_native"):
        cases.insert(
            2,
            BenchmarkCase(
                "native-duffing-rk4-20k-steps",
                lambda: integrate_fixed_native(duffing, 0.0, 20.0, state, 0.001, "RK4"),
                "ode",
                iterations,
                warmup,
            ),
        )
        cases.insert(
            4,
            BenchmarkCase(
                "native-rpe-dimensionless-rk4-25k-steps",
                lambda: integrate_fixed_native(
                    rpe, 0.0, 25.0, rpe_state, 0.001, "RK4", coordinate_system="dimensionless"
                ),
                "ode",
                iterations,
                warmup,
            ),
        )
        module = __import__("_dynsyslab_native")
        if module.capabilities().get("krom_rollout", False):
            cases.append(BenchmarkCase(
                "native-krom-rollout-4001x64x64",
                lambda: module.krom_rollout(krom_eigenvalues, krom_modes, krom_initial, krom_time),
                "krom",
                iterations,
                warmup,
            ))
        if module.capabilities().get("ising2d", False):
            cases.append(BenchmarkCase(
                "native-ising2d-32x32-20-sweeps",
                lambda: module.ising2d_run(ising_initial, 1.0, 0.0, 2.269, 20260819, 20),
                "lattice",
                iterations,
                warmup,
            ))
    if sundials_backend().available:
        cases.append(BenchmarkCase(
            "sundials-duffing-adams-20s",
            lambda: integrate_adaptive_sundials(duffing.rhs, 0.0, 20.0, state, method="RK45", rtol=1e-8, atol=1e-10, max_step=None),
            "ode-adaptive",
            iterations,
            warmup,
        ))
    return tuple(cases)


def _select_cases(cases: tuple[BenchmarkCase,...], names: list[str]) -> tuple[BenchmarkCase,...]:
    if not names: return cases
    requested=set(names); available={case.name for case in cases}; unknown=sorted(requested-available)
    if unknown: raise ValueError(f"Unknown benchmark case(s): {', '.join(unknown)}.")
    return tuple(case for case in cases if case.name in requested)


def main(argv: list[str] | None=None) -> int:
    """Run selected benchmarks and optional versioned performance gates."""
    parser = argparse.ArgumentParser(description="Run reproducible DynSysLab numerical benchmarks.")
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output", help="Optional JSON report path.")
    parser.add_argument("--case",action="append",default=[],help="Run only this named case; may be repeated.")
    parser.add_argument("--list-cases",action="store_true",help="List available case names without running them.")
    parser.add_argument("--baseline", help="Optional prior JSON report to compare on the same runtime.")
    parser.add_argument(
        "--require-speedup",
        action="append",
        default=[],
        metavar="CASE=FACTOR",
        help="Required median speedup against --baseline; may be repeated.",
    )
    parser.add_argument(
        "--require-relative-speedup",
        action="append",
        default=[],
        metavar="CANDIDATE/REFERENCE=FACTOR",
        help="Required speedup between two cases in this run; may be repeated.",
    )
    args = parser.parse_args(argv)
    cases=_cases(args.iterations,args.warmup)
    if args.list_cases:
        print(json.dumps({"format":"dynsyslab.benchmark-cases/v1","cases":[case.name for case in cases]},indent=2))
        return 0
    try: selected=_select_cases(cases,args.case)
    except ValueError as error: parser.error(str(error))
    report = run_benchmarks(selected)
    payload = report.payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        target = Path(args.output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)
    exit_status = 0
    if args.require_speedup and not args.baseline:
        parser.error("--require-speedup requires --baseline.")
    if args.baseline:
        requirements = {}
        for expression in args.require_speedup:
            try:
                name, factor = expression.rsplit("=", 1)
                requirements[name] = float(factor)
            except ValueError as error:
                parser.error(f"Invalid speedup expression {expression!r}; use CASE=FACTOR.")
        if not requirements:
            requirements = {sample.name: 0.95 for sample in report.samples}
        baseline = json.loads(Path(args.baseline).read_text(encoding="utf-8"))
        gate = evaluate_speedups(baseline, payload, requirements)
        print(json.dumps(gate.payload(), indent=2, sort_keys=True))
        if not gate.passed:
            exit_status = 2
    if args.require_relative_speedup:
        requirements = {}
        for expression in args.require_relative_speedup:
            try:
                pair, factor = expression.rsplit("=", 1)
                candidate, reference = pair.split("/", 1)
                requirements[(candidate, reference)] = float(factor)
            except ValueError:
                parser.error(f"Invalid relative speedup {expression!r}; use CANDIDATE/REFERENCE=FACTOR.")
        gate = evaluate_relative_speedups(payload, requirements)
        print(json.dumps(gate.payload(), indent=2, sort_keys=True))
        if not gate.passed:
            exit_status = 2
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
