"""Locate the asymptotic viscous RPE Koopman regime from echo history only.

This is an isolated experiment.  It fits one pooled Hankel basis from five
long, unforced viscous RPE training trajectories and never recomputes that
basis on held-out data.  Polynomial EDMD fits are repeated on progressively
later suffixes of the same fixed coordinates.  Direct (r, rdot) fits are
included only as a coordinate-system audit.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eig, eigh
from scipy.signal import fftconvolve
from scipy.sparse.linalg import svds

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_echo_sampled_koopman import (
    PAIR_STRIDE,
    PolynomialEDMD,
    eigenfunction_metrics,
    orient_coefficient,
    phase_regression,
    polynomial_exponents,
    polynomial_values,
)
from analyze_rpe_inviscid_echo_hankel import far_field_echo, load_native_source
from compare_rpe_echo_pooled_training_injectivity import PooledHankelOperator
from project_rpe_inviscid_echo_basis_out_of_sample import (
    COLORS,
    EMBEDDING_DECIMATION,
    EXPECTED_STEP,
)
from sweep_rpe_echo_delay_injectivity import hankel_total_energy


ROOT = Path(__file__).resolve().parents[1]
TRAINING_GLOB = (
    "rpe-viscous-echo-long-training-five-ics-trajectory-*"
    "__20261002T001350Z.npz"
)
VALIDATION_GLOB = (
    "rpe-viscous-echo-long-validation-five-ics-trajectory-*"
    "__20261002T001403Z.npz"
)
DELAY_SAMPLES = 8644
MAX_LATENT_RANK = 4
START_TIMES = (5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
END_TIME = 80.0
DEGREES = {"latent2": 5, "latent4": 3, "state": 5}
RIDGE_RELATIVE = 1.0e-10
METRICS = OUTPUT_ROOT / "rpe-viscous-echo-asymptotic-koopman.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-viscous-echo-asymptotic-koopman.npz"
SUMMARY_PLOT = "rpe-viscous-echo-asymptotic-koopman-summary"
EVOLUTION_PLOT = "rpe-viscous-echo-asymptotic-koopman-evolution"


def load_records(pattern: str, role: str) -> list[dict]:
    """Load the five expected free RPE records and prepare their echo histories."""
    sources = sorted((ROOT / "archives/solutions/rpe").glob(pattern))
    if len(sources) != 5:
        raise RuntimeError(f"expected five {role} records, found {len(sources)}")
    records = []
    for source in sources:
        time, state, metadata = load_native_source(source)
        if metadata.get("forcing_type") != "free":
            raise RuntimeError("representation training and validation must be unforced")
        if metadata.get("controller_kind") != "none":
            raise RuntimeError("unexpected controller in unforced data")
        if metadata.get("data_role") != role:
            raise RuntimeError(f"unexpected data role {metadata.get('data_role')}")
        if not np.isclose(float(metadata["resolved_system_parameters"]["mu_pa_s"]), 0.001):
            raise RuntimeError("comparison requires standard air-water viscosity")
        records.append(
            {
                "source": str(source),
                "time": time,
                "state": state,
                "echo": far_field_echo(state, metadata),
                "initial_state": np.asarray(metadata["initial_state"], dtype=float),
                "metadata": metadata,
            }
        )
    return records


def theoretical_eigenvalue(metadata: dict) -> complex:
    """Compute the positive-frequency eigenvalue of the RPE equilibrium linearization."""
    groups = metadata["dimensionless_groups"]
    parameters = metadata["resolved_system_parameters"]
    stiffness = (
        -3.0 * float(parameters["kappa"]) * float(groups["gas_pressure_ratio"])
        + float(groups["capillarity"])
    )
    jacobian = np.asarray(
        [[0.0, 1.0], [stiffness, -float(groups["viscosity"])]], dtype=float
    )
    values = np.linalg.eigvals(jacobian)
    return complex(values[np.argmax(values.imag)])


def fit_hankel_basis(records: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # The equilibrium echo is exactly zero for the model.  Centering on the
    # transient sample mean would create a spurious nonzero late-time history.
    """Fit a pooled RMS-scaled echo basis without subtracting the transient mean."""
    signals = [record["echo"][::EMBEDDING_DECIMATION] for record in records]
    global_scale = float(np.sqrt(np.mean(np.concatenate(signals) ** 2)))
    normalized = [signal / global_scale for signal in signals]
    operator = PooledHankelOperator(normalized, DELAY_SAMPLES)
    left, singular, right = svds(
        operator,
        k=MAX_LATENT_RANK,
        which="LM",
        solver="arpack",
        tol=1.0e-9,
        maxiter=10000,
        return_singular_vectors=True,
        rng=np.random.default_rng(20261002),
    )
    order = np.argsort(singular)[::-1]
    left = left[:, order]
    singular = singular[order]
    right = right[order]
    for mode in range(MAX_LATENT_RANK):
        pivot = int(np.argmax(np.abs(left[:, mode])))
        if left[pivot, mode] < 0.0:
            left[:, mode] *= -1.0
            right[mode] *= -1.0
    total_energy = sum(hankel_total_energy(signal, DELAY_SAMPLES) for signal in normalized)
    cumulative = np.cumsum(singular * singular) / total_energy
    return left, singular, cumulative


def project_echo(record: dict, modes: np.ndarray, singular: np.ndarray) -> dict:
    """Apply the frozen Hankel modes to causal echo windows and align the results."""
    signal = record["echo"][::EMBEDDING_DECIMATION]
    # Scaling of the raw signal cancels against the fitted singular values.
    coordinates = np.column_stack(
        [
            fftconvolve(signal, modes[:, mode][::-1], mode="valid")
            / singular[mode]
            for mode in range(MAX_LATENT_RANK)
        ]
    )
    return {
        "source": record["source"],
        "initial_state": record["initial_state"],
        "time": record["time"][::EMBEDDING_DECIMATION][DELAY_SAMPLES - 1 :],
        "state": record["state"][::EMBEDDING_DECIMATION][DELAY_SAMPLES - 1 :],
        "coordinates": coordinates,
    }


def fit_edmd(coordinate_runs: list[np.ndarray], degree: int) -> PolynomialEDMD:
    """Fit a regularized polynomial sampled operator and its eigensystem from separate runs."""
    combined = np.vstack(coordinate_runs)
    center = np.mean(combined, axis=0)
    scale = np.std(combined, axis=0)
    if np.any(scale <= np.finfo(float).eps):
        raise RuntimeError("coordinate normalization is singular")
    exponents = polynomial_exponents(combined.shape[1], degree)
    normalized = [(run - center) / scale for run in coordinate_runs]
    raw_features = [polynomial_values(run, exponents) for run in normalized]
    feature_scale = np.sqrt(np.mean(np.vstack(raw_features) ** 2, axis=0))
    feature_scale[feature_scale <= np.finfo(float).eps] = 1.0
    features = [values / feature_scale for values in raw_features]
    minus = np.vstack([values[:-PAIR_STRIDE] for values in features])
    plus = np.vstack([values[PAIR_STRIDE:] for values in features])
    gram = minus.T @ minus / minus.shape[0]
    cross = minus.T @ plus / minus.shape[0]
    ridge = RIDGE_RELATIVE * float(np.trace(gram) / gram.shape[0])
    operator = np.linalg.solve(gram + ridge * np.eye(gram.shape[0]), cross)
    values, vectors = eig(operator)
    return PolynomialEDMD(
        center=center,
        scale=scale,
        feature_scale=feature_scale,
        exponents=exponents,
        operator=operator,
        eigenvalues=values,
        eigenvectors=vectors,
        pair_step_tau=EXPECTED_STEP * EMBEDDING_DECIMATION * PAIR_STRIDE,
        training_runs=len(coordinate_runs),
        ridge=ridge,
    )


def prescribed_refit(
    model: PolynomialEDMD,
    coordinate_runs: list[np.ndarray],
    multiplier: complex,
) -> np.ndarray:
    """Fit polynomial coefficients at a fixed multiplier using the normalized sampled defect."""
    features = [model.features(run) for run in coordinate_runs]
    minus = np.vstack([values[:-PAIR_STRIDE] for values in features])
    plus = np.vstack([values[PAIR_STRIDE:] for values in features])
    residual = plus - multiplier * minus
    numerator = residual.conj().T @ residual / residual.shape[0]
    denominator = minus.conj().T @ minus / minus.shape[0]
    ridge = 1.0e-12 * float(np.trace(denominator).real / denominator.shape[0])
    values, vectors = eigh(
        numerator,
        denominator + ridge * np.eye(denominator.shape[0]),
        check_finite=True,
    )
    coefficient = vectors[:, int(np.argmin(values))]
    rms = float(
        np.sqrt(
            np.mean(
                np.abs(
                    np.concatenate(
                        [model.evaluate(run, coefficient) for run in coordinate_runs]
                    )
                )
                ** 2
            )
        )
    )
    coefficient /= rms
    return orient_coefficient(model, coefficient, coordinate_runs[0])


def discovered_candidate(
    model: PolynomialEDMD, coordinate_runs: list[np.ndarray]
) -> tuple[complex, np.ndarray, list[dict]]:
    """Select an admissible oscillatory EDMD mode using the training diagnostics."""
    rows = []
    for index, (value, vector) in enumerate(
        zip(model.eigenvalues, model.eigenvectors.T, strict=True)
    ):
        if value.imag <= 0.0 or not (0.90 <= abs(value) <= 1.005):
            continue
        continuous = np.log(value) / model.pair_step_tau
        if not (4.0 <= continuous.imag <= 7.0 and -1.0 <= continuous.real <= 0.1):
            continue
        values_by_run = [model.evaluate(run, vector) for run in coordinate_runs]
        norm = float(np.sqrt(np.mean(np.abs(np.concatenate(values_by_run)) ** 2)))
        if norm <= np.finfo(float).tiny:
            continue
        vector = vector / norm
        residuals = []
        phase_r2 = []
        for run, values in zip(coordinate_runs, values_by_run, strict=True):
            values = values / norm
            residuals.append(
                float(
                    np.linalg.norm(values[PAIR_STRIDE:] - value * values[:-PAIR_STRIDE])
                    / max(np.linalg.norm(values[PAIR_STRIDE:]), np.finfo(float).tiny)
                )
            )
            time = np.arange(values.size) * EXPECTED_STEP * EMBEDDING_DECIMATION
            phase_r2.append(phase_regression(time, values)["r_squared"])
        rows.append(
            {
                "index": index,
                "multiplier": value,
                "continuous_eigenvalue": continuous,
                "coefficient": vector,
                "mean_residual": float(np.mean(residuals)),
                "minimum_phase_r_squared": float(np.min(phase_r2)),
            }
        )
    if not rows:
        raise RuntimeError("no admissible oscillatory EDMD candidate")
    rows.sort(key=lambda row: (row["mean_residual"], -row["minimum_phase_r_squared"]))
    best = rows[0]
    coefficient = orient_coefficient(
        model, np.asarray(best["coefficient"], dtype=complex), coordinate_runs[0]
    )
    return complex(best["multiplier"]), coefficient, rows


def suffix_runs(
    projected: list[dict], start: float, coordinate_key: str, rank: int | None
) -> tuple[list[np.ndarray], list[dict]]:
    """Select aligned training or validation suffixes in the requested coordinates."""
    coordinates = []
    metric_runs = []
    for index, record in enumerate(projected, start=1):
        mask = (record["time"] >= start) & (record["time"] <= END_TIME)
        if coordinate_key == "coordinates":
            values = record[coordinate_key][mask, :rank]
        else:
            values = record[coordinate_key][mask]
        time = record["time"][mask]
        coordinates.append(values)
        metric_runs.append(
            {
                "label": f"trajectory #{index}",
                "initial_state": record["initial_state"],
                "coordinates": values,
                "latent_time": time,
                "echo_period_tau": None,
            }
        )
    return coordinates, metric_runs


def anchored_error(
    model: PolynomialEDMD,
    coefficient: np.ndarray,
    multiplier: complex,
    run: dict,
) -> float:
    """Measure relative L2 departure from the exponential anchored at the first sample."""
    values = model.evaluate(run["coordinates"], coefficient)
    elapsed = run["latent_time"] - run["latent_time"][0]
    continuous = np.log(multiplier) / model.pair_step_tau
    predicted = values[0] * np.exp(continuous * elapsed)
    return float(
        np.linalg.norm(values - predicted)
        / max(np.linalg.norm(values), np.finfo(float).tiny)
    )


def summarize_metrics(rows: list[dict], anchors: list[float]) -> dict:
    """Aggregate one-step, horizon, anchored, and phase diagnostics across trajectories."""
    return {
        "median_one_step": float(np.median([row["one_step_relative_l2"] for row in rows])),
        "median_one_tau": float(np.median([row["horizon_relative_l2"]["1"] for row in rows])),
        "median_four_tau": float(np.median([row["horizon_relative_l2"]["4"] for row in rows])),
        "median_anchored": float(np.median(anchors)),
        "median_phase_r_squared": float(np.median([row["phase_linearity_r_squared"] for row in rows])),
        "median_observed_frequency": float(
            np.median([row["observed_eigenfunction_frequency_rad_per_tau"] for row in rows])
        ),
    }


def fit_window(
    name: str,
    start: float,
    training_coordinates: list[np.ndarray],
    training_runs: list[dict],
    validation_coordinates: list[np.ndarray],
    validation_runs: list[dict],
    degree: int,
    theoretical_lambda: complex,
) -> dict:
    """Compare discovered and prescribed eigenfunctions on one fitting time window."""
    model = fit_edmd(training_coordinates, degree)
    discovered_multiplier, discovered_coefficient, candidates = discovered_candidate(
        model, training_coordinates
    )
    prescribed_multiplier = complex(np.exp(theoretical_lambda * model.pair_step_tau))
    prescribed_coefficient = prescribed_refit(
        model, training_coordinates, prescribed_multiplier
    )
    result = {"name": name, "start_tau": start, "degree": degree, "cases": {}}
    for case_name, multiplier, coefficient in (
        ("discovered", discovered_multiplier, discovered_coefficient),
        ("prescribed", prescribed_multiplier, prescribed_coefficient),
    ):
        train = eigenfunction_metrics(model, coefficient, multiplier, training_runs)
        validation = eigenfunction_metrics(model, coefficient, multiplier, validation_runs)
        train_anchor = [anchored_error(model, coefficient, multiplier, run) for run in training_runs]
        validation_anchor = [
            anchored_error(model, coefficient, multiplier, run) for run in validation_runs
        ]
        result["cases"][case_name] = {
            "multiplier": multiplier,
            "continuous_eigenvalue": np.log(multiplier) / model.pair_step_tau,
            "coefficient": coefficient,
            "training_summary": summarize_metrics(train, train_anchor),
            "validation_summary": summarize_metrics(validation, validation_anchor),
            "training_raw": train,
            "validation_raw": validation,
        }
    result["candidate_count"] = len(candidates)
    return {"model": model, **result}


def clean_case(case: dict) -> dict:
    """Serialize the spectrum and accuracy summary without large sample arrays."""
    return {
        "multiplier": {
            "real": float(case["multiplier"].real),
            "imaginary": float(case["multiplier"].imag),
        },
        "continuous_eigenvalue": {
            "real": float(case["continuous_eigenvalue"].real),
            "imaginary": float(case["continuous_eigenvalue"].imag),
        },
        "training_summary": case["training_summary"],
        "validation_summary": case["validation_summary"],
    }


def clean_window(window: dict) -> dict:
    """Serialize the discovered-versus-prescribed summaries for one time window."""
    return {
        "name": window["name"],
        "start_tau": window["start_tau"],
        "degree": window["degree"],
        "candidate_count": window["candidate_count"],
        "discovered": clean_case(window["cases"]["discovered"]),
        "prescribed": clean_case(window["cases"]["prescribed"]),
    }


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    training = load_records(TRAINING_GLOB, "free_training")
    validation = load_records(VALIDATION_GLOB, "out_of_sample_validation")
    modes, singular, cumulative_energy = fit_hankel_basis(training)
    training_projected = [project_echo(record, modes, singular) for record in training]
    validation_projected = [project_echo(record, modes, singular) for record in validation]
    theoretical_lambda = theoretical_eigenvalue(training[0]["metadata"])

    configurations = (
        ("latent2", "coordinates", 2),
        ("latent4", "coordinates", 4),
        ("state", "state", None),
    )
    sweeps: dict[str, list[dict]] = {name: [] for name, _, _ in configurations}
    for start in START_TIMES:
        for name, key, rank in configurations:
            train_coordinates, train_runs = suffix_runs(
                training_projected, start, key, rank
            )
            validation_coordinates, validation_runs = suffix_runs(
                validation_projected, start, key, rank
            )
            sweeps[name].append(
                fit_window(
                    name,
                    start,
                    train_coordinates,
                    train_runs,
                    validation_coordinates,
                    validation_runs,
                    DEGREES[name],
                    theoretical_lambda,
                )
            )
            print(f"fit {name} from tau={start:g}")

    output = {
        "schema": "dynsyslab.rpe-viscous-echo-asymptotic-koopman/v1",
        "status": "isolated research experiment; not promoted to a KROM",
        "training_policy": "unforced only",
        "held_out_basis_policy": "fixed training U; no held-out SVD or refit",
        "viscosity_pa_s": 0.001,
        "delay_samples": DELAY_SAMPLES,
        "delay_window_tau": (DELAY_SAMPLES - 1) * EXPECTED_STEP * EMBEDDING_DECIMATION,
        "basis_cumulative_energy": cumulative_energy.tolist(),
        "theoretical_continuous_eigenvalue": {
            "real": float(theoretical_lambda.real),
            "imaginary": float(theoretical_lambda.imag),
        },
        "fit_start_times_tau": list(START_TIMES),
        "sweeps": {
            name: [clean_window(window) for window in windows]
            for name, windows in sweeps.items()
        },
    }
    METRICS.write_text(json.dumps(output, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        delay_modes=modes,
        delay_singular=singular,
        cumulative_energy=cumulative_energy,
        theoretical_continuous_eigenvalue=np.asarray(theoretical_lambda),
        metadata_json=np.asarray(json.dumps(output, separators=(",", ":"))),
    )

    styles = {
        "latent2": ("echo latent rank 2", "#0891b2"),
        "latent4": ("echo latent rank 4", "#7c3aed"),
        "state": ("state-space audit", "#f59e0b"),
    }
    starts = np.asarray(START_TIMES)
    summary = [
        figure(
            "line",
            "Discovered decay rate versus fit start",
            "fit starts at tau",
            "Re lambda",
            [
                trace(
                    label,
                    starts,
                    np.asarray(
                        [window["cases"]["discovered"]["continuous_eigenvalue"].real for window in sweeps[name]]
                    ),
                    color=color,
                    marker=True,
                    width=2.0,
                )
                for name, (label, color) in styles.items()
            ]
            + [
                trace(
                    "equilibrium linearization",
                    starts,
                    np.full(starts.shape, theoretical_lambda.real),
                    color="#111827",
                    width=1.5,
                )
            ],
        ),
        figure(
            "line",
            "Discovered angular frequency versus fit start",
            "fit starts at tau",
            "Im lambda",
            [
                trace(
                    label,
                    starts,
                    np.asarray(
                        [window["cases"]["discovered"]["continuous_eigenvalue"].imag for window in sweeps[name]]
                    ),
                    color=color,
                    marker=True,
                    width=2.0,
                )
                for name, (label, color) in styles.items()
            ]
            + [
                trace(
                    "equilibrium linearization",
                    starts,
                    np.full(starts.shape, theoretical_lambda.imag),
                    color="#111827",
                    width=1.5,
                )
            ],
        ),
        figure(
            "line",
            "Held-out four-tau error (discovered eigenfunction)",
            "fit starts at tau",
            "median relative L2 error",
            [
                trace(
                    label,
                    starts,
                    np.asarray(
                        [window["cases"]["discovered"]["validation_summary"]["median_four_tau"] for window in sweeps[name]]
                    ),
                    color=color,
                    marker=True,
                    width=2.0,
                )
                for name, (label, color) in styles.items()
            ],
        ),
        figure(
            "line",
            "Held-out four-tau error (prescribed asymptotic eigenvalue)",
            "fit starts at tau",
            "median relative L2 error",
            [
                trace(
                    label,
                    starts,
                    np.asarray(
                        [window["cases"]["prescribed"]["validation_summary"]["median_four_tau"] for window in sweeps[name]]
                    ),
                    color=color,
                    marker=True,
                    width=2.0,
                )
                for name, (label, color) in styles.items()
            ],
        ),
    ]
    render_viewer(
        Path(training[0]["source"]),
        summary,
        name=SUMMARY_PLOT,
        label="Viscous RPE asymptotic Koopman sweep",
        columns=2,
        width=1900,
        height=1900,
    )

    selected = sweeps["latent2"][-2]  # Fit from tau=50: nonlinear remnants remain.
    case = selected["cases"]["discovered"]
    continuous = case["continuous_eigenvalue"]
    evolution = []
    for index, (row, color) in enumerate(
        zip(case["validation_raw"][:4], COLORS[:4], strict=True), start=1
    ):
        sample = np.arange(0, row["time"].size, 20)
        prediction = row["values"][0] * np.exp(
            continuous * (row["time"] - row["time"][0])
        )
        evolution.append(
            figure(
                "line",
                f"Held-out IC #{index}: rank-2 psi, fit from tau=50",
                "tau",
                "Re psi",
                [
                    trace("observed", row["time"][sample], row["values"][sample].real, color=color, width=1.8),
                    trace("exponential prediction", row["time"][sample], prediction[sample].real, color="#111827", width=1.4),
                ],
            )
        )
    render_viewer(
        Path(training[0]["source"]),
        evolution,
        name=EVOLUTION_PLOT,
        label="Viscous echo-only asymptotic eigenfunction evolution",
        columns=2,
        width=1900,
        height=1900,
    )

    print(f"basis cumulative energy: {cumulative_energy}")
    print(f"theoretical lambda: {theoretical_lambda:.9g}")
    for name, windows in sweeps.items():
        print(name)
        for window in windows:
            case = window["cases"]["discovered"]
            metrics = case["validation_summary"]
            print(
                f"  start={window['start_tau']:g}: lambda={case['continuous_eigenvalue']:.7g}, "
                f"heldout four={metrics['median_four_tau']:.5g}, "
                f"anchor={metrics['median_anchored']:.5g}"
            )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
