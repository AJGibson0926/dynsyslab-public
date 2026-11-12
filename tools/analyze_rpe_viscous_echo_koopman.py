"""Compare viscous RPE Koopman fits in echo-latent and physical coordinates."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.signal import fftconvolve
from scipy.sparse.linalg import svds

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_echo_sampled_koopman import (
    PAIR_STRIDE,
    candidate_table,
    echo_period,
    eigenfunction_metrics,
    fit_edmd,
    orient_coefficient,
)
from analyze_rpe_inviscid_echo_hankel import far_field_echo, load_native_source
from compare_rpe_echo_pooled_training_injectivity import PooledHankelOperator
from project_rpe_inviscid_echo_basis_out_of_sample import (
    COLORS,
    EMBEDDING_DECIMATION,
    EXPECTED_STEP,
)


ROOT = Path(__file__).resolve().parents[1]
TRAINING_GLOB = (
    "rpe-viscous-echo-pooled-training-five-ics-trajectory-*"
    "__20261002T000717Z.npz"
)
VALIDATION_GLOB = (
    "rpe-viscous-echo-validation-five-ics-trajectory-*"
    "__20261002T000717Z.npz"
)
DELAY_SAMPLES = 8644
LATENT_RANK = 2
METRICS = OUTPUT_ROOT / "rpe-viscous-echo-sampled-koopman.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-viscous-echo-sampled-koopman.npz"
SUMMARY_PLOT = "rpe-viscous-echo-sampled-koopman-summary"
EVOLUTION_PLOT = "rpe-viscous-echo-sampled-koopman-evolution"
CONSERVATIVE_METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-sampled-koopman.metrics.json"


def validate_source(metadata: dict, role: str) -> None:
    """Require unforced uncontrolled data with viscosity exactly 0.001 Pa s."""
    if metadata.get("forcing_type") != "free":
        raise RuntimeError("viscous Koopman data must be unforced")
    if metadata.get("controller_kind") != "none":
        raise RuntimeError("viscous Koopman data contain control")
    if float(metadata["resolved_system_parameters"]["mu_pa_s"]) != 0.001:
        raise RuntimeError("viscous comparison requires mu_pa_s=0.001 exactly")
    if metadata.get("data_role") != role:
        raise RuntimeError(f"unexpected data role: {metadata.get('data_role')}")


def load_records(pattern: str, role: str) -> list[dict]:
    """Load the five expected free RPE records and prepare their echo histories."""
    sources = sorted((ROOT / "archives/solutions/rpe").glob(pattern))
    if len(sources) != 5:
        raise RuntimeError(f"expected five {role} sources, found {len(sources)}")
    records = []
    for source in sources:
        time, state, metadata = load_native_source(source)
        validate_source(metadata, role)
        records.append(
            {
                "source": str(source),
                "time": time,
                "state": state,
                "initial_state": np.asarray(metadata["initial_state"], dtype=float),
                "echo": far_field_echo(state, metadata),
                "metadata": metadata,
            }
        )
    return records


def fit_hankel_basis(
    records: list[dict], mean: float, scale: float
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a pooled Hankel basis using the supplied frozen echo mean and scale."""
    signals = [
        (record["echo"][::EMBEDDING_DECIMATION] - mean) / scale
        for record in records
    ]
    operator = PooledHankelOperator(signals, DELAY_SAMPLES)
    left, singular, right = svds(
        operator,
        k=LATENT_RANK,
        which="LM",
        solver="arpack",
        tol=1.0e-10,
        maxiter=10000,
        return_singular_vectors=True,
        rng=np.random.default_rng(20261002),
    )
    order = np.argsort(singular)[::-1]
    left = left[:, order]
    singular = singular[order]
    right = right[order]
    for mode in range(LATENT_RANK):
        pivot = int(np.argmax(np.abs(left[:, mode])))
        if left[pivot, mode] < 0.0:
            left[:, mode] *= -1.0
            right[mode] *= -1.0
    return left, singular


def project_echo(
    echo: np.ndarray,
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
) -> np.ndarray:
    """Apply the frozen Hankel modes to causal echo windows and align the results."""
    normalized = (echo[::EMBEDDING_DECIMATION] - mean) / scale
    return np.column_stack(
        [
            fftconvolve(normalized, modes[:, mode][::-1], mode="valid")
            / singular[mode]
            for mode in range(LATENT_RANK)
        ]
    )


def prepare_runs(
    records: list[dict],
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
    prefix: str,
) -> tuple[list[dict], list[np.ndarray]]:
    """Project each echo history and attach aligned timestamps and diagnostic labels."""
    latent_runs = []
    state_runs = []
    for index, record in enumerate(records, start=1):
        coordinates = project_echo(record["echo"], modes, singular, mean, scale)
        time = record["time"][::EMBEDDING_DECIMATION][DELAY_SAMPLES - 1 :]
        state = record["state"][::EMBEDDING_DECIMATION][DELAY_SAMPLES - 1 :]
        latent_runs.append(
            {
                "label": f"{prefix} #{index}",
                "initial_state": record["initial_state"],
                "coordinates": coordinates,
                "latent_time": time,
                "echo_period_tau": echo_period(record["time"], record["echo"]),
            }
        )
        state_runs.append(state)
    return latent_runs, state_runs


def theoretical_eigenvalue(metadata: dict) -> complex:
    """Compute the positive-frequency eigenvalue of the RPE equilibrium linearization."""
    groups = metadata["dimensionless_groups"]
    parameters = metadata["resolved_system_parameters"]
    derivative = (
        -3.0 * float(parameters["kappa"]) * float(groups["gas_pressure_ratio"])
        + float(groups["capillarity"])
    )
    jacobian = np.asarray(
        [[0.0, 1.0], [derivative, -float(groups["viscosity"])]], dtype=float
    )
    values = np.linalg.eigvals(jacobian)
    positive = values[values.imag > 0.0]
    if positive.size != 1:
        raise RuntimeError("viscous equilibrium is not the expected stable focus")
    return complex(positive[0])


def prescribed_refit(
    model,
    coordinate_runs: list[np.ndarray],
    multiplier: complex,
) -> np.ndarray:
    """Fit polynomial coefficients at a fixed multiplier using the normalized sampled defect."""
    feature_runs = [model.features(run) for run in coordinate_runs]
    minus = np.vstack([values[:-PAIR_STRIDE] for values in feature_runs])
    plus = np.vstack([values[PAIR_STRIDE:] for values in feature_runs])
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
    all_coordinates = np.vstack(coordinate_runs)
    rms = float(np.sqrt(np.mean(np.abs(model.evaluate(all_coordinates, coefficient)) ** 2)))
    coefficient /= rms
    return orient_coefficient(model, coefficient, coordinate_runs[0])


def state_metric_runs(template: list[dict], states: list[np.ndarray]) -> list[dict]:
    """Substitute aligned physical-state coordinates for the coordinate audit."""
    result = []
    for source, state in zip(template, states, strict=True):
        result.append({**source, "coordinates": state})
    return result


def anchored_error(model, coefficient: np.ndarray, multiplier: complex, run: dict) -> float:
    """Measure relative L2 departure from the exponential anchored at the first sample."""
    values = model.evaluate(run["coordinates"], coefficient)
    continuous = np.log(multiplier) / model.pair_step_tau
    elapsed = run["latent_time"] - run["latent_time"][0]
    prediction = values[0] * np.exp(continuous * elapsed)
    return float(
        np.linalg.norm(values - prediction)
        / max(np.linalg.norm(values), np.finfo(float).tiny)
    )


def serializable_metrics(rows: list[dict], anchored: list[float]) -> list[dict]:
    """Remove sample arrays from per-trajectory error summaries."""
    result = []
    for row, anchor in zip(rows, anchored, strict=True):
        result.append(
            {
                **{
                    key: value
                    for key, value in row.items()
                    if key not in {"values", "phase", "fitted_phase", "time"}
                },
                "anchored_full_record_relative_l2": anchor,
            }
        )
    return result


def fit_case(
    training_coordinates: list[np.ndarray],
    training_runs: list[dict],
    validation_coordinates: list[np.ndarray],
    validation_runs: list[dict],
    theoretical_lambda: complex,
) -> dict:
    """Fit and compare eigenfunction candidates for a specified coordinate representation."""
    model = fit_edmd(training_coordinates)
    candidates = candidate_table(model, training_coordinates)
    if not candidates:
        raise RuntimeError("no fundamental viscous EDMD mode was found")
    discovered_multiplier = complex(candidates[0]["multiplier"])
    discovered_coefficient = orient_coefficient(
        model,
        np.asarray(candidates[0]["coefficient"], dtype=complex),
        training_coordinates[0],
    )
    prescribed_multiplier = complex(
        np.exp(theoretical_lambda * model.pair_step_tau)
    )
    prescribed_coefficient = prescribed_refit(
        model, training_coordinates, prescribed_multiplier
    )

    cases = {}
    for name, multiplier, coefficient in (
        ("discovered", discovered_multiplier, discovered_coefficient),
        ("prescribed", prescribed_multiplier, prescribed_coefficient),
    ):
        train = eigenfunction_metrics(
            model, coefficient, multiplier, training_runs
        )
        validation = eigenfunction_metrics(
            model, coefficient, multiplier, validation_runs
        )
        train_anchor = [
            anchored_error(model, coefficient, multiplier, run) for run in training_runs
        ]
        validation_anchor = [
            anchored_error(model, coefficient, multiplier, run)
            for run in validation_runs
        ]
        continuous = np.log(multiplier) / model.pair_step_tau
        cases[name] = {
            "multiplier": multiplier,
            "continuous_eigenvalue": continuous,
            "coefficient": coefficient,
            "training_raw": train,
            "validation_raw": validation,
            "training": serializable_metrics(train, train_anchor),
            "validation": serializable_metrics(validation, validation_anchor),
        }
    return {"model": model, "candidates": candidates, **cases}


def complex_value(value: complex) -> dict:
    """Convert a complex scalar to real and imaginary JSON fields."""
    return {"real": float(value.real), "imaginary": float(value.imag)}


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    training_records = load_records(TRAINING_GLOB, "free_training")
    validation_records = load_records(VALIDATION_GLOB, "out_of_sample_validation")
    sampled_training_echo = np.concatenate(
        [record["echo"][::EMBEDDING_DECIMATION] for record in training_records]
    )
    echo_mean = float(np.mean(sampled_training_echo))
    echo_scale = float(np.std(sampled_training_echo))
    modes, singular = fit_hankel_basis(training_records, echo_mean, echo_scale)
    latent_training, state_training = prepare_runs(
        training_records, modes, singular, echo_mean, echo_scale, "training"
    )
    latent_validation, state_validation = prepare_runs(
        validation_records, modes, singular, echo_mean, echo_scale, "validation"
    )
    theoretical_lambda = theoretical_eigenvalue(training_records[0]["metadata"])

    latent = fit_case(
        [run["coordinates"] for run in latent_training],
        latent_training,
        [run["coordinates"] for run in latent_validation],
        latent_validation,
        theoretical_lambda,
    )
    state_training_runs = state_metric_runs(latent_training, state_training)
    state_validation_runs = state_metric_runs(latent_validation, state_validation)
    state = fit_case(
        state_training,
        state_training_runs,
        state_validation,
        state_validation_runs,
        theoretical_lambda,
    )

    with CONSERVATIVE_METRICS.open("r", encoding="utf-8") as stream:
        conservative = json.load(stream)
    conservative_validation = conservative["pooled_orbit"]["eigenfunction"][
        "per_trajectory"
    ][5:]

    def case_summary(case: dict) -> dict:
        return {
            "multiplier": complex_value(case["multiplier"]),
            "continuous_eigenvalue": complex_value(case["continuous_eigenvalue"]),
            "coefficient_real": case["coefficient"].real.tolist(),
            "coefficient_imaginary": case["coefficient"].imag.tolist(),
            "training": case["training"],
            "validation": case["validation"],
        }

    output = {
        "schema": "dynsyslab.rpe-viscous-echo-koopman-comparison/v1",
        "status": "research experiment; not promoted to a KROM",
        "training_policy": "unforced only",
        "viscosity_pa_s": 0.001,
        "representation": {
            "training_trajectories": 5,
            "validation_trajectories": 5,
            "delay_samples": DELAY_SAMPLES,
            "delay_window_tau": (DELAY_SAMPLES - 1)
            * EXPECTED_STEP
            * EMBEDDING_DECIMATION,
            "latent_rank": LATENT_RANK,
            "echo_mean": echo_mean,
            "echo_scale": echo_scale,
        },
        "equilibrium_linearization_eigenvalue": complex_value(theoretical_lambda),
        "echo_latent": {
            "discovered": case_summary(latent["discovered"]),
            "prescribed": case_summary(latent["prescribed"]),
        },
        "state_space_audit": {
            "hidden_state_used_for_echo_fit": False,
            "discovered": case_summary(state["discovered"]),
            "prescribed": case_summary(state["prescribed"]),
        },
        "conservative_reference": {
            "validation_median_one_step": float(
                np.median([row["one_step_relative_l2"] for row in conservative_validation])
            ),
            "validation_median_four_tau": float(
                np.median([row["horizon_relative_l2"]["4"] for row in conservative_validation])
            ),
            "validation_median_anchored_full_record": None,
        },
    }
    METRICS.write_text(json.dumps(output, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        delay_modes=modes,
        delay_singular=singular,
        echo_mean=np.asarray(echo_mean),
        echo_scale=np.asarray(echo_scale),
        theoretical_continuous_eigenvalue=np.asarray(theoretical_lambda),
        latent_discovered_multiplier=np.asarray(latent["discovered"]["multiplier"]),
        latent_discovered_coefficient=latent["discovered"]["coefficient"],
        latent_prescribed_multiplier=np.asarray(latent["prescribed"]["multiplier"]),
        latent_prescribed_coefficient=latent["prescribed"]["coefficient"],
        state_discovered_multiplier=np.asarray(state["discovered"]["multiplier"]),
        state_discovered_coefficient=state["discovered"]["coefficient"],
        state_prescribed_multiplier=np.asarray(state["prescribed"]["multiplier"]),
        state_prescribed_coefficient=state["prescribed"]["coefficient"],
        metadata_json=np.asarray(json.dumps(output, separators=(",", ":"))),
    )

    labels = ["echo EDMD", "echo prescribed", "state EDMD", "state prescribed"]
    cases = [latent["discovered"], latent["prescribed"], state["discovered"], state["prescribed"]]
    eigenvalues = np.asarray([case["continuous_eigenvalue"] for case in cases])
    validation_x = np.arange(1, 6)
    summary = [
        figure(
            "phase",
            "Continuous eigenvalue comparison",
            "Re lambda",
            "Im lambda",
            [
                trace(label, np.asarray([value.real]), np.asarray([value.imag]), color=color, marker=True, width=0.0, marker_size=7.0)
                for label, value, color in zip(labels, eigenvalues, ("#0891b2", "#2563eb", "#f59e0b", "#dc2626"), strict=True)
            ]
            + [
                trace("equilibrium linearization", np.asarray([theoretical_lambda.real]), np.asarray([theoretical_lambda.imag]), color="#111827", marker=True, width=0.0, marker_size=8.0)
            ],
        ),
        figure(
            "line",
            "Held-out four-tau semigroup error",
            "validation trajectory index",
            "relative L2 error",
            [
                trace(
                    label,
                    validation_x,
                    np.asarray([row["horizon_relative_l2"]["4"] for row in case["validation"]]),
                    color=color,
                    marker=True,
                    width=2.0,
                )
                for label, case, color in zip(labels, cases, ("#0891b2", "#2563eb", "#f59e0b", "#dc2626"), strict=True)
            ],
        ),
        figure(
            "line",
            "Held-out anchored full-record error",
            "validation trajectory index",
            "relative L2 error",
            [
                trace(
                    label,
                    validation_x,
                    np.asarray([row["anchored_full_record_relative_l2"] for row in case["validation"]]),
                    color=color,
                    marker=True,
                    width=2.0,
                )
                for label, case, color in zip(labels, cases, ("#0891b2", "#2563eb", "#f59e0b", "#dc2626"), strict=True)
            ],
        ),
        figure(
            "line",
            "Conservative versus viscous echo-latent error",
            "comparison",
            "median held-out relative L2 error",
            [
                trace(
                    "one sampled step",
                    np.asarray([1, 2]),
                    np.asarray(
                        [
                            output["conservative_reference"]["validation_median_one_step"],
                            np.median([row["one_step_relative_l2"] for row in latent["discovered"]["validation"]]),
                        ]
                    ),
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "four tau",
                    np.asarray([1, 2]),
                    np.asarray(
                        [
                            output["conservative_reference"]["validation_median_four_tau"],
                            np.median([row["horizon_relative_l2"]["4"] for row in latent["discovered"]["validation"]]),
                        ]
                    ),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
    ]
    summary[0]["square"] = True
    render_viewer(
        Path(training_records[0]["source"]),
        summary,
        name=SUMMARY_PLOT,
        label="Viscous echo-only Koopman comparison",
        columns=2,
        width=1900,
        height=1900,
    )

    evolution = []
    case = latent["discovered"]
    continuous = case["continuous_eigenvalue"]
    for index, (run, row, color) in enumerate(
        zip(latent_validation[:4], case["validation_raw"][:4], COLORS[:4], strict=True), start=1
    ):
        sample = np.arange(0, row["time"].size, 8)
        prediction = row["values"][0] * np.exp(
            continuous * (row["time"] - row["time"][0])
        )
        evolution.append(
            figure(
                "line",
                f"Viscous held-out IC #{index}: eigenfunction evolution",
                "tau",
                "Re psi",
                [
                    trace("observed", row["time"][sample], row["values"][sample].real, color=color, width=1.8),
                    trace("exponential prediction", row["time"][sample], prediction[sample].real, color="#111827", width=1.4),
                ],
            )
        )
    render_viewer(
        Path(training_records[0]["source"]),
        evolution,
        name=EVOLUTION_PLOT,
        label="Viscous echo-only eigenfunction evolution",
        columns=2,
        width=1900,
        height=1900,
    )

    for coordinate, result in (("echo latent", latent), ("state audit", state)):
        print(coordinate)
        for name in ("discovered", "prescribed"):
            case = result[name]
            print(
                f"  {name}: lambda={case['continuous_eigenvalue']:.8g}, "
                f"validation median step/four/anchored="
                f"{np.median([row['one_step_relative_l2'] for row in case['validation']]):.5g}/"
                f"{np.median([row['horizon_relative_l2']['4'] for row in case['validation']]):.5g}/"
                f"{np.median([row['anchored_full_record_relative_l2'] for row in case['validation']]):.5g}"
            )
    print(f"theoretical lambda={theoretical_lambda:.8g}")
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
