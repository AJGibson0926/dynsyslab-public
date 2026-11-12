"""Validate and sparsify the successful orbit-local echo-only eigenfunction."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_echo_sampled_koopman import (
    DELAY_SAMPLES,
    PAIR_STRIDE,
    build_runs,
    candidate_table,
    coefficient_summary,
    eigenfunction_metrics,
    fit_edmd,
    orient_coefficient,
)
from analyze_rpe_inviscid_echo_hankel import far_field_echo, load_native_source
from compare_rpe_echo_pooled_training_injectivity import (
    ARTIFACT as HANKEL_ARTIFACT,
    METRICS as HANKEL_METRICS,
    load_training_records,
)
from project_rpe_inviscid_echo_basis_out_of_sample import (
    COLORS,
    EXPECTED_STEP,
    TRAINING,
    validate_validation_source,
)


ROOT = Path(__file__).resolve().parents[1]
PHASE_GLOB = (
    "rpe-inviscid-echo-local-eigenfunction-phase-validation-trajectory-*"
    "__20261001T235951Z.npz"
)
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-local-koopman-validation.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-local-koopman-validation.npz"
PLOT_NAME = "rpe-inviscid-echo-local-koopman-validation"
EVOLUTION_PLOT = "rpe-inviscid-echo-local-koopman-phase-replay"
THRESHOLDS = (0.001, 0.003, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2)


def load_phase_records() -> list[dict]:
    """Load the four expected phase-shift validation trajectories."""
    sources = sorted((ROOT / "archives/solutions/rpe").glob(PHASE_GLOB))
    if len(sources) != 4:
        raise RuntimeError(f"expected four phase-shift validation records, found {len(sources)}")
    records = []
    for source in sources:
        time, state, metadata = load_native_source(source)
        validate_validation_source(metadata)
        records.append(
            {
                "source": str(source),
                "time": time,
                "state": state,
                "initial_state": np.asarray(metadata["initial_state"], dtype=float),
                "echo": far_field_echo(state, metadata),
            }
        )
    return records


def fixed_multiplier_refit(
    model,
    training_coordinates: np.ndarray,
    multiplier: complex,
    active: np.ndarray,
) -> np.ndarray:
    """Refit active polynomial terms at a prescribed multiplier and normalize their values."""
    features = model.features(training_coordinates)
    minus = features[:-PAIR_STRIDE, active]
    plus = features[PAIR_STRIDE:, active]
    residual = plus - multiplier * minus
    numerator = residual.conj().T @ residual / residual.shape[0]
    denominator = minus.conj().T @ minus / minus.shape[0]
    ridge = 1.0e-12 * float(np.trace(denominator).real / denominator.shape[0])
    values, vectors = eigh(
        numerator,
        denominator + ridge * np.eye(denominator.shape[0]),
        check_finite=True,
    )
    local = vectors[:, int(np.argmin(values))]
    coefficient = np.zeros(features.shape[1], dtype=complex)
    coefficient[active] = local
    rms = float(np.sqrt(np.mean(np.abs(features @ coefficient) ** 2)))
    if rms <= np.finfo(float).tiny:
        raise RuntimeError("sparse refit produced a numerically zero eigenfunction")
    return orient_coefficient(
        model, coefficient / rms, training_coordinates
    )


def compact_metrics(rows: list[dict]) -> list[dict]:
    """Strip numerical sample arrays from the serializable trajectory metrics."""
    return [
        {
            key: value
            for key, value in row.items()
            if key not in {"values", "phase", "fitted_phase", "time"}
        }
        for row in rows
    ]


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    with np.load(HANKEL_ARTIFACT, allow_pickle=False) as archive:
        modes = np.asarray(archive[f"delay_{DELAY_SAMPLES}_modes"], dtype=float)
        singular = np.asarray(
            archive[f"delay_{DELAY_SAMPLES}_singular"], dtype=float
        )
    with HANKEL_METRICS.open("r", encoding="utf-8") as stream:
        hankel_metrics = json.load(stream)
    echo_mean = float(hankel_metrics["pooling"]["mean"])
    echo_scale = float(hankel_metrics["pooling"]["scale"])

    training_records, _ = load_training_records()
    training_runs = build_runs(
        [training_records[0]], modes, singular, echo_mean, echo_scale, "training"
    )
    phase_records = load_phase_records()
    phase_runs = build_runs(
        phase_records, modes, singular, echo_mean, echo_scale, "phase shift"
    )
    training_coordinates = training_runs[0]["coordinates"]
    model = fit_edmd([training_coordinates])
    candidates = candidate_table(model, [training_coordinates])
    if not candidates:
        raise RuntimeError("the retained local oscillatory mode was not recovered")
    multiplier = complex(candidates[0]["multiplier"])
    dense = orient_coefficient(
        model,
        np.asarray(candidates[0]["coefficient"], dtype=complex),
        training_coordinates,
    )
    dense_training = eigenfunction_metrics(
        model, dense, multiplier, training_runs
    )
    dense_validation = eigenfunction_metrics(
        model, dense, multiplier, phase_runs
    )

    variants = []
    maximum = float(np.max(np.abs(dense)))
    for threshold in THRESHOLDS:
        active = np.flatnonzero(np.abs(dense) >= threshold * maximum)
        if active.size < 2:
            continue
        coefficient = fixed_multiplier_refit(
            model, training_coordinates, multiplier, active
        )
        training_result = eigenfunction_metrics(
            model, coefficient, multiplier, training_runs
        )[0]
        variants.append(
            {
                "threshold": threshold,
                "active": active,
                "coefficient": coefficient,
                "training": training_result,
            }
        )

    dense_train = dense_training[0]
    acceptable = [
        item
        for item in variants
        if item["training"]["one_step_relative_l2"]
        <= max(0.005, 2.0 * dense_train["one_step_relative_l2"])
        and item["training"]["horizon_relative_l2"]["4"]
        <= max(0.01, 2.0 * dense_train["horizon_relative_l2"]["4"])
        and item["training"]["phase_linearity_r_squared"] >= 0.9999
    ]
    if not acceptable:
        raise RuntimeError("no sparse candidate preserved the training semigroup test")
    selected = min(
        acceptable,
        key=lambda item: (
            item["active"].size,
            item["training"]["horizon_relative_l2"]["4"],
        ),
    )
    sparse = np.asarray(selected["coefficient"], dtype=complex)
    dense_training_values = model.evaluate(training_coordinates, dense)
    sparse_training_values = model.evaluate(training_coordinates, sparse)
    alignment = np.vdot(sparse_training_values, dense_training_values) / np.vdot(
        sparse_training_values, sparse_training_values
    )
    sparse *= alignment
    selected["training"] = eigenfunction_metrics(
        model, sparse, multiplier, training_runs
    )[0]
    sparse_validation = eigenfunction_metrics(
        model, sparse, multiplier, phase_runs
    )

    measured_periods = np.asarray([run["echo_period_tau"] for run in phase_runs])
    period_relative_error = np.abs(
        measured_periods - training_runs[0]["echo_period_tau"]
    ) / training_runs[0]["echo_period_tau"]
    output = {
        "schema": "dynsyslab.echo-local-koopman-validation/v1",
        "status": "research experiment; not promoted to a KROM",
        "training_policy": "one unforced echo trajectory",
        "validation_policy": {
            "source_count": len(phase_runs),
            "forcing": "free",
            "construction": "fresh solves initialized at four phase points on the training orbit",
            "hidden_state_use": "validation-data construction and audit only; not fit or evaluation",
            "maximum_echo_period_relative_difference": float(np.max(period_relative_error)),
        },
        "representation": {
            "parent_hankel_training_sources": 5,
            "delay_samples": DELAY_SAMPLES,
            "latent_dimension": 2,
            "online_input": "echo history only",
        },
        "eigenfunction": {
            "multiplier": {
                "real": float(multiplier.real),
                "imaginary": float(multiplier.imag),
            },
            "pair_step_tau": model.pair_step_tau,
            "implied_period_tau": float(
                2.0 * np.pi / (np.angle(multiplier) / model.pair_step_tau)
            ),
            "dense_term_count": int(dense.size),
            "sparse_term_count": int(selected["active"].size),
            "sparse_threshold_relative_to_dense_maximum": selected["threshold"],
            "dense_coefficients": {
                "real": dense.real.tolist(),
                "imaginary": dense.imag.tolist(),
                "terms": coefficient_summary(model, dense, relative_floor=0.0),
            },
            "sparse_coefficients": {
                "real": sparse.real.tolist(),
                "imaginary": sparse.imag.tolist(),
                "terms": coefficient_summary(model, sparse, relative_floor=1.0e-14),
            },
        },
        "dense_training": compact_metrics(dense_training),
        "dense_phase_shift_validation": compact_metrics(dense_validation),
        "sparse_training": compact_metrics([selected["training"]]),
        "sparse_phase_shift_validation": compact_metrics(sparse_validation),
        "sparsity_sweep": [
            {
                "threshold": item["threshold"],
                "term_count": int(item["active"].size),
                "training_one_step_relative_l2": item["training"][
                    "one_step_relative_l2"
                ],
                "training_four_tau_relative_l2": item["training"][
                    "horizon_relative_l2"
                ]["4"],
            }
            for item in variants
        ],
    }
    METRICS.write_text(json.dumps(output, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        multiplier=np.asarray(multiplier),
        dense_coefficient=dense,
        sparse_coefficient=sparse,
        polynomial_exponents=np.asarray(model.exponents, dtype=int),
        coordinate_center=model.center,
        coordinate_scale=model.scale,
        feature_scale=model.feature_scale,
        phase_initial_states=np.vstack([run["initial_state"] for run in phase_runs]),
        metadata_json=np.asarray(json.dumps(output, separators=(",", ":"))),
    )

    horizons = np.asarray([0.01, 0.1, 0.5, 1.0, 2.0, 4.0])
    dense_horizon_median = np.asarray(
        [
            np.median([row["horizon_relative_l2"][f"{value:g}"] for row in dense_validation])
            for value in horizons
        ]
    )
    sparse_horizon_median = np.asarray(
        [
            np.median([row["horizon_relative_l2"][f"{value:g}"] for row in sparse_validation])
            for value in horizons
        ]
    )
    figures = [
        figure(
            "line",
            "Independent same-orbit semigroup error",
            "prediction horizon (tau)",
            "median relative L2 error",
            [
                trace("dense", horizons, dense_horizon_median, color="#2563eb", marker=True, width=2.0),
                trace("sparse", horizons, sparse_horizon_median, color="#dc2626", marker=True, width=2.0),
            ],
        ),
        figure(
            "line",
            "Four-tau error by phase shift",
            "phase-shift validation index",
            "relative L2 error",
            [
                trace(
                    "dense",
                    np.arange(1, 5),
                    np.asarray([row["horizon_relative_l2"]["4"] for row in dense_validation]),
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "sparse",
                    np.arange(1, 5),
                    np.asarray([row["horizon_relative_l2"]["4"] for row in sparse_validation]),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "line",
            "Training sparsity tradeoff",
            "retained polynomial terms",
            "four-tau relative L2 error",
            [
                trace(
                    "fixed-multiplier refit",
                    np.asarray([item["active"].size for item in variants]),
                    np.asarray([item["training"]["horizon_relative_l2"]["4"] for item in variants]),
                    color="#7c3aed",
                    marker=True,
                    width=2.0,
                )
            ],
        ),
        figure(
            "line",
            "Independent same-orbit periods",
            "phase-shift validation index",
            "echo period (tau)",
            [
                trace(
                    "measured",
                    np.arange(1, 5),
                    measured_periods,
                    color="#0891b2",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "eigenfunction implied",
                    np.arange(1, 5),
                    np.full(4, output["eigenfunction"]["implied_period_tau"]),
                    color="#111827",
                    width=1.8,
                ),
            ],
        ),
    ]
    render_viewer(
        TRAINING,
        figures,
        name=PLOT_NAME,
        label="Orbit-local echo eigenfunction validation",
        columns=2,
        width=1900,
        height=1900,
    )

    replay_figures = []
    continuous = np.log(multiplier) / model.pair_step_tau
    for index, (dense_row, sparse_row, color) in enumerate(
        zip(dense_validation, sparse_validation, COLORS[:4], strict=True), start=1
    ):
        sample = np.arange(0, dense_row["time"].size, 8)
        predicted = dense_row["values"][0] * np.exp(
            continuous * (dense_row["time"] - dense_row["time"][0])
        )
        sparse_predicted = sparse_row["values"][0] * np.exp(
            continuous * (sparse_row["time"] - sparse_row["time"][0])
        )
        replay_figures.append(
            figure(
                "line",
                f"Phase shift #{index}: local eigenfunction replay",
                "tau",
                "Re psi",
                [
                    trace(
                        "dense observed",
                        dense_row["time"][sample],
                        dense_row["values"][sample].real,
                        color=color,
                        width=1.8,
                    ),
                    trace(
                        "sparse observed",
                        sparse_row["time"][sample],
                        sparse_row["values"][sample].real,
                        color="#dc2626",
                        width=1.4,
                    ),
                    trace(
                        "dense exponential prediction",
                        dense_row["time"][sample],
                        predicted[sample].real,
                        color="#111827",
                        width=1.2,
                    ),
                    trace(
                        "sparse exponential prediction",
                        sparse_row["time"][sample],
                        sparse_predicted[sample].real,
                        color="#f59e0b",
                        width=1.0,
                    ),
                ],
            )
        )
    render_viewer(
        TRAINING,
        replay_figures,
        name=EVOLUTION_PLOT,
        label="Orbit-local echo eigenfunction phase replay",
        columns=2,
        width=1900,
        height=1900,
    )

    print(
        f"dense validation median step/four-tau="
        f"{np.median([row['one_step_relative_l2'] for row in dense_validation]):.6g}/"
        f"{np.median([row['horizon_relative_l2']['4'] for row in dense_validation]):.6g}"
    )
    print(
        f"sparse terms={selected['active'].size}/{dense.size}, threshold={selected['threshold']}, "
        f"validation median step/four-tau="
        f"{np.median([row['one_step_relative_l2'] for row in sparse_validation]):.6g}/"
        f"{np.median([row['horizon_relative_l2']['4'] for row in sparse_validation]):.6g}"
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
