"""Formally validate the late-time echo-Hankel RPE eigenfunction.

The experimental fit uses canonical eigenfunction/error contracts but is not
promoted to a KROM because its coordinates require a causal echo ring buffer.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_viscous_echo_asymptotic_koopman import (
    VALIDATION_GLOB,
    TRAINING_GLOB,
    discovered_candidate,
    fit_edmd,
    load_records,
    project_echo,
    suffix_runs,
)
from spectral.contracts import (
    DynamicalSystemReference,
    KoopmanEigenfunction,
    KoopmanEigenvalue,
)
from spectral.errors import FreeTrajectoryDataset
from time_series import TimeSeries
from rpe_echo_eigenfunction_viewer import (
    eigenfunction_history,
    heat_trace,
    trajectory_overlays,
)


PARENT_ARTIFACT = OUTPUT_ROOT / "rpe-viscous-echo-asymptotic-koopman.npz"
FIT_START = 50.0
FULL_AUDIT_START = 5.0
FIT_END = 80.0
LATENT_RANK = 2
POLYNOMIAL_DEGREE = 5
FIELD_POINTS = 121
METRICS = OUTPUT_ROOT / "rpe-viscous-echo-asymptotic-eigenfunction-domain.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-viscous-echo-asymptotic-eigenfunction-domain.npz"
FIELD_PLOT = "rpe-viscous-echo-asymptotic-eigenfunction-field"
ACCURACY_PLOT = "rpe-viscous-echo-asymptotic-eigenfunction-accuracy"
def polynomial_gradient(model, coefficient: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    """Differentiate the normalized polynomial with respect to its input coordinates."""
    points = np.asarray(coordinates, dtype=float)
    if points.ndim == 1:
        points = points[None, :]
    normalized = (points - model.center) / model.scale
    result = np.zeros((points.shape[0], points.shape[1]), dtype=complex)
    for column, exponent in enumerate(model.exponents):
        weight = coefficient[column] / model.feature_scale[column]
        for axis, power in enumerate(exponent):
            if not power:
                continue
            term = np.full(points.shape[0], weight * power / model.scale[axis], dtype=complex)
            for factor_axis, factor_power in enumerate(exponent):
                effective = factor_power - (1 if factor_axis == axis else 0)
                if effective:
                    term *= normalized[:, factor_axis] ** effective
            result[:, axis] += term
    return result
def canonical_eigenfunction(model, coefficient: np.ndarray, multiplier: complex, metadata: dict):
    """Wrap fitted coefficients and spectrum in the standard reference eigenfunction object."""
    reference = DynamicalSystemReference(
        metadata.get("source_archive"),
        metadata.get("system", "rpe"),
        dict(metadata.get("resolved_system_parameters") or {}),
    )
    eigenvalue = KoopmanEigenvalue.from_discrete(
        multiplier,
        model.pair_step_tau,
        reference,
        {"method": "sampled polynomial EDMD in fixed echo-Hankel coordinates"},
    )

    def psi(value: np.ndarray) -> complex:
        return complex(model.evaluate(np.asarray(value, dtype=float).reshape(1, -1), coefficient)[0])

    def gradient(value: np.ndarray) -> np.ndarray:
        return polynomial_gradient(model, coefficient, value)[0]

    def residual(points: np.ndarray, velocities: np.ndarray) -> np.ndarray:
        x = np.asarray(points, dtype=float)
        velocity = np.asarray(velocities, dtype=float)
        values = model.evaluate(x, coefficient)
        return np.einsum("ij,ij->i", polynomial_gradient(model, coefficient, x), velocity) - eigenvalue.continuous_time * values

    return KoopmanEigenfunction(
        eigenvalue,
        "degree-5 polynomial in causal echo-Hankel coordinates",
        coefficients=np.asarray(coefficient, dtype=complex),
        state_callable=psi,
        gradient_callable=gradient,
        residual_operator=residual,
        representation_kind="analytic",
        metadata={
            "method": "isolated sampled-EDMD experiment",
            "coordinate_names": ["xi1", "xi2"],
            "coordinate_source": "fixed causal echo ring-buffer projection",
            "history_dependent_coordinates": True,
            "fit_start_tau": FIT_START,
            "fit_end_tau": FIT_END,
            "polynomial_degree": POLYNOMIAL_DEGREE,
            "canonical_krom_promotion": False,
        },
    )


def free_dataset(projected: list[dict], start: float, *, label: str, role: str) -> FreeTrajectoryDataset:
    """Build a canonical free-trajectory validation dataset over the chosen time interval."""
    trajectories = []
    for record in projected:
        mask = (record["time"] >= start) & (record["time"] <= FIT_END)
        trajectories.append(
            TimeSeries(record["time"][mask], record["coordinates"][mask, :LATENT_RANK])
        )
    return FreeTrajectoryDataset(
        trajectories,
        label=label,
        role=role,
        coordinate_system="model",
        provenance={
            "coordinates": "xi1,xi2 from fixed training-only echo Hankel basis",
            "fit_start_tau": FIT_START,
            "source_records": [record["source"] for record in projected],
        },
    )


def binned_defect(
    model,
    coefficient: np.ndarray,
    multiplier: complex,
    projected: list[dict],
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    horizon: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Bin logarithmic normalized finite-horizon defects over occupied coordinate cells."""
    step = float(projected[0]["time"][1] - projected[0]["time"][0])
    shift = int(round(horizon / step))
    continuous = np.log(multiplier) / model.pair_step_tau
    all_points = []
    all_errors = []
    all_weights = []
    for record in projected:
        mask = (record["time"] >= FIT_START) & (record["time"] <= FIT_END)
        coordinates = record["coordinates"][mask, :LATENT_RANK]
        values = model.evaluate(coordinates, coefficient)
        if shift >= values.size:
            continue
        actual = values[shift:]
        predicted = np.exp(continuous * horizon) * values[:-shift]
        scale = max(float(np.max(np.abs(actual))), np.finfo(float).tiny)
        floor = np.sqrt(np.finfo(float).eps) * scale
        error = np.abs(actual - predicted) / np.maximum(np.abs(actual), floor)
        all_points.append(coordinates[:-shift])
        all_errors.append(error)
        all_weights.append(np.abs(actual))
    points = np.vstack(all_points)
    errors = np.concatenate(all_errors)
    weights = np.concatenate(all_weights)
    columns = x_edges.size - 1
    rows = y_edges.size - 1
    x_index = np.clip(np.digitize(points[:, 0], x_edges) - 1, 0, columns - 1)
    y_index = np.clip(np.digitize(points[:, 1], y_edges) - 1, 0, rows - 1)
    flat = y_index * columns + x_index
    numerator = np.zeros(rows * columns)
    denominator = np.zeros(rows * columns)
    np.add.at(numerator, flat, errors * weights)
    np.add.at(denominator, flat, weights)
    field = np.full(rows * columns, np.nan)
    occupied = denominator > 0.0
    field[occupied] = np.log10(np.maximum(numerator[occupied] / denominator[occupied], 1.0e-12))
    return field.reshape(rows, columns), occupied.reshape(rows, columns)


def segment_error(
    model,
    coefficient: np.ndarray,
    multiplier: complex,
    projected: list[dict],
    starts: np.ndarray,
    horizon: float,
) -> np.ndarray:
    """Compute relative L2 exponential-prediction errors within successive time segments."""
    step = float(projected[0]["time"][1] - projected[0]["time"][0])
    shift = int(round(horizon / step))
    continuous = np.log(multiplier) / model.pair_step_tau
    errors = []
    for lower, upper in zip(starts[:-1], starts[1:], strict=True):
        actual_parts = []
        predicted_parts = []
        for record in projected:
            mask = (record["time"] >= lower) & (record["time"] <= upper)
            coordinates = record["coordinates"][mask, :LATENT_RANK]
            values = model.evaluate(coordinates, coefficient)
            if shift >= values.size:
                continue
            actual_parts.append(values[shift:])
            predicted_parts.append(np.exp(continuous * horizon) * values[:-shift])
        actual = np.concatenate(actual_parts)
        predicted = np.concatenate(predicted_parts)
        errors.append(float(np.linalg.norm(actual - predicted) / np.linalg.norm(actual)))
    return np.asarray(errors)


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    with np.load(PARENT_ARTIFACT, allow_pickle=False) as parent:
        modes = np.asarray(parent["delay_modes"], dtype=float)
        singular = np.asarray(parent["delay_singular"], dtype=float)
    training_records = load_records(TRAINING_GLOB, "free_training")
    validation_records = load_records(VALIDATION_GLOB, "out_of_sample_validation")
    training = [project_echo(record, modes, singular) for record in training_records]
    validation = [project_echo(record, modes, singular) for record in validation_records]

    train_coordinates, train_runs = suffix_runs(training, FIT_START, "coordinates", LATENT_RANK)
    validation_coordinates, validation_runs = suffix_runs(validation, FIT_START, "coordinates", LATENT_RANK)
    model = fit_edmd(train_coordinates, POLYNOMIAL_DEGREE)
    multiplier, coefficient, candidates = discovered_candidate(model, train_coordinates)
    function = canonical_eigenfunction(
        model,
        coefficient,
        multiplier,
        training_records[0]["metadata"],
    )

    combined_training = np.vstack(train_coordinates)
    lower = np.min(combined_training, axis=0)
    upper = np.max(combined_training, axis=0)
    span = upper - lower
    field_lower = lower - 0.04 * span
    field_upper = upper + 0.04 * span
    training_domain = {
        "encoding": "axis-aligned-observable-coordinate-domain/v1",
        "minimum": lower.tolist(),
        "maximum": upper.tolist(),
        "sample_count": int(combined_training.shape[0]),
        "coordinate_names": ["xi1", "xi2"],
        "source": "five unforced training suffixes, tau 50 through 80",
    }
    training_evaluation = function.evaluate_errors(
        free_dataset(training, FIT_START, label="viscous late training", role="training"),
        horizons=(0.01, 0.1, 0.5, 1.0, 2.0, 4.0),
        training_domain=training_domain,
        metadata={"history_basis_refit": False},
    )
    validation_evaluation = function.evaluate_errors(
        free_dataset(validation, FIT_START, label="viscous late held-out", role="test"),
        horizons=(0.01, 0.1, 0.5, 1.0, 2.0, 4.0),
        training_domain=training_domain,
        metadata={"history_basis_refit": False, "eigenfunction_refit": False},
    )

    x = np.linspace(field_lower[0], field_upper[0], FIELD_POINTS)
    y = np.linspace(field_lower[1], field_upper[1], FIELD_POINTS)
    xx, yy = np.meshgrid(x, y)
    grid_coordinates = np.column_stack((xx.reshape(-1), yy.reshape(-1)))
    psi = model.evaluate(grid_coordinates, coefficient).reshape(xx.shape)
    overlays = trajectory_overlays(training, FIT_START, FIT_END)
    field_figures = []
    for title, values, color_label, history_component in (
        (
            "Real part of late-time eigenfunction",
            psi.real,
            r"\operatorname{Re}(\psi)",
            "real",
        ),
        (
            "Imaginary part of late-time eigenfunction",
            psi.imag,
            r"\operatorname{Im}(\psi)",
            "imaginary",
        ),
    ):
        item = figure(
            "krom_heatmap",
            title,
            r"\xi_{1}",
            r"\xi_{2}",
            [heat_trace(color_label, x, y, values), *overlays],
            color_label=color_label,
        )
        item["square"] = True
        item["value_scale"] = "linear"
        field_figures.append(item)
        field_figures.append(
            eigenfunction_history(
                model,
                coefficient,
                training,
                FIT_START,
                FIT_END,
                component=history_component,
            )
        )
    normalized_magnitude = np.abs(psi) / max(float(np.max(np.abs(psi))), np.finfo(float).tiny)
    phase_item = figure(
        "krom_heatmap",
        "Late-time eigenfunction: phase (hue), magnitude (brightness)",
        r"\xi_{1}",
        r"\xi_{2}",
        [
            heat_trace(
                r"\operatorname{arg}(\psi)",
                x,
                y,
                np.angle(psi),
                encoding="domain-color",
                magnitudes=normalized_magnitude,
            ),
            *overlays,
        ],
        color_label=r"\operatorname{arg}(\psi)\;\mathrm{(hue)},\;|\psi|\;\mathrm{(brightness)}",
    )
    phase_item["square"] = True
    phase_item["value_scale"] = "linear"
    phase_item["row_span"] = 2
    field_figures.append(phase_item)
    field_figures.append(
        eigenfunction_history(
            model,
            coefficient,
            training,
            FIT_START,
            FIT_END,
            component="magnitude",
        )
    )
    field_figures.append(
        eigenfunction_history(
            model,
            coefficient,
            training,
            FIT_START,
            FIT_END,
            component="phase",
        )
    )
    render_viewer(
        Path(training_records[0]["source"]),
        field_figures,
        name=FIELD_PLOT,
        label="Viscous echo-Hankel eigenfunction field",
        columns=2,
        width=1900,
        height=3300,
    )

    x_edges = np.linspace(lower[0], upper[0], 82)
    y_edges = np.linspace(lower[1], upper[1], 82)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    defect_step, occupied_step = binned_defect(
        model, coefficient, multiplier, training, x_edges, y_edges, 0.01
    )
    defect_four, occupied_four = binned_defect(
        model, coefficient, multiplier, training, x_edges, y_edges, 4.0
    )
    full_traces = []
    for index, record in enumerate(training, start=1):
        full_mask = (record["time"] >= FULL_AUDIT_START) & (record["time"] <= FIT_END)
        late_mask = (record["time"] >= FIT_START) & (record["time"] <= FIT_END)
        full_traces.append(
            trace(
                "earlier transient" if index == 1 else "",
                record["coordinates"][full_mask, 0][::80],
                record["coordinates"][full_mask, 1][::80],
                color="#94a3b8",
                width=1.0,
            )
        )
        full_traces.append(
            trace(
                f"fit suffix #{index}",
                record["coordinates"][late_mask, 0][::30],
                record["coordinates"][late_mask, 1][::30],
                color="#26313d",
                line_style="dash",
                opacity=max(0.48, 0.88 - 0.10 * (index - 1)),
            )
        )
    support = figure(
        "phase",
        "Full observed latent region versus fitted suffix",
        "xi_1",
        "xi_2",
        full_traces,
    )
    support["square"] = True
    step_map = figure(
        "krom_heatmap",
        "Occupied-cell 0.01 tau Koopman defect",
        "xi_1",
        "xi_2",
        [heat_trace("log10 relative defect", x_centers, y_centers, defect_step), *overlays],
        color_label="log10 relative defect",
    )
    step_map["square"] = True
    step_map["value_scale"] = "linear"
    four_map = figure(
        "krom_heatmap",
        "Occupied-cell 4 tau Koopman defect",
        "xi_1",
        "xi_2",
        [heat_trace("log10 relative defect", x_centers, y_centers, defect_four), *overlays],
        color_label="log10 relative defect",
    )
    four_map["square"] = True
    four_map["value_scale"] = "linear"
    segment_edges = np.asarray((5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0))
    segment_midpoints = 0.5 * (segment_edges[:-1] + segment_edges[1:])
    error_history = figure(
        "line",
        "Fixed late-time eigenfunction outside its fit interval",
        "segment midpoint tau",
        "4 tau relative L2 error",
        [
            trace(
                "training trajectories",
                segment_midpoints,
                segment_error(model, coefficient, multiplier, training, segment_edges, 4.0),
                color="#0891b2",
                marker=True,
                width=2.0,
            ),
            trace(
                "held-out trajectories",
                segment_midpoints,
                segment_error(model, coefficient, multiplier, validation, segment_edges, 4.0),
                color="#dc2626",
                marker=True,
                width=2.0,
            ),
        ],
    )
    render_viewer(
        Path(training_records[0]["source"]),
        [support, step_map, four_map, error_history],
        name=ACCURACY_PLOT,
        label="Viscous echo-Hankel eigenfunction support and accuracy",
        columns=2,
        width=1900,
        height=1900,
    )

    output = {
        "schema": "dynsyslab.rpe-viscous-echo-eigenfunction-domain/v1",
        "status": "formal eigenfunction evaluation; not promoted to a KROM",
        "coordinate_contract": {
            "polynomial_variables": ["xi1", "xi2"],
            "history_dependent": True,
            "basis": str(PARENT_ARTIFACT),
            "online_input": "causal echo ring buffer",
            "core_history_coordinate_provider_available": False,
        },
        "fit": {
            "start_tau": FIT_START,
            "end_tau": FIT_END,
            "training_trajectories": len(training),
            "polynomial_degree": POLYNOMIAL_DEGREE,
            "feature_count": len(model.exponents),
            "candidate_count": len(candidates),
            "multiplier": {"real": float(multiplier.real), "imaginary": float(multiplier.imag)},
            "continuous_eigenvalue": {
                "real": float(function.eigenvalue.continuous_time.real),
                "imaginary": float(function.eigenvalue.continuous_time.imag),
            },
            "coordinate_domain": training_domain,
        },
        "canonical_eigenfunction": function.summary(),
        "canonical_training_evaluation": training_evaluation.summary(),
        "canonical_validation_evaluation": validation_evaluation.summary(),
        "occupied_grid_fraction": {
            "0.01_tau": float(np.mean(occupied_step)),
            "4_tau": float(np.mean(occupied_four)),
        },
    }
    METRICS.write_text(json.dumps(output, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        coefficient=coefficient,
        exponents=np.asarray(model.exponents, dtype=int),
        coordinate_center=model.center,
        coordinate_scale=model.scale,
        feature_scale=model.feature_scale,
        multiplier=np.asarray(multiplier),
        continuous_eigenvalue=np.asarray(function.eigenvalue.continuous_time),
        training_domain_minimum=lower,
        training_domain_maximum=upper,
        metadata_json=np.asarray(json.dumps(output, separators=(",", ":"))),
    )
    print(f"lambda={function.eigenvalue.continuous_time:.10g}")
    print("training", training_evaluation.summary())
    print("validation", validation_evaluation.summary())
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
