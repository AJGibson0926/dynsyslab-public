"""Preserved numerical diagnostics for the historical sampled-EDMD experiment.

Research reference only; production Koopman fitting lives in the native core.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations_with_replacement
import numpy as np
from scipy.linalg import eig
from scipy.signal import find_peaks
from project_rpe_inviscid_echo_basis_out_of_sample import EMBEDDING_DECIMATION, EXPECTED_STEP

POLYNOMIAL_DEGREE = 5
PAIR_STRIDE = 20
RIDGE_RELATIVE = 1.0e-10
HORIZONS_TAU = (0.01, 0.1, 0.5, 1.0, 2.0, 4.0)


def polynomial_exponents(dimension: int, degree: int) -> list[tuple[int, ...]]:
    """Enumerate total-degree monomials including the constant in canonical order."""
    exponents = [tuple([0] * dimension)]
    for order in range(1, degree + 1):
        for indices in combinations_with_replacement(range(dimension), order):
            value = [0] * dimension
            for index in indices:
                value[index] += 1
            exponents.append(tuple(value))
    return exponents


def polynomial_values(
    coordinates: np.ndarray, exponents: list[tuple[int, ...]]
) -> np.ndarray:
    """Evaluate the supplied monomial exponent library on a batch of coordinates."""
    coordinates = np.asarray(coordinates, dtype=float)
    values = np.ones((coordinates.shape[0], len(exponents)), dtype=float)
    for column, exponent in enumerate(exponents[1:], start=1):
        for axis, power in enumerate(exponent):
            if power:
                values[:, column] *= coordinates[:, axis] ** power
    return values


def exponent_label(exponent: tuple[int, ...]) -> str:
    """Format a monomial exponent tuple as a readable normalized-coordinate term."""
    if not any(exponent):
        return "1"
    pieces = []
    for axis, power in enumerate(exponent, start=1):
        if power == 1:
            pieces.append(f"z{axis}")
        elif power:
            pieces.append(f"z{axis}^{power}")
    return "*".join(pieces)


@dataclass
class PolynomialEDMD:
    """Retain the research polynomial dictionary, normalization, operator, and spectrum."""
    center: np.ndarray
    scale: np.ndarray
    feature_scale: np.ndarray
    exponents: list[tuple[int, ...]]
    operator: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    pair_step_tau: float
    training_runs: int
    ridge: float

    def features(self, coordinates: np.ndarray) -> np.ndarray:
        normalized = (np.asarray(coordinates) - self.center) / self.scale
        return polynomial_values(normalized, self.exponents) / self.feature_scale

    def evaluate(self, coordinates: np.ndarray, coefficient: np.ndarray) -> np.ndarray:
        return self.features(coordinates) @ coefficient


def fit_edmd(coordinate_runs: list[np.ndarray]) -> PolynomialEDMD:
    """Fit a regularized polynomial sampled operator and its eigensystem from separate runs."""
    combined = np.vstack(coordinate_runs)
    center = np.mean(combined, axis=0)
    scale = np.std(combined, axis=0)
    if np.any(scale <= np.finfo(float).eps):
        raise RuntimeError("latent-coordinate normalization is singular")
    exponents = polynomial_exponents(combined.shape[1], POLYNOMIAL_DEGREE)
    normalized_runs = [(run - center) / scale for run in coordinate_runs]
    all_features = np.vstack(
        [polynomial_values(run, exponents) for run in normalized_runs]
    )
    feature_scale = np.sqrt(np.mean(all_features * all_features, axis=0))
    feature_scale[feature_scale <= np.finfo(float).eps] = 1.0
    feature_runs = [
        polynomial_values(run, exponents) / feature_scale for run in normalized_runs
    ]
    minus = np.vstack([values[:-PAIR_STRIDE] for values in feature_runs])
    plus = np.vstack([values[PAIR_STRIDE:] for values in feature_runs])
    gram = minus.T @ minus / minus.shape[0]
    cross = minus.T @ plus / minus.shape[0]
    ridge = RIDGE_RELATIVE * float(np.trace(gram) / gram.shape[0])
    operator = np.linalg.solve(gram + ridge * np.eye(gram.shape[0]), cross)
    eigenvalues, eigenvectors = eig(operator)
    return PolynomialEDMD(
        center=center,
        scale=scale,
        feature_scale=feature_scale,
        exponents=exponents,
        operator=operator,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        pair_step_tau=EXPECTED_STEP * EMBEDDING_DECIMATION * PAIR_STRIDE,
        training_runs=len(coordinate_runs),
        ridge=ridge,
    )


def echo_period(time: np.ndarray, echo: np.ndarray) -> float:
    """Estimate the scalar echo period from detected prominent peaks."""
    sample_step = float(np.median(np.diff(time)))
    peaks, _ = find_peaks(
        echo,
        distance=max(1, int(round(0.5 / sample_step))),
        prominence=0.15 * float(np.ptp(echo)),
    )
    if peaks.size < 3:
        raise RuntimeError("too few echo peaks for a period estimate")
    return float(np.median(np.diff(time[peaks])))


def phase_regression(time: np.ndarray, values: np.ndarray) -> dict:
    """Fit an unwrapped phase line and report its frequency and residual diagnostics."""
    phase = np.unwrap(np.angle(values))
    centered_time = time - np.mean(time)
    centered_phase = phase - np.mean(phase)
    denominator = float(centered_time @ centered_time)
    slope = float(centered_time @ centered_phase / denominator)
    intercept = float(np.mean(phase) - slope * np.mean(time))
    fitted = intercept + slope * time
    residual = phase - fitted
    total = phase - np.mean(phase)
    r_squared = 1.0 - float(residual @ residual) / max(
        float(total @ total), np.finfo(float).tiny
    )
    return {
        "frequency_rad_per_tau": slope,
        "period_tau": float(2.0 * np.pi / abs(slope)) if slope else None,
        "r_squared": r_squared,
        "phase": phase,
        "fitted_phase": fitted,
    }


def eigenfunction_metrics(
    model: PolynomialEDMD,
    coefficient: np.ndarray,
    multiplier: complex,
    runs: list[dict],
) -> list[dict]:
    """Measure sampled, finite-horizon, amplitude, and phase defects along separate runs."""
    continuous_value = np.log(multiplier) / model.pair_step_tau
    results = []
    for run in runs:
        values = model.evaluate(run["coordinates"], coefficient)
        time = run["latent_time"]
        phase_fit = phase_regression(time, values)
        pair_residual = float(
            np.linalg.norm(values[PAIR_STRIDE:] - multiplier * values[:-PAIR_STRIDE])
            / max(np.linalg.norm(values[PAIR_STRIDE:]), np.finfo(float).tiny)
        )
        horizon_errors = {}
        for horizon in HORIZONS_TAU:
            shift = int(round(horizon / (EXPECTED_STEP * EMBEDDING_DECIMATION)))
            if shift >= values.size:
                continue
            propagated = np.exp(continuous_value * horizon) * values[:-shift]
            horizon_errors[f"{horizon:g}"] = float(
                np.linalg.norm(values[shift:] - propagated)
                / max(np.linalg.norm(values[shift:]), np.finfo(float).tiny)
            )
        magnitude = np.abs(values)
        results.append(
            {
                "label": run["label"],
                "initial_state": run["initial_state"].tolist(),
                "echo_period_tau": run["echo_period_tau"],
                "one_step_relative_l2": pair_residual,
                "horizon_relative_l2": horizon_errors,
                "observed_eigenfunction_frequency_rad_per_tau": phase_fit[
                    "frequency_rad_per_tau"
                ],
                "observed_eigenfunction_period_tau": phase_fit["period_tau"],
                "phase_linearity_r_squared": phase_fit["r_squared"],
                "magnitude_coefficient_of_variation": float(
                    np.std(magnitude) / max(np.mean(magnitude), np.finfo(float).tiny)
                ),
                "values": values,
                "phase": phase_fit["phase"],
                "fitted_phase": phase_fit["fitted_phase"],
                "time": time,
            }
        )
    return results


def candidate_table(
    model: PolynomialEDMD, coordinate_runs: list[np.ndarray]
) -> list[dict]:
    """Rank admissible oscillatory modes using training residual and phase diagnostics."""
    rows = []
    for index, (value, vector) in enumerate(
        zip(model.eigenvalues, model.eigenvectors.T, strict=True)
    ):
        if value.imag <= 0.0 or not (0.9 <= abs(value) <= 1.02):
            continue
        frequency = float(np.angle(value) / model.pair_step_tau)
        if not 2.5 <= frequency <= 6.5:
            continue
        norm = np.sqrt(
            np.mean(
                np.abs(
                    np.concatenate(
                        [model.evaluate(run, vector) for run in coordinate_runs]
                    )
                )
                ** 2
            )
        )
        if norm <= np.finfo(float).tiny:
            continue
        vector = vector / norm
        residuals = []
        phase_r2 = []
        phase_frequencies = []
        for run in coordinate_runs:
            values = model.evaluate(run, vector)
            residuals.append(
                float(
                    np.linalg.norm(
                        values[PAIR_STRIDE:] - value * values[:-PAIR_STRIDE]
                    )
                    / max(
                        np.linalg.norm(values[PAIR_STRIDE:]), np.finfo(float).tiny
                    )
                )
            )
            time = np.arange(values.size) * EXPECTED_STEP * EMBEDDING_DECIMATION
            fit = phase_regression(time, values)
            phase_r2.append(fit["r_squared"])
            phase_frequencies.append(fit["frequency_rad_per_tau"])
        rows.append(
            {
                "index": index,
                "multiplier": value,
                "coefficient": vector,
                "frequency_rad_per_tau": frequency,
                "period_tau": float(2.0 * np.pi / frequency),
                "mean_one_step_residual": float(np.mean(residuals)),
                "maximum_one_step_residual": float(np.max(residuals)),
                "minimum_phase_r_squared": float(np.min(phase_r2)),
                "observed_frequency_spread": float(np.ptp(phase_frequencies)),
            }
        )
    rows.sort(
        key=lambda row: (
            row["mean_one_step_residual"],
            -row["minimum_phase_r_squared"],
        )
    )
    return rows


def orient_coefficient(
    model: PolynomialEDMD, coefficient: np.ndarray, reference: np.ndarray
) -> np.ndarray:
    """Fix the eigenfunction phase using the selected training reference."""
    values = model.evaluate(reference, coefficient)
    pivot = values[np.argmax(np.abs(values))]
    if abs(pivot) > np.finfo(float).tiny:
        coefficient = coefficient * np.exp(-1j * np.angle(pivot))
    return coefficient


def coefficient_summary(
    model: PolynomialEDMD, coefficient: np.ndarray, relative_floor: float = 1.0e-3
) -> list[dict]:
    """Rank polynomial terms by coefficient magnitude for inspection."""
    maximum = float(np.max(np.abs(coefficient)))
    rows = []
    for exponent, value in zip(model.exponents, coefficient, strict=True):
        if abs(value) < relative_floor * maximum:
            continue
        rows.append(
            {
                "term": exponent_label(exponent),
                "real": float(value.real),
                "imaginary": float(value.imag),
                "magnitude": float(abs(value)),
            }
        )
    rows.sort(key=lambda row: row["magnitude"], reverse=True)
    return rows


def fit_action_invariant(
    model: PolynomialEDMD, coordinate_runs: list[np.ndarray]
) -> tuple[np.ndarray, dict]:
    """Fit a nonconstant mu=1 eigenfunction using run identity only.

    This is a Fisher-style generalized eigenproblem: maximize variation of
    the polynomial's mean between unforced trajectories while minimizing its
    variation within each trajectory.  The constant feature is excluded from
    the optimization and restored only to center the final invariant.
    """
    features = [model.features(run) for run in coordinate_runs]
    reduced = [values[:, 1:] for values in features]
    means = [np.mean(values, axis=0) for values in reduced]
    total_samples = int(sum(values.shape[0] for values in reduced))
    global_mean = sum(
        values.shape[0] * mean for values, mean in zip(reduced, means, strict=True)
    ) / total_samples
    within = np.zeros((reduced[0].shape[1], reduced[0].shape[1]))
    between = np.zeros_like(within)
    for values, mean in zip(reduced, means, strict=True):
        centered = values - mean
        within += centered.T @ centered
        offset = mean - global_mean
        between += values.shape[0] * np.outer(offset, offset)
    within /= total_samples
    between /= total_samples
    ridge = 1.0e-8 * float(np.trace(within) / within.shape[0])
    values, vectors = eig(between, within + ridge * np.eye(within.shape[0]))
    admissible = np.flatnonzero(
        np.isfinite(values.real)
        & (np.abs(values.imag) < 1.0e-8)
        & (values.real > 0.0)
    )
    if admissible.size == 0:
        raise RuntimeError("no nonconstant invariant candidate was found")
    index = int(admissible[np.argmax(values.real[admissible])])
    reduced_coefficient = np.asarray(vectors[:, index].real, dtype=float)
    pivot = int(np.argmax(np.abs(reduced_coefficient)))
    if reduced_coefficient[pivot] < 0.0:
        reduced_coefficient *= -1.0
    coefficient = np.zeros(len(model.exponents), dtype=float)
    coefficient[1:] = reduced_coefficient
    coefficient[0] = -float(global_mean @ reduced_coefficient)
    run_means = np.asarray(
        [np.mean(values @ coefficient) for values in features], dtype=float
    )
    mean_scale = float(np.std(run_means))
    if mean_scale <= np.finfo(float).eps:
        raise RuntimeError("action invariant failed to distinguish training runs")
    coefficient /= mean_scale
    return coefficient, {
        "generalized_eigenvalue": float(values[index].real),
        "within_regularization": ridge,
        "constraint": "maximize between-run variance / within-run variance",
        "constant_feature_in_optimization": False,
        "hidden_state_used": False,
    }


def invariant_metrics(
    model: PolynomialEDMD,
    coefficient: np.ndarray,
    runs: list[dict],
) -> list[dict]:
    """Measure within-orbit variation and sampled defects of a candidate invariant."""
    results = []
    for run in runs:
        values = np.asarray(
            model.evaluate(run["coordinates"], coefficient), dtype=float
        )
        difference = values[PAIR_STRIDE:] - values[:-PAIR_STRIDE]
        results.append(
            {
                "label": run["label"],
                "initial_state": run["initial_state"].tolist(),
                "echo_period_tau": run["echo_period_tau"],
                "echo_frequency_rad_per_tau": float(
                    2.0 * np.pi / run["echo_period_tau"]
                ),
                "mean": float(np.mean(values)),
                "standard_deviation": float(np.std(values)),
                "peak_to_peak": float(np.ptp(values)),
                "sampled_invariance_relative_l2": float(
                    np.linalg.norm(difference)
                    / max(np.linalg.norm(values[PAIR_STRIDE:]), np.finfo(float).tiny)
                ),
                "time": run["latent_time"],
                "values": values,
            }
        )
    return results


def invariant_serializable(
    model: PolynomialEDMD,
    coefficient: np.ndarray,
    fit_diagnostics: dict,
    training_metrics: list[dict],
    validation_metrics: list[dict],
) -> dict:
    """Serialize invariant accuracy and the diagnostic frequency-versus-invariant fit."""
    training_action = np.asarray([row["mean"] for row in training_metrics])
    training_frequency = np.asarray(
        [row["echo_frequency_rad_per_tau"] for row in training_metrics]
    )
    frequency_coefficients = np.polyfit(training_action, training_frequency, 2)
    validation_action = np.asarray([row["mean"] for row in validation_metrics])
    validation_frequency = np.asarray(
        [row["echo_frequency_rad_per_tau"] for row in validation_metrics]
    )
    validation_prediction = np.polyval(frequency_coefficients, validation_action)
    frequency_relative_l2 = float(
        np.linalg.norm(validation_frequency - validation_prediction)
        / np.linalg.norm(validation_frequency)
    )
    return {
        "sampled_multiplier": 1.0,
        "continuous_eigenvalue": 0.0,
        "interpretation": "data-driven action-like invariant",
        "fit_diagnostics": fit_diagnostics,
        "coefficient_real": coefficient.tolist(),
        "coefficient_terms_above_relative_1e_3": coefficient_summary(
            model, coefficient.astype(complex)
        ),
        "training": [
            {key: value for key, value in row.items() if key not in {"time", "values"}}
            for row in training_metrics
        ],
        "validation": [
            {key: value for key, value in row.items() if key not in {"time", "values"}}
            for row in validation_metrics
        ],
        "frequency_conditioning_diagnostic": {
            "model": "omega(I) quadratic least squares",
            "coefficients_descending_power": frequency_coefficients.tolist(),
            "held_out_relative_l2": frequency_relative_l2,
            "note": "omega(I) is a state-dependent phase-rate model, not a fixed-eigenvalue Koopman eigenfunction",
        },
    }


def serializable_eigenfunction(
    model: PolynomialEDMD,
    multiplier: complex,
    coefficient: np.ndarray,
    candidates: list[dict],
    metrics: list[dict],
) -> dict:
    """Serialize fitted normalization, coefficients, eigenvalues, and error summaries."""
    continuous = np.log(multiplier) / model.pair_step_tau
    return {
        "sampled_multiplier": {
            "real": float(multiplier.real),
            "imaginary": float(multiplier.imag),
            "magnitude": float(abs(multiplier)),
            "angle_rad": float(np.angle(multiplier)),
        },
        "continuous_log_eigenvalue": {
            "real": float(continuous.real),
            "imaginary": float(continuous.imag),
        },
        "implied_period_tau": float(2.0 * np.pi / continuous.imag),
        "pair_step_tau": model.pair_step_tau,
        "polynomial_degree": POLYNOMIAL_DEGREE,
        "feature_count": len(model.exponents),
        "ridge": model.ridge,
        "coordinate_center": model.center.tolist(),
        "coordinate_scale": model.scale.tolist(),
        "feature_scale": model.feature_scale.tolist(),
        "coefficient_real": coefficient.real.tolist(),
        "coefficient_imaginary": coefficient.imag.tolist(),
        "coefficient_terms_above_relative_1e_3": coefficient_summary(
            model, coefficient
        ),
        "candidate_count": len(candidates),
        "leading_candidates": [
            {
                key: (
                    {
                        "real": float(value.real),
                        "imaginary": float(value.imag),
                    }
                    if key == "multiplier"
                    else value
                )
                for key, value in row.items()
                if key not in {"coefficient"}
            }
            for row in candidates[:8]
        ],
        "per_trajectory": [
            {key: value for key, value in row.items() if key not in {"values", "phase", "fitted_phase", "time"}}
            for row in metrics
        ],
    }
