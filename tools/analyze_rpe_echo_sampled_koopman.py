"""Test sampled Koopman eigenfunctions on echo-history latent coordinates.

This is an isolated research experiment.  It uses the retained three-period,
five-trajectory Hankel basis and only its online-accessible top-two coordinates
to fit polynomial EDMD eigenfunctions.  Hidden RPE state is never used in the
fit.  A single-orbit fit is compared with a pooled-orbit fit because an
inviscid nonlinear oscillator need not possess one global oscillatory Koopman
multiplier across energy levels with amplitude-dependent frequency.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.signal import find_peaks

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from compare_rpe_echo_pooled_training_injectivity import (
    ARTIFACT as HANKEL_ARTIFACT,
    METRICS as HANKEL_METRICS,
    load_training_records,
    load_validations,
)
from project_rpe_inviscid_echo_basis_out_of_sample import (
    COLORS,
    EMBEDDING_DECIMATION,
    EXPECTED_STEP,
    TRAINING,
)
from sweep_rpe_echo_delay_injectivity import fixed_projection


ROOT = Path(__file__).resolve().parents[1]
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-sampled-koopman.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-sampled-koopman.npz"
SUMMARY_PLOT = "rpe-inviscid-echo-sampled-koopman-summary"
EVOLUTION_PLOT = "rpe-inviscid-echo-sampled-koopman-evolution"
INVARIANT_PLOT = "rpe-inviscid-echo-sampled-koopman-invariant"

DELAY_SAMPLES = 8644
POLYNOMIAL_DEGREE = 5
PAIR_STRIDE = 20
RIDGE_RELATIVE = 1.0e-10
HORIZONS_TAU = (0.01, 0.1, 0.5, 1.0, 2.0, 4.0)


from echo_sampled_koopman_plots import render_experiment
from echo_sampled_koopman_diagnostics import (
    polynomial_exponents,
    polynomial_values,
    exponent_label,
    PolynomialEDMD,
    fit_edmd,
    echo_period,
    phase_regression,
    eigenfunction_metrics,
    candidate_table,
    orient_coefficient,
    coefficient_summary,
    fit_action_invariant,
    invariant_metrics,
    invariant_serializable,
    serializable_eigenfunction,
)


def build_runs(
    records: list[dict],
    modes: np.ndarray,
    singular: np.ndarray,
    echo_mean: float,
    echo_scale: float,
    prefix: str,
) -> list[dict]:
    """Apply the retained echo basis and align latent windows to their final sample."""
    runs = []
    for index, record in enumerate(records, start=1):
        coordinates = fixed_projection(
            record["echo"], modes, singular, echo_mean, echo_scale
        )
        sampled_time = record["time"][::EMBEDDING_DECIMATION]
        latent_time = sampled_time[DELAY_SAMPLES - 1 :]
        runs.append(
            {
                "label": f"{prefix} #{index}",
                "initial_state": np.asarray(record["initial_state"], dtype=float),
                "coordinates": coordinates,
                "latent_time": latent_time,
                "echo_period_tau": echo_period(record["time"], record["echo"]),
            }
        )
    return runs


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
    validation_records = load_validations()
    for record in validation_records:
        record["initial_state"] = np.asarray(record["state"][0], dtype=float)
    training_runs = build_runs(
        training_records, modes, singular, echo_mean, echo_scale, "training"
    )
    validation_runs = build_runs(
        validation_records, modes, singular, echo_mean, echo_scale, "validation"
    )

    single_model = fit_edmd([training_runs[0]["coordinates"]])
    pooled_model = fit_edmd([run["coordinates"] for run in training_runs])
    single_candidates = candidate_table(
        single_model, [training_runs[0]["coordinates"]]
    )
    pooled_candidates = candidate_table(
        pooled_model, [run["coordinates"] for run in training_runs]
    )
    if not single_candidates or not pooled_candidates:
        raise RuntimeError("no admissible fundamental oscillatory EDMD mode was found")

    single_multiplier = complex(single_candidates[0]["multiplier"])
    pooled_multiplier = complex(pooled_candidates[0]["multiplier"])
    single_coefficient = orient_coefficient(
        single_model,
        np.asarray(single_candidates[0]["coefficient"], dtype=complex),
        training_runs[0]["coordinates"],
    )
    pooled_coefficient = orient_coefficient(
        pooled_model,
        np.asarray(pooled_candidates[0]["coefficient"], dtype=complex),
        training_runs[0]["coordinates"],
    )

    single_training_metrics = eigenfunction_metrics(
        single_model,
        single_coefficient,
        single_multiplier,
        [training_runs[0]],
    )
    single_validation_metrics = eigenfunction_metrics(
        single_model,
        single_coefficient,
        single_multiplier,
        validation_runs,
    )
    pooled_training_metrics = eigenfunction_metrics(
        pooled_model,
        pooled_coefficient,
        pooled_multiplier,
        training_runs,
    )
    pooled_validation_metrics = eigenfunction_metrics(
        pooled_model,
        pooled_coefficient,
        pooled_multiplier,
        validation_runs,
    )
    invariant_coefficient, invariant_fit = fit_action_invariant(
        pooled_model, [run["coordinates"] for run in training_runs]
    )
    invariant_training_metrics = invariant_metrics(
        pooled_model, invariant_coefficient, training_runs
    )
    invariant_validation_metrics = invariant_metrics(
        pooled_model, invariant_coefficient, validation_runs
    )
    invariant_result = invariant_serializable(
        pooled_model,
        invariant_coefficient,
        invariant_fit,
        invariant_training_metrics,
        invariant_validation_metrics,
    )

    output = {
        "schema": "dynsyslab.echo-sampled-koopman-experiment/v1",
        "status": "research experiment; not promoted to a KROM",
        "training_policy": "unforced echo only",
        "representation": {
            "hankel_training_sources": 5,
            "delay_samples": DELAY_SAMPLES,
            "delay_window_tau": (DELAY_SAMPLES - 1)
            * EXPECTED_STEP
            * EMBEDDING_DECIMATION,
            "latent_dimension": 2,
            "online_inputs": ["echo history"],
            "hidden_state_used_for_fit": False,
        },
        "fit": {
            "method": "polynomial sampled EDMD",
            "pair_stride": PAIR_STRIDE,
            "pair_step_tau": pooled_model.pair_step_tau,
            "polynomial_degree": POLYNOMIAL_DEGREE,
            "feature_count": len(pooled_model.exponents),
        },
        "single_orbit": {
            "fit_training_source_count": 1,
            "eigenfunction": serializable_eigenfunction(
                single_model,
                single_multiplier,
                single_coefficient,
                single_candidates,
                single_training_metrics + single_validation_metrics,
            ),
        },
        "pooled_orbit": {
            "fit_training_source_count": len(training_runs),
            "eigenfunction": serializable_eigenfunction(
                pooled_model,
                pooled_multiplier,
                pooled_coefficient,
                pooled_candidates,
                pooled_training_metrics + pooled_validation_metrics,
            ),
        },
        "pooled_invariant": invariant_result,
    }
    METRICS.write_text(json.dumps(output, indent=2), encoding="utf-8")

    saved: dict[str, np.ndarray] = {
        "delay_modes": modes,
        "delay_singular": singular,
        "single_multiplier": np.asarray(single_multiplier),
        "single_coefficient": single_coefficient,
        "single_coordinate_center": single_model.center,
        "single_coordinate_scale": single_model.scale,
        "single_feature_scale": single_model.feature_scale,
        "pooled_multiplier": np.asarray(pooled_multiplier),
        "pooled_coefficient": pooled_coefficient,
        "invariant_multiplier": np.asarray(1.0),
        "invariant_coefficient": invariant_coefficient,
        "polynomial_exponents": np.asarray(pooled_model.exponents, dtype=int),
        "pooled_coordinate_center": pooled_model.center,
        "pooled_coordinate_scale": pooled_model.scale,
        "pooled_feature_scale": pooled_model.feature_scale,
        "metadata_json": np.asarray(json.dumps(output, separators=(",", ":"))),
    }
    for prefix, metrics in (
        ("single_validation", single_validation_metrics),
        ("pooled_training", pooled_training_metrics),
        ("pooled_validation", pooled_validation_metrics),
    ):
        for index, row in enumerate(metrics, start=1):
            saved[f"{prefix}_{index}_time"] = row["time"]
            saved[f"{prefix}_{index}_values"] = row["values"]
    for prefix, metrics in (
        ("invariant_training", invariant_training_metrics),
        ("invariant_validation", invariant_validation_metrics),
    ):
        for index, row in enumerate(metrics, start=1):
            saved[f"{prefix}_{index}_time"] = row["time"]
            saved[f"{prefix}_{index}_values"] = row["values"]
    np.savez_compressed(ARTIFACT, **saved)

    render_experiment(
        training_runs, validation_runs, pooled_training_metrics, pooled_validation_metrics, single_validation_metrics, pooled_multiplier, pooled_model, invariant_training_metrics, invariant_validation_metrics, invariant_result
    )

    print("single-orbit candidates")
    for row in single_candidates[:6]:
        print(
            f"  mu={row['multiplier']:.8f}, period={row['period_tau']:.6f}, "
            f"one-step={row['mean_one_step_residual']:.4g}, "
            f"phase-R2={row['minimum_phase_r_squared']:.6f}"
        )
    print("pooled candidates")
    for row in pooled_candidates[:6]:
        print(
            f"  mu={row['multiplier']:.8f}, period={row['period_tau']:.6f}, "
            f"one-step={row['mean_one_step_residual']:.4g}, "
            f"phase-R2={row['minimum_phase_r_squared']:.6f}, "
            f"frequency-spread={row['observed_frequency_spread']:.4g}"
        )
    print("action-like invariant")
    print(
        "  training within-orbit std: "
        + ", ".join(
            f"{row['standard_deviation']:.4g}" for row in invariant_training_metrics
        )
    )
    print(
        "  held-out within-orbit std: "
        + ", ".join(
            f"{row['standard_deviation']:.4g}" for row in invariant_validation_metrics
        )
    )
    print(
        "  held-out omega(I) relative L2: "
        f"{invariant_result['frequency_conditioning_diagnostic']['held_out_relative_l2']:.4g}"
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
