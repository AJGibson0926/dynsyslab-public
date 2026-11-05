"""Escalate echo-only Koopman tests from two to four latent coordinates."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve
from scipy.sparse.linalg import svds

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_echo_sampled_koopman import (
    DELAY_SAMPLES,
    candidate_table,
    eigenfunction_metrics,
    fit_action_invariant,
    fit_edmd,
    invariant_metrics,
    invariant_serializable,
    orient_coefficient,
)
from compare_rpe_echo_pooled_training_injectivity import (
    METRICS as HANKEL_METRICS,
    PooledHankelOperator,
    load_training_records,
    load_validations,
)
from project_rpe_inviscid_echo_basis_out_of_sample import (
    EMBEDDING_DECIMATION,
    TRAINING,
)


ROOT = Path(__file__).resolve().parents[1]
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-latent-rank-koopman.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-latent-rank-koopman.npz"
PLOT_NAME = "rpe-inviscid-echo-latent-rank-koopman"
RANKS = (2, 3, 4)


def fit_rank_four(signals: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Fit the leading pooled Hankel modes with deterministic singular-vector orientation."""
    operator = PooledHankelOperator(signals, DELAY_SAMPLES)
    left, singular, right = svds(
        operator,
        k=max(RANKS),
        which="LM",
        solver="arpack",
        tol=1.0e-10,
        maxiter=10000,
        return_singular_vectors=True,
        rng=np.random.default_rng(20261001),
    )
    order = np.argsort(singular)[::-1]
    left = left[:, order]
    singular = singular[order]
    right = right[order]
    for mode in range(left.shape[1]):
        pivot = int(np.argmax(np.abs(left[:, mode])))
        if left[pivot, mode] < 0.0:
            left[:, mode] *= -1.0
            right[mode] *= -1.0
    return left, singular


def project(
    echo: np.ndarray,
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
) -> np.ndarray:
    """Project normalized echo windows through fixed modes and singular-value scaling."""
    normalized = (echo[::EMBEDDING_DECIMATION] - mean) / scale
    return np.column_stack(
        [
            fftconvolve(normalized, modes[:, mode][::-1], mode="valid")
            / singular[mode]
            for mode in range(modes.shape[1])
        ]
    )


def prepare_runs(
    records: list[dict],
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
    prefix: str,
) -> list[dict]:
    """Project each echo history and attach aligned timestamps and diagnostic labels."""
    from analyze_rpe_echo_sampled_koopman import echo_period

    runs = []
    for index, record in enumerate(records, start=1):
        coordinates = project(record["echo"], modes, singular, mean, scale)
        runs.append(
            {
                "label": f"{prefix} #{index}",
                "initial_state": np.asarray(record["initial_state"], dtype=float),
                "coordinates": coordinates,
                "latent_time": record["time"][::EMBEDDING_DECIMATION][
                    DELAY_SAMPLES - 1 :
                ],
                "echo_period_tau": echo_period(record["time"], record["echo"]),
            }
        )
    return runs


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    with HANKEL_METRICS.open("r", encoding="utf-8") as stream:
        hankel_metrics = json.load(stream)
    mean = float(hankel_metrics["pooling"]["mean"])
    scale = float(hankel_metrics["pooling"]["scale"])
    training_records, _ = load_training_records()
    validation_records = load_validations()
    for record in validation_records:
        record["initial_state"] = np.asarray(record["state"][0], dtype=float)
    normalized_echoes = [
        (record["echo"][::EMBEDDING_DECIMATION] - mean) / scale
        for record in training_records
    ]
    modes, singular = fit_rank_four(normalized_echoes)

    rows = []
    saved: dict[str, np.ndarray] = {
        "delay_modes_rank4": modes,
        "delay_singular_rank4": singular,
    }
    for rank in RANKS:
        training_runs = prepare_runs(
            training_records, modes[:, :rank], singular[:rank], mean, scale, "training"
        )
        validation_runs = prepare_runs(
            validation_records,
            modes[:, :rank],
            singular[:rank],
            mean,
            scale,
            "validation",
        )
        model = fit_edmd([run["coordinates"] for run in training_runs])
        candidates = candidate_table(
            model, [run["coordinates"] for run in training_runs]
        )
        if not candidates:
            raise RuntimeError(f"rank {rank} has no admissible oscillatory candidate")
        multiplier = complex(candidates[0]["multiplier"])
        coefficient = orient_coefficient(
            model,
            np.asarray(candidates[0]["coefficient"], dtype=complex),
            training_runs[0]["coordinates"],
        )
        phase_validation = eigenfunction_metrics(
            model, coefficient, multiplier, validation_runs
        )
        invariant_coefficient, invariant_fit = fit_action_invariant(
            model, [run["coordinates"] for run in training_runs]
        )
        invariant_training = invariant_metrics(
            model, invariant_coefficient, training_runs
        )
        invariant_validation = invariant_metrics(
            model, invariant_coefficient, validation_runs
        )
        invariant = invariant_serializable(
            model,
            invariant_coefficient,
            invariant_fit,
            invariant_training,
            invariant_validation,
        )
        row = {
            "rank": rank,
            "feature_count": len(model.exponents),
            "oscillatory_multiplier": {
                "real": float(multiplier.real),
                "imaginary": float(multiplier.imag),
            },
            "oscillatory_implied_period_tau": float(
                2.0 * np.pi / (np.angle(multiplier) / model.pair_step_tau)
            ),
            "oscillatory_validation_median_one_step": float(
                np.median([item["one_step_relative_l2"] for item in phase_validation])
            ),
            "oscillatory_validation_median_four_tau": float(
                np.median(
                    [item["horizon_relative_l2"]["4"] for item in phase_validation]
                )
            ),
            "invariant_validation_median_std": float(
                np.median([item["standard_deviation"] for item in invariant_validation])
            ),
            "invariant_validation_maximum_std": float(
                np.max([item["standard_deviation"] for item in invariant_validation])
            ),
            "invariant_validation_median_step_residual": float(
                np.median(
                    [item["sampled_invariance_relative_l2"] for item in invariant_validation]
                )
            ),
            "frequency_conditioning_validation_relative_l2": invariant[
                "frequency_conditioning_diagnostic"
            ]["held_out_relative_l2"],
            "invariant": invariant,
        }
        rows.append(row)
        saved[f"rank_{rank}_oscillatory_coefficient"] = coefficient
        saved[f"rank_{rank}_invariant_coefficient"] = invariant_coefficient
        print(
            f"rank={rank}, features={row['feature_count']}, "
            f"phase step/four-tau={row['oscillatory_validation_median_one_step']:.4g}/"
            f"{row['oscillatory_validation_median_four_tau']:.4g}, "
            f"invariant median/max std={row['invariant_validation_median_std']:.4g}/"
            f"{row['invariant_validation_maximum_std']:.4g}, "
            f"omega(I)={row['frequency_conditioning_validation_relative_l2']:.4g}",
            flush=True,
        )

    output = {
        "schema": "dynsyslab.echo-latent-rank-koopman/v1",
        "training_policy": "unforced echo only",
        "controlled_comparison": {
            "ranks": list(RANKS),
            "same_training_runs": True,
            "same_validation_runs": True,
            "same_delay_samples": DELAY_SAMPLES,
            "same_normalization": True,
            "same_parent_svd": True,
            "same_polynomial_degree": True,
        },
        "rows": rows,
    }
    METRICS.write_text(json.dumps(output, indent=2), encoding="utf-8")
    saved["metadata_json"] = np.asarray(json.dumps(output, separators=(",", ":")))
    np.savez_compressed(ARTIFACT, **saved)

    x = np.asarray(RANKS)
    figures = [
        figure(
            "line",
            "Held-out oscillatory eigenfunction error",
            "retained latent rank",
            "relative L2 error",
            [
                trace(
                    "one sampled step",
                    x,
                    np.asarray([row["oscillatory_validation_median_one_step"] for row in rows]),
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "four tau",
                    x,
                    np.asarray([row["oscillatory_validation_median_four_tau"] for row in rows]),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "line",
            "Held-out invariant variation",
            "retained latent rank",
            "standard deviation of I",
            [
                trace(
                    "median",
                    x,
                    np.asarray([row["invariant_validation_median_std"] for row in rows]),
                    color="#0891b2",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "maximum",
                    x,
                    np.asarray([row["invariant_validation_maximum_std"] for row in rows]),
                    color="#f59e0b",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "line",
            "Invariant sampled-step residual",
            "retained latent rank",
            "median relative L2 residual",
            [
                trace(
                    "mu = 1",
                    x,
                    np.asarray(
                        [row["invariant_validation_median_step_residual"] for row in rows]
                    ),
                    color="#7c3aed",
                    marker=True,
                    width=2.0,
                )
            ],
        ),
        figure(
            "line",
            "Held-out frequency error from omega(I)",
            "retained latent rank",
            "relative L2 error",
            [
                trace(
                    "quadratic frequency law",
                    x,
                    np.asarray(
                        [row["frequency_conditioning_validation_relative_l2"] for row in rows]
                    ),
                    color="#059669",
                    marker=True,
                    width=2.0,
                )
            ],
        ),
    ]
    render_viewer(
        TRAINING,
        figures,
        name=PLOT_NAME,
        label="Echo-only Koopman latent-rank comparison",
        columns=2,
        width=1900,
        height=1900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
