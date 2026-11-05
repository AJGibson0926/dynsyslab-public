"""Viewer composition for the historical echo sampled-Koopman experiment."""
from __future__ import annotations
import numpy as np
from analyze_rpe_echo_hankel_latents import figure, render_viewer, trace
from project_rpe_inviscid_echo_basis_out_of_sample import COLORS, TRAINING

SUMMARY_PLOT = "rpe-inviscid-echo-sampled-koopman-summary"
EVOLUTION_PLOT = "rpe-inviscid-echo-sampled-koopman-evolution"
INVARIANT_PLOT = "rpe-inviscid-echo-sampled-koopman-invariant"


def render_experiment(
    training_runs, validation_runs, pooled_training_metrics, pooled_validation_metrics,
    single_validation_metrics, pooled_multiplier, pooled_model, invariant_training_metrics,
    invariant_validation_metrics, invariant_result,
):
    """Render the historical period, eigenfunction-evolution, and invariant comparisons."""
    all_period_runs = training_runs + validation_runs
    period_x = np.arange(1, len(all_period_runs) + 1)
    pooled_train_residual = np.asarray(
        [row["one_step_relative_l2"] for row in pooled_training_metrics]
    )
    pooled_test_residual = np.asarray(
        [row["one_step_relative_l2"] for row in pooled_validation_metrics]
    )
    pooled_test_long = np.asarray(
        [row["horizon_relative_l2"]["4"] for row in pooled_validation_metrics]
    )
    single_test_long = np.asarray(
        [row["horizon_relative_l2"]["4"] for row in single_validation_metrics]
    )
    summary_figures = [
        figure(
            "line",
            "Amplitude-dependent periods in the unforced data",
            "trajectory index",
            "echo period (tau)",
            [
                trace(
                    "training",
                    period_x[: len(training_runs)],
                    np.asarray([run["echo_period_tau"] for run in training_runs]),
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "held out",
                    period_x[len(training_runs) :],
                    np.asarray([run["echo_period_tau"] for run in validation_runs]),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "line",
            "Pooled eigenfunction one-step residual",
            "trajectory index",
            "relative L2 residual",
            [
                trace(
                    "training",
                    np.arange(1, 6),
                    pooled_train_residual,
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "held out",
                    np.arange(1, 6),
                    pooled_test_residual,
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "line",
            "Four-tau semigroup error on held-out trajectories",
            "held-out trajectory index",
            "relative L2 error",
            [
                trace(
                    "fit on one orbit",
                    np.arange(1, 6),
                    single_test_long,
                    color="#6b7280",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "fit on five orbits",
                    np.arange(1, 6),
                    pooled_test_long,
                    color="#7c3aed",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "phase",
            "Pooled eigenfunction image on held-out trajectories",
            "Re psi",
            "Im psi",
            [
                trace(
                    f"IC #{index}",
                    row["values"][::8].real,
                    row["values"][::8].imag,
                    color=color,
                    width=1.8,
                )
                for index, (row, color) in enumerate(
                    zip(pooled_validation_metrics, COLORS, strict=True), start=1
                )
            ],
        ),
    ]
    summary_figures[-1]["square"] = True
    render_viewer(
        TRAINING,
        summary_figures,
        name=SUMMARY_PLOT,
        label="Echo-only sampled Koopman test",
        columns=2,
        width=1900,
        height=1900,
    )

    evolution_figures = []
    for index, (row, color) in enumerate(
        zip(pooled_validation_metrics[:4], COLORS[:4], strict=True), start=1
    ):
        sample = np.arange(0, row["time"].size, 8)
        continuous = np.log(pooled_multiplier) / pooled_model.pair_step_tau
        predicted = (
            row["values"][0]
            * np.exp(continuous * (row["time"] - row["time"][0]))
        )
        evolution_figures.append(
            figure(
                "line",
                f"Held-out IC #{index}: eigenfunction evolution",
                "tau",
                "Re psi",
                [
                    trace(
                        "observed",
                        row["time"][sample],
                        row["values"][sample].real,
                        color=color,
                        width=1.8,
                    ),
                    trace(
                        "global exponential prediction",
                        row["time"][sample],
                        predicted[sample].real,
                        color="#111827",
                        width=1.5,
                    ),
                ],
            )
        )
    render_viewer(
        TRAINING,
        evolution_figures,
        name=EVOLUTION_PLOT,
        label="Echo-only eigenfunction evolution",
        columns=2,
        width=1900,
        height=1900,
    )

    training_action = np.asarray(
        [row["mean"] for row in invariant_training_metrics]
    )
    validation_action = np.asarray(
        [row["mean"] for row in invariant_validation_metrics]
    )
    training_frequency = np.asarray(
        [row["echo_frequency_rad_per_tau"] for row in invariant_training_metrics]
    )
    validation_frequency = np.asarray(
        [row["echo_frequency_rad_per_tau"] for row in invariant_validation_metrics]
    )
    frequency_coefficients = np.asarray(
        invariant_result["frequency_conditioning_diagnostic"][
            "coefficients_descending_power"
        ]
    )
    action_grid = np.linspace(
        min(float(np.min(training_action)), float(np.min(validation_action))),
        max(float(np.max(training_action)), float(np.max(validation_action))),
        300,
    )
    invariant_figures = [
        figure(
            "line",
            "Action-like invariant on training trajectories",
            "tau",
            "I(xi_1, xi_2)",
            [
                trace(
                    f"training #{index}",
                    row["time"][::8],
                    row["values"][::8],
                    color=color,
                    width=1.8,
                )
                for index, (row, color) in enumerate(
                    zip(invariant_training_metrics, COLORS, strict=True), start=1
                )
            ],
        ),
        figure(
            "line",
            "Action-like invariant on held-out trajectories",
            "tau",
            "I(xi_1, xi_2)",
            [
                trace(
                    f"IC #{index}",
                    row["time"][::8],
                    row["values"][::8],
                    color=color,
                    width=1.8,
                )
                for index, (row, color) in enumerate(
                    zip(invariant_validation_metrics, COLORS, strict=True), start=1
                )
            ],
        ),
        figure(
            "line",
            "Frequency conditioned on the invariant",
            "mean I",
            "measured echo frequency (rad/tau)",
            [
                trace(
                    "quadratic fit on training",
                    action_grid,
                    np.polyval(frequency_coefficients, action_grid),
                    color="#111827",
                    width=2.0,
                ),
                trace(
                    "training",
                    training_action,
                    training_frequency,
                    color="#2563eb",
                    marker=True,
                    width=0.0,
                    marker_size=6.0,
                ),
                trace(
                    "held out",
                    validation_action,
                    validation_frequency,
                    color="#dc2626",
                    marker=True,
                    width=0.0,
                    marker_size=6.0,
                ),
            ],
        ),
        figure(
            "line",
            "Invariant within-orbit variation",
            "trajectory index",
            "standard deviation of I",
            [
                trace(
                    "training",
                    np.arange(1, 6),
                    np.asarray(
                        [row["standard_deviation"] for row in invariant_training_metrics]
                    ),
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "held out",
                    np.arange(1, 6),
                    np.asarray(
                        [row["standard_deviation"] for row in invariant_validation_metrics]
                    ),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
    ]
    render_viewer(
        TRAINING,
        invariant_figures,
        name=INVARIANT_PLOT,
        label="Echo-only Koopman invariant test",
        columns=2,
        width=1900,
        height=1900,
    )
