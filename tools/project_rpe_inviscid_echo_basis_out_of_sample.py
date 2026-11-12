"""Project five unseen inviscid RPE echoes through one fixed Hankel basis.

This is deliberately not a decomposition or refit.  The training mean, scale,
left singular vectors U, and singular values are loaded from the retained
high-resolution training artifact.  Each new echo is sampled causally into a
ring buffer and evaluated online as xi_i(k) = u_i.T d_k / sigma_i.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

from analyze_rpe_echo_hankel_latents import (
    OUTPUT_ROOT,
    figure,
    reference,
    render_viewer,
    trace,
)
from analyze_rpe_inviscid_echo_hankel import (
    far_field_echo,
    load_native_source,
)


ROOT = Path(__file__).resolve().parents[1]
TRAINING = OUTPUT_ROOT / "rpe-inviscid-echo-hankel-hires-latents.npz"
SOURCE_GLOB = (
    "rpe-inviscid-echo-out-of-sample-five-ics-trajectory-*"
    "__20261001T220042Z.npz"
)
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-fixed-basis-five-ics.npz"
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-fixed-basis-five-ics.metrics.json"
PLOT_NAME = "rpe-inviscid-echo-fixed-basis-five-ics"

EMBEDDING_DECIMATION = 5
EXPECTED_STEP = 1.0e-4
COLORS = ["#0891b2", "#2563eb", "#7c3aed", "#f59e0b", "#dc2626"]
EXPECTED_ICS = np.asarray(
    [[1.1, 0.3], [1.5, -2.0], [2.0, 0.0], [2.45, 0.0], [2.6, 0.0]]
)


def causal_fixed_basis_projection(
    echo: np.ndarray,
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
) -> tuple[np.ndarray, float]:
    """Evaluate two fixed latent coordinates from a live echo ring buffer."""
    sampled = (echo[::EMBEDDING_DECIMATION] - mean) / scale
    delay = int(modes.shape[0])
    if sampled.size < delay:
        raise ValueError("echo history is shorter than the retained delay basis")

    ring = np.empty(2 * delay, dtype=float)
    coordinates = np.empty((sampled.size - delay + 1, modes.shape[1]))
    captured_energy = 0.0
    total_energy = 0.0
    output = 0
    for sample_index, value in enumerate(sampled):
        position = sample_index % delay
        ring[position] = value
        ring[position + delay] = value
        if sample_index < delay - 1:
            continue
        history_start = (sample_index + 1) % delay
        history = ring[history_start : history_start + delay]
        coefficients = history @ modes
        coordinates[output] = coefficients / singular
        captured_energy += float(coefficients @ coefficients)
        total_energy += float(history @ history)
        output += 1
    return coordinates, captured_energy / total_energy


def phase_tile(
    title: str,
    x_label: str,
    y_label: str,
    traces: list[dict],
    *,
    references: list[dict] | None = None,
) -> dict:
    """Compose a square Viewer phase-plane tile from the supplied curves."""
    tile = figure(
        "phase",
        title,
        x_label,
        y_label,
        traces,
        references=references,
    )
    tile["square"] = True
    return tile


def validate_validation_source(metadata: dict) -> None:
    """Require the out-of-sample data to be unforced, uncontrolled, and inviscid."""
    if metadata.get("forcing_type") != "free":
        raise RuntimeError("out-of-sample trajectory is not unforced")
    if metadata.get("controller_kind") != "none":
        raise RuntimeError("out-of-sample trajectory contains control")
    if float(metadata["resolved_system_parameters"]["mu_pa_s"]) != 0.0:
        raise RuntimeError("out-of-sample trajectory is not inviscid")


def closure_metric(
    echo: np.ndarray,
    coordinates: np.ndarray,
    delay: int,
    sample_step: float,
) -> tuple[float | None, int]:
    """Measure latent orbit closure at echo peaks relative to the orbit scale."""
    sampled_echo = echo[::EMBEDDING_DECIMATION]
    peaks, _ = find_peaks(
        sampled_echo,
        distance=max(1, int(round(0.5 / sample_step))),
        prominence=0.15 * float(np.ptp(sampled_echo)),
    )
    coordinate_indices = peaks - (delay - 1)
    coordinate_indices = coordinate_indices[coordinate_indices >= 0]
    if coordinate_indices.size < 3:
        return None, int(coordinate_indices.size)
    peak_coordinates = coordinates[coordinate_indices]
    closure = np.linalg.norm(np.diff(peak_coordinates, axis=0), axis=1)
    center = np.mean(coordinates, axis=0)
    orbit_scale = float(
        np.sqrt(np.mean(np.sum((coordinates - center) ** 2, axis=1)))
    )
    return float(np.median(closure) / orbit_scale), int(coordinate_indices.size)


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    with np.load(TRAINING, allow_pickle=False) as training:
        modes = np.asarray(training["delay_modes"][:, :2], dtype=float)
        singular = np.asarray(
            training["leading_singular_values"][:2], dtype=float
        )
        mean = float(training["echo_mean"])
        scale = float(training["echo_scale"])
        training_echo = np.asarray(training["echo"], dtype=float)
        training_latent = np.asarray(
            training["latent_coordinates"][:, :2], dtype=float
        )
        training_radius = np.asarray(training["radius"], dtype=float)
        training_velocity = np.asarray(
            training["radial_velocity"], dtype=float
        )
        training_metadata = json.loads(str(training["metadata_json"]))

    delay = int(modes.shape[0])
    sample_step = EXPECTED_STEP * EMBEDDING_DECIMATION
    if delay != int(training_metadata["delay_samples"]):
        raise RuntimeError("retained U is inconsistent with its metadata")

    projected_training, _ = causal_fixed_basis_projection(
        training_echo, modes, singular, mean, scale
    )
    replay_error = float(np.max(np.abs(projected_training - training_latent)))
    if replay_error > 5.0e-13:
        raise RuntimeError("fixed-U online projection does not replay training V")

    sources = sorted(
        (ROOT / "archives/solutions/rpe").glob(SOURCE_GLOB)
    )
    if len(sources) != 5:
        raise RuntimeError(f"expected five retained validation sources, got {len(sources)}")

    state_traces: list[dict] = []
    echo_traces: list[dict] = []
    latent_traces: list[dict] = []
    state_references = [reference(1.0, 0.0, "EC", True)]
    echo_references = [reference(0.0, 0.0, "EC", True)]
    results: list[dict] = []
    saved: dict[str, np.ndarray] = {}

    training_envelope = {
        "radius_min": float(np.min(training_radius)),
        "radius_max": float(np.max(training_radius)),
        "radial_velocity_abs_max": float(np.max(np.abs(training_velocity))),
    }

    for index, (source, expected_ic, color) in enumerate(
        zip(sources, EXPECTED_ICS, COLORS, strict=True), start=1
    ):
        time, state, metadata = load_native_source(source)
        validate_validation_source(metadata)
        initial = np.asarray(metadata["initial_state"], dtype=float)
        if not np.array_equal(initial, expected_ic):
            raise RuntimeError(f"unexpected IC ordering in {source.name}")
        step = float(np.median(np.diff(time)))
        if not np.isclose(step, EXPECTED_STEP, rtol=0.0, atol=1e-12):
            raise RuntimeError("validation trajectory resolution changed")

        echo = far_field_echo(state, metadata)
        echo_dot = np.gradient(echo, time, edge_order=2)
        coordinates, captured_fraction = causal_fixed_basis_projection(
            echo, modes, singular, mean, scale
        )
        latent_time = time[::EMBEDDING_DECIMATION][delay - 1 :]
        closure, cycles = closure_metric(
            echo, coordinates, delay, sample_step
        )

        inside = bool(
            np.min(state[:, 0]) >= training_envelope["radius_min"]
            and np.max(state[:, 0]) <= training_envelope["radius_max"]
            and np.max(np.abs(state[:, 1]))
            <= training_envelope["radial_velocity_abs_max"]
        )
        regime = "inside" if inside else "outside"
        label = (
            f"IC #{index}: ({initial[0]:.2f}, {initial[1]:+.2f}) "
            f"({regime})"
        )

        state_display = np.arange(0, time.size, 20)
        echo_display = np.arange(0, time.size, 4)
        latent_display = np.arange(0, coordinates.shape[0], 2)
        state_traces.append(
            trace(
                label,
                state[state_display, 0],
                state[state_display, 1],
                color=color,
                width=1.6,
            )
        )
        echo_traces.append(
            trace(
                label,
                echo[echo_display],
                echo_dot[echo_display],
                color=color,
                width=1.4,
            )
        )
        latent_traces.append(
            trace(
                label,
                coordinates[latent_display, 0],
                coordinates[latent_display, 1],
                color=color,
                width=1.8,
            )
        )
        state_references.append(
            reference(float(initial[0]), float(initial[1]), f"IC #{index}", False)
        )
        echo_references.append(
            reference(
                float(echo[0]), float(echo_dot[0]), f"IC #{index}", False
            )
        )

        results.append(
            {
                "trajectory": index,
                "source": str(source),
                "initial_state": initial.tolist(),
                "regime_relative_to_training_state_envelope": regime,
                "radius_range": [
                    float(np.min(state[:, 0])),
                    float(np.max(state[:, 0])),
                ],
                "radial_velocity_abs_max": float(
                    np.max(np.abs(state[:, 1]))
                ),
                "echo_range": [float(np.min(echo)), float(np.max(echo))],
                "fixed_U_first_two_window_energy_fraction": captured_fraction,
                "latent_cycle_closure_relative": closure,
                "completed_latent_peak_cycles": cycles,
                "latent_range": [
                    [float(np.min(coordinates[:, 0])), float(np.max(coordinates[:, 0]))],
                    [float(np.min(coordinates[:, 1])), float(np.max(coordinates[:, 1]))],
                ],
            }
        )
        prefix = f"trajectory_{index}"
        saved[f"{prefix}_time"] = time
        saved[f"{prefix}_state"] = state
        saved[f"{prefix}_echo"] = echo
        saved[f"{prefix}_echo_dot"] = echo_dot
        saved[f"{prefix}_latent_time"] = latent_time
        saved[f"{prefix}_latent"] = coordinates

    metrics = {
        "schema": "dynsyslab.echo-fixed-basis-oos/v1",
        "method": "causal fixed-U ring-buffer projection",
        "new_decomposition_performed": False,
        "basis_adaptation_performed": False,
        "training_artifact": str(TRAINING),
        "formula": "xi_i(k) = u_i^T d_k / sigma_i",
        "normalization": {
            "mean": mean,
            "scale": scale,
            "source": "training trajectory only",
        },
        "delay_samples": delay,
        "delay_window_tau": (delay - 1) * sample_step,
        "sample_step_tau": sample_step,
        "training_projection_replay_max_abs_error": replay_error,
        "training_state_envelope": training_envelope,
        "trajectories": results,
    }
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        fixed_delay_modes=modes,
        fixed_singular_values=singular,
        fixed_echo_mean=np.asarray(mean),
        fixed_echo_scale=np.asarray(scale),
        metadata_json=np.asarray(json.dumps(metrics, separators=(",", ":"))),
        **saved,
    )

    figures = [
        phase_tile(
            "RPE state phase space — five unseen ICs",
            "r",
            "dr/dtau",
            state_traces,
            references=state_references,
        ),
        phase_tile(
            "Echo phase space — five unseen ICs",
            "e",
            "de/dtau",
            echo_traces,
            references=echo_references,
        ),
        phase_tile(
            "Fixed-U latent phase space — causal echo histories",
            "xi_1",
            "xi_2",
            latent_traces,
        ),
    ]
    render_viewer(
        TRAINING,
        figures,
        name=PLOT_NAME,
        label="Out-of-sample fixed-U echo projection",
        columns=3,
        width=2400,
        height=900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
