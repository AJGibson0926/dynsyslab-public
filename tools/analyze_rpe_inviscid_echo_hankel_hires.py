"""Resolve the nonlinear inviscid RPE echo at tenfold time resolution.

The solve output step is 1e-4 tau.  Hankel samples are retained every fifth
output, giving a 5e-4 tau latent sampling interval (ten times finer than the
earlier 5e-3 tau analysis), while the physical 1.275 tau delay window is held
fixed.  Radius and velocity remain comparison-only truth channels.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.sparse.linalg import svds

from analyze_rpe_echo_hankel_latents import (
    OUTPUT_ROOT,
    figure,
    phase_figure,
    reference,
    render_viewer,
    stem_traces,
    trace,
)
from analyze_rpe_inviscid_echo_hankel import (
    EXPECTED_INITIAL_STATE,
    far_field_echo,
    load_native_source,
    periodicity_metrics,
    validate_source,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "archives/solutions/rpe/"
    "rpe-inviscid-periodic-tau25-dtau0p0001__20261001T212708Z.npz"
)
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-hankel-hires-latents.npz"
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-hankel-hires.metrics.json"
OVERVIEW_NAME = "rpe-inviscid-echo-hankel-hires-overview"
RESOLUTION_NAME = "rpe-inviscid-echo-hires-resolution"

SOURCE_OUTPUT_STEP = 1.0e-4
EMBEDDING_DECIMATION = 5
EMBEDDING_STEP = SOURCE_OUTPUT_STEP * EMBEDDING_DECIMATION
DELAY_WINDOW_TAU = 1.275
DELAY_SAMPLES = int(round(DELAY_WINDOW_TAU / EMBEDDING_STEP)) + 1
LEADING_MODE_COUNT = 48
LATENT_COUNT = 6


def leading_hankel_svd(
    echo: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Fit leading high-resolution echo Hankel modes and retain normalization and energy."""
    mean = float(np.mean(echo))
    scale = float(np.std(echo))
    normalized = (echo - mean) / scale
    windows = np.lib.stride_tricks.sliding_window_view(
        normalized, DELAY_SAMPLES
    )
    hankel = np.ascontiguousarray(windows.T)
    total_energy = float(np.sum(hankel * hankel))
    print(
        f"high-resolution Hankel {hankel.shape}, "
        f"{hankel.nbytes / 1.0e9:.3f} GB; computing leading modes",
        flush=True,
    )
    left, singular, right_transpose = svds(
        hankel,
        k=LEADING_MODE_COUNT,
        which="LM",
        solver="arpack",
        tol=1.0e-9,
        maxiter=10000,
        return_singular_vectors=True,
        rng=np.random.default_rng(20261001),
    )
    order = np.argsort(singular)[::-1]
    left = left[:, order]
    singular = singular[order]
    right_transpose = right_transpose[order, :]
    for mode in range(singular.size):
        pivot = int(np.argmax(np.abs(left[:, mode])))
        if left[pivot, mode] < 0.0:
            left[:, mode] *= -1.0
            right_transpose[mode, :] *= -1.0
    energy_fraction = singular**2 / total_energy
    latent = right_transpose[:LATENT_COUNT, :].T
    return left, singular, energy_fraction, latent, mean, scale


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    time, state, metadata = load_native_source(SOURCE)
    validate_source(metadata)
    output_step = float(np.median(np.diff(time)))
    if not np.isclose(output_step, SOURCE_OUTPUT_STEP, rtol=0.0, atol=1e-12):
        raise RuntimeError("source is not the requested tenfold-resolution solve")

    echo = far_field_echo(state, metadata)
    echo_dot = np.gradient(echo, time, edge_order=2)
    periodic = periodicity_metrics(time, state, echo, metadata)

    peak_indices, _ = find_peaks(
        echo,
        distance=max(1, int(round(0.5 / output_step))),
        prominence=0.15 * float(np.ptp(echo)),
    )
    widths = peak_widths(echo, peak_indices, rel_height=0.5)[0] * output_step
    median_width = float(np.median(widths))

    sample_indices = np.arange(0, time.size, EMBEDDING_DECIMATION)
    sampled_time = time[sample_indices]
    sampled_echo = echo[sample_indices]
    (
        delay_modes,
        singular,
        leading_energy,
        latent,
        echo_mean,
        echo_scale,
    ) = leading_hankel_svd(sampled_echo)
    latent_time = sampled_time[DELAY_SAMPLES - 1 :]
    aligned_indices = sample_indices[DELAY_SAMPLES - 1 :]

    metrics = {
        "schema": "dynsyslab.echo-hankel-inviscid-hires/v1",
        "source": str(SOURCE),
        "training_policy": "unforced-only",
        "training_input_channels": ["far_field_acoustic_echo"],
        "comparison_only_channels": ["r", "dr/dtau"],
        "comparison_channels_used_for_fit": False,
        "initial_state": EXPECTED_INITIAL_STATE.tolist(),
        "mu_pa_s": 0.0,
        "source_output_step_tau": output_step,
        "resolution_increase_over_previous_source": 10.0,
        "embedding_decimation": EMBEDDING_DECIMATION,
        "embedding_sample_step_tau": EMBEDDING_STEP,
        "resolution_increase_over_previous_embedding": 10.0,
        "delay_samples": DELAY_SAMPLES,
        "delay_window_tau": float((DELAY_SAMPLES - 1) * EMBEDDING_STEP),
        "leading_modes_computed": LEADING_MODE_COUNT,
        "first_mode_energy": float(leading_energy[0]),
        "first_two_energy": float(np.sum(leading_energy[:2])),
        "first_six_energy": float(np.sum(leading_energy[:6])),
        "first_48_energy": float(np.sum(leading_energy)),
        "echo_peak_fwhm_tau": median_width,
        "source_samples_per_echo_fwhm": median_width / output_step,
        "embedding_samples_per_echo_fwhm": median_width / EMBEDDING_STEP,
        "periodicity": periodic,
    }

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        ARTIFACT,
        time=time,
        radius=state[:, 0],
        radial_velocity=state[:, 1],
        echo=echo,
        echo_dot_diagnostic=echo_dot,
        embedding_time=latent_time,
        latent_coordinates=latent,
        leading_singular_values=singular,
        leading_energy_fraction=leading_energy,
        delay_modes=delay_modes,
        echo_mean=np.asarray(echo_mean),
        echo_scale=np.asarray(echo_scale),
        aligned_source_indices=aligned_indices,
        metadata_json=np.asarray(json.dumps(metrics, separators=(",", ":"))),
    )
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    energy = figure(
        "stem",
        "Leading echo Hankel energy — 1.275 tau window",
        "mode j",
        "sigma_j^2 / sum sigma_i^2",
        stem_traces(np.arange(1, LEADING_MODE_COUNT + 1), leading_energy),
    )
    energy["value_scale"] = "log"
    energy["square"] = False

    state_display = np.arange(0, time.size, 40)
    echo_display = np.arange(0, time.size, 4)
    overview = [
        energy,
        phase_figure(
            "RPE state phase space",
            "r",
            "dr/dtau",
            state[state_display, 0],
            state[state_display, 1],
            time[state_display],
            references=[
                reference(float(state[0, 0]), float(state[0, 1]), "IC", False),
                reference(1.0, 0.0, "EC", True),
            ],
        ),
        phase_figure(
            "Echo phase space — 10x resolution",
            "e",
            "de/dtau",
            echo[echo_display],
            echo_dot[echo_display],
            time[echo_display],
            references=[
                reference(float(echo[0]), float(echo_dot[0]), "IC", False),
                reference(0.0, 0.0, "EC", True),
            ],
        ),
        phase_figure(
            "Latent phase space 1–2 — 10x sampling",
            "xi_1",
            "xi_2",
            latent[:, 0],
            latent[:, 1],
            latent_time,
        ),
    ]
    render_viewer(
        SOURCE,
        overview,
        name=OVERVIEW_NAME,
        label="Inviscid RPE echo-only Hankel — 10x resolution",
        columns=2,
        width=1900,
        height=1900,
    )

    peak = int(peak_indices[0])
    zoom = np.flatnonzero(np.abs(time - time[peak]) <= 0.02)
    old_grid = zoom[::10]
    cycle = np.flatnonzero(
        (time >= time[peak_indices[0]]) & (time <= time[peak_indices[1]])
    )
    resolution_figures = [
        figure(
            "line",
            "Resolved collapse echo pulse",
            "tau",
            "e",
            [
                trace(
                    "Delta tau = 0.0001",
                    time[zoom],
                    echo[zoom],
                    color="#0ea5e9",
                    width=2.0,
                ),
                trace(
                    "old 0.001 grid",
                    time[old_grid],
                    echo[old_grid],
                    color="#f97316",
                    marker=True,
                    marker_size=4.0,
                    width=1.0,
                ),
            ],
        ),
        phase_figure(
            "One resolved echo cycle",
            "e",
            "de/dtau",
            echo[cycle],
            echo_dot[cycle],
            time[cycle],
        ),
    ]
    render_viewer(
        SOURCE,
        resolution_figures,
        name=RESOLUTION_NAME,
        label="Inviscid RPE echo resolution check",
        columns=2,
        width=1900,
        height=900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
