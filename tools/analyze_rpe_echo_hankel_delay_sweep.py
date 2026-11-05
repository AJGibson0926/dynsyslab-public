"""Sweep causal Hankel delay length on one unforced RPE echo record.

The fit uses only the scalar echo. Radius and radial velocity never enter the
Hankel decomposition. The experiment is deliberately separate from the core
training and control pipelines and writes only ignored Viewer outputs.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import svds

from analyze_rpe_echo_hankel_latents import (
    DECIMATION,
    OUTPUT_ROOT,
    SOURCE,
    dimensionless_history,
    far_field_echo,
    figure,
    load_source,
    phase_figure,
    render_viewer,
    stem_traces,
    trace,
)


ARTIFACT = OUTPUT_ROOT / "rpe-unforced-echo-hankel-delay-sweep.npz"
METRICS = OUTPUT_ROOT / "rpe-unforced-echo-hankel-delay-sweep.metrics.json"
SUMMARY_PLOT_NAME = "rpe-unforced-echo-hankel-delay-sweep-summary"
ENERGY_PLOT_NAME = "rpe-unforced-echo-hankel-delay-sweep-energy"
PHASE_PLOT_NAME = "rpe-unforced-echo-hankel-delay-sweep-phase"

DELAY_SAMPLES = (
    2,
    3,
    4,
    5,
    6,
    8,
    12,
    16,
    24,
    32,
    48,
    64,
    96,
    128,
    160,
    192,
    224,
    256,
    320,
    384,
    448,
    512,
    640,
    768,
    896,
    1024,
    1280,
    1536,
    1792,
    2048,
)
SPECTRUM_MODES = 48


def dominant_period(signal: np.ndarray, sample_interval: float) -> float:
    """Return the strongest nonzero tapered-FFT period in tau units."""
    centered = signal - np.mean(signal)
    tapered = centered * np.hanning(centered.size)
    frequencies = np.fft.rfftfreq(centered.size, sample_interval)
    power = np.abs(np.fft.rfft(tapered)) ** 2
    usable = (frequencies >= 0.1) & (frequencies <= 8.0)
    peak = int(np.flatnonzero(usable)[np.argmax(power[usable])])
    return float(1.0 / frequencies[peak])


def hankel_matrix(normalized_echo: np.ndarray, delay_samples: int) -> np.ndarray:
    """Build contiguous delay-by-window columns from the normalized scalar echo."""
    windows = np.lib.stride_tricks.sliding_window_view(
        normalized_echo, delay_samples
    )
    return np.ascontiguousarray(windows.T)


def leading_svd(
    hankel: np.ndarray,
    count: int,
    *,
    vectors: bool,
) -> tuple[np.ndarray | None, np.ndarray, np.ndarray | None]:
    """Return descending leading singular triplets deterministically."""
    available = min(hankel.shape)
    requested = min(count, available)
    if requested == available:
        if vectors:
            left, singular, right_transpose = np.linalg.svd(
                hankel, full_matrices=False
            )
            return left, singular, right_transpose
        singular = np.linalg.svd(hankel, compute_uv=False)
        return None, singular, None

    result = svds(
        hankel,
        k=requested,
        which="LM",
        return_singular_vectors=vectors,
        solver="propack",
        rng=np.random.default_rng(20261001),
    )
    if vectors:
        left, singular, right_transpose = result
        order = np.argsort(singular)[::-1]
        left = left[:, order]
        singular = singular[order]
        right_transpose = right_transpose[order, :]
        for mode in range(singular.size):
            pivot = int(np.argmax(np.abs(left[:, mode])))
            if left[pivot, mode] < 0.0:
                left[:, mode] *= -1.0
                right_transpose[mode, :] *= -1.0
        return left, singular, right_transpose

    singular = np.sort(np.asarray(result, dtype=float))[::-1]
    return None, singular, None


def energy_summary(hankel: np.ndarray) -> tuple[float, float, float]:
    """Summarize leading Hankel singular-value energy relative to the full matrix."""
    total_energy = float(np.sum(hankel * hankel))
    _, singular, _ = leading_svd(hankel, 2, vectors=False)
    modal_energy = singular**2 / total_energy
    return (
        float(modal_energy[0]),
        float(np.sum(modal_energy[:2])),
        float(modal_energy[1] / modal_energy[0]),
    )


def full_energy_spectrum(hankel: np.ndarray) -> np.ndarray:
    """Return the exact modal-energy spectrum through the smaller Gramian."""
    gramian = hankel @ hankel.T
    eigenvalues = np.linalg.eigvalsh(gramian)[::-1]
    energy = np.maximum(eigenvalues, 0.0)
    return energy / np.sum(energy)


def energy_rank(energy_fraction: np.ndarray, threshold: float) -> int:
    """Return the first cumulative-energy rank reaching the requested fraction."""
    return int(np.searchsorted(np.cumsum(energy_fraction), threshold) + 1)


def latent_pair(
    hankel: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return the leading right-singular coordinates, energy fractions, and total energy."""
    total_energy = float(np.sum(hankel * hankel))
    _, singular, right_transpose = leading_svd(hankel, 2, vectors=True)
    if right_transpose is None:
        raise RuntimeError("right singular vectors were not returned")
    modal_energy = singular**2 / total_energy
    return right_transpose.T, modal_energy, total_energy


def nearest_delay(target: float, sample_interval: float) -> int:
    """Choose the tested delay span nearest the requested physical duration."""
    return min(
        DELAY_SAMPLES,
        key=lambda delay: abs((delay - 1) * sample_interval - target),
    )


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    physical_time, physical_state, source_metadata = load_source(SOURCE)
    forcing = source_metadata.get("forcing", {})
    if source_metadata.get("forcing_type") != "free" or forcing.get(
        "forcing_present", False
    ):
        raise RuntimeError("training-policy violation: source trajectory is forced")
    if source_metadata.get("control_history", {}).get("present", False):
        raise RuntimeError("training-policy violation: control data were found")

    tau, hidden_state, radius0, _ = dimensionless_history(
        physical_time, physical_state, source_metadata
    )
    echo = far_field_echo(hidden_state, source_metadata, radius0)
    sample = np.arange(0, echo.size, DECIMATION)
    sampled_tau = tau[sample]
    sampled_echo = echo[sample]
    sample_interval = float(np.median(np.diff(sampled_tau)))
    echo_mean = float(np.mean(sampled_echo))
    echo_scale = float(np.std(sampled_echo))
    normalized_echo = (sampled_echo - echo_mean) / echo_scale
    period = dominant_period(normalized_echo, sample_interval)

    rows: list[dict] = []
    for delay in DELAY_SAMPLES:
        hankel = hankel_matrix(normalized_echo, delay)
        first, first_two, balance = energy_summary(hankel)
        rows.append(
            {
                "delay_samples": delay,
                "delay_window_tau": float((delay - 1) * sample_interval),
                "delay_periods": float(
                    (delay - 1) * sample_interval / period
                ),
                "hankel_columns": int(hankel.shape[1]),
                "first_mode_energy": first,
                "first_two_energy": first_two,
                "tail_after_two": 1.0 - first_two,
                "second_to_first_energy": balance,
            }
        )
        print(
            f"delay={delay:4d} window={rows[-1]['delay_window_tau']:.3f} "
            f"periods={rows[-1]['delay_periods']:.3f} "
            f"E12={first_two:.8f}"
        )

    selected_delays = sorted(
        {
            2,
            8,
            nearest_delay(0.25 * period, sample_interval),
            256,
            nearest_delay(2.0 * period, sample_interval),
            nearest_delay(4.0 * period, sample_interval),
        }
    )

    spectra: dict[int, np.ndarray] = {}
    latent_pairs: dict[int, np.ndarray] = {}
    latent_times: dict[int, np.ndarray] = {}
    selected_details: list[dict] = []
    for delay in selected_delays:
        hankel = hankel_matrix(normalized_echo, delay)
        latent, pair_energy, _ = latent_pair(hankel)
        full_spectrum = full_energy_spectrum(hankel)
        spectrum = full_spectrum[:SPECTRUM_MODES]
        spectra[delay] = spectrum
        latent_pairs[delay] = latent
        latent_times[delay] = sampled_tau[delay - 1 :]

        positive = full_spectrum[full_spectrum > 0.0]
        entropy_dimension = float(
            np.exp(-np.sum(positive * np.log(positive)))
        )

        source_row = next(row for row in rows if row["delay_samples"] == delay)
        selected_details.append(
            {
                **source_row,
                "first_48_energy": float(np.sum(spectrum)),
                "first_two_energy_recomputed": float(np.sum(pair_energy)),
                "energy_rank": {
                    "90_percent": energy_rank(full_spectrum, 0.90),
                    "95_percent": energy_rank(full_spectrum, 0.95),
                    "99_percent": energy_rank(full_spectrum, 0.99),
                },
                "spectral_entropy_dimension": entropy_dimension,
            }
        )

    meaningful = [
        row
        for row in rows
        if 0.75 <= row["delay_periods"] <= 4.0
    ]
    best_meaningful = max(meaningful, key=lambda row: row["first_two_energy"])

    metadata = {
        "schema": "dynsyslab.echo-hankel-delay-sweep/v1",
        "label": "RPE free decay — echo-only Hankel delay sweep",
        "source": str(SOURCE),
        "source_forcing_type": "free",
        "training_policy": "unforced-only",
        "training_input_channels": ["far_field_acoustic_echo"],
        "comparison_channels_used_for_fit": False,
        "decimation": DECIMATION,
        "sample_interval_tau": sample_interval,
        "dominant_echo_period_tau": period,
        "dominant_echo_period_samples": period / sample_interval,
        "normalization": "global echo mean and standard deviation",
        "energy_definition": "sigma_j^2 / Frobenius(H)^2",
        "phase_coordinate_convention": "first two columns of V, sign stabilized",
        "sweep": rows,
        "selected_delays": selected_details,
        "best_meaningful_window": best_meaningful,
        "meaningful_window_rule": "0.75 to 4 dominant echo periods",
    }

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        ARTIFACT,
        delay_samples=np.asarray(DELAY_SAMPLES, dtype=int),
        delay_window_tau=np.asarray(
            [row["delay_window_tau"] for row in rows], dtype=float
        ),
        delay_periods=np.asarray(
            [row["delay_periods"] for row in rows], dtype=float
        ),
        first_mode_energy=np.asarray(
            [row["first_mode_energy"] for row in rows], dtype=float
        ),
        first_two_energy=np.asarray(
            [row["first_two_energy"] for row in rows], dtype=float
        ),
        tail_after_two=np.asarray(
            [row["tail_after_two"] for row in rows], dtype=float
        ),
        selected_delay_samples=np.asarray(selected_delays, dtype=int),
        metadata_json=np.asarray(json.dumps(metadata, separators=(",", ":"))),
    )
    METRICS.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    delay_periods = np.asarray([row["delay_periods"] for row in rows])
    first_mode = np.asarray([row["first_mode_energy"] for row in rows])
    first_two = np.asarray([row["first_two_energy"] for row in rows])
    balance = np.asarray([row["second_to_first_energy"] for row in rows])
    summary_figures = [
        figure(
            "line",
            "Energy captured by the first two latent variables",
            "delay window / dominant echo period",
            "cumulative energy fraction",
            [
                trace(
                    "mode 1",
                    delay_periods,
                    first_mode,
                    color="#f59e0b",
                    marker=True,
                ),
                trace(
                    "modes 1–2",
                    delay_periods,
                    first_two,
                    color="#0ea5e9",
                    marker=True,
                ),
            ],
        ),
        figure(
            "line",
            "Balance within the leading oscillatory pair",
            "delay window / dominant echo period",
            "E_2 / E_1",
            [
                trace(
                    "E_2 / E_1",
                    delay_periods,
                    balance,
                    color="#8b5cf6",
                    marker=True,
                )
            ],
        ),
    ]

    energy_figures: list[dict] = []
    phase_figures: list[dict] = []
    for delay in selected_delays:
        row = next(row for row in rows if row["delay_samples"] == delay)
        spectrum = spectra[delay]
        energy_plot = figure(
            "stem",
            (
                f"Delay {row['delay_periods']:.2f} periods "
                f"({row['delay_window_tau']:.3f} tau)"
            ),
            "mode j",
            "sigma_j^2 / sum sigma_i^2",
            stem_traces(np.arange(1, spectrum.size + 1), spectrum),
        )
        energy_plot["value_scale"] = "log"
        energy_figures.append(energy_plot)

        latent = latent_pairs[delay]
        phase_figures.append(
            phase_figure(
                (
                    f"Latent phase space — {row['delay_periods']:.2f} "
                    "periods"
                ),
                "xi_1",
                "xi_2",
                latent[:, 0],
                latent[:, 1],
                latent_times[delay],
            )
        )

    render_viewer(
        ARTIFACT,
        summary_figures,
        name=SUMMARY_PLOT_NAME,
        label="RPE free decay — delay sweep summary",
        columns=2,
        width=1900,
        height=900,
    )
    render_viewer(
        ARTIFACT,
        energy_figures,
        name=ENERGY_PLOT_NAME,
        label="RPE free decay — delay-dependent Hankel energy",
        columns=2,
        width=1900,
        height=2600,
    )
    render_viewer(
        ARTIFACT,
        phase_figures,
        name=PHASE_PLOT_NAME,
        label="RPE free decay — delay-dependent latent phase spaces",
        columns=2,
        width=1900,
        height=2700,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
