"""Legacy reproduction only: this archived study uses the pre-correction RPE echo.
For new studies use the registered far_field_acoustic_echo and a freshly fitted basis.

Analyze an unforced nonlinear inviscid RPE echo with causal Hankel POD.

The retained solve differs from the viscous reference only by ``mu_pa_s=0``.
Hankel fitting uses the scalar echo alone; radius and radial velocity are
comparison-only truth channels for the requested Viewer phase-space panels.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from scipy.sparse.linalg import svds

from analyze_rpe_echo_hankel_delay_sweep import (
    DELAY_SAMPLES,
    SPECTRUM_MODES,
    energy_rank,
    full_energy_spectrum,
    hankel_matrix,
    nearest_delay,
)
from analyze_rpe_echo_hankel_latents import (
    DECIMATION,
    OUTPUT_ROOT,
    figure,
    fit_echo_hankel,
    phase_figure,
    reference,
    render_viewer,
    stem_traces,
    trace,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "archives/solutions/rpe/"
    "rpe-inviscid-periodic-tau25-dtau0p001__20261001T211832Z.npz"
)
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-hankel-latents.npz"
SWEEP_ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-hankel-delay-sweep.npz"
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-hankel-analysis.metrics.json"
OVERVIEW_PLOT_NAME = "rpe-inviscid-echo-hankel-baseline-overview"
SUMMARY_PLOT_NAME = "rpe-inviscid-echo-hankel-delay-sweep-summary"
ENERGY_PLOT_NAME = "rpe-inviscid-echo-hankel-delay-sweep-energy"
PHASE_PLOT_NAME = "rpe-inviscid-echo-hankel-delay-sweep-phase"

DELAY_SAMPLES_BASELINE = 256
EXPECTED_INITIAL_STATE = np.asarray([2.333386598601735, 0.0])


def energy_summary(hankel: np.ndarray) -> tuple[float, float, float]:
    """Return the leading energy fractions with a robust ARPACK solve.

    The inviscid echo contains narrow collapse spikes and a rich harmonic
    spectrum.  PROPACK's default 20-step budget can fail to converge for some
    delay lengths, whereas the symmetric ARPACK iteration remains reliable.
    """
    total_energy = float(np.sum(hankel * hankel))
    if min(hankel.shape) <= 2:
        singular = np.linalg.svd(hankel, compute_uv=False)
        modal_energy = singular**2 / total_energy
        first = float(modal_energy[0])
        first_two = float(np.sum(modal_energy[:2]))
        balance = float(modal_energy[1] / modal_energy[0])
        return first, first_two, balance
    _, singular, _ = svds(
        hankel,
        k=2,
        which="LM",
        solver="arpack",
        tol=1.0e-12,
        maxiter=10000,
        return_singular_vectors=True,
        rng=np.random.default_rng(20261001),
    )
    singular = np.sort(np.asarray(singular, dtype=float))[::-1]
    modal_energy = singular**2 / total_energy
    return (
        float(modal_energy[0]),
        float(np.sum(modal_energy)),
        float(modal_energy[1] / modal_energy[0]),
    )


def latent_pair(
    hankel: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return the leading two right-singular coordinates robustly."""
    total_energy = float(np.sum(hankel * hankel))
    if min(hankel.shape) <= 2:
        left, singular, right_transpose = np.linalg.svd(
            hankel, full_matrices=False
        )
    else:
        left, singular, right_transpose = svds(
            hankel,
            k=2,
            which="LM",
            solver="arpack",
            tol=1.0e-12,
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
    modal_energy = singular**2 / total_energy
    return right_transpose.T, modal_energy, total_energy


def load_native_source(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load native saved time, state samples, and metadata without pickle."""
    with np.load(path, allow_pickle=False) as archive:
        time = np.asarray(archive["time"], dtype=float)
        state = np.asarray(archive["values"], dtype=float)
        metadata = json.loads(str(archive["metadata_json"]))
    return time, state, metadata


def validate_source(metadata: dict) -> None:
    """Require free, uncontrolled, zero-viscosity data with the retained initial state."""
    parameters = metadata["resolved_system_parameters"]
    if metadata.get("forcing_type") != "free":
        raise RuntimeError("training-policy violation: source trajectory is forced")
    if metadata.get("controller_kind") != "none":
        raise RuntimeError("training-policy violation: controller data were found")
    if float(parameters["mu_pa_s"]) != 0.0:
        raise RuntimeError("inviscid analysis requires mu_pa_s=0 exactly")
    initial = np.asarray(metadata["initial_state"], dtype=float)
    if not np.array_equal(initial, EXPECTED_INITIAL_STATE):
        raise RuntimeError("the retained initial condition does not match the reference")


def far_field_echo(state: np.ndarray, metadata: dict) -> np.ndarray:
    """Evaluate the model far-field echo from archived states and nondimensional groups."""
    groups = metadata["dimensionless_groups"]
    parameters = metadata["resolved_system_parameters"]
    radius = state[:, 0]
    velocity = state[:, 1]
    momentum = (
        float(groups["gas_pressure_ratio"])
        * radius ** (-3.0 * float(parameters["kappa"]))
        - 1.5 * velocity**2
        + float(groups["vapor_pressure_ratio"])
        - float(groups["ambient_pressure_ratio"])
        - float(groups["capillarity"]) / radius
        - float(groups["viscosity"]) * velocity / radius
    )
    return (
        float(parameters["equilibrium_radius_m"])
        * float(groups["pressure_scale_pa"])
        * (2.0 * radius * velocity**2 + momentum)
    )


def echo_period(
    time: np.ndarray, echo: np.ndarray
) -> tuple[float, np.ndarray]:
    """Estimate the echo period from prominent peaks and return their sample indices."""
    sample_interval = float(np.median(np.diff(time)))
    peaks, _ = find_peaks(
        echo,
        distance=max(1, int(round(0.5 / sample_interval))),
        prominence=0.15 * float(np.ptp(echo)),
    )
    if peaks.size < 4:
        raise RuntimeError("too few echo peaks to establish periodicity")
    return float(np.mean(np.diff(time[peaks]))), peaks


def periodicity_metrics(
    time: np.ndarray,
    state: np.ndarray,
    echo: np.ndarray,
    metadata: dict,
) -> dict:
    """Summarize period, orbit variation, and nonlinear-periodicity diagnostics."""
    period, echo_peaks = echo_period(time, echo)
    peak_periods = np.diff(time[echo_peaks])
    radius = state[:, 0]
    radius_peaks, _ = find_peaks(
        radius,
        distance=max(1, int(round(0.5 / np.median(np.diff(time))))),
        prominence=0.1,
    )
    radius_peak_drift = float(
        (radius[radius_peaks[-1]] - radius[radius_peaks[0]])
        / radius[radius_peaks[0]]
    )
    analytic = metadata.get("analytic_solution", {})
    return {
        "period_tau": period,
        "period_relative_to_linear_natural_period": period,
        "completed_cycles": int(echo_peaks.size),
        "period_coefficient_of_variation": float(
            np.std(peak_periods) / np.mean(peak_periods)
        ),
        "radius_peak_relative_drift": radius_peak_drift,
        "radius_min": float(np.min(radius)),
        "radius_max": float(np.max(radius)),
        "radial_velocity_abs_max": float(np.max(np.abs(state[:, 1]))),
        "echo_min": float(np.min(echo)),
        "echo_max": float(np.max(echo)),
        "excursion_asymmetry": float(
            (np.max(radius) - 1.0) / (1.0 - np.min(radius))
        ),
        "analytic_solution_id": analytic.get("solution", {}).get("id"),
        "analytic_solution_status": analytic.get("status"),
        "analytic_maximum_absolute_error": analytic.get(
            "maximum_absolute_error"
        ),
        "nonlinear_regime": bool(
            np.min(radius) < 0.5
            and np.max(radius) > 2.0
            and period > 1.2
        ),
    }


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    time, state, source_metadata = load_native_source(SOURCE)
    validate_source(source_metadata)
    echo = far_field_echo(state, source_metadata)
    echo_dot = np.gradient(echo, time, edge_order=2)
    periodic = periodicity_metrics(time, state, echo, source_metadata)
    if not periodic["nonlinear_regime"]:
        raise RuntimeError("selected inviscid trajectory is not sufficiently nonlinear")

    sample = np.arange(0, time.size, DECIMATION)
    sampled_time = time[sample]
    sampled_echo = echo[sample]
    sample_interval = float(np.median(np.diff(sampled_time)))
    echo_mean = float(np.mean(sampled_echo))
    echo_scale = float(np.std(sampled_echo))
    normalized_echo = (sampled_echo - echo_mean) / echo_scale

    left, singular, energy_fraction, latent, _, _ = fit_echo_hankel(
        sampled_echo, DELAY_SAMPLES_BASELINE
    )
    latent_time = sampled_time[DELAY_SAMPLES_BASELINE - 1 :]
    aligned_sample = sample[DELAY_SAMPLES_BASELINE - 1 :]

    rows: list[dict] = []
    for delay in DELAY_SAMPLES:
        hankel = hankel_matrix(normalized_echo, delay)
        first, first_two, balance = energy_summary(hankel)
        rows.append(
            {
                "delay_samples": delay,
                "delay_window_tau": float((delay - 1) * sample_interval),
                "delay_periods": float(
                    (delay - 1) * sample_interval / periodic["period_tau"]
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
            nearest_delay(0.25 * periodic["period_tau"], sample_interval),
            DELAY_SAMPLES_BASELINE,
            nearest_delay(2.0 * periodic["period_tau"], sample_interval),
            nearest_delay(4.0 * periodic["period_tau"], sample_interval),
        }
    )
    spectra: dict[int, np.ndarray] = {}
    latent_pairs: dict[int, np.ndarray] = {}
    latent_times: dict[int, np.ndarray] = {}
    selected_details: list[dict] = []
    for delay in selected_delays:
        hankel = hankel_matrix(normalized_echo, delay)
        latent_pair_values, pair_energy, _ = latent_pair(hankel)
        full_spectrum = full_energy_spectrum(hankel)
        spectra[delay] = full_spectrum[:SPECTRUM_MODES]
        latent_pairs[delay] = latent_pair_values
        latent_times[delay] = sampled_time[delay - 1 :]
        positive = full_spectrum[full_spectrum > 0.0]
        row = next(item for item in rows if item["delay_samples"] == delay)
        selected_details.append(
            {
                **row,
                "first_48_energy": float(np.sum(spectra[delay])),
                "first_two_energy_recomputed": float(np.sum(pair_energy)),
                "energy_rank": {
                    "90_percent": energy_rank(full_spectrum, 0.90),
                    "95_percent": energy_rank(full_spectrum, 0.95),
                    "99_percent": energy_rank(full_spectrum, 0.99),
                },
                "spectral_entropy_dimension": float(
                    np.exp(-np.sum(positive * np.log(positive)))
                ),
            }
        )

    meaningful = [
        row for row in rows if 0.75 <= row["delay_periods"] <= 4.0
    ]
    best_meaningful = max(meaningful, key=lambda row: row["first_two_energy"])
    analysis_metadata = {
        "schema": "dynsyslab.echo-hankel-inviscid-analysis/v1",
        "label": "RPE inviscid nonlinear periodic echo-only Hankel analysis",
        "source": str(SOURCE),
        "source_forcing_type": "free",
        "training_policy": "unforced-only",
        "training_input_channels": ["far_field_acoustic_echo"],
        "comparison_only_channels": ["r", "dr/dtau"],
        "comparison_channels_used_for_fit": False,
        "initial_state": EXPECTED_INITIAL_STATE.tolist(),
        "mu_pa_s": 0.0,
        "decimation": DECIMATION,
        "sample_interval_tau": sample_interval,
        "baseline_delay_samples": DELAY_SAMPLES_BASELINE,
        "baseline_delay_window_tau": float(
            (DELAY_SAMPLES_BASELINE - 1) * sample_interval
        ),
        "periodicity": periodic,
        "normalization": "global echo mean and standard deviation",
        "latent_coordinate_convention": "rows of V*, sign stabilized",
        "sweep": rows,
        "selected_delays": selected_details,
        "best_meaningful_window": best_meaningful,
        "meaningful_window_rule": "0.75 to 4 echo periods",
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
        singular_values=singular,
        energy_fraction=energy_fraction,
        cumulative_energy=np.cumsum(energy_fraction),
        delay_modes=left,
        echo_mean=np.asarray(echo_mean),
        echo_scale=np.asarray(echo_scale),
        aligned_source_indices=aligned_sample,
        metadata_json=np.asarray(
            json.dumps(analysis_metadata, separators=(",", ":"))
        ),
    )
    np.savez_compressed(
        SWEEP_ARTIFACT,
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
        metadata_json=np.asarray(
            json.dumps(analysis_metadata, separators=(",", ":"))
        ),
    )
    METRICS.write_text(json.dumps(analysis_metadata, indent=2), encoding="utf-8")

    baseline_spectrum = figure(
        "stem",
        "Echo Hankel energy — 256 delays",
        "mode j",
        "sigma_j^2 / sum sigma_i^2",
        stem_traces(
            np.arange(1, energy_fraction.size + 1), energy_fraction
        ),
    )
    baseline_spectrum["value_scale"] = "log"
    baseline_spectrum["square"] = False
    display = np.arange(0, time.size, 4)
    overview_figures = [
        baseline_spectrum,
        phase_figure(
            "RPE state phase space",
            "r",
            "dr/dtau",
            state[display, 0],
            state[display, 1],
            time[display],
            references=[
                reference(float(state[0, 0]), float(state[0, 1]), "IC", False),
                reference(1.0, 0.0, "EC", True),
            ],
        ),
        phase_figure(
            "Echo phase space",
            "e",
            "de/dtau",
            echo[display],
            echo_dot[display],
            time[display],
            references=[
                reference(float(echo[0]), float(echo_dot[0]), "IC", False),
                reference(0.0, 0.0, "EC", True),
            ],
        ),
        phase_figure(
            "Latent phase space 1–2",
            "xi_1",
            "xi_2",
            latent[:, 0],
            latent[:, 1],
            latent_time,
        ),
    ]

    delay_periods = np.asarray([row["delay_periods"] for row in rows])
    first_mode = np.asarray([row["first_mode_energy"] for row in rows])
    first_two = np.asarray([row["first_two_energy"] for row in rows])
    balance = np.asarray([row["second_to_first_energy"] for row in rows])
    summary_figures = [
        figure(
            "line",
            "Energy captured by the first two latent variables",
            "delay window / echo period",
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
            "delay window / echo period",
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
        row = next(item for item in rows if item["delay_samples"] == delay)
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
        pair_values = latent_pairs[delay]
        phase_figures.append(
            phase_figure(
                f"Latent phase space — {row['delay_periods']:.2f} periods",
                "xi_1",
                "xi_2",
                pair_values[:, 0],
                pair_values[:, 1],
                latent_times[delay],
            )
        )

    render_viewer(
        ARTIFACT,
        overview_figures,
        name=OVERVIEW_PLOT_NAME,
        label="RPE inviscid nonlinear periodic — 256-delay overview",
        columns=2,
        width=1900,
        height=1900,
    )
    render_viewer(
        SWEEP_ARTIFACT,
        summary_figures,
        name=SUMMARY_PLOT_NAME,
        label="RPE inviscid nonlinear periodic — delay sweep summary",
        columns=2,
        width=1900,
        height=900,
    )
    render_viewer(
        SWEEP_ARTIFACT,
        energy_figures,
        name=ENERGY_PLOT_NAME,
        label="RPE inviscid nonlinear periodic — delay-dependent energy",
        columns=2,
        width=1900,
        height=2600,
    )
    render_viewer(
        SWEEP_ARTIFACT,
        phase_figures,
        name=PHASE_PLOT_NAME,
        label="RPE inviscid nonlinear periodic — latent phase spaces",
        columns=2,
        width=1900,
        height=2700,
    )
    print(METRICS)
    print(ARTIFACT)
    print(SWEEP_ARTIFACT)


if __name__ == "__main__":
    main()
