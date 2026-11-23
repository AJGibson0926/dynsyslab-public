"""Thin native adapters for analytic periodic sampling and orbit-family filtering."""
from __future__ import annotations

from collections.abc import Mapping
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any

from .workflow import installation_root


def analytic_histories(request: Mapping[str, Any], *, executable: str | Path | None = None,
                       timeout: float = 1800) -> dict[str, Any]:
    """Sample named observables on analytic flows from a batch of initial states.

    Supply system, initial_states, time_offsets, observables and optionally an
    immutable output_archive. Handles complex and vector observables. This uses
    the actual native/physical clock and does not require periodic trajectories.
    """
    return _query("analytic", {"operation": "histories", **dict(request)},
                  "dynsyslab.analytic-query/v1", "dynsyslab.analytic-result/v1", executable, timeout)


def spectral_sweep(request: Mapping[str, Any], *, executable: str | Path | None = None,
                   timeout: float = 1800) -> dict[str, Any]:
    """Fit matched periodic mpEDMD delay dictionaries and compare spectral packets.

    Native code owns quadrature, fitting, frozen-dictionary evaluation and errors.
    See PERIODIC_SPECTRAL_COMPARISON.md for the versioned request contract.
    """
    return _query("periodic-analysis", {"operation": "spectral_sweep", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def fit_delay_model(request: Mapping[str, Any], *, executable: str | Path | None = None,
                    timeout: float = 1800) -> dict[str, Any]:
    """Fit EDMD, ResDMD or mpEDMD on weighted periodic Fourier history data.

    Supply training, channels, time_offsets, phase_count, sample_interval and
    method. Dense/streaming QR preserve the same linear history problem. No
    spectral smoothing or physical equations are needed. Archives are immutable.
    """
    return _query("periodic-analysis", {"operation": "delay_fit", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def audit_delay_forecast(request: Mapping[str, Any], *, executable: str | Path | None = None,
                         timeout: float = 1800) -> dict[str, Any]:
    """Audit a frozen history model by native Fourier translation and Parseval.

    Supply archive/model, validation Fourier orbits and steps_per_orbit, or
    evaluation_steps (one strictly increasing integer list starting at zero per orbit).
    Sparse horizons use native matrix powers, not temporal interpolation. Inputs
    end at zero; the latest coordinate is forecast at the model's discrete clock.
    Returned phase means are for finite series, without a phase sampling grid.
    """
    return _query("periodic-analysis", {"operation": "delay_forecast_audit", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def spectral_compare(request: Mapping[str, Any], *, executable: str | Path | None = None,
                     timeout: float = 1800) -> dict[str, Any]:
    """Compare frozen periodic spectral KROMs without refitting their dictionaries."""
    return _query("periodic-analysis", {"operation": "spectral_compare", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def packet_features(request: Mapping[str, Any], *, executable: str | Path | None = None,
                    timeout: float = 1800) -> dict[str, Any]:
    """Export frozen newest-observable packet coefficients for past-only histories.

    Supply archive/model, channel, packet_centers, bandwidth and order. Returns
    a real/imaginary linear map and its model/clock provenance. Does not refit.
    """
    return _query("periodic-analysis", {"operation": "packet_features", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def identify_history_features(request: Mapping[str, Any], *, executable: str | Path | None = None,
                              timeout: float = 1800) -> dict[str, Any]:
    """Infer orbit label/phase from real linear transforms of measured histories.

    Supply knots, features, histories and independent input noise scales. Full
    transformed covariance whitening defaults on. Coarse training-knot/phase
    seeds precede local fits; truth is not an initializer. Forecasts use the
    fitted atlas frequency, not autonomous linear packet evolution.
    """
    return _query("periodic-analysis", {"operation": "feature_identify", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def quadrature_audit(request: Mapping[str, Any], *, executable: str | Path | None = None,
                     timeout: float = 120) -> dict[str, Any]:
    """Compare explicit product-quadrature moments with finite-Fourier Parseval."""
    return _query("periodic-analysis", {"operation": "quadrature_audit", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1", executable, timeout)


def evaluate_quartic_orbit(request: Mapping[str, Any], *, executable: str | Path | None = None,
                           timeout: float = 120) -> dict[str, Any]:
    """Evaluate a quartic first integral through the native elliptic representation.

    Required: five ascending polynomial coefficients, initial, velocity, elapsed.
    The returned parameter_period is an auxiliary-time period, not a declaration
    of physical time. Use a supplied clock for a time-reparameterized system.
    """
    return _query("special-function", {"function": "quartic_orbit", **dict(request)},
                  "dynsyslab.special-function-query/v1", "dynsyslab.special-function-result/v1",
                  executable, timeout)


def _query(command: str, request: Mapping[str, Any], schema: str,
           result_schema: str, executable: str | Path | None,
           timeout: float) -> dict[str, Any]:
    cli = (Path(executable).expanduser().resolve() if executable else
           installation_root() / "bin" / ("dynsyslab-cli.exe" if os.name == "nt"
                                          else "dynsyslab-cli"))
    if not cli.is_file():
        raise FileNotFoundError(f"Native CLI is missing: {cli}")
    with tempfile.TemporaryDirectory(prefix="dynsyslab-periodic-") as temporary:
        path = Path(temporary) / "query.json"
        path.write_text(json.dumps({"schema": schema, **dict(request)},
                                   allow_nan=False), encoding="utf-8")
        result = subprocess.run([str(cli), command, str(path)], capture_output=True,
                                text=True, encoding="utf-8", timeout=timeout)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    response = json.loads(result.stdout)
    if response.get("schema") != result_schema:
        raise RuntimeError("Unexpected native periodic response schema")
    return response


def evaluate_trajectory(request: Mapping[str, Any], *, executable: str | Path | None = None,
                        timeout: float = 120) -> dict[str, Any]:
    """Evaluate a registered analytic trajectory at explicit native model times.

    Required: system, initial_state, time. No periodicity is assumed here.
    A family may return not_applicable with a reason for the configured inputs.
    """
    return _query("analytic", {"operation": "evaluate", **dict(request)},
                  "dynsyslab.analytic-query/v1", "dynsyslab.analytic-result/v1",
                  executable, timeout)


def sample(request: Mapping[str, Any], *, executable: str | Path | None = None,
           timeout: float = 120) -> dict[str, Any]:
    """Sample a registered analytic orbit uniformly in model time, optionally FFT it.

    Required: system, initial_state, samples, observables. A known native period
    or explicit period is required. No repeated endpoint or Python numerics.
    """
    return _query("analytic", {"operation": "periodic", **dict(request)},
                  "dynsyslab.analytic-query/v1", "dynsyslab.analytic-result/v1",
                  executable, timeout)


def check_coefficient_shift(request: Mapping[str, Any], *, executable: str | Path | None = None,
                            timeout: float = 120) -> dict[str, Any]:
    """Check Fourier coefficient covariance after restarting at shifted initial states.

    Required: system, initial_state, samples, maximum_harmonic, observables, shifts.
    Shifts use native model time. This finite numerical test assumes an autonomous
    periodic flow and does not certify a global Koopman eigenfunction expansion.
    """
    return _query("analytic", {"operation": "coefficient_shift", **dict(request)},
                  "dynsyslab.analytic-query/v1", "dynsyslab.analytic-result/v1",
                  executable, timeout)


def snapshot_pairs(request: Mapping[str, Any], *, executable: str | Path | None = None,
                   timeout: float = 120) -> dict[str, Any]:
    """Write canonical pairs from Fourier orbit channels and explicit orbit masses.

    Required: orbits, phase_count, time_offsets, sample_interval, output_archive.
    Negative offsets are past delays; positive offsets use future values.
    The destination must not exist. No Python numerical construction is used.
    """
    return _query("periodic-analysis", {"operation": "snapshot_pairs", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def filter_orbits(request: Mapping[str, Any], *, executable: str | Path | None = None,
                  timeout: float = 120) -> dict[str, Any]:
    """Causally filter measured channels using a frozen periodic Fourier atlas.

    Uses deterministic parameter/phase particles, with no process noise or
    resampling. The current backend applies only to unforced periodic families.
    """
    return _query("periodic-analysis", {"operation": "filter", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def predict_controlled(request: Mapping[str, Any], *, executable: str | Path | None = None,
                       timeout: float = 120) -> dict[str, Any]:
    """Advance a periodic chart under a scalar input, with optional physical reference.

    Required: system, knots, parameter, initial_phase, final_time, output_step.
    The native backend precomputes the input map using the registered system's B.
    input accepts bias, sines [amplitude,angular_frequency,phase], state_gain,
    state_reference and limit, all in native model units. Feedback uses decoded
    physical state. Only [parameter,phase] evolve in the reduced model. This is
    model-based nonlinear transport, not KLQR or an echo-only estimator.
    """
    return _query("periodic-analysis", {"operation": "controlled_predict", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def identify_controlled(request: Mapping[str, Any], *, executable: str | Path | None = None,
                        timeout: float = 120) -> dict[str, Any]:
    """Fit an initial condition from past measurements under a known scalar input.

    Native single shooting compares model='atlas' or 'physical' using the same
    parameter/phase initial candidates. Input accepts bias and analytic sines,
    or a recorded ZOH applied-input history, not state feedback. Feedthrough is explicit. forecast_final_time extends
    the free rollout without including future measurements in the objective.
    See PERIODIC_CONTROL_TRANSPORT.md for the complete contract and limitations.
    """
    return _query("periodic-analysis", {"operation": "controlled_identify", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def simulate_feedback(request: Mapping[str, Any], *, executable: str | Path | None = None,
                      timeout: float = 120) -> dict[str, Any]:
    """Simulate a noisy sensor -> online shooting observer -> held control loop.

    Requires model/system/atlas, initial_parameters for the plant only, time
    endpoints, sample_time, first_fit_index, fit_every, measurement mapping and
    variances, plus controller state_gain and limit. Initial candidates or a
    coarse_initialization grid initialize the observer from its acquired prefix.
    prediction_substeps controls dense truth/forecast plots. feedback_source is
    estimate (default), truth (explicit full-state baseline), or zero.
    This nonlinear simulation omits process noise and computational delays.
    """
    return _query("periodic-analysis", {"operation": "controlled_feedback", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def replay_controlled(request: Mapping[str, Any], *, executable: str | Path | None = None,
                      timeout: float = 120) -> dict[str, Any]:
    """Forward solve physical/chart dynamics with known analytic or recorded input.

    Requires system, knots, initial_parameters [parameter, phase], time endpoints,
    output_step, measurement_channels. Physical model also needs observables.
    Recorded input: time, values, before_first, measurement_timing (before_command
    or after_command). Applied values are right-held; no extrapolation is allowed.
    This forward solve does not estimate the supplied initial parameters.
    """
    return _query("periodic-analysis", {"operation": "controlled_replay", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def identify_controlled_sequential(request: Mapping[str, Any], *, executable: str | Path | None = None,
                                   timeout: float = 120) -> dict[str, Any]:
    """Causally refit from successive past prefixes and predict between updates.

    Supply the controlled_identify data plus update_indices and optionally
    prediction_substeps. Returned information_time records each estimate's
    measurement cutoff. No past estimates are revised. This is expanding-window
    deterministic identification, not Kalman filtering or finite-window MHE.
    """
    return _query("periodic-analysis", {"operation": "controlled_sequential_identify", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def analyze_torus(request: Mapping[str, Any], *, executable: str | Path | None = None,
                  timeout: float = 120) -> dict[str, Any]:
    """Qualify and sample a registered analytic torus; return multiindex Fourier data.

    Requires system and initial_state. Optional shape, observables, and
    maximum_harmonics select a product angle mesh. Numerical checks are not
    certificates of global integrability or of a KAM theorem.
    """
    return _query("analytic", {"operation": "torus", **dict(request)},
                  "dynsyslab.analytic-query/v1", "dynsyslab.analytic-result/v1",
                  executable, timeout)


def predict_orbit(request: Mapping[str, Any], *, executable: str | Path | None = None,
                  timeout: float = 120) -> dict[str, Any]:
    """Evaluate a frozen atlas at arbitrary times without refitting or observations.

    Required: knots, parameter, initial_phase, time. The phase is anchored at
    initial_time (default zero); table_samples and frozen_frequency are optional.
    """
    return _query("periodic-analysis", {"operation": "predict", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)


def identify_orbit(request: Mapping[str, Any], *, executable: str | Path | None = None,
                   timeout: float = 120) -> dict[str, Any]:
    """Refine orbit parameter and initial phase from a supplied past measurement window.

    Requires initial_candidates. Optional prediction_times produce a free forecast
    after fitting; no observations at those additional times enter the objective.
    """
    return _query("periodic-analysis", {"operation": "identify", **dict(request)},
                  "dynsyslab.periodic-query/v1", "dynsyslab.periodic-result/v1",
                  executable, timeout)
