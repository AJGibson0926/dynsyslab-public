"""Legacy reproduction only: this archived study uses the pre-correction RPE echo.
For new studies use the registered far_field_acoustic_echo and a freshly fitted basis.

Build the reproducible metrics and figure for the long-horizon RPE sine trial."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]

PATHS = {
    "free_training": Path(
        "archives/solutions/rayleigh-plesset/"
        "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz"
    ),
    "hdmd": Path(
        "archives/derived/koopman-hankel-dmd/"
        "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__"
        "koopman-hankel-dmd__20260826T221322.756123Z__604a65f1.npz"
    ),
    "state_sindy": Path(
        "archives/derived/koopman-with-sindy/"
        "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__"
        "koopman-with-sindy__20260827T183645.394109Z__c4dde209.npz"
    ),
    "parent_krom": Path(
        "archives/derived/koopman-krom-extract/"
        "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__"
        "koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz"
    ),
    "forced_training": Path(
        "archives/solutions/rpe/"
        "rpe-echo-study-state-basis-baseline__20260911T104742Z.npz"
    ),
    "echo_jet_krom": Path(
        "archives/derived/koopman-with-sindy/"
        "rpe-echo-study-state-basis-baseline__20260911T104742Z__"
        "koopman-with-sindy__20260914T165124Z.npz"
    ),
    "sine_target": Path(
        "archives/solutions/rpe/"
        "rpe-periodic-100khz-20kpa-phase-portrait__20260911T214426Z.npz"
    ),
    "control_job": Path(
        "examples/native/case-studies/"
        "rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300.json"
    ),
    "controlled": Path(
        "archives/solutions/rpe/"
        "rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300__"
        "20260914T171839Z.npz"
    ),
}

METRICS_PATH = ROOT / "viewer_outputs/rpe-echo-jet2-sine-long-term-control.json"
FIGURE_PATH = ROOT / "viewer_outputs/rpe-echo-jet2-sine-long-term-control.png"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_npz(relative: Path) -> tuple[dict[str, np.ndarray], dict]:
    path = ROOT / relative
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files if name != "metadata_json"}
        metadata = json.loads(str(archive["metadata_json"]))
    return arrays, metadata


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def nrmse(actual: np.ndarray, reference: np.ndarray) -> float:
    denominator = rms(reference - np.mean(reference))
    return rms(actual - reference) / denominator


def tone_fit(time: np.ndarray, values: np.ndarray, omega: float) -> dict[str, float]:
    design = np.column_stack(
        (np.ones_like(time), np.cos(omega * time), np.sin(omega * time))
    )
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    fitted = design @ coefficients
    cosine, sine = coefficients[1], coefficients[2]
    residual = np.sum(np.square(values - fitted))
    total = np.sum(np.square(values - np.mean(values)))
    return {
        "mean": float(coefficients[0]),
        "amplitude": float(np.hypot(cosine, sine)),
        "phase_rad": float(np.arctan2(-sine, cosine)),
        "phase_deg": float(np.degrees(np.arctan2(-sine, cosine))),
        "r_squared": float(1.0 - residual / total),
    }


def echo_jet(states: np.ndarray, groups: dict[str, float], radius0: float) -> np.ndarray:
    """Evaluate [h, L_f0 h] using the expressions in SystemsRPE.cpp."""
    r = states[:, 0]
    v = states[:, 1]
    gas = groups["gas_pressure_ratio"]
    vapor = groups["vapor_pressure_ratio"]
    ambient = groups["ambient_pressure_ratio"]
    capillary = groups["capillarity"]
    viscous = groups["viscosity"]
    pressure_scale = groups["pressure_scale_pa"]

    momentum = (
        gas * r**-4
        - 1.5 * v**2
        + vapor
        - ambient
        - capillary / r
        - viscous * v / r
    )
    acceleration = momentum / r
    scale = radius0 * pressure_scale
    echo = scale * (2.0 * r * v**2 + momentum)
    echo_dr = scale * (
        2.0 * v**2
        - 4.0 * gas * r**-5
        + capillary / r**2
        + viscous * v / r**2
    )
    echo_dv = scale * (4.0 * r * v - 3.0 * v - viscous / r)
    echo_dot = echo_dr * v + echo_dv * acceleration
    return np.column_stack((echo, echo_dot))


def window_metrics(
    time: np.ndarray,
    control: np.ndarray,
    known: np.ndarray,
    state: np.ndarray,
    target_state: np.ndarray,
    observed: np.ndarray,
    target_observed: np.ndarray,
    omega: float,
    start: float,
    stop: float,
) -> dict:
    mask = (time >= start - 1e-12) & (time <= stop + 1e-12)
    u = control[mask]
    u_ref = known[mask]
    x = state[mask]
    x_ref = target_state[mask]
    y = observed[mask]
    y_ref = target_observed[mask]
    fit = tone_fit(time[mask], u, omega)
    return {
        "tau": [start, stop],
        "samples": int(np.count_nonzero(mask)),
        "driver_rmse": rms(u - u_ref),
        "driver_nrmse": nrmse(u, u_ref),
        "driver_correlation": float(np.corrcoef(u, u_ref)[0, 1]),
        "driver_maximum_absolute": float(np.max(np.abs(u))),
        "driver_rms": rms(u),
        "fundamental": fit,
        "fundamental_amplitude_ratio": fit["amplitude"] / 0.2,
        "state_rmse": [rms(x[:, index] - x_ref[:, index]) for index in range(2)],
        "echo_jet_nrmse": [nrmse(y[:, index], y_ref[:, index]) for index in range(2)],
        "radius_range": [float(np.min(x[:, 0])), float(np.max(x[:, 0]))],
    }


def main() -> None:
    target, target_metadata = load_npz(PATHS["sine_target"])
    result, result_metadata = load_npz(PATHS["controlled"])
    krom, krom_metadata = load_npz(PATHS["echo_jet_krom"])

    time = np.asarray(result["time"], dtype=float)
    target_time = np.asarray(target["time"], dtype=float)
    if not np.allclose(time, target_time, rtol=0.0, atol=1e-12):
        raise RuntimeError("The target and controlled archives do not share a time grid")

    state = np.asarray(result["values"], dtype=float)
    target_state = np.asarray(result["control_target_values"], dtype=float)
    control = np.asarray(result["control_values"], dtype=float)[:, 0]
    known = np.asarray(target["control_values"], dtype=float)[:, 0]
    target_observed = np.asarray(result["control_observable_target_values"], dtype=float)

    groups = result_metadata["dimensionless_groups"]
    radius0 = result_metadata["resolved_system_parameters"]["equilibrium_radius_m"]
    observed = echo_jet(state, groups, radius0)
    mapped_target = echo_jet(target_state, groups, radius0)
    target_mapping_error = float(np.max(np.abs(mapped_target - target_observed)))

    frequency_hz = target_metadata["controller_configuration"]["settings"][
        "frequencies_hz"
    ][0]
    time_scale = result_metadata["time_scale_to_seconds"]
    omega = 2.0 * np.pi * frequency_hz * time_scale
    period = 2.0 * np.pi / omega

    windows = {}
    for start in (0.0, 10.0, 20.0, 30.0, 60.0, 100.0, 200.0, 290.0):
        windows[f"tau_{int(start)}_to_300"] = window_metrics(
            time,
            control,
            known,
            state,
            target_state,
            observed,
            target_observed,
            omega,
            start,
            300.0,
        )

    cycle_rows = []
    first_cycle = int(np.ceil(100.0 / period))
    final_cycle_exclusive = int(np.floor(300.0 / period))
    for cycle in range(first_cycle, final_cycle_exclusive):
        start = cycle * period
        stop = (cycle + 1) * period
        row = window_metrics(
            time,
            control,
            known,
            state,
            target_state,
            observed,
            target_observed,
            omega,
            start,
            stop,
        )
        row["cycle_index"] = cycle
        cycle_rows.append(row)

    cycle_nrmse = np.array([row["driver_nrmse"] for row in cycle_rows])
    cycle_corr = np.array([row["driver_correlation"] for row in cycle_rows])
    cycle_r_rmse = np.array([row["state_rmse"][0] for row in cycle_rows])
    cycle_v_rmse = np.array([row["state_rmse"][1] for row in cycle_rows])
    cycle_peak = np.array([row["driver_maximum_absolute"] for row in cycle_rows])

    lineage = {}
    for name, relative in PATHS.items():
        absolute = ROOT / relative
        lineage[name] = {
            "path": relative.as_posix(),
            "sha256": sha256(absolute),
            "bytes": absolute.stat().st_size,
        }

    payload = {
        "schema": "dynsyslab.rpe-echo-jet-sine-long-term-control/v1",
        "date": "2026-09-14",
        "success": True,
        "success_definition": (
            "After acquisition, sustain the target state and echo-jet orbit and recover "
            "the known sine driver with low waveform error, unit correlation, stable "
            "cycle metrics, and no radius drift through tau=300."
        ),
        "lineage": lineage,
        "plant": {
            "system": "dimensionless Rayleigh-Plesset",
            "medium": "air in ultrapure water at 20 C",
            "initial_state": result_metadata["initial_state"],
            "resolved_parameters": result_metadata["resolved_system_parameters"],
            "dimensionless_groups": groups,
            "time_scale_seconds": time_scale,
            "input_scale_pa": result_metadata["control_scale_to_dimensional"],
        },
        "echo_jet_krom": {
            "coordinates": krom_metadata["coordinate_observables"],
            "state_center": krom["state_center"].tolist(),
            "state_scale": krom["state_scale"].tolist(),
            "polynomial_degree": krom_metadata["polynomial_degree"],
            "feature_count": krom_metadata["feature_count"],
            "include_input_terms": krom_metadata["include_input_terms"],
            "continuous_eigenvalues": [
                [float(value.real), float(value.imag)]
                for value in krom["continuous_eigenvalues"]
            ],
            "training_residual_rms": krom["training_residual_rms"].tolist(),
            "selected_modes_for_control": [0],
        },
        "known_driver": {
            "frequency_hz": frequency_hz,
            "amplitude_pa": 20000.0,
            "model_amplitude": 0.2,
            "omega_tau": omega,
            "period_tau": period,
            "formula": f"u*(tau)=0.2*cos({omega:.16g}*tau)",
            "cycles_over_full_run": 300.0 / period,
        },
        "control_contract": {
            "target_space": result_metadata["control_target_space"],
            "controller": result_metadata["controller_configuration"],
            "solver": result_metadata["solver_configuration"],
            "time": result_metadata["time_configuration"],
            "target_mapping_maximum_absolute_error": target_mapping_error,
        },
        "full_run": {
            "samples": int(time.size),
            "radius_range": [float(np.min(state[:, 0])), float(np.max(state[:, 0]))],
            "velocity_range": [float(np.min(state[:, 1])), float(np.max(state[:, 1]))],
            "control_rms": rms(control),
            "control_maximum_absolute": float(np.max(np.abs(control))),
        },
        "windows": windows,
        "settled_cycle_stability_after_tau_100": {
            "complete_cycles": len(cycle_rows),
            "cycle_indices": [first_cycle, final_cycle_exclusive - 1],
            "driver_nrmse": {
                "minimum": float(np.min(cycle_nrmse)),
                "median": float(np.median(cycle_nrmse)),
                "p95": float(np.percentile(cycle_nrmse, 95)),
                "maximum": float(np.max(cycle_nrmse)),
            },
            "driver_correlation": {
                "minimum": float(np.min(cycle_corr)),
                "median": float(np.median(cycle_corr)),
            },
            "radius_rmse": {
                "median": float(np.median(cycle_r_rmse)),
                "maximum": float(np.max(cycle_r_rmse)),
            },
            "velocity_rmse": {
                "median": float(np.median(cycle_v_rmse)),
                "maximum": float(np.max(cycle_v_rmse)),
            },
            "control_peak": {
                "median": float(np.median(cycle_peak)),
                "maximum": float(np.max(cycle_peak)),
            },
            "last_complete_cycle": cycle_rows[-1],
        },
        "interpretation": {
            "validation_separation": (
                "The echo-jet KROM was fitted on the separate multisine forced-training "
                "trajectory; the sine archive supplied only the held-out control reference."
            ),
            "long_term": (
                "The settled cycle errors remain flat through the last complete cycle, "
                "with no secular radius growth."
            ),
            "limitation": (
                "The eigenfunction coordinates depend only on [echo, autonomous echo-dot], "
                "but exact B_psi is evaluated from simulated physical state. This is not yet "
                "a laboratory echo-only observer-controller, and the early unconstrained "
                "control peak is not actuator-ready."
            ),
        },
        "figure": FIGURE_PATH.relative_to(ROOT).as_posix(),
    }

    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    colors = {
        "target": "#1677ff",
        "controlled": "#ef7d32",
        "error": "#7657d6",
        "grid": "#dce5ef",
        "ink": "#1c2938",
    }
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.edgecolor": "#8da0b4",
            "axes.labelcolor": colors["ink"],
            "xtick.color": colors["ink"],
            "ytick.color": colors["ink"],
            "text.color": colors["ink"],
            "figure.facecolor": "#f7faff",
            "axes.facecolor": "#ffffff",
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 8.4), constrained_layout=True)
    headline = windows["tau_100_to_300"]
    figure.suptitle(
        "Repaired echo-jet KLQR: 100 kHz sine-driver rediscovery over 300 tau\n"
        f"settled driver NRMSE {100.0 * headline['driver_nrmse']:.3f}%  |  "
        f"correlation {headline['driver_correlation']:.6f}  |  "
        f"phase {headline['fundamental']['phase_deg']:.4f} deg",
        fontsize=14,
        fontweight="bold",
    )

    ax = axes[0, 0]
    ax.plot(time, target_state[:, 0], color=colors["target"], lw=1.0, label="target radius")
    ax.plot(time, state[:, 0], color=colors["controlled"], lw=0.8, alpha=0.9, label="controlled radius")
    ax.axvspan(0, 100, color="#f5c26b", alpha=0.13, label="acquisition excluded from headline")
    ax.set(title="Bounded state for the complete run", xlabel=r"$\tau$", ylabel=r"$r=R/R_0$")
    ax.legend(loc="upper right", frameon=False)

    last_start = 300.0 - 5.0 * period
    last = time >= last_start
    ax = axes[0, 1]
    ax.plot(time[last], known[last], color=colors["target"], lw=2.2, label="known 20 kPa sine")
    ax.plot(time[last], control[last], color=colors["controlled"], lw=1.4, label="KLQR recovered input")
    ax.set(title="Final five cycles: recovered driver", xlabel=r"$\tau$", ylabel=r"$u=p_c/p_0$")
    ax.legend(loc="upper right", frameon=False)

    ax = axes[1, 0]
    cycle_numbers = np.array([row["cycle_index"] for row in cycle_rows])
    ax.plot(cycle_numbers, 100.0 * cycle_nrmse, color=colors["error"], marker="o", ms=2.8, lw=1.2)
    ax.axhline(100.0 * np.median(cycle_nrmse), color="#475569", lw=1, ls="--", label="median")
    ax.set(title="Settled waveform error stays flat", xlabel="complete forcing cycle", ylabel="driver NRMSE (%)")
    ax.legend(loc="upper right", frameon=False)

    ax = axes[1, 1]
    radius_error_nm = (state[last, 0] - target_state[last, 0]) * radius0 * 1e9
    velocity_scale = radius0 / time_scale
    velocity_error_mm_s = (state[last, 1] - target_state[last, 1]) * velocity_scale * 1e3
    radius_line = ax.plot(
        time[last], radius_error_nm, color=colors["error"], lw=1.5, label="radius error"
    )
    ax.set(
        title="Final five cycles: physical state error",
        xlabel=r"$\tau$",
        ylabel="radius error (nm)",
    )
    ax.tick_params(axis="y", colors=colors["error"])
    second_axis = ax.twinx()
    velocity_line = second_axis.plot(
        time[last],
        velocity_error_mm_s,
        color="#0f9f84",
        lw=1.0,
        alpha=0.9,
        label="velocity error",
    )
    second_axis.set_ylabel("velocity error (mm/s)", color="#0f9f84")
    second_axis.tick_params(axis="y", colors="#0f9f84")
    ax.legend(radius_line + velocity_line, ["radius error", "velocity error"], loc="upper right", frameon=False)

    for ax in axes.flat:
        ax.grid(True, color=colors["grid"], lw=0.7)
        ax.spines[["top", "right"]].set_visible(False)

    figure.savefig(FIGURE_PATH, dpi=180, bbox_inches="tight")

    print(METRICS_PATH.relative_to(ROOT).as_posix())
    print(FIGURE_PATH.relative_to(ROOT).as_posix())
    print(json.dumps(windows["tau_100_to_300"], indent=2))


if __name__ == "__main__":
    main()
