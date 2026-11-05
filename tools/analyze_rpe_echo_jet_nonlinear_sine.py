"""Reproduce the nonlinear RPE echo-jet KLQR stress-test evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOLUTION = Path("archives/solutions/rpe")
PATHS = {
    "training": SOLUTION / "rpe-echo-study-state-basis-baseline__20260911T104742Z.npz",
    "krom": Path("archives/derived/koopman-with-sindy/rpe-echo-study-state-basis-baseline__20260911T104742Z__koopman-with-sindy__20260914T165124Z.npz"),
    "target_20": SOLUTION / "rpe-periodic-100khz-20kpa-phase-portrait__20260911T214426Z.npz",
    "echo_20": SOLUTION / "rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300__20260914T171839Z.npz",
    "target_28": SOLUTION / "rpe-periodic-100khz-28kpa-nonlinear-tau300__20260914T184520Z.npz",
    "echo_28": SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-28kpa-tau300__20260914T184551Z.npz",
    "target_40": SOLUTION / "rpe-periodic-100khz-40kpa-settled-target__20260914T183110Z.npz",
    "echo_40_default": SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-settled-tau300__20260914T183133Z.npz",
    "echo_40_tuned": SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3-tau300__20260914T192234Z.npz",
    "echo_40_branch_hopping": SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r3e-3-tau300__20260914T190937Z.npz",
    "echo_40_soft_comparator": SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r1e-2-tau300__20260914T192234Z.npz",
    "state_40": SOLUTION / "rpe-state-klqr-rediscover-100khz-40kpa-settled-tau300__20260914T183611Z.npz",
}

GAIN_PROBES = {
    0.0001: PATHS["echo_40_default"],
    0.0003: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r3e-4__20260914T185707Z.npz",
    0.001: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r1e-3__20260914T185707Z.npz",
    0.002: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r2e-3__20260914T190321Z.npz",
    0.003: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r3e-3__20260914T185707Z.npz",
    0.004: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r4e-3__20260914T190321Z.npz",
    0.006: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3__20260914T190321Z.npz",
    0.01: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r1e-2__20260914T190027Z.npz",
    0.03: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r3e-2__20260914T190027Z.npz",
    0.1: SOLUTION / "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r1e-1__20260914T190028Z.npz",
}

OUTPUT_JSON = ROOT / "viewer_outputs/rpe-echo-jet2-nonlinear-sine-stress.json"
OUTPUT_PNG = ROOT / "viewer_outputs/rpe-echo-jet2-nonlinear-sine-stress.png"


def load(relative: Path) -> tuple[dict[str, np.ndarray], dict]:
    with np.load(ROOT / relative, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files if key != "metadata_json"}
        metadata = json.loads(str(archive["metadata_json"]))
    return arrays, metadata


def digest(relative: Path) -> str:
    value = hashlib.sha256()
    with (ROOT / relative).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def nrmse(actual: np.ndarray, reference: np.ndarray) -> float:
    return rms(actual - reference) / rms(reference - np.mean(reference))


def tone_fit(time: np.ndarray, values: np.ndarray, omega: float) -> dict[str, float]:
    design = np.column_stack((np.ones_like(time), np.cos(omega * time), np.sin(omega * time)))
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    fitted = design @ coefficients
    cosine, sine = coefficients[1:]
    total = np.sum(np.square(values - np.mean(values)))
    return {
        "mean": float(coefficients[0]),
        "amplitude": float(np.hypot(cosine, sine)),
        "phase_deg": float(np.degrees(np.arctan2(-sine, cosine))),
        "r_squared": float(1.0 - np.sum(np.square(values - fitted)) / total),
    }


def harmonic_content(time: np.ndarray, values: np.ndarray, omega: float, count: int = 5) -> dict:
    columns = [np.ones_like(time)]
    for order in range(1, count + 1):
        columns.extend((np.cos(order * omega * time), np.sin(order * omega * time)))
    coefficients, *_ = np.linalg.lstsq(np.column_stack(columns), values, rcond=None)
    amplitudes = [float(np.hypot(coefficients[2 * order - 1], coefficients[2 * order])) for order in range(1, count + 1)]
    return {
        "amplitudes": amplitudes,
        "second_to_fundamental": amplitudes[1] / amplitudes[0],
        "higher_rms_to_fundamental": float(np.linalg.norm(amplitudes[1:]) / amplitudes[0]),
    }


def echo_jet(states: np.ndarray, groups: dict[str, float], radius0: float) -> np.ndarray:
    r = states[:, 0]
    v = states[:, 1]
    gas = groups["gas_pressure_ratio"]
    vapor = groups["vapor_pressure_ratio"]
    ambient = groups["ambient_pressure_ratio"]
    capillary = groups["capillarity"]
    viscous = groups["viscosity"]
    pressure_scale = groups["pressure_scale_pa"]
    momentum = gas * r**-4 - 1.5 * v**2 + vapor - ambient - capillary / r - viscous * v / r
    acceleration = momentum / r
    scale = radius0 * pressure_scale
    echo = scale * (2.0 * r * v**2 + momentum)
    echo_dr = scale * (2.0 * v**2 - 4.0 * gas * r**-5 + capillary / r**2 + viscous * v / r**2)
    echo_dv = scale * (4.0 * r * v - 3.0 * v - viscous / r)
    return np.column_stack((echo, echo_dr * v + echo_dv * acceleration))


def compare(target_path: Path, result_path: Path, start: float, stop: float | None = None) -> dict:
    target, _ = load(target_path)
    result, result_metadata = load(result_path)
    time = np.asarray(result["time"], dtype=float)
    target_time = np.asarray(target["time"], dtype=float)[: len(time)]
    if not np.allclose(time, target_time, rtol=0.0, atol=1e-12):
        raise RuntimeError(f"time grids differ for {result_path}")
    stop = float(time[-1]) if stop is None else stop
    mask = (time >= start) & (time <= stop)
    state = np.asarray(result["values"], dtype=float)
    target_state = np.asarray(result["control_target_values"], dtype=float)
    control = np.asarray(result["control_values"], dtype=float)[:, 0]
    known = np.asarray(target["control_values"], dtype=float)[: len(time), 0]
    frequency_hz = 100000.0
    omega = 2.0 * np.pi * frequency_hz * result_metadata["time_scale_to_seconds"]
    fit = tone_fit(time[mask], control[mask], omega)
    known_fit = tone_fit(time[mask], known[mask], omega)
    row = {
        "tau": [start, stop],
        "driver_nrmse": nrmse(control[mask], known[mask]),
        "driver_correlation": float(np.corrcoef(control[mask], known[mask])[0, 1]),
        "driver_peak_absolute": float(np.max(np.abs(control[mask]))),
        "fundamental": fit,
        "known_fundamental": known_fit,
        "fundamental_amplitude_ratio": fit["amplitude"] / known_fit["amplitude"],
        "fundamental_phase_error_deg": float((fit["phase_deg"] - known_fit["phase_deg"] + 180.0) % 360.0 - 180.0),
        "state_rmse": [rms(state[mask, index] - target_state[mask, index]) for index in range(2)],
        "state_range": [[float(np.min(state[mask, index])), float(np.max(state[mask, index]))] for index in range(2)],
    }
    if "control_observable_target_values" in result:
        target_observed = np.asarray(result["control_observable_target_values"], dtype=float)
        groups = result_metadata["dimensionless_groups"]
        radius0 = result_metadata["resolved_system_parameters"]["equilibrium_radius_m"]
        observed = echo_jet(state, groups, radius0)
        row["echo_jet_nrmse"] = [nrmse(observed[mask, index], target_observed[mask, index]) for index in range(2)]

    period = 2.0 * np.pi / omega
    cycle_rows = []
    for cycle in range(int(np.ceil(start / period)), int(np.floor(stop / period))):
        cycle_mask = (time >= cycle * period) & (time < (cycle + 1) * period)
        if np.count_nonzero(cycle_mask) > 2:
            cycle_rows.append({
                "cycle": cycle,
                "midpoint_tau": (cycle + 0.5) * period,
                "driver_nrmse": nrmse(control[cycle_mask], known[cycle_mask]),
            })
    row["cycles"] = cycle_rows
    values = np.array([entry["driver_nrmse"] for entry in cycle_rows])
    row["cycle_driver_nrmse_summary"] = {
        "count": int(values.size),
        "minimum": float(np.min(values)),
        "median": float(np.median(values)),
        "maximum": float(np.max(values)),
        "last": float(values[-1]),
    }
    return row


def target_metrics(path: Path, start: float) -> dict:
    arrays, metadata = load(path)
    time = np.asarray(arrays["time"], dtype=float)
    state = np.asarray(arrays["values"], dtype=float)
    mask = time >= start
    omega = 2.0 * np.pi * 100000.0 * metadata["time_scale_to_seconds"]
    harmonics = harmonic_content(time[mask], state[mask, 0], omega)
    return {
        "tau": [start, float(time[-1])],
        "radius_range": [float(np.min(state[mask, 0])), float(np.max(state[mask, 0]))],
        "velocity_range": [float(np.min(state[mask, 1])), float(np.max(state[mask, 1]))],
        "radius_peak_to_peak": float(np.ptp(state[mask, 0])),
        "radius_harmonics": harmonics,
        "omega_dimensionless": omega,
        "period_dimensionless": 2.0 * np.pi / omega,
    }


def main() -> None:
    metrics_20 = compare(PATHS["target_20"], PATHS["echo_20"], 100.0, 300.0)
    metrics_28 = compare(PATHS["target_28"], PATHS["echo_28"], 100.0, 300.0)
    metrics_40_default = compare(PATHS["target_40"], PATHS["echo_40_default"], 100.0, 300.0)
    metrics_40_tuned = compare(PATHS["target_40"], PATHS["echo_40_tuned"], 100.0, 300.0)
    metrics_40_branch_hopping = compare(PATHS["target_40"], PATHS["echo_40_branch_hopping"], 100.0, 300.0)
    metrics_40_soft_comparator = compare(PATHS["target_40"], PATHS["echo_40_soft_comparator"], 100.0, 300.0)
    metrics_40_state = compare(PATHS["target_40"], PATHS["state_40"], 100.0, 300.0)
    targets = {
        "20_kpa": target_metrics(PATHS["target_20"], 100.0),
        "28_kpa": target_metrics(PATHS["target_28"], 100.0),
        "40_kpa": target_metrics(PATHS["target_40"], 100.0),
    }

    gain_sweep = []
    for penalty, path in GAIN_PROBES.items():
        start = 100.0 if penalty == 0.0001 else 60.0
        stop = 300.0 if penalty == 0.0001 else 100.0
        row = compare(PATHS["target_40"], path, start, stop)
        gain_sweep.append({"R": penalty, "Q": 5000.0, **{key: value for key, value in row.items() if key != "cycles"}})

    lineage = {}
    for name, path in PATHS.items():
        lineage[name] = {"path": path.as_posix(), "sha256": digest(path), "bytes": (ROOT / path).stat().st_size}

    payload = {
        "schema": "dynsyslab.rpe-echo-jet-nonlinear-sine-stress/v1",
        "date": "2026-09-14",
        "conclusion": (
            "The unchanged echo/echo-dot eigenfunction robustly rediscovers a 28 kPa sine driver. "
            "At 40 kPa it requires a much softer KLQR gain (R=0.006); that controller remains "
            "on the correct branch over all 59 complete cycles after tau=100, while the original "
            "R=0.0001 controller fails and the state-coordinate control supplies the representation baseline."
        ),
        "targets": targets,
        "controls": {
            "echo_20_default": metrics_20,
            "echo_28_default": metrics_28,
            "echo_40_default": metrics_40_default,
            "echo_40_tuned": metrics_40_tuned,
            "echo_40_branch_hopping_R_0p003": metrics_40_branch_hopping,
            "echo_40_soft_comparator_R_0p01": metrics_40_soft_comparator,
            "state_40_default": metrics_40_state,
        },
        "gain_sweep_40_kpa": gain_sweep,
        "lineage": lineage,
    }
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2), constrained_layout=True)
    fig.patch.set_facecolor("#f7fafc")
    for axis in axes.flat:
        axis.set_facecolor("#ffffff")
        axis.grid(True, color="#dbe4ee", linewidth=0.7, alpha=0.8)

    colors = {20: "#4f7cac", 28: "#2a9d8f", 40: "#e76f51"}
    for amplitude, key in ((20, "target_20"), (28, "target_28"), (40, "target_40")):
        arrays, metadata = load(PATHS[key])
        time = arrays["time"]
        omega = 2.0 * np.pi * 100000.0 * metadata["time_scale_to_seconds"]
        period = 2.0 * np.pi / omega
        mask = time >= time[-1] - period
        phase = ((time[mask] - time[mask][0]) / period) % 1.0
        axes[0, 0].plot(phase, arrays["values"][mask, 0], color=colors[amplitude], linewidth=2.0, label=f"{amplitude} kPa")
    axes[0, 0].set(title="Nonlinear target radius over one forcing cycle", xlabel="cycle phase", ylabel=r"$r/R_0$")
    axes[0, 0].legend(frameon=False)

    target_40, _ = load(PATHS["target_40"])
    tuned_40, _ = load(PATHS["echo_40_tuned"])
    default_40, _ = load(PATHS["echo_40_default"])
    state_40, _ = load(PATHS["state_40"])
    time = tuned_40["time"]
    mask = time >= 290.0
    axes[0, 1].plot(time[mask], target_40["control_values"][mask, 0], color="#172b4d", linewidth=2.1, label="known sine")
    axes[0, 1].plot(time[mask], tuned_40["control_values"][mask, 0], color="#2a9d8f", linewidth=1.5, label="echo jet, tuned")
    axes[0, 1].plot(time[mask], state_40["control_values"][mask, 0], color="#4f7cac", linewidth=1.2, linestyle="--", label="state baseline")
    axes[0, 1].plot(time[mask], default_40["control_values"][mask, 0], color="#e76f51", linewidth=1.0, alpha=0.8, label="echo jet, old gain")
    axes[0, 1].set(title="40 kPa driver recovery, late time", xlabel=r"dimensionless time $\tau$", ylabel=r"$u=p_c/100\,\mathrm{kPa}$")
    axes[0, 1].legend(frameon=False, ncol=2)

    penalties = np.array([row["R"] for row in gain_sweep])
    errors = 100.0 * np.array([row["driver_nrmse"] for row in gain_sweep])
    axes[1, 0].semilogx(penalties, errors, "o-", color="#6b5ca5", linewidth=1.5)
    axes[1, 0].axhline(20.0, color="#2a9d8f", linestyle="--", linewidth=1.0, label="20% screen")
    axes[1, 0].scatter([0.006], [100.0 * gain_sweep[[row["R"] for row in gain_sweep].index(0.006)]["driver_nrmse"]], color="#2a9d8f", s=60, zorder=5, label="selected for long run")
    axes[1, 0].set(title="40 kPa echo-jet gain sweep", xlabel="input penalty R (Q = 5000)", ylabel="driver NRMSE (%)")
    axes[1, 0].legend(frameon=False)

    for label, row, color in (
        ("28 kPa echo jet", metrics_28, "#2a9d8f"),
        ("40 kPa echo jet, tuned", metrics_40_tuned, "#e9c46a"),
        ("40 kPa echo jet, old gain", metrics_40_default, "#e76f51"),
        ("40 kPa state baseline", metrics_40_state, "#4f7cac"),
    ):
        midpoints = [entry["midpoint_tau"] for entry in row["cycles"]]
        cycle_errors = 100.0 * np.array([entry["driver_nrmse"] for entry in row["cycles"]])
        axes[1, 1].plot(midpoints, cycle_errors, color=color, linewidth=1.4, label=label)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(title="Long-term cycle-by-cycle driver error", xlabel=r"dimensionless time $\tau$", ylabel="cycle driver NRMSE (%)")
    axes[1, 1].legend(frameon=False)

    fig.suptitle("RPE echo/echo-dot KLQR: nonlinear sine-driver stress test", fontsize=15, color="#172b4d")
    fig.savefig(OUTPUT_PNG, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(OUTPUT_JSON)
    print(OUTPUT_PNG)


if __name__ == "__main__":
    main()
