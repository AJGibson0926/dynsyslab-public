"""Legacy reproduction only: this archived study uses the pre-correction RPE echo.
For new studies use the registered far_field_acoustic_echo and a freshly fitted basis.

Extract causal Hankel/HAVOK coordinates from one unforced RPE echo record.

The delay basis is fitted from the scalar far-field echo only.  Radius and
radial velocity are loaded after the fit and are retained solely as hidden
truth for Viewer comparison panels.
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess

import numpy as np

from plot_rpe_echo_jet_pipeline_evidence import figure, reference, trace


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "archives/solutions/rayleigh-plesset/"
    "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz"
)
OUTPUT_ROOT = ROOT / "viewer_outputs"
ARTIFACT = OUTPUT_ROOT / "rpe-unforced-echo-hankel-latents.npz"
METRICS = OUTPUT_ROOT / "rpe-unforced-echo-hankel-latents.metrics.json"
PLOT_NAME = "rpe-unforced-echo-hankel-latents"
ENERGY_PLOT_NAME = "rpe-unforced-echo-hankel-energy"
HISTORY_PLOT_NAME = "rpe-unforced-echo-hankel-histories"
PHASE_PLOT_NAME = "rpe-unforced-echo-hankel-phase-spaces"

DECIMATION = 5
DELAY_SAMPLES = 256
LATENT_COUNT = 6


def viewer_cli() -> Path:
    """Locate the compiled Viewer renderer or fail before generating plots."""
    built = ROOT / "build/native-suite-fast/bin/DynSysLabViewerCli.exe"
    if built.is_file():
        return built
    raise RuntimeError(f"built Viewer CLI is unavailable: {built}")


def render_viewer(
    source: Path,
    figures: list[dict],
    *,
    name: str,
    label: str,
    columns: int,
    width: int,
    height: int,
) -> None:
    """Write a local Viewer bundle, render its PNG, and verify PlotCanvas provenance."""
    bundle = OUTPUT_ROOT / f"{name}-viewer"
    payload_path = bundle / "payload.json"
    catalog_path = bundle / "catalog.json"
    output_path = OUTPUT_ROOT / f"{name}.png"
    report_path = OUTPUT_ROOT / f"{name}.render.json"
    payload = {
        "schema": "dynsyslab.viewer.payload/v1",
        "source": str(source),
        "view": "data_overview",
        "tile": {
            "name": label,
            "columns": columns,
        },
        "page_labels": ["Unforced echo-only latent extraction"],
        "warnings": [],
        "layout": "",
        "complex_display": "",
        "selected_modes": [],
        "figures": figures,
    }
    catalog = {
        "schema": "dynsyslab.viewer.catalog/v1",
        "items": [
            {
                "id": name,
                "label": label,
                "kind": "trajectory",
                "source": payload_path.name,
                "views": ["data_overview"],
            }
        ],
    }
    (bundle / "archives").mkdir(parents=True, exist_ok=True)
    (bundle / "cache").mkdir(parents=True, exist_ok=True)
    payload_path.write_text(
        json.dumps(payload, separators=(",", ":")), encoding="utf-8"
    )
    catalog_path.write_text(json.dumps(catalog, indent=2), encoding="utf-8")
    completed = subprocess.run(
        [
            str(viewer_cli()),
            "render",
            "--id",
            name,
            "--view",
            "data_overview",
            "--project",
            str(bundle),
            "--archives",
            str(bundle / "archives"),
            "--cache",
            str(bundle / "cache"),
            "--catalog",
            str(catalog_path),
            "--output",
            str(output_path),
            "--width",
            str(width),
            "--height",
            str(height),
            "--page",
            "0",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())
    report = json.loads(completed.stdout)
    if report.get("renderer") != "Viewer PlotCanvas":
        raise RuntimeError("plot did not use Viewer PlotCanvas")
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(output_path)


def load_source(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load physical trajectory samples and metadata from the local archive."""
    with np.load(path, allow_pickle=False) as archive:
        time = np.asarray(archive["time"], dtype=float)
        physical_state = np.asarray(archive["values"], dtype=float)
        metadata = json.loads(str(archive["metadata_json"]))
    return time, physical_state, metadata


def dimensionless_history(
    time: np.ndarray, physical_state: np.ndarray, metadata: dict
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Convert saved physical time and states with the archive radius and time scales."""
    parameters = metadata["system_parameters"]
    nondimensionalization = metadata["nondimensionalization"]
    radius0 = float(parameters["equilibrium_radius_m"])
    time_scale = float(nondimensionalization["time_scale_s"])
    state = np.column_stack(
        (
            physical_state[:, 0] / radius0,
            physical_state[:, 1] * time_scale / radius0,
        )
    )
    return time / time_scale, state, radius0, time_scale


def far_field_echo(
    state: np.ndarray, metadata: dict, radius0: float
) -> np.ndarray:
    """Evaluate the archived scalar echo measurement convention."""
    scales = metadata["nondimensionalization"]
    radius = state[:, 0]
    velocity = state[:, 1]
    momentum = (
        float(scales["gas_pressure"]) * radius**-4
        - 1.5 * velocity**2
        + float(scales["vapor_pressure"])
        - float(scales["ambient_pressure"])
        - float(scales["capillarity"]) / radius
        - float(scales["viscosity"]) * velocity / radius
    )
    return (
        radius0
        * float(scales["pressure_scale_pa"])
        * (2.0 * radius * velocity**2 + momentum)
    )


def fit_echo_hankel(
    echo: np.ndarray, delay_samples: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Fit a Hankel POD from one scalar channel and return HAVOK coordinates."""
    if echo.ndim != 1:
        raise ValueError("echo-only fitting requires exactly one scalar channel")
    if echo.size <= delay_samples:
        raise ValueError("the echo record is shorter than the delay window")
    mean = float(np.mean(echo))
    scale = float(np.std(echo))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        raise ValueError("the echo record has no usable variation")
    normalized_echo = (echo - mean) / scale
    windows = np.lib.stride_tricks.sliding_window_view(
        normalized_echo, delay_samples
    )
    hankel = np.ascontiguousarray(windows.T)
    left, singular, right_transpose = np.linalg.svd(
        hankel, full_matrices=False
    )

    # Remove the arbitrary SVD sign so regenerated coordinates are stable.
    for mode in range(singular.size):
        pivot = int(np.argmax(np.abs(left[:, mode])))
        if left[pivot, mode] < 0.0:
            left[:, mode] *= -1.0
            right_transpose[mode, :] *= -1.0

    energy = singular**2
    energy_fraction = energy / np.sum(energy)
    latent = right_transpose[:LATENT_COUNT, :].T
    projected = (left[:, :LATENT_COUNT].T @ hankel).T / singular[:LATENT_COUNT]
    if not np.allclose(latent, projected, rtol=2e-11, atol=2e-12):
        raise RuntimeError("online delay projection and SVD coordinates differ")
    return left, singular, energy_fraction, latent, mean, scale


def energy_rank(energy_fraction: np.ndarray, threshold: float) -> int:
    """Return the first cumulative-energy rank reaching the requested fraction."""
    return int(np.searchsorted(np.cumsum(energy_fraction), threshold) + 1)


def stem_traces(
    modes: np.ndarray, energy_fraction: np.ndarray
) -> list[dict]:
    """Represent a stem chart with ordinary Viewer line and marker traces."""
    positive = energy_fraction[energy_fraction > 0.0]
    floor = max(float(np.min(positive)) * 0.55, 1e-32)
    traces: list[dict] = []
    for mode, energy in zip(modes, energy_fraction, strict=True):
        shown = max(float(energy), floor)
        traces.append(
            trace(
                "",
                np.asarray([mode, mode]),
                np.asarray([floor, shown]),
                color="#0ea5e9",
            )
        )
        traces.append(
            trace(
                "",
                np.asarray([mode]),
                np.asarray([shown]),
                color="#0ea5e9",
                marker=True,
            )
        )
    return traces


def phase_figure(
    title: str,
    x_label: str,
    y_label: str,
    x: np.ndarray,
    y: np.ndarray,
    time: np.ndarray,
    *,
    references: list[dict] | None = None,
) -> dict:
    """Compose a time-colored square phase portrait with optional reference markers."""
    result = figure(
        "phase",
        title,
        x_label,
        y_label,
        [
            trace(
                "trajectory",
                x,
                y,
                encoding="time",
                color_values=time,
            )
        ],
        color_label="tau",
        references=references,
    )
    result["square"] = True
    return result


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

    tau, hidden_state, radius0, time_scale = dimensionless_history(
        physical_time, physical_state, source_metadata
    )
    echo = far_field_echo(hidden_state, source_metadata, radius0)
    echo_dot = np.gradient(echo, tau, edge_order=2)

    sample = np.arange(0, echo.size, DECIMATION)
    sampled_tau = tau[sample]
    sampled_echo = echo[sample]
    left, singular, energy_fraction, latent, echo_mean, echo_scale = (
        fit_echo_hankel(sampled_echo, DELAY_SAMPLES)
    )
    latent_time = sampled_tau[DELAY_SAMPLES - 1 :]
    cumulative_energy = np.cumsum(energy_fraction)
    aligned_sample = sample[DELAY_SAMPLES - 1 :]

    analysis_metadata = {
        "schema": "dynsyslab.echo-hankel-latent-analysis/v1",
        "label": "RPE free decay — echo-only Hankel coordinates",
        "source": str(SOURCE),
        "source_forcing_type": "free",
        "training_policy": "unforced-only",
        "training_input_channels": ["far_field_acoustic_echo"],
        "comparison_only_channels": ["r", "dr/dtau"],
        "comparison_channels_used_for_fit": False,
        "decimation": DECIMATION,
        "sample_interval_tau": float(
            np.median(np.diff(sampled_tau))
        ),
        "delay_samples": DELAY_SAMPLES,
        "delay_window_tau": float(
            sampled_tau[DELAY_SAMPLES - 1] - sampled_tau[0]
        ),
        "normalization": "global echo mean and standard deviation",
        "latent_coordinate_convention": (
            "xi_j(k)=sigma_j^-1 u_j^T dbar_k, equivalently row j of V*"
        ),
        "online_evaluation": "causal ring buffer ending at the current sample",
        "radius0_m": radius0,
        "time_scale_s": time_scale,
    }
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        ARTIFACT,
        time=tau,
        radius=hidden_state[:, 0],
        radial_velocity=hidden_state[:, 1],
        echo=echo,
        echo_dot_diagnostic=echo_dot,
        embedding_time=latent_time,
        latent_coordinates=latent,
        singular_values=singular,
        energy_fraction=energy_fraction,
        cumulative_energy=cumulative_energy,
        delay_modes=left,
        echo_mean=np.asarray(echo_mean),
        echo_scale=np.asarray(echo_scale),
        aligned_source_indices=aligned_sample,
        metadata_json=np.asarray(json.dumps(analysis_metadata, separators=(",", ":"))),
    )

    positive_energy = energy_fraction[energy_fraction > 0.0]
    entropy_dimension = float(
        np.exp(-np.sum(positive_energy * np.log(positive_energy)))
    )
    metrics = {
        **analysis_metadata,
        "source_samples": int(echo.size),
        "fitted_echo_samples": int(sampled_echo.size),
        "hankel_shape": [DELAY_SAMPLES, int(latent_time.size)],
        "latent_coordinates_shown": LATENT_COUNT,
        "energy_rank": {
            "90_percent": energy_rank(energy_fraction, 0.90),
            "95_percent": energy_rank(energy_fraction, 0.95),
            "99_percent": energy_rank(energy_fraction, 0.99),
            "99_9_percent": energy_rank(energy_fraction, 0.999),
        },
        "spectral_entropy_dimension": entropy_dimension,
        "first_10_energy_fraction": energy_fraction[:10].tolist(),
        "artifact": str(ARTIFACT),
    }
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    spectrum = figure(
        "stem",
        "Echo Hankel energy spectrum",
        "mode j",
        "sigma_j^2 / sum sigma_i^2",
        stem_traces(np.arange(1, singular.size + 1), energy_fraction),
    )
    spectrum["value_scale"] = "log"
    spectrum["square"] = False

    display = np.arange(0, tau.size, 4)
    colors = ["#0ea5e9", "#f59e0b", "#8b5cf6", "#10b981", "#ef4444", "#64748b"]
    figures = [
        spectrum,
        figure(
            "time",
            "Unforced radius history (comparison only)",
            "tau",
            "r",
            [trace("r", tau[display], hidden_state[display, 0], color="#0ea5e9")],
        ),
        figure(
            "time",
            "Unforced radial-velocity history (comparison only)",
            "tau",
            "dr/dtau",
            [
                trace(
                    "dr/dtau",
                    tau[display],
                    hidden_state[display, 1],
                    color="#f59e0b",
                )
            ],
        ),
        figure(
            "observable_time",
            "Echo history supplied to the delay embedding",
            "tau",
            "e",
            [trace("echo", tau[display], echo[display], color="#10b981")],
        ),
        figure(
            "time",
            "Echo-history latent variables 1–3",
            "tau",
            "xi_j",
            [
                trace(f"xi_{{{mode + 1}}}", latent_time, latent[:, mode], color=colors[mode])
                for mode in range(3)
            ],
        ),
        figure(
            "time",
            "Echo-history latent variables 4–6",
            "tau",
            "xi_j",
            [
                trace(f"xi_{{{mode + 1}}}", latent_time, latent[:, mode], color=colors[mode])
                for mode in range(3, LATENT_COUNT)
            ],
        ),
    ]
    phase_display = np.arange(0, tau.size, 4)
    latent_display = np.arange(0, latent_time.size, 1)
    latent_3d_trace = trace(
        "trajectory",
        latent[latent_display, 0],
        latent[latent_display, 1],
        encoding="time",
        color_values=latent_time[latent_display],
    )
    latent_3d_trace["z"] = latent[latent_display, 2].tolist()
    latent_3d = figure(
        "three_d",
        "First three echo-history latent variables",
        "xi_1",
        "xi_2",
        [latent_3d_trace],
        color_label="tau",
    )
    latent_3d["z_label"] = "xi_3"
    latent_3d["square"] = True
    phase_figures = [
        phase_figure(
            "Unforced RPE state phase space (comparison only)",
            "r",
            "dr/dtau",
            hidden_state[phase_display, 0],
            hidden_state[phase_display, 1],
            tau[phase_display],
            references=[
                reference(
                    float(hidden_state[0, 0]),
                    float(hidden_state[0, 1]),
                    "IC",
                    False,
                ),
                reference(1.0, 0.0, "EC", True),
            ],
        ),
        phase_figure(
            "Echo phase space (numerical derivative diagnostic)",
            "e",
            "de/dtau",
            echo[phase_display],
            echo_dot[phase_display],
            tau[phase_display],
            references=[
                reference(float(echo[0]), float(echo_dot[0]), "IC", False),
                reference(0.0, 0.0, "EC", True),
            ],
        ),
        phase_figure(
            "Latent phase space 1–2",
            "xi_1",
            "xi_2",
            latent[latent_display, 0],
            latent[latent_display, 1],
            latent_time[latent_display],
        ),
        latent_3d,
        phase_figure(
            "Latent phase space 3–4",
            "xi_3",
            "xi_4",
            latent[latent_display, 2],
            latent[latent_display, 3],
            latent_time[latent_display],
        ),
        phase_figure(
            "Latent phase space 5–6",
            "xi_5",
            "xi_6",
            latent[latent_display, 4],
            latent[latent_display, 5],
            latent_time[latent_display],
        ),
    ]
    history_figures = figures[1:]
    history_figures[-1]["grid_span"] = 2
    render_viewer(
        ARTIFACT,
        [spectrum],
        name=ENERGY_PLOT_NAME,
        label="RPE free decay — echo Hankel energy",
        columns=1,
        width=1600,
        height=1000,
    )
    render_viewer(
        ARTIFACT,
        history_figures,
        name=HISTORY_PLOT_NAME,
        label="RPE free decay — echo-only latent histories",
        columns=2,
        width=1900,
        height=1850,
    )
    render_viewer(
        ARTIFACT,
        phase_figures,
        name=PHASE_PLOT_NAME,
        label="RPE free decay — echo-only latent phase spaces",
        columns=2,
        width=1900,
        height=1900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
