"""Render the training and KROM evidence used by the echo-jet KLQR report.

The output figures use the installed Viewer PlotCanvas so the report's
identification evidence follows the same presentation conventions as the
final controlled trajectories.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import zipfile

import numpy as np

from plot_rpe_echo_jet_success_detail import echo_jet, load


ROOT = Path(__file__).resolve().parents[1]
FREE_PATH = ROOT / (
    "archives/solutions/rayleigh-plesset/"
    "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz"
)
BRANCH_PATH = ROOT / (
    "archives/solutions/rpe/"
    "rpe-echo-study-state-basis-baseline__20260911T104742Z.npz"
)
HDMD_PATH = ROOT / (
    "archives/derived/koopman-hankel-dmd/"
    "rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__"
    "koopman-hankel-dmd__20260826T221322.756123Z__604a65f1.npz"
)
KROM_PATH = ROOT / (
    "archives/derived/koopman-with-sindy/"
    "rpe-echo-study-state-basis-baseline__20260911T104742Z__"
    "koopman-with-sindy__20260914T165124Z.npz"
)
OUTPUT_ROOT = ROOT / "viewer_outputs"


def numbers(array: np.ndarray) -> list[float]:
    return np.asarray(array, dtype=float).tolist()


def trace(
    label: str,
    x: np.ndarray,
    y: np.ndarray,
    *,
    role: str = "data",
    color: str = "",
    line_style: str = "solid",
    encoding: str = "",
    color_values: np.ndarray | None = None,
    marker: bool = False,
    marker_size: float = 3.0,
    width: float = 1.7,
    opacity: float = 1.0,
) -> dict:
    return {
        "label": label,
        "role": role,
        "color": color,
        "line_style": line_style,
        "encoding": encoding,
        "selection_id": -1,
        "selected": False,
        "marker": marker,
        "marker_size": marker_size,
        "width": width,
        "opacity": opacity,
        "heatmap": False,
        "x": numbers(x),
        "y": numbers(y),
        "z": [],
        "color_values": [] if color_values is None else numbers(color_values),
        "magnitudes": [],
        "z_shape": [],
    }


def figure(
    kind: str,
    title: str,
    x_label: str,
    y_label: str,
    traces: list[dict],
    *,
    color_label: str = "",
    references: list[dict] | None = None,
) -> dict:
    return {
        "kind": kind,
        "value_scale": "auto",
        "observable_identity": "",
        "title": title,
        "x_label": x_label,
        "y_label": y_label,
        "square": kind == "koopman_spectrum",
        "title_color": "",
        "dashboard_page": 0,
        "grid_span": 1,
        "row_span": 1,
        "coordinates": "dimensionless",
        "color_label": color_label,
        "traces": traces,
        "references": references or [],
    }


def reference(x: float, y: float, label: str, equilibrium: bool) -> dict:
    return {"x": x, "y": y, "label": label, "equilibrium": equilibrium}


def viewer_cli() -> Path:
    local = os.environ.get("LOCALAPPDATA")
    if not local:
        raise RuntimeError("LOCALAPPDATA is unavailable")
    executable = Path(local) / "Programs/DynSysLab/bin/DynSysLabViewerCli.exe"
    if not executable.is_file():
        raise RuntimeError(f"installed Viewer CLI is unavailable: {executable}")
    return executable


def render(
    *,
    name: str,
    label: str,
    page_label: str,
    source: Path,
    figures: list[dict],
    columns: int,
    height: int,
) -> None:
    bundle = OUTPUT_ROOT / f"{name}-viewer"
    payload_path = bundle / "payload.json"
    catalog_path = bundle / "catalog.json"
    output_path = OUTPUT_ROOT / f"{name}.png"
    report_path = OUTPUT_ROOT / f"{name}.render.json"
    catalog_id = name

    payload = {
        "schema": "dynsyslab.viewer.payload/v1",
        "source": str(source),
        "view": "data_overview",
        "tile": {"name": label, "columns": columns},
        "page_labels": [page_label],
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
                "id": catalog_id,
                "label": label,
                "kind": "trajectory",
                "source": payload_path.name,
                "views": ["data_overview"],
            }
        ],
    }
    viewer_archives = bundle / "archives"
    viewer_cache = bundle / "cache"
    viewer_archives.mkdir(parents=True, exist_ok=True)
    viewer_cache.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    catalog_path.write_text(json.dumps(catalog, indent=2), encoding="utf-8")

    completed = subprocess.run(
        [
            str(viewer_cli()),
            "render",
            "--id",
            catalog_id,
            "--view",
            "data_overview",
            "--project",
            str(bundle),
            "--archives",
            str(viewer_archives),
            "--cache",
            str(viewer_cache),
            "--catalog",
            str(catalog_path),
            "--output",
            str(output_path),
            "--width",
            "1900",
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
        raise RuntimeError(
            f"Viewer render failed for {name}:\n"
            + (completed.stderr.strip() or completed.stdout.strip())
        )
    render_report = json.loads(completed.stdout)
    if render_report.get("renderer") != "Viewer PlotCanvas":
        raise RuntimeError(f"{name} did not use Viewer PlotCanvas")
    report_path.write_text(json.dumps(render_report, indent=2), encoding="utf-8")
    print(output_path)


def free_training_figures(groups: dict, radius0: float) -> list[dict]:
    arrays, metadata = load(FREE_PATH)
    nondimensionalization = metadata["nondimensionalization"]
    time_scale = float(nondimensionalization["time_scale_s"])
    physical = np.asarray(arrays["values"], dtype=float)
    state = np.column_stack(
        (
            physical[:, 0] / radius0,
            physical[:, 1] * time_scale / radius0,
        )
    )
    time = np.asarray(arrays["time"], dtype=float) / time_scale
    observed = echo_jet(state, groups, radius0)
    sample = np.arange(0, time.size, 4)
    t = time[sample]
    x = state[sample]
    y = observed[sample]
    return [
        figure("time", "Radius across the free transient", "\\tau", "r", [trace("r", t, x[:, 0])]),
        figure(
            "time",
            "Radial velocity across the free transient",
            "\\tau",
            "\\frac{dr}{d\\tau}",
            [trace("radial velocity", t, x[:, 1])],
        ),
        figure(
            "phase",
            "Free-training state coverage",
            "r",
            "\\frac{dr}{d\\tau}",
            [trace("free transient", x[:, 0], x[:, 1], encoding="time", color_values=t)],
            color_label="\\tau",
            references=[
                reference(float(x[0, 0]), float(x[0, 1]), "IC", False),
                reference(1.0, 0.0, "EC", True),
            ],
        ),
        figure(
            "phase",
            "The same transient folds in echo coordinates",
            "e",
            "L_{f_0}e",
            [trace("echo jet", y[:, 0], y[:, 1], encoding="time", color_values=t)],
            color_label="\\tau",
            references=[
                reference(float(y[0, 0]), float(y[0, 1]), "IC", False),
                reference(0.0, 0.0, "EC", True),
            ],
        ),
    ]


def branch_training_figures(groups: dict, radius0: float) -> list[dict]:
    arrays, _ = load(BRANCH_PATH)
    time = np.asarray(arrays["time"], dtype=float)
    state = np.asarray(arrays["values"], dtype=float)
    target = np.asarray(arrays["control_target_values"], dtype=float)
    control = np.asarray(arrays["control_values"], dtype=float)[:, 0]
    observed = echo_jet(state, groups, radius0)
    target_observed = echo_jet(target, groups, radius0)
    sample = np.arange(0, time.size, 3)
    t = time[sample]
    x = state[sample]
    xt = target[sample]
    y = observed[sample]
    yt = target_observed[sample]
    target_style = {"role": "control_target", "color": "#ef4444", "line_style": "dash"}
    return [
        figure("time", "Branch-covering control", "\\tau", "u", [trace("control", t, control[sample], role="control")]),
        figure(
            "time",
            "Radius on the training branch",
            "\\tau",
            "r",
            [trace("r", t, x[:, 0]), trace("target", t, xt[:, 0], **target_style)],
        ),
        figure(
            "time",
            "Radial velocity on the training branch",
            "\\tau",
            "\\frac{dr}{d\\tau}",
            [trace("radial velocity", t, x[:, 1]), trace("target", t, xt[:, 1], **target_style)],
        ),
        figure(
            "phase",
            "State-space branch coverage",
            "r",
            "\\frac{dr}{d\\tau}",
            [
                trace("training trajectory", x[:, 0], x[:, 1], encoding="time", color_values=t),
                trace("reference", xt[:, 0], xt[:, 1], **target_style),
            ],
            color_label="\\tau",
            references=[
                reference(float(x[0, 0]), float(x[0, 1]), "IC", False),
                reference(1.0, 0.0, "EC", True),
            ],
        ),
        figure(
            "observable_time",
            "Echo on the training branch",
            "\\tau",
            "e",
            [trace("echo", t, y[:, 0]), trace("reference", t, yt[:, 0], **target_style)],
        ),
        figure(
            "phase",
            "Echo-jet branch used for the refit",
            "e",
            "L_{f_0}e",
            [
                trace("training trajectory", y[:, 0], y[:, 1], encoding="time", color_values=t),
                trace("reference", yt[:, 0], yt[:, 1], **target_style),
            ],
            color_label="\\tau",
            references=[
                reference(float(y[0, 0]), float(y[0, 1]), "IC", False),
                reference(0.0, 0.0, "EC", True),
            ],
        ),
    ]


def polynomial_eigenfunction(
    coordinates: np.ndarray,
    coordinate_velocities: np.ndarray,
    center: np.ndarray,
    scale: np.ndarray,
    exponents: np.ndarray,
    coefficients: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    eta = (coordinates - center) / scale
    powers = exponents.astype(int)
    theta = np.prod(eta[:, None, :] ** powers[None, :, :], axis=2)
    psi = theta @ coefficients
    basis_gradient = np.zeros((eta.shape[0], powers.shape[0], eta.shape[1]))
    for term, alpha in enumerate(powers):
        for coordinate, degree in enumerate(alpha):
            if not degree:
                continue
            derivative = np.full(eta.shape[0], float(degree))
            for factor, factor_degree in enumerate(alpha):
                derivative *= eta[:, factor] ** (
                    factor_degree - (1 if factor == coordinate else 0)
                )
            basis_gradient[:, term, coordinate] = derivative
    eta_dot = coordinate_velocities / scale
    psi_dot = np.einsum("ntj,nj,t->n", basis_gradient, eta_dot, coefficients)
    return psi, psi_dot


def krom_figures() -> list[dict]:
    with zipfile.ZipFile(HDMD_PATH) as archive:
        with archive.open("eigenvalues.npy") as stream:
            spectrum = np.asarray(np.load(stream), dtype=complex)
    with np.load(KROM_PATH, allow_pickle=False) as archive:
        eigenvalues = np.asarray(archive["discrete_eigenvalues"], dtype=complex)
        lambdas = np.asarray(archive["continuous_eigenvalues"], dtype=complex)
        coordinates = np.asarray(archive["training_coordinates"], dtype=float)
        coordinate_velocities = np.asarray(
            archive["training_coordinate_velocities"], dtype=float
        )
        time = np.asarray(archive["training_time"], dtype=float)
        center = np.asarray(archive["state_center"], dtype=float)
        scale = np.asarray(archive["state_scale"], dtype=float)
        exponents = np.asarray(archive["polynomial_exponents"], dtype=int)
        coefficients = np.asarray(archive["coefficients"], dtype=complex)[0]

    psi, actual_generator = polynomial_eigenfunction(
        coordinates,
        coordinate_velocities,
        center,
        scale,
        exponents,
        coefficients,
    )
    predicted_generator = lambdas[0] * psi
    circle = np.exp(1j * np.linspace(0.0, 2.0 * np.pi, 361))
    finite = spectrum[np.isfinite(spectrum)]
    spectrum_sample = finite[:: max(1, finite.size // 180)]
    sample = np.arange(0, time.size, 3)
    closure = np.arange(0, time.size, 5)
    lo = min(actual_generator.real.min(), predicted_generator.real.min())
    hi = max(actual_generator.real.max(), predicted_generator.real.max())
    identity = np.asarray([lo, hi])
    lo_i = min(actual_generator.imag.min(), predicted_generator.imag.min())
    hi_i = max(actual_generator.imag.max(), predicted_generator.imag.max())
    identity_i = np.asarray([lo_i, hi_i])
    return [
        figure(
            "koopman_spectrum",
            "Hankel-DMD candidates and the retained pair",
            "Re(\\mu)",
            "Im(\\mu)",
            [
                trace("unit circle", circle.real, circle.imag, color="#808793", width=1.0),
                *[
                    trace(
                        "candidate spectrum" if index == 0 else "",
                        np.asarray([value.real]),
                        np.asarray([value.imag]),
                        color="#8792a2",
                        marker=True,
                        marker_size=2.8,
                        opacity=0.62,
                    )
                    for index, value in enumerate(spectrum_sample)
                ],
                *[
                    trace(
                        "retained pair" if index == 0 else "",
                        np.asarray([value.real]),
                        np.asarray([value.imag]),
                        color="#ef4444",
                        marker=True,
                        marker_size=11.0,
                    )
                    for index, value in enumerate(eigenvalues)
                ],
            ],
        ),
        figure(
            "phase",
            "Normalized echo-jet chart",
            "(e-c_e)/s_e",
            "(L_{f_0}e-c_{e_1})/s_{e_1}",
            [
                trace(
                    "training branch",
                    ((coordinates - center) / scale)[sample, 0],
                    ((coordinates - center) / scale)[sample, 1],
                    encoding="time",
                    color_values=time[sample],
                )
            ],
            color_label="\\tau",
        ),
        figure(
            "phase",
            "Learned Koopman eigenfunction on the branch",
            "Re(\\psi)",
            "Im(\\psi)",
            [trace("eigenfunction", psi.real[sample], psi.imag[sample], encoding="time", color_values=time[sample])],
            color_label="\\tau",
        ),
        figure(
            "phase",
            "Generator closure: learned derivative versus \\lambda\\psi",
            "measured \\dot{\\psi}",
            "predicted \\lambda\\psi",
            [
                trace("real part", actual_generator.real[closure], predicted_generator.real[closure], color="#3478c8", marker=True, marker_size=2.5),
                trace("real identity", identity, identity, color="#3478c8", line_style="dash", width=1.0),
                trace("imaginary part", actual_generator.imag[closure], predicted_generator.imag[closure], color="#ef4444", marker=True, marker_size=2.5),
                trace("imaginary identity", identity_i, identity_i, color="#ef4444", line_style="dash", width=1.0),
            ],
        ),
    ]


def main() -> None:
    _, branch_metadata = load(BRANCH_PATH)
    groups = branch_metadata["dimensionless_groups"]
    radius0 = float(
        branch_metadata["resolved_system_parameters"]["equilibrium_radius_m"]
    )
    render(
        name="rpe-echo-jet2-free-training-evidence",
        label="Free spectral training",
        page_label="Broad free-transient identification data",
        source=FREE_PATH,
        figures=free_training_figures(groups, radius0),
        columns=2,
        height=1200,
    )
    render(
        name="rpe-echo-jet2-branch-training-evidence",
        label="Forced branch training",
        page_label="Target-local echo-jet identification data",
        source=BRANCH_PATH,
        figures=branch_training_figures(groups, radius0),
        columns=2,
        height=1700,
    )
    render(
        name="rpe-echo-jet2-spectrum-fit-evidence",
        label="Spectrum selection and echo-jet refit",
        page_label="KROM extraction evidence",
        source=KROM_PATH,
        figures=krom_figures(),
        columns=2,
        height=1200,
    )


if __name__ == "__main__":
    main()
