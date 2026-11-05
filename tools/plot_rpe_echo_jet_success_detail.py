"""Render the nonlinear echo-jet success case through Viewer PlotCanvas."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TARGET_PATH = ROOT / (
    "archives/solutions/rpe/"
    "rpe-periodic-100khz-40kpa-settled-target__20260914T183110Z.npz"
)
CONTROL_PATH = ROOT / (
    "archives/solutions/rpe/"
    "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3-tau300__"
    "20260914T192234Z.npz"
)
VIEWER_BUNDLE = ROOT / "viewer_outputs/rpe-echo-jet2-40kpa-success-viewer"
PAYLOAD_PATH = VIEWER_BUNDLE / "payload.json"
CATALOG_PATH = VIEWER_BUNDLE / "catalog.json"
RENDER_REPORT_PATH = ROOT / "viewer_outputs/rpe-echo-jet2-40kpa-success-detail.render.json"
OUTPUT_PATH = ROOT / "viewer_outputs/rpe-echo-jet2-40kpa-success-detail.png"
CATALOG_ID = "rpe-echo-jet2-40kpa-success-detail"
TILE_NAME = "RPE 40 kPa echo-jet KLQR — final six cycles"
CATALOG_LABEL = "40 kPa echo-jet KLQR"
PAGE_LABEL = "Successful nonlinear control"


def load(path: Path) -> tuple[dict[str, np.ndarray], dict]:
    with np.load(path, allow_pickle=False) as archive:
        arrays = {
            name: archive[name]
            for name in archive.files
            if name != "metadata_json"
        }
        metadata = json.loads(str(archive["metadata_json"]))
    return arrays, metadata


def echo_jet(
    states: np.ndarray, groups: dict[str, float], radius0: float
) -> np.ndarray:
    """Evaluate [echo, L_f0 echo] using the archived RPE convention."""
    radius = states[:, 0]
    velocity = states[:, 1]
    gas = groups["gas_pressure_ratio"]
    vapor = groups["vapor_pressure_ratio"]
    ambient = groups["ambient_pressure_ratio"]
    capillary = groups["capillarity"]
    viscous = groups["viscosity"]
    pressure_scale = groups["pressure_scale_pa"]

    momentum = (
        gas * radius**-4
        - 1.5 * velocity**2
        + vapor
        - ambient
        - capillary / radius
        - viscous * velocity / radius
    )
    acceleration = momentum / radius
    scale = radius0 * pressure_scale
    echo = scale * (2.0 * radius * velocity**2 + momentum)
    echo_dr = scale * (
        2.0 * velocity**2
        - 4.0 * gas * radius**-5
        + capillary / radius**2
        + viscous * velocity / radius**2
    )
    echo_dv = scale * (
        4.0 * radius * velocity - 3.0 * velocity - viscous / radius
    )
    echo_dot = echo_dr * velocity + echo_dv * acceleration
    return np.column_stack((echo, echo_dot))


def values(array: np.ndarray) -> list[float]:
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
) -> dict:
    return {
        "label": label,
        "role": role,
        "color": color,
        "line_style": line_style,
        "encoding": encoding,
        "selection_id": -1,
        "selected": False,
        "marker": False,
        "heatmap": False,
        "x": values(x),
        "y": values(y),
        "z": [],
        "color_values": [] if color_values is None else values(color_values),
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
        "square": False,
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
    path = Path(local) / "Programs/DynSysLab/bin/DynSysLabViewerCli.exe"
    if not path.is_file():
        raise RuntimeError(f"installed Viewer CLI is unavailable: {path}")
    return path


def main() -> None:
    target, _ = load(TARGET_PATH)
    controlled, metadata = load(CONTROL_PATH)

    time = np.asarray(controlled["time"], dtype=float)
    if not np.allclose(time, target["time"], rtol=0.0, atol=1e-12):
        raise RuntimeError("controlled and target archives do not share a time grid")

    state = np.asarray(controlled["values"], dtype=float)
    target_state = np.asarray(controlled["control_target_values"], dtype=float)
    if not np.allclose(target_state, target["values"], rtol=0.0, atol=1e-12):
        raise RuntimeError("embedded and source state targets differ")

    target_echo_jet = np.asarray(
        controlled["control_observable_target_values"], dtype=float
    )
    actual_echo_jet = echo_jet(
        state,
        metadata["dimensionless_groups"],
        metadata["resolved_system_parameters"]["equilibrium_radius_m"],
    )
    control = np.asarray(controlled["control_values"], dtype=float)[:, 0]
    known = np.asarray(target["control_values"], dtype=float)[:, 0]

    omega = (
        2.0
        * np.pi
        * 100000.0
        * float(metadata["time_scale_to_seconds"])
    )
    period = 2.0 * np.pi / omega
    cycles_to_show = 6
    mask = time >= time[-1] - cycles_to_show * period
    t = time[mask]
    x = state[mask]
    x_target = target_state[mask]
    observed = actual_echo_jet[mask]
    observed_target = target_echo_jet[mask]

    target_style = {
        "role": "control_target",
        "color": "#ef4444",
        "line_style": "dash",
    }
    figures = [
        figure(
            "time",
            "Control history",
            "\\tau",
            "u",
            [
                trace("u_{KLQR}", t, control[mask], role="control"),
                trace("known sine driver", t, known[mask], **target_style),
            ],
        ),
        figure(
            "time",
            "r versus time",
            "\\tau",
            "r",
            [
                trace("r", t, x[:, 0]),
                trace("target", t, x_target[:, 0], **target_style),
            ],
        ),
        figure(
            "time",
            "\\frac{dr}{d\\tau} versus time",
            "\\tau",
            "\\frac{dr}{d\\tau}",
            [
                trace("\\frac{dr}{d\\tau}", t, x[:, 1]),
                trace("target", t, x_target[:, 1], **target_style),
            ],
        ),
        figure(
            "observable_time",
            "Far-field acoustic echo",
            "\\tau",
            "e",
            [
                trace("echo", t, observed[:, 0]),
                trace("target", t, observed_target[:, 0], **target_style),
            ],
        ),
        figure(
            "phase",
            "Phase space",
            "r",
            "\\frac{dr}{d\\tau}",
            [
                trace(
                    "solution",
                    x[:, 0],
                    x[:, 1],
                    encoding="time",
                    color_values=t,
                ),
                trace("target", x_target[:, 0], x_target[:, 1], **target_style),
            ],
            color_label="\\tau",
            references=[
                reference(float(x[0, 0]), float(x[0, 1]), "IC", False),
                reference(1.0, 0.0, "EC", True),
            ],
        ),
        figure(
            "phase",
            "Echo phase space",
            "e",
            "L_{f_0}e",
            [
                trace(
                    "echo jet",
                    observed[:, 0],
                    observed[:, 1],
                    encoding="time",
                    color_values=t,
                ),
                trace(
                    "target",
                    observed_target[:, 0],
                    observed_target[:, 1],
                    **target_style,
                ),
            ],
            color_label="\\tau",
            references=[
                reference(
                    float(observed[0, 0]), float(observed[0, 1]), "IC", False
                ),
                reference(0.0, 0.0, "EC", True),
            ],
        ),
    ]

    payload = {
        "schema": "dynsyslab.viewer.payload/v1",
        "source": str(CONTROL_PATH),
        "view": "control",
        "tile": {
            "name": TILE_NAME,
            "columns": 2,
        },
        "page_labels": [PAGE_LABEL],
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
                "id": CATALOG_ID,
                "label": CATALOG_LABEL,
                "kind": "trajectory",
                "source": PAYLOAD_PATH.name,
                "views": ["control"],
            }
        ],
    }
    PAYLOAD_PATH.parent.mkdir(parents=True, exist_ok=True)
    viewer_archives = VIEWER_BUNDLE / "archives"
    viewer_cache = VIEWER_BUNDLE / "cache"
    viewer_archives.mkdir(parents=True, exist_ok=True)
    viewer_cache.mkdir(parents=True, exist_ok=True)
    PAYLOAD_PATH.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    CATALOG_PATH.write_text(json.dumps(catalog, indent=2), encoding="utf-8")

    completed = subprocess.run(
        [
            str(viewer_cli()),
            "render",
            "--id",
            CATALOG_ID,
            "--view",
            "control",
            "--project",
            str(VIEWER_BUNDLE),
            "--archives",
            str(viewer_archives),
            "--cache",
            str(viewer_cache),
            "--catalog",
            str(CATALOG_PATH),
            "--output",
            str(OUTPUT_PATH),
            "--width",
            "1900",
            "--height",
            "1700",
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
            "Viewer render failed:\n"
            + (completed.stderr.strip() or completed.stdout.strip())
        )
    render_report = json.loads(completed.stdout)
    if render_report.get("renderer") != "Viewer PlotCanvas":
        raise RuntimeError("plot did not use Viewer PlotCanvas")
    RENDER_REPORT_PATH.write_text(
        json.dumps(render_report, indent=2), encoding="utf-8"
    )
    print(OUTPUT_PATH)
    print(RENDER_REPORT_PATH)


if __name__ == "__main__":
    main()
