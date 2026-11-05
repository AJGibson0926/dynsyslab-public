"""Viewer presentation helpers for echo-Hankel eigenfunction experiments."""

from __future__ import annotations

import numpy as np

from plot_rpe_echo_jet_pipeline_evidence import figure, trace


TRAJECTORY_STYLES = (
    ("#111827", 0.88),
    ("#26313d", 0.78),
    ("#3f4b57", 0.68),
    ("#5f6b76", 0.58),
    ("#7f8993", 0.48),
)


def heat_trace(
    label: str,
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    *,
    encoding: str = "",
    magnitudes: np.ndarray | None = None,
) -> dict:
    """Serialize a heat-map trace, including optional domain-coloring magnitudes."""
    flattened = np.asarray(values, dtype=float).reshape(-1)
    return {
        "label": label,
        "role": "data",
        "color": "",
        "line_style": "solid",
        "encoding": encoding,
        "selection_id": -1,
        "selected": False,
        "marker": False,
        "marker_size": 3.0,
        "width": 1.0,
        "opacity": 1.0,
        "heatmap": True,
        "x": np.asarray(x, dtype=float).tolist(),
        "y": np.asarray(y, dtype=float).tolist(),
        "z": [float(value) if np.isfinite(value) else None for value in flattened],
        "color_values": [],
        "magnitudes": [] if magnitudes is None else np.asarray(magnitudes).reshape(-1).tolist(),
        "z_shape": list(np.asarray(values).shape),
    }


def trajectory_overlays(
    projected: list[dict], start: float, end: float, *, stride: int = 30
) -> list[dict]:
    """Build dashed grayscale overlays for aligned training trajectories."""
    overlays = []
    for index, (record, (color, opacity)) in enumerate(
        zip(projected, TRAJECTORY_STYLES, strict=True), start=1
    ):
        mask = (record["time"] >= start) & (record["time"] <= end)
        coordinates = record["coordinates"][mask, :2]
        overlays.append(
            trace(
                f"training #{index}",
                coordinates[::stride, 0],
                coordinates[::stride, 1],
                color=color,
                line_style="dash",
                opacity=opacity,
            )
        )
    return overlays


def eigenfunction_history(
    model,
    coefficient: np.ndarray,
    projected: list[dict],
    start: float,
    end: float,
    *,
    component: str,
) -> dict:
    """Compose histories of the requested real, imaginary, magnitude, or phase component."""
    components = {
        "real": (np.real, "Real part along training trajectories", r"\operatorname{Re}(\psi)"),
        "imaginary": (
            np.imag,
            "Imaginary part along training trajectories",
            r"\operatorname{Im}(\psi)",
        ),
        "magnitude": (np.abs, "Magnitude along training trajectories", r"|\psi|"),
        "phase": (np.angle, "Phase along training trajectories", r"\operatorname{arg}(\psi)"),
    }
    if component not in components:
        raise ValueError(f"unsupported eigenfunction history component: {component}")
    transform, title, y_label = components[component]
    traces = []
    for index, (record, (color, opacity)) in enumerate(
        zip(projected, TRAJECTORY_STYLES, strict=True), start=1
    ):
        mask = (record["time"] >= start) & (record["time"] <= end)
        time = record["time"][mask]
        values = model.evaluate(record["coordinates"][mask, :2], coefficient)
        traces.append(
            trace(
                f"training #{index}",
                time,
                transform(values),
                color=color,
                opacity=max(opacity, 0.58),
            )
        )
    return figure(
        "krom_time",
        title,
        r"\tau",
        y_label,
        traces,
    )
