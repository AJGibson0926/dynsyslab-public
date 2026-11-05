"""Lazy numerical payload generation for the C++ Viewer workbench."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .contracts import PAYLOAD_SCHEMA


def _plain(value: Any) -> Any:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, complex):
        return {"real": float(value.real), "imaginary": float(value.imag)}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _flatten(value: Any) -> tuple[list[float], list[int]]:
    array = np.asarray(value)
    if np.iscomplexobj(array):
        array = np.abs(array)
    return np.asarray(array, dtype=float).reshape(-1).tolist(), list(array.shape)


def _trace_payload(trace: dict[str, Any]) -> dict[str, Any]:
    result = {
        "label": str(trace.get("label", "series")),
        "role": str(trace.get("role", "data")),
        "color": str(trace.get("color", "")),
        "line_style": str(trace.get("line_style", "solid")),
        "marker": bool(trace.get("marker", False)),
        "heatmap": bool(trace.get("heatmap", False)),
        "encoding": str(trace.get("encoding", "")),
    }
    for key in ("x", "y", "z"):
        if key in trace:
            values, shape = _flatten(trace[key])
            result[key], result[f"{key}_shape"] = values, shape
    if "value" in trace:
        values = np.asarray(trace["value"])
        result["color_values"] = np.angle(values).reshape(-1).tolist() if np.iscomplexobj(values) else values.astype(float).reshape(-1).tolist()
        result["magnitudes"] = np.abs(values).reshape(-1).tolist()
    return result


def _report_payload(report: Any, source: Path, view: str) -> dict[str, Any]:
    figures = []
    for figure in report.figures:
        metadata = dict(figure.metadata or {})
        axes = metadata.get("axes") or {}
        figures.append({
            "kind": figure.kind,
            "title": figure.title,
            "x_label": str(axes.get("x", "")),
            "y_label": str(axes.get("y", "")),
            "square": bool(metadata.get("square") or figure.kind in {"phase", "koopman_spectrum", "control_error", "observable_heatmap"}),
            "detail_inset": _plain(metadata.get("detail_inset")),
            "metadata": _plain({key: value for key, value in metadata.items()
                               if key not in {"detail_inset", "archive_metadata"}}),
            "traces": [_trace_payload(trace) for trace in figure.traces],
        })
    return {
        "schema": PAYLOAD_SCHEMA,
        "source": str(source),
        "view": view,
        "tile": {"name": report.tile_plan.name, "columns": report.tile_plan.columns},
        "warnings": list(report.warnings),
        "figures": figures,
    }


def _viewer_payload(source: Path, view: str, *, observable: str = "state_0") -> dict[str, Any]:
    from viewer import (ArchiveSource, CampaignSource, ControlErrorSource, KLQESource, KoopmanSource,
                        ObservabilitySource, PlotRequest, TilePlan, Viewer)

    if source.suffix.lower() == ".json" and view in {"time", "phase", "fourier", "bifurcation"}:
        viewer = Viewer(tile_plan=TilePlan("campaign", columns=2)).add(
            PlotRequest(view, CampaignSource(source), max_display_samples=1800))
    elif view in {"overview", "control"}:
        viewer = Viewer(tile_plan=TilePlan("workbench", columns=2))
        viewer.add(PlotRequest("phase", ArchiveSource(source), max_display_samples=1800))
        viewer.add(PlotRequest("time", ArchiveSource(source), max_display_samples=1800))
    elif view in {"time", "phase", "three_d", "fourier"}:
        components = (0, 1, 2) if view == "three_d" else (0, 1)
        viewer = Viewer(tile_plan=TilePlan("workbench", columns=2)).add(
            PlotRequest(view, ArchiveSource(source), state_components=components,
                        max_display_samples=2200))
    elif view in {"observable_time", "observable_heatmap"}:
        viewer = Viewer().add(PlotRequest(view, ArchiveSource(source), observable=observable,
                                          domain_density=80, max_display_samples=2200))
    elif view == "observable_spectrum_comparison":
        from viewer.observable_spectra import observable_spectrum_comparison_report
        return _report_payload(observable_spectrum_comparison_report(source),source,view)
    elif view in {"koopman_spectrum", "krom"}:
        viewer = Viewer(tile_plan=TilePlan("krom_workbench", columns=2)).add(
            PlotRequest(view, KoopmanSource(source), krom_mode_limit=5,
                        max_display_samples=1600))
    elif view == "control_error":
        viewer = Viewer(tile_plan=TilePlan("errors", columns=2)).add(
            PlotRequest("control_error", ControlErrorSource(source),
                        control_error_measure="relative"))
    elif view == "observability":
        from rom import load_observability_analysis
        viewer = Viewer(tile_plan=TilePlan("observability", columns=2)).add(
            PlotRequest("observability", ObservabilitySource(load_observability_analysis(source))))
    elif view == "klqe":
        from estimation import load_klqe_result
        viewer = Viewer(tile_plan=TilePlan("klqe", columns=2)).add(
            PlotRequest("klqe", KLQESource(load_klqe_result(source))))
    else:
        raise ValueError(f"View {view!r} is not backed by the Viewer report API.")
    return _report_payload(viewer.render(), source, view)


def _select_array(data: Any, names: Iterable[str], *, dimensions: int | None = None) -> np.ndarray | None:
    for name in names:
        if name not in data.files:
            continue
        value = np.asarray(data[name])
        if dimensions is None or value.ndim == dimensions:
            return value
    return None


def _generic_payload(source: Path) -> dict[str, Any]:
    figures = []
    with np.load(source, allow_pickle=False) as data:
        time = _select_array(data, ("time", "controlled_time", "control_time", "free_time"), dimensions=1)
        states = _select_array(data, ("states", "controlled_states", "values", "free_training_values"), dimensions=2)
        controls = _select_array(data, ("controls", "controlled_controls", "control_values", "controlled_controls_sample_aligned"))
        if states is not None:
            shown = states[: min(states.shape[0], 2200)]
            axis = np.arange(shown.shape[0]) if time is None or time.size < shown.shape[0] else time[:shown.shape[0]]
            traces = []
            for index in range(min(4, shown.shape[1])):
                traces.append({"label": f"state {index}", "x": axis.tolist(),
                               "y": shown[:, index].astype(float).tolist(), "color": ""})
            figures.append({"kind": "time", "title": "Dataset state components",
                            "x_label": "sample" if time is None else "time", "y_label": "state",
                            "square": False, "traces": traces, "metadata": {}})
            if shown.shape[1] >= 2:
                figures.append({"kind": "phase", "title": "Dataset phase coverage",
                                "x_label": "state 0", "y_label": "state 1", "square": True,
                                "traces": [{"label": "samples", "x": shown[:, 0].tolist(),
                                            "y": shown[:, 1].tolist(), "marker": True, "color": ""}],
                                "metadata": {}})
        if controls is not None:
            values = np.asarray(controls).reshape(np.asarray(controls).shape[0], -1)
            values = values[: min(values.shape[0], 2200)]
            figures.append({"kind": "time", "title": "Control excitation",
                            "x_label": "sample", "y_label": "control", "square": False,
                            "traces": [{"label": f"u{index}", "x": list(range(values.shape[0])),
                                        "y": values[:, index].astype(float).tolist(), "color": ""}
                                       for index in range(min(4, values.shape[1]))], "metadata": {}})
        if not figures:
            numeric = [(name, np.asarray(data[name])) for name in data.files
                       if name != "metadata_json" and np.asarray(data[name]).dtype.kind in "iufc"]
            for name, values in numeric[:4]:
                flat = np.abs(values).reshape(-1)[:2200].astype(float)
                figures.append({"kind": "time", "title": name, "x_label": "index",
                                "y_label": name, "square": False,
                                "traces": [{"label": name, "x": list(range(flat.size)),
                                            "y": flat.tolist(), "color": ""}], "metadata": {}})
    return {"schema": PAYLOAD_SCHEMA, "source": str(source), "view": "data_overview",
            "tile": {"name": "dataset", "columns": 2}, "warnings": [], "figures": figures}


def build_payload(source: str | Path, view: str, *, observable: str = "state_0") -> dict[str, Any]:
    path = Path(source).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if view == "data_overview":
        return _generic_payload(path)
    return _viewer_payload(path, view, observable=observable)


def write_payload(source: str | Path, view: str, output: str | Path, *, observable: str = "state_0") -> Path:
    target = Path(output).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(build_payload(source, view, observable=observable),
                                 separators=(",", ":")), encoding="utf-8")
    return target
