"""Audit injectivity of the fixed echo-history coordinates on unseen orbits.

The retained U and singular values are reused without any SVD, HAVOK fit, or
basis adaptation.  The audit detects proper intersections of representative
closed curves in (xi_1, xi_2), measures the physical-state ambiguity at each
crossing, and tests whether already-retained coordinates xi_3...xi_6 separate
the colliding histories.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

from analyze_rpe_echo_hankel_latents import (
    OUTPUT_ROOT,
    figure,
    reference,
    render_viewer,
    trace,
)
from analyze_rpe_inviscid_echo_hankel import far_field_echo, load_native_source
from project_rpe_inviscid_echo_basis_out_of_sample import (
    COLORS,
    EMBEDDING_DECIMATION,
    EXPECTED_ICS,
    EXPECTED_STEP,
    SOURCE_GLOB,
    TRAINING,
    causal_fixed_basis_projection,
    validate_validation_source,
)


ROOT = Path(__file__).resolve().parents[1]
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-latent-injectivity.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-latent-injectivity.npz"
PLOT_NAME = "rpe-inviscid-echo-latent-injectivity"
AUDIT_DIMENSIONS = (2, 3, 4, 6)
PHASE_SAMPLES = 1201


def cross2(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Return the signed two-dimensional vector cross product."""
    return left[..., 0] * right[..., 1] - left[..., 1] * right[..., 0]


def representative_cycle(
    echo: np.ndarray,
    coordinates: np.ndarray,
    state: np.ndarray,
    delay: int,
    sample_step: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a representative late orbit segment for geometric comparisons."""
    sampled_echo = echo[::EMBEDDING_DECIMATION]
    latent_state = state[::EMBEDDING_DECIMATION][delay - 1 :]
    peaks, _ = find_peaks(
        sampled_echo,
        distance=max(1, int(round(0.5 / sample_step))),
        prominence=0.15 * float(np.ptp(sampled_echo)),
    )
    latent_peaks = peaks - (delay - 1)
    latent_peaks = latent_peaks[
        (latent_peaks >= 0) & (latent_peaks < coordinates.shape[0])
    ]
    if latent_peaks.size < 3:
        raise RuntimeError("too few complete latent cycles for injectivity audit")
    start, stop = (int(value) for value in latent_peaks[-3:-1])
    source_index = np.arange(start, stop + 1, dtype=float)
    target_index = np.linspace(start, stop, PHASE_SAMPLES)
    cycle_coordinates = np.column_stack(
        [
            np.interp(target_index, source_index, coordinates[start : stop + 1, mode])
            for mode in range(coordinates.shape[1])
        ]
    )
    cycle_state = np.column_stack(
        [
            np.interp(target_index, source_index, latent_state[start : stop + 1, mode])
            for mode in range(latent_state.shape[1])
        ]
    )
    return cycle_coordinates, cycle_state


def proper_intersections(
    left: np.ndarray,
    left_state: np.ndarray,
    right: np.ndarray,
    right_state: np.ndarray,
    latent_scale: np.ndarray,
    state_scale: np.ndarray,
) -> list[dict]:
    """Return deduplicated proper crossings of two closed projected curves."""
    left_xy = left[:, :2] / latent_scale[:2]
    right_xy = right[:, :2] / latent_scale[:2]
    q = right_xy[:-1]
    s = np.diff(right_xy, axis=0)
    q_min = np.minimum(q, q + s)
    q_max = np.maximum(q, q + s)
    found: list[dict] = []

    for left_index, (p, r) in enumerate(
        zip(left_xy[:-1], np.diff(left_xy, axis=0), strict=True)
    ):
        p_min = np.minimum(p, p + r)
        p_max = np.maximum(p, p + r)
        candidates = np.flatnonzero(
            (q_max[:, 0] >= p_min[0])
            & (q_min[:, 0] <= p_max[0])
            & (q_max[:, 1] >= p_min[1])
            & (q_min[:, 1] <= p_max[1])
        )
        if candidates.size == 0:
            continue
        q_candidate = q[candidates]
        s_candidate = s[candidates]
        denominator = cross2(np.broadcast_to(r, s_candidate.shape), s_candidate)
        usable = np.abs(denominator) > 1.0e-13
        if not np.any(usable):
            continue
        candidates = candidates[usable]
        q_candidate = q_candidate[usable]
        s_candidate = s_candidate[usable]
        denominator = denominator[usable]
        offset = q_candidate - p
        left_fraction = cross2(offset, s_candidate) / denominator
        right_fraction = cross2(
            offset, np.broadcast_to(r, offset.shape)
        ) / denominator
        crossing = (
            (left_fraction > 1.0e-9)
            & (left_fraction < 1.0 - 1.0e-9)
            & (right_fraction > 1.0e-9)
            & (right_fraction < 1.0 - 1.0e-9)
        )
        for right_index, alpha, beta in zip(
            candidates[crossing],
            left_fraction[crossing],
            right_fraction[crossing],
            strict=True,
        ):
            left_latent = left[left_index] + alpha * (
                left[left_index + 1] - left[left_index]
            )
            right_latent = right[right_index] + beta * (
                right[right_index + 1] - right[right_index]
            )
            left_physical = left_state[left_index] + alpha * (
                left_state[left_index + 1] - left_state[left_index]
            )
            right_physical = right_state[right_index] + beta * (
                right_state[right_index + 1] - right_state[right_index]
            )
            separation = {
                str(dimension): float(
                    np.linalg.norm(
                        (left_latent[:dimension] - right_latent[:dimension])
                        / latent_scale[:dimension]
                    )
                )
                for dimension in AUDIT_DIMENSIONS
            }
            found.append(
                {
                    "point": ((left_latent[:2] + right_latent[:2]) / 2.0).tolist(),
                    "left_segment": left_index,
                    "right_segment": int(right_index),
                    "left_fraction": float(alpha),
                    "right_fraction": float(beta),
                    "left_state": left_physical.tolist(),
                    "right_state": right_physical.tolist(),
                    "normalized_state_separation": float(
                        np.linalg.norm(
                            (left_physical - right_physical) / state_scale
                        )
                    ),
                    "normalized_latent_separation": separation,
                }
            )

    deduplicated: list[dict] = []
    for candidate in found:
        point = np.asarray(candidate["point"]) / latent_scale[:2]
        if any(
            np.linalg.norm(
                point - np.asarray(existing["point"]) / latent_scale[:2]
            )
            < 2.0e-4
            for existing in deduplicated
        ):
            continue
        deduplicated.append(candidate)
    return deduplicated


def polygon_area(curve: np.ndarray) -> float:
    """Compute the planar polygon area of a sampled closed curve."""
    x = curve[:, 0]
    y = curve[:, 1]
    return float(0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def pairwise_minimum_separation(
    left: np.ndarray,
    right: np.ndarray,
    latent_scale: np.ndarray,
    dimension: int,
) -> float:
    """Measure nearest pairwise separations between the sampled orbit curves."""
    right_normalized = right[:, :dimension] / latent_scale[:dimension]
    tree = cKDTree(right_normalized)
    distances, _ = tree.query(left[:, :dimension] / latent_scale[:dimension], k=1)
    return float(np.min(distances))


def phase_tile(title: str, x_label: str, y_label: str, traces: list[dict]) -> dict:
    """Compose a square Viewer phase-plane tile from the supplied curves."""
    tile = figure("phase", title, x_label, y_label, traces)
    tile["square"] = True
    return tile


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    with np.load(TRAINING, allow_pickle=False) as training:
        modes = np.asarray(training["delay_modes"][:, :6], dtype=float)
        singular = np.asarray(training["leading_singular_values"][:6], dtype=float)
        mean = float(training["echo_mean"])
        scale = float(training["echo_scale"])
        training_radius = np.asarray(training["radius"], dtype=float)
        training_velocity = np.asarray(training["radial_velocity"], dtype=float)

    sources = sorted((ROOT / "archives/solutions/rpe").glob(SOURCE_GLOB))
    if len(sources) != 5:
        raise RuntimeError("the five-trajectory validation set is unavailable")

    delay = modes.shape[0]
    sample_step = EXPECTED_STEP * EMBEDDING_DECIMATION
    trajectories: list[dict] = []
    saved: dict[str, np.ndarray] = {}
    for index, (source, expected_ic) in enumerate(
        zip(sources, EXPECTED_ICS, strict=True), start=1
    ):
        time, state, metadata = load_native_source(source)
        validate_validation_source(metadata)
        initial = np.asarray(metadata["initial_state"], dtype=float)
        if not np.array_equal(initial, expected_ic):
            raise RuntimeError("validation IC order changed")
        echo = far_field_echo(state, metadata)
        coordinates, captured = causal_fixed_basis_projection(
            echo, modes, singular, mean, scale
        )
        cycle, cycle_state = representative_cycle(
            echo, coordinates, state, delay, sample_step
        )
        trajectories.append(
            {
                "index": index,
                "initial": initial,
                "state": state,
                "coordinates": coordinates,
                "cycle": cycle,
                "cycle_state": cycle_state,
                "captured_six": captured,
            }
        )
        saved[f"trajectory_{index}_latent"] = coordinates
        saved[f"trajectory_{index}_cycle_latent"] = cycle
        saved[f"trajectory_{index}_cycle_state"] = cycle_state

    all_cycle_latent = np.vstack([item["cycle"] for item in trajectories])
    latent_scale = np.ptp(all_cycle_latent, axis=0)
    latent_scale[latent_scale <= np.finfo(float).eps] = 1.0
    state_scale = np.asarray(
        [
            float(np.ptp(training_radius)),
            2.0 * float(np.max(np.abs(training_velocity))),
        ]
    )

    intersections: list[dict] = []
    pair_metrics: list[dict] = []
    for left, right in itertools.combinations(trajectories, 2):
        pair_crossings = proper_intersections(
            left["cycle"],
            left["cycle_state"],
            right["cycle"],
            right["cycle_state"],
            latent_scale,
            state_scale,
        )
        for crossing in pair_crossings:
            crossing["trajectory_pair"] = [left["index"], right["index"]]
        intersections.extend(pair_crossings)
        pair_metrics.append(
            {
                "trajectory_pair": [left["index"], right["index"]],
                "proper_intersections_xi1_xi2": len(pair_crossings),
                "minimum_normalized_latent_separation": {
                    str(dimension): pairwise_minimum_separation(
                        left["cycle"], right["cycle"], latent_scale, dimension
                    )
                    for dimension in AUDIT_DIMENSIONS
                },
            }
        )

    areas = [polygon_area(item["cycle"][:, :2]) for item in trajectories]
    area_monotone = bool(np.all(np.diff(areas) > 0.0))
    if intersections:
        dimension_minimum = {
            str(dimension): float(
                np.min(
                    [
                        crossing["normalized_latent_separation"][str(dimension)]
                        for crossing in intersections
                    ]
                )
            )
            for dimension in AUDIT_DIMENSIONS
        }
        dimension_median = {
            str(dimension): float(
                np.median(
                    [
                        crossing["normalized_latent_separation"][str(dimension)]
                        for crossing in intersections
                    ]
                )
            )
            for dimension in AUDIT_DIMENSIONS
        }
        dimension_maximum = {
            str(dimension): float(
                np.max(
                    [
                        crossing["normalized_latent_separation"][str(dimension)]
                        for crossing in intersections
                    ]
                )
            )
            for dimension in AUDIT_DIMENSIONS
        }
    else:
        dimension_minimum = {
            str(dimension): float(
                np.min(
                    [
                        pair["minimum_normalized_latent_separation"][str(dimension)]
                        for pair in pair_metrics
                    ]
                )
            )
            for dimension in AUDIT_DIMENSIONS
        }
        dimension_median = {
            str(dimension): float(
                np.median(
                    [
                        pair["minimum_normalized_latent_separation"][str(dimension)]
                        for pair in pair_metrics
                    ]
                )
            )
            for dimension in AUDIT_DIMENSIONS
        }
        dimension_maximum = dict(dimension_median)

    metrics = {
        "schema": "dynsyslab.echo-latent-injectivity-audit/v1",
        "method": "fixed-U causal projection plus representative-cycle intersection audit",
        "new_decomposition_performed": False,
        "basis_adaptation_performed": False,
        "training_artifact": str(TRAINING),
        "tested_dimensions": list(AUDIT_DIMENSIONS),
        "delay_samples": int(delay),
        "delay_window_tau": float((delay - 1) * sample_step),
        "representative_cycle_samples": PHASE_SAMPLES,
        "latent_scale_for_distance": latent_scale.tolist(),
        "state_scale_for_distance": state_scale.tolist(),
        "xi1_xi2_signed_area_magnitude_by_trajectory": areas,
        "areas_strictly_increase_with_orbit_energy": area_monotone,
        "proper_intersection_count_xi1_xi2": len(intersections),
        "intersection_physical_ambiguity": {
            "median_normalized_state_separation": (
                float(
                    np.median(
                        [item["normalized_state_separation"] for item in intersections]
                    )
                )
                if intersections
                else 0.0
            ),
            "maximum_normalized_state_separation": (
                float(
                    np.max(
                        [item["normalized_state_separation"] for item in intersections]
                    )
                )
                if intersections
                else 0.0
            ),
        },
        "intersection_separation_by_retained_dimension": {
            "minimum": dimension_minimum,
            "median": dimension_median,
            "maximum": dimension_maximum,
        },
        "pairwise": pair_metrics,
        "intersections": intersections,
        "interpretation": (
            "xi1-xi2 has proper cross-orbit intersections and is not globally injective"
            if intersections
            else "no proper cross-orbit intersections were detected on representative cycles"
        ),
    }
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        fixed_delay_modes=modes,
        fixed_singular_values=singular,
        metadata_json=np.asarray(json.dumps(metrics, separators=(",", ":"))),
        **saved,
    )

    labels = [
        f"IC #{index}: ({initial[0]:.2f}, {initial[1]:+.2f})"
        for index, initial in enumerate(EXPECTED_ICS, start=1)
    ]
    latent_12_traces = [
        trace(
            label,
            item["cycle"][:, 0],
            item["cycle"][:, 1],
            color=color,
            width=2.0,
        )
        for item, label, color in zip(trajectories, labels, COLORS, strict=True)
    ]
    if intersections:
        points = np.asarray([item["point"] for item in intersections])
        for point_index, point in enumerate(points):
            latent_12_traces.append(
                trace(
                    "proper intersections" if point_index == 0 else "",
                    point[0:1],
                    point[1:2],
                    color="#111827",
                    marker=True,
                    marker_size=7.0,
                    width=0.0,
                )
            )
    latent_123_traces = []
    for item, label, color in zip(trajectories, labels, COLORS, strict=True):
        plotted = trace(
            label,
            item["cycle"][:, 0],
            item["cycle"][:, 1],
            color=color,
            width=2.0,
        )
        plotted["z"] = [float(value) for value in item["cycle"][:, 2]]
        latent_123_traces.append(plotted)
    state_traces = [
        trace(
            label,
            item["cycle_state"][:, 0],
            item["cycle_state"][:, 1],
            color=color,
            width=2.0,
        )
        for item, label, color in zip(trajectories, labels, COLORS, strict=True)
    ]
    dimensions = np.asarray(AUDIT_DIMENSIONS, dtype=float)
    separation_tile = figure(
        "line",
        "Separation of xi1-xi2 collisions in added coordinates",
        "retained latent dimension",
        "normalized latent separation",
        [
            trace(
                "minimum",
                dimensions,
                np.maximum(
                    np.asarray([dimension_minimum[str(int(d))] for d in dimensions]),
                    1.0e-14,
                ),
                color="#f59e0b",
                marker=True,
                marker_size=5.0,
                width=2.0,
            ),
            trace(
                "median",
                dimensions,
                np.maximum(
                    np.asarray([dimension_median[str(int(d))] for d in dimensions]),
                    1.0e-14,
                ),
                color="#2563eb",
                marker=True,
                marker_size=5.0,
                width=2.0,
            ),
            trace(
                "maximum",
                dimensions,
                np.maximum(
                    np.asarray([dimension_maximum[str(int(d))] for d in dimensions]),
                    1.0e-14,
                ),
                color="#dc2626",
                marker=True,
                marker_size=5.0,
                width=2.0,
            ),
        ],
    )
    separation_tile["value_scale"] = "log"

    latent_123_tile = figure(
        "three_d",
        "Fixed coordinates xi1-xi2-xi3",
        "xi_1",
        "xi_2",
        latent_123_traces,
    )
    latent_123_tile["z_label"] = "xi_3"
    latent_123_tile["square"] = True

    figures = [
        phase_tile(
            "Fixed-U phase map xi1-xi2 — representative cycles",
            "xi_1",
            "xi_2",
            latent_12_traces,
        ),
        latent_123_tile,
        phase_tile(
            "Physical state for audited cycles",
            "r",
            "dr/dtau",
            state_traces,
        ),
        separation_tile,
    ]
    figures[2]["references"] = [reference(1.0, 0.0, "EC", True)] + [
        reference(float(ic[0]), float(ic[1]), f"IC #{index}", False)
        for index, ic in enumerate(EXPECTED_ICS, start=1)
    ]
    render_viewer(
        TRAINING,
        figures,
        name=PLOT_NAME,
        label="Fixed echo-history coordinate injectivity audit",
        columns=2,
        width=1900,
        height=1900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
