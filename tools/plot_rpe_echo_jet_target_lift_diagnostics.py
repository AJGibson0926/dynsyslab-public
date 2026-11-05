"""Render how the echo-jet target becomes the lifted KLQR target."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from plot_rpe_echo_jet_pipeline_evidence import figure, render, trace
from plot_rpe_echo_jet_success_detail import echo_jet, load


ROOT = Path(__file__).resolve().parents[1]
KROM_PATH = ROOT / (
    "archives/derived/koopman-with-sindy/"
    "rpe-echo-study-state-basis-baseline__20260911T104742Z__"
    "koopman-with-sindy__20260914T165124Z.npz"
)
CASES = (
    {
        "name": "rpe-echo-jet2-20kpa-target-lift-diagnostic",
        "label": "20 kPa target-to-lift map",
        "page": "Last six cycles of the low-amplitude case",
        "path": ROOT / (
            "archives/solutions/rpe/"
            "rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300__"
            "20260914T171839Z.npz"
        ),
    },
    {
        "name": "rpe-echo-jet2-40kpa-target-lift-diagnostic",
        "label": "40 kPa target-to-lift map",
        "page": "Last six cycles of the nonlinear case",
        "path": ROOT / (
            "archives/solutions/rpe/"
            "rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3-tau300__"
            "20260914T192234Z.npz"
        ),
    },
)
METRICS_PATH = ROOT / "viewer_outputs/rpe-echo-jet2-target-lift-diagnostics.json"


def polynomial_values(
    coordinates: np.ndarray,
    center: np.ndarray,
    scale: np.ndarray,
    exponents: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    eta = (coordinates - center) / scale
    theta = np.prod(
        eta[:, None, :] ** exponents.astype(int)[None, :, :], axis=2
    )
    return theta @ coefficients


def relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    centered = reference - np.mean(reference)
    denominator = float(np.sqrt(np.mean(np.abs(centered) ** 2)))
    return float(np.sqrt(np.mean(np.abs(reference - candidate) ** 2)) / denominator)


def case_figures(case: dict, model: dict) -> tuple[list[dict], dict]:
    arrays, metadata = load(case["path"])
    time = np.asarray(arrays["time"], dtype=float)
    state = np.asarray(arrays["values"], dtype=float)
    target_state = np.asarray(arrays["control_target_values"], dtype=float)
    target_coordinates = np.asarray(
        arrays["control_observable_target_values"], dtype=float
    )
    groups = metadata["dimensionless_groups"]
    radius0 = float(
        metadata["resolved_system_parameters"]["equilibrium_radius_m"]
    )
    coordinates = echo_jet(state, groups, radius0)
    recomputed_target = echo_jet(target_state, groups, radius0)
    coordinate_match = float(np.max(np.abs(target_coordinates - recomputed_target)))

    psi = polynomial_values(coordinates, **model)
    psi_target = polynomial_values(target_coordinates, **model)
    centered = np.broadcast_to(model["center"], target_coordinates.shape)
    psi_echo_only = polynomial_values(
        np.column_stack((target_coordinates[:, 0], centered[:, 1])), **model
    )
    psi_echo_dot_only = polynomial_values(
        np.column_stack((centered[:, 0], target_coordinates[:, 1])), **model
    )

    omega = 2.0 * np.pi * 100000.0 * float(metadata["time_scale_to_seconds"])
    period = 2.0 * np.pi / omega
    mask = time >= time[-1] - 6.0 * period
    t = time[mask]
    y = coordinates[mask]
    yt = target_coordinates[mask]
    z = psi[mask]
    zt = psi_target[mask]
    target_style = {
        "role": "control_target",
        "color": "#ef4444",
        "line_style": "dash",
    }
    figures = [
        figure(
            "observable_time",
            "Echo coordinate supplied to the target lift",
            "\\tau",
            "e",
            [trace("actual", t, y[:, 0]), trace("target", t, yt[:, 0], **target_style)],
        ),
        figure(
            "observable_time",
            "Autonomous echo-dot coordinate supplied to the target lift",
            "\\tau",
            "L_{f_0}e",
            [trace("actual", t, y[:, 1]), trace("target", t, yt[:, 1], **target_style)],
        ),
        figure(
            "time",
            "Real lifted coordinate",
            "\\tau",
            "Re(\\psi)",
            [trace("actual", t, z.real), trace("target", t, zt.real, **target_style)],
        ),
        figure(
            "time",
            "Imaginary lifted coordinate",
            "\\tau",
            "Im(\\psi)",
            [trace("actual", t, z.imag), trace("target", t, zt.imag, **target_style)],
        ),
        figure(
            "phase",
            "Target path through the two-component echo chart",
            "e",
            "L_{f_0}e",
            [
                trace(
                    "actual",
                    y[:, 0],
                    y[:, 1],
                    encoding="time",
                    color_values=t,
                ),
                trace("target", yt[:, 0], yt[:, 1], **target_style),
            ],
            color_label="\\tau",
        ),
        figure(
            "time",
            "The errors actually penalized by Q",
            "\\tau",
            "z-z_{ref}",
            [
                trace("Re error", t, (z - zt).real, color="#3478c8"),
                trace("Im error", t, (z - zt).imag, color="#ef4444"),
            ],
        ),
    ]
    settled = time >= 100.0
    metrics = {
        "archive": str(case["path"].relative_to(ROOT)).replace("\\", "/"),
        "sample_count": int(time.size),
        "observable_target_shape": list(target_coordinates.shape),
        "observable_target_names": str(arrays["control_observable_target_names"]).splitlines(),
        "max_archived_vs_recomputed_target_coordinate_error": coordinate_match,
        "settled_target_ranges": {
            "echo": [
                float(target_coordinates[settled, 0].min()),
                float(target_coordinates[settled, 0].max()),
            ],
            "autonomous_echo_dot": [
                float(target_coordinates[settled, 1].min()),
                float(target_coordinates[settled, 1].max()),
            ],
            "psi_real": [
                float(psi_target[settled].real.min()),
                float(psi_target[settled].real.max()),
            ],
            "psi_imag": [
                float(psi_target[settled].imag.min()),
                float(psi_target[settled].imag.max()),
            ],
        },
        "settled_full_lift_difference_if_coordinate_frozen_at_training_center": {
            "freeze_autonomous_echo_dot": relative_rms(
                psi_target[settled], psi_echo_only[settled]
            ),
            "freeze_echo": relative_rms(
                psi_target[settled], psi_echo_dot_only[settled]
            ),
        },
        "settled_lift_tracking_rms": {
            "real": float(np.sqrt(np.mean((psi[settled].real - psi_target[settled].real) ** 2))),
            "imaginary": float(np.sqrt(np.mean((psi[settled].imag - psi_target[settled].imag) ** 2))),
        },
    }
    return figures, metrics


def main() -> None:
    with np.load(KROM_PATH, allow_pickle=False) as archive:
        model = {
            "center": np.asarray(archive["state_center"], dtype=float),
            "scale": np.asarray(archive["state_scale"], dtype=float),
            "exponents": np.asarray(archive["polynomial_exponents"], dtype=int),
            "coefficients": np.asarray(archive["coefficients"], dtype=complex)[0],
        }
    evidence = {"schema": "dynsyslab.echo-jet-target-lift-diagnostics/v1", "cases": {}}
    for case in CASES:
        figures, metrics = case_figures(case, model)
        render(
            name=case["name"],
            label=case["label"],
            page_label=case["page"],
            source=case["path"],
            figures=figures,
            columns=2,
            height=1700,
        )
        evidence["cases"][case["name"]] = metrics
    METRICS_PATH.write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    print(METRICS_PATH)


if __name__ == "__main__":
    main()
