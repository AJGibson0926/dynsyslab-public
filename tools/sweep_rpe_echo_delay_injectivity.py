"""Measure top-two echo-coordinate injectivity versus delay-window length."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve
from scipy.sparse.linalg import LinearOperator, svds

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_inviscid_echo_hankel import far_field_echo, load_native_source
from audit_rpe_echo_latent_injectivity import (
    proper_intersections,
    representative_cycle,
)
from project_rpe_inviscid_echo_basis_out_of_sample import (
    COLORS,
    EMBEDDING_DECIMATION,
    EXPECTED_ICS,
    EXPECTED_STEP,
    SOURCE_GLOB,
    TRAINING,
    closure_metric,
    validate_validation_source,
)


ROOT = Path(__file__).resolve().parents[1]
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-delay-injectivity-sweep.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-delay-injectivity-sweep.npz"
SUMMARY_PLOT = "rpe-inviscid-echo-delay-injectivity-sweep"
PHASE_PLOT = "rpe-inviscid-echo-delay-injectivity-phase-comparison"
PERIOD_MULTIPLIERS = (0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0)
BASELINE_WINDOW_TAU = 1.275


class HankelOperator(LinearOperator):
    """FFT-backed matrix-free Hankel operator for one scalar record."""

    def __init__(self, signal: np.ndarray, delay: int):
        self.signal = np.asarray(signal, dtype=float)
        self.delay = int(delay)
        columns = self.signal.size - self.delay + 1
        super().__init__(dtype=np.dtype(float), shape=(self.delay, columns))

    def _matvec(self, vector: np.ndarray) -> np.ndarray:
        return fftconvolve(
            self.signal, np.asarray(vector).reshape(-1)[::-1], mode="valid"
        )

    def _rmatvec(self, vector: np.ndarray) -> np.ndarray:
        return fftconvolve(
            self.signal, np.asarray(vector).reshape(-1)[::-1], mode="valid"
        )

    def _matmat(self, matrix: np.ndarray) -> np.ndarray:
        return np.column_stack(
            [self._matvec(matrix[:, column]) for column in range(matrix.shape[1])]
        )

    def _rmatmat(self, matrix: np.ndarray) -> np.ndarray:
        return np.column_stack(
            [self._rmatvec(matrix[:, column]) for column in range(matrix.shape[1])]
        )


def hankel_total_energy(signal: np.ndarray, delay: int) -> float:
    """Compute Hankel Frobenius energy using window multiplicities without forming the matrix."""
    columns = signal.size - delay + 1
    index = np.arange(signal.size)
    weights = np.minimum.reduce(
        (
            index + 1,
            np.full(signal.size, delay),
            np.full(signal.size, columns),
            signal.size - index,
        )
    )
    return float(np.dot(weights, signal * signal))


def fit_top_two(signal: np.ndarray, delay: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit two Hankel singular modes and their fraction of the full window energy."""
    operator = HankelOperator(signal, delay)
    left, singular, right = svds(
        operator,
        k=2,
        which="LM",
        solver="arpack",
        tol=1.0e-10,
        maxiter=10000,
        return_singular_vectors=True,
        rng=np.random.default_rng(20261001),
    )
    order = np.argsort(singular)[::-1]
    left = left[:, order]
    singular = singular[order]
    right = right[order]
    for mode in range(2):
        pivot = int(np.argmax(np.abs(left[:, mode])))
        if left[pivot, mode] < 0.0:
            left[:, mode] *= -1.0
            right[mode] *= -1.0
    energy = float(np.sum(singular * singular) / hankel_total_energy(signal, delay))
    return left, singular, energy


def fixed_projection(
    echo: np.ndarray,
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
) -> np.ndarray:
    """Project echo windows with the frozen normalization, two modes, and singular scales."""
    normalized = (echo[::EMBEDDING_DECIMATION] - mean) / scale
    coefficients = np.column_stack(
        [
            fftconvolve(normalized, modes[:, mode][::-1], mode="valid")
            for mode in range(2)
        ]
    )
    return coefficients / singular


def phase_tile(title: str, traces: list[dict]) -> dict:
    """Compose a square Viewer phase-plane tile from the supplied curves."""
    tile = figure("phase", title, "xi_1", "xi_2", traces)
    tile["square"] = True
    return tile


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    with np.load(TRAINING, allow_pickle=False) as training:
        training_echo = np.asarray(training["echo"], dtype=float)
        training_radius = np.asarray(training["radius"], dtype=float)
        training_velocity = np.asarray(training["radial_velocity"], dtype=float)
        training_metadata = json.loads(str(training["metadata_json"]))
    mean = float(np.mean(training_echo[::EMBEDDING_DECIMATION]))
    scale = float(np.std(training_echo[::EMBEDDING_DECIMATION]))
    normalized_training = (
        training_echo[::EMBEDDING_DECIMATION] - mean
    ) / scale
    period = float(training_metadata["periodicity"]["period_tau"])
    sample_step = EXPECTED_STEP * EMBEDDING_DECIMATION

    windows = sorted(
        set([BASELINE_WINDOW_TAU] + [period * value for value in PERIOD_MULTIPLIERS])
    )
    delays = [int(round(window / sample_step)) + 1 for window in windows]

    sources = sorted((ROOT / "archives/solutions/rpe").glob(SOURCE_GLOB))
    validations: list[dict] = []
    for source, expected_ic in zip(sources, EXPECTED_ICS, strict=True):
        time, state, metadata = load_native_source(source)
        validate_validation_source(metadata)
        if not np.array_equal(np.asarray(metadata["initial_state"]), expected_ic):
            raise RuntimeError("validation IC order changed")
        validations.append(
            {"time": time, "state": state, "echo": far_field_echo(state, metadata)}
        )

    state_scale = np.asarray(
        [
            float(np.ptp(training_radius)),
            2.0 * float(np.max(np.abs(training_velocity))),
        ]
    )
    regimes = {1: "inside", 2: "inside", 3: "inside", 4: "outside", 5: "outside"}
    rows: list[dict] = []
    cycles_by_delay: dict[int, list[np.ndarray]] = {}
    saved: dict[str, np.ndarray] = {}

    for delay, requested_window in zip(delays, windows, strict=True):
        actual_window = (delay - 1) * sample_step
        modes, singular, training_energy = fit_top_two(normalized_training, delay)
        cycles: list[np.ndarray] = []
        cycle_states: list[np.ndarray] = []
        closures: dict[str, list[float]] = {"inside": [], "outside": []}
        for validation in validations:
            coordinates = fixed_projection(
                validation["echo"], modes, singular, mean, scale
            )
            cycle, cycle_state = representative_cycle(
                validation["echo"],
                coordinates,
                validation["state"],
                delay,
                sample_step,
            )
            cycles.append(cycle)
            cycle_states.append(cycle_state)
            closure, _ = closure_metric(
                validation["echo"], coordinates, delay, sample_step
            )
            regime = regimes[len(cycles)]
            if closure is not None:
                closures[regime].append(closure)

        all_latent = np.vstack(cycles)
        latent_scale = np.ones(6)
        latent_scale[:2] = np.ptp(all_latent, axis=0)
        latent_scale[latent_scale <= np.finfo(float).eps] = 1.0
        padded_cycles = [
            np.column_stack((cycle, np.zeros((cycle.shape[0], 4))))
            for cycle in cycles
        ]
        counts = {"inside_inside": 0, "inside_outside": 0, "outside_outside": 0}
        ambiguities: list[float] = []
        for left_index, right_index in itertools.combinations(range(5), 2):
            crossings = proper_intersections(
                padded_cycles[left_index],
                cycle_states[left_index],
                padded_cycles[right_index],
                cycle_states[right_index],
                latent_scale,
                state_scale,
            )
            left_regime = regimes[left_index + 1]
            right_regime = regimes[right_index + 1]
            category = (
                "inside_inside"
                if left_regime == right_regime == "inside"
                else "outside_outside"
                if left_regime == right_regime == "outside"
                else "inside_outside"
            )
            counts[category] += len(crossings)
            ambiguities.extend(
                crossing["normalized_state_separation"] for crossing in crossings
            )

        row = {
            "requested_window_tau": requested_window,
            "delay_samples": delay,
            "delay_window_tau": actual_window,
            "delay_periods": actual_window / period,
            "training_first_two_energy": training_energy,
            "proper_intersections": {**counts, "total": int(sum(counts.values()))},
            "inside_median_cycle_closure": float(np.median(closures["inside"])),
            "outside_median_cycle_closure": float(np.median(closures["outside"])),
            "intersection_state_ambiguity_median": (
                float(np.median(ambiguities)) if ambiguities else 0.0
            ),
            "intersection_state_ambiguity_maximum": (
                float(np.max(ambiguities)) if ambiguities else 0.0
            ),
        }
        rows.append(row)
        cycles_by_delay[delay] = cycles
        saved[f"delay_{delay}_modes"] = modes
        saved[f"delay_{delay}_singular"] = singular
        print(
            f"{row['delay_periods']:.3f} periods, delay={delay}, "
            f"intersections={row['proper_intersections']['total']}, "
            f"closure in/out={row['inside_median_cycle_closure']:.4g}/"
            f"{row['outside_median_cycle_closure']:.4g}",
            flush=True,
        )

    best = min(
        rows,
        key=lambda item: (
            item["proper_intersections"]["total"],
            item["outside_median_cycle_closure"],
            item["delay_samples"],
        ),
    )
    baseline = min(rows, key=lambda item: abs(item["delay_window_tau"] - BASELINE_WINDOW_TAU))
    metrics = {
        "schema": "dynsyslab.echo-delay-injectivity-sweep/v1",
        "training_source_count": 1,
        "training_policy": "unforced-only",
        "validation_source_count": 5,
        "embedding_sample_step_tau": sample_step,
        "training_period_tau": period,
        "normalization": {"mean": mean, "scale": scale, "source": "training only"},
        "rows": rows,
        "baseline": baseline,
        "best_by_intersection_then_outside_closure": best,
    }
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        delay_samples=np.asarray(delays),
        delay_window_tau=np.asarray([row["delay_window_tau"] for row in rows]),
        delay_periods=np.asarray([row["delay_periods"] for row in rows]),
        intersection_total=np.asarray(
            [row["proper_intersections"]["total"] for row in rows]
        ),
        metadata_json=np.asarray(json.dumps(metrics, separators=(",", ":"))),
        **saved,
    )

    x = np.asarray([row["delay_periods"] for row in rows])
    summary = [
        figure(
            "line",
            "Top-two coordinate intersections versus delay length",
            "delay window / training period",
            "proper intersections",
            [
                trace(
                    label,
                    x,
                    np.asarray([row["proper_intersections"][key] for row in rows]),
                    color=color,
                    marker=True,
                    marker_size=4.0,
                    width=2.0,
                )
                for label, key, color in (
                    ("inside-inside", "inside_inside", "#2563eb"),
                    ("inside-outside", "inside_outside", "#7c3aed"),
                    ("outside-outside", "outside_outside", "#dc2626"),
                    ("total", "total", "#111827"),
                )
            ],
        ),
        figure(
            "line",
            "Latent cycle closure versus delay length",
            "delay window / training period",
            "median relative closure error",
            [
                trace(
                    "inside",
                    x,
                    np.asarray([row["inside_median_cycle_closure"] for row in rows]),
                    color="#2563eb",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "outside",
                    x,
                    np.asarray([row["outside_median_cycle_closure"] for row in rows]),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
        figure(
            "line",
            "Training energy in the top two coordinates",
            "delay window / training period",
            "energy fraction",
            [
                trace(
                    "top two",
                    x,
                    np.asarray([row["training_first_two_energy"] for row in rows]),
                    color="#0891b2",
                    marker=True,
                    width=2.0,
                )
            ],
        ),
        figure(
            "line",
            "Physical ambiguity at detected intersections",
            "delay window / training period",
            "normalized state separation",
            [
                trace(
                    "median",
                    x,
                    np.asarray(
                        [row["intersection_state_ambiguity_median"] for row in rows]
                    ),
                    color="#f59e0b",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "maximum",
                    x,
                    np.asarray(
                        [row["intersection_state_ambiguity_maximum"] for row in rows]
                    ),
                    color="#dc2626",
                    marker=True,
                    width=2.0,
                ),
            ],
        ),
    ]
    render_viewer(
        TRAINING,
        summary,
        name=SUMMARY_PLOT,
        label="Delay-window injectivity sweep",
        columns=2,
        width=1900,
        height=1900,
    )

    selected_rows = []
    for target in (0.25, baseline["delay_periods"], 2.0, best["delay_periods"]):
        candidate = min(rows, key=lambda item: abs(item["delay_periods"] - target))
        if candidate not in selected_rows:
            selected_rows.append(candidate)
    phase_figures: list[dict] = []
    for row in selected_rows:
        traces = [
            trace(
                f"IC #{index}",
                cycle[:, 0],
                cycle[:, 1],
                color=color,
                width=2.0,
            )
            for index, (cycle, color) in enumerate(
                zip(cycles_by_delay[row["delay_samples"]], COLORS, strict=True),
                start=1,
            )
        ]
        phase_figures.append(
            phase_tile(
                f"{row['delay_periods']:.2f} periods — "
                f"{row['proper_intersections']['total']} intersections",
                traces,
            )
        )
    render_viewer(
        TRAINING,
        phase_figures,
        name=PHASE_PLOT,
        label="Delay-window latent phase comparison",
        columns=2,
        width=1900,
        height=1900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
