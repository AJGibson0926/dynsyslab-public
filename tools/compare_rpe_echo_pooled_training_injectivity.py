"""Compare one-orbit and pooled-orbit echo Hankel bases for injectivity.

The pooled decomposition is formed as D = [D_1 D_2 ... D_m], where every
D_i is built independently from one unforced echo record.  A delay window is
never permitted to cross an initial-condition boundary.  The resulting fixed
U basis is then evaluated on the same five held-out trajectories used by the
single-trajectory delay sweep.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import LinearOperator, svds

from analyze_rpe_echo_hankel_latents import OUTPUT_ROOT, figure, render_viewer, trace
from analyze_rpe_inviscid_echo_hankel import far_field_echo, load_native_source
from audit_rpe_echo_latent_injectivity import proper_intersections, representative_cycle
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
from sweep_rpe_echo_delay_injectivity import (
    HankelOperator,
    fixed_projection,
    hankel_total_energy,
)


ROOT = Path(__file__).resolve().parents[1]
NEW_TRAINING_GLOB = (
    "rpe-inviscid-echo-pooled-training-four-ics-trajectory-*"
    "__20261001T232136Z.npz"
)
SINGLE_METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-delay-injectivity-sweep.metrics.json"
METRICS = OUTPUT_ROOT / "rpe-inviscid-echo-pooled-training-comparison.metrics.json"
ARTIFACT = OUTPUT_ROOT / "rpe-inviscid-echo-pooled-training-comparison.npz"
SUMMARY_PLOT = "rpe-inviscid-echo-pooled-training-comparison"
PHASE_PLOT = "rpe-inviscid-echo-pooled-training-phase-comparison"
TRAINING_DURATION_TAU = 12.0


class PooledHankelOperator(LinearOperator):
    """Horizontal concatenation of per-trajectory Hankel operators."""

    def __init__(self, signals: list[np.ndarray], delay: int):
        self.blocks = [HankelOperator(signal, delay) for signal in signals]
        self.column_counts = [block.shape[1] for block in self.blocks]
        super().__init__(
            dtype=np.dtype(float),
            shape=(delay, int(sum(self.column_counts))),
        )

    def _matvec(self, vector: np.ndarray) -> np.ndarray:
        vector = np.asarray(vector).reshape(-1)
        result = np.zeros(self.shape[0], dtype=float)
        start = 0
        for block, count in zip(self.blocks, self.column_counts, strict=True):
            result += block @ vector[start : start + count]
            start += count
        return result

    def _rmatvec(self, vector: np.ndarray) -> np.ndarray:
        vector = np.asarray(vector).reshape(-1)
        return np.concatenate([block.T @ vector for block in self.blocks])

    def _matmat(self, matrix: np.ndarray) -> np.ndarray:
        return np.column_stack(
            [self._matvec(matrix[:, column]) for column in range(matrix.shape[1])]
        )

    def _rmatmat(self, matrix: np.ndarray) -> np.ndarray:
        return np.column_stack(
            [self._rmatvec(matrix[:, column]) for column in range(matrix.shape[1])]
        )


def fit_pooled_top_two(
    signals: list[np.ndarray], delay: int
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit two modes of the concatenated trajectory Hankel operator and report their energy."""
    operator = PooledHankelOperator(signals, delay)
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
    total_energy = sum(hankel_total_energy(signal, delay) for signal in signals)
    energy = float(np.sum(singular * singular) / total_energy)
    return left, singular, energy


def load_training_records() -> tuple[list[dict], float]:
    """Load the retained free training trajectories and their reference period."""
    records: list[dict] = []
    with np.load(TRAINING, allow_pickle=False) as archive:
        time = np.asarray(archive["time"], dtype=float)
        mask = time <= TRAINING_DURATION_TAU + 1.0e-12
        metadata = json.loads(str(archive["metadata_json"]))
        records.append(
            {
                "source": str(metadata["source"]),
                "initial_state": np.asarray(metadata["initial_state"], dtype=float),
                "time": time[mask],
                "state": np.column_stack(
                    (
                        np.asarray(archive["radius"], dtype=float)[mask],
                        np.asarray(archive["radial_velocity"], dtype=float)[mask],
                    )
                ),
                "echo": np.asarray(archive["echo"], dtype=float)[mask],
                "period_tau": float(metadata["periodicity"]["period_tau"]),
            }
        )

    sources = sorted((ROOT / "archives/solutions/rpe").glob(NEW_TRAINING_GLOB))
    if len(sources) != 4:
        raise RuntimeError(f"expected four new training sources, found {len(sources)}")
    for source in sources:
        time, state, metadata = load_native_source(source)
        validate_validation_source(metadata)
        if not np.isclose(time[-1], TRAINING_DURATION_TAU, atol=1.0e-12):
            raise RuntimeError("pooled training records must have equal duration")
        records.append(
            {
                "source": str(source),
                "initial_state": np.asarray(metadata["initial_state"], dtype=float),
                "time": time,
                "state": state,
                "echo": far_field_echo(state, metadata),
                "period_tau": None,
            }
        )
    return records, float(records[0]["period_tau"])


def load_validations() -> list[dict]:
    """Load and check the five held-out trajectories and expected initial conditions."""
    sources = sorted((ROOT / "archives/solutions/rpe").glob(SOURCE_GLOB))
    if len(sources) != 5:
        raise RuntimeError(f"expected five validation sources, found {len(sources)}")
    validations: list[dict] = []
    for source, expected_ic in zip(sources, EXPECTED_ICS, strict=True):
        time, state, metadata = load_native_source(source)
        validate_validation_source(metadata)
        if not np.array_equal(np.asarray(metadata["initial_state"]), expected_ic):
            raise RuntimeError("validation IC order changed")
        validations.append(
            {
                "source": str(source),
                "time": time,
                "state": state,
                "echo": far_field_echo(state, metadata),
            }
        )
    return validations


def measure_basis(
    validations: list[dict],
    modes: np.ndarray,
    singular: np.ndarray,
    mean: float,
    scale: float,
    state_scale: np.ndarray,
) -> tuple[dict, list[np.ndarray]]:
    """Compare frozen-basis orbit closure and geometric separation inside and outside training."""
    delay = modes.shape[0]
    sample_step = EXPECTED_STEP * EMBEDDING_DECIMATION
    regimes = {1: "inside", 2: "inside", 3: "inside", 4: "outside", 5: "outside"}
    cycles: list[np.ndarray] = []
    cycle_states: list[np.ndarray] = []
    closures: dict[str, list[float]] = {"inside": [], "outside": []}
    for index, validation in enumerate(validations, start=1):
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
        if closure is not None:
            closures[regimes[index]].append(closure)

    all_latent = np.vstack(cycles)
    latent_scale = np.ones(6)
    latent_scale[:2] = np.ptp(all_latent, axis=0)
    latent_scale[latent_scale <= np.finfo(float).eps] = 1.0
    padded_cycles = [
        np.column_stack((cycle, np.zeros((cycle.shape[0], 4)))) for cycle in cycles
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
    return (
        {
            "proper_intersections": {**counts, "total": int(sum(counts.values()))},
            "inside_median_cycle_closure": float(np.median(closures["inside"])),
            "outside_median_cycle_closure": float(np.median(closures["outside"])),
            "intersection_state_ambiguity_median": (
                float(np.median(ambiguities)) if ambiguities else 0.0
            ),
            "intersection_state_ambiguity_maximum": (
                float(np.max(ambiguities)) if ambiguities else 0.0
            ),
        },
        cycles,
    )


def phase_tile(title: str, cycles: list[np.ndarray]) -> dict:
    """Compose a square Viewer phase-plane tile from the supplied curves."""
    tile = figure(
        "phase",
        title,
        "xi_1",
        "xi_2",
        [
            trace(
                f"IC #{index}",
                cycle[:, 0],
                cycle[:, 1],
                color=color,
                width=2.0,
            )
            for index, (cycle, color) in enumerate(
                zip(cycles, COLORS, strict=True), start=1
            )
        ],
    )
    tile["square"] = True
    return tile


def main() -> None:
    """Run this isolated research experiment and save its local metrics and Viewer figures."""
    training_records, reference_period = load_training_records()
    validations = load_validations()
    sampled_echoes = [
        np.asarray(record["echo"], dtype=float)[::EMBEDDING_DECIMATION]
        for record in training_records
    ]
    if len({echo.size for echo in sampled_echoes}) != 1:
        raise RuntimeError("equal-record pooling is required for balanced training")
    combined_echo = np.concatenate(sampled_echoes)
    mean = float(np.mean(combined_echo))
    scale = float(np.std(combined_echo))
    normalized = [(echo - mean) / scale for echo in sampled_echoes]
    sample_step = EXPECTED_STEP * EMBEDDING_DECIMATION

    with SINGLE_METRICS.open("r", encoding="utf-8") as stream:
        single_metrics = json.load(stream)
    single_rows = single_metrics["rows"]
    delays = [int(row["delay_samples"]) for row in single_rows]
    state_stack = np.vstack([record["state"] for record in training_records])
    state_scale = np.asarray(
        [
            float(np.ptp(state_stack[:, 0])),
            2.0 * float(np.max(np.abs(state_stack[:, 1]))),
        ]
    )

    rows: list[dict] = []
    cycles_by_delay: dict[int, list[np.ndarray]] = {}
    saved: dict[str, np.ndarray] = {}
    for delay in delays:
        modes, singular, energy = fit_pooled_top_two(normalized, delay)
        measured, cycles = measure_basis(
            validations, modes, singular, mean, scale, state_scale
        )
        row = {
            "delay_samples": delay,
            "delay_window_tau": (delay - 1) * sample_step,
            "delay_periods": ((delay - 1) * sample_step) / reference_period,
            "training_first_two_energy": energy,
            **measured,
        }
        rows.append(row)
        cycles_by_delay[delay] = cycles
        saved[f"delay_{delay}_modes"] = modes
        saved[f"delay_{delay}_singular"] = singular
        print(
            f"pooled: {row['delay_periods']:.3f} periods, delay={delay}, "
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
    baseline = min(rows, key=lambda item: abs(item["delay_window_tau"] - 1.275))
    training_sources = [
        {
            "source": record["source"],
            "initial_state": record["initial_state"].tolist(),
            "duration_tau": float(record["time"][-1] - record["time"][0]),
        }
        for record in training_records
    ]
    metrics = {
        "schema": "dynsyslab.echo-pooled-training-injectivity/v1",
        "training_policy": "unforced-only",
        "training_source_count": len(training_records),
        "training_sources": training_sources,
        "pooling": {
            "matrix": "horizontal concatenation of independent per-run Hankel blocks",
            "cross_boundary_windows": False,
            "equal_duration_per_source_tau": TRAINING_DURATION_TAU,
            "normalization": "one mean and standard deviation from all training echoes",
            "mean": mean,
            "scale": scale,
        },
        "validation_source_count": len(validations),
        "embedding_sample_step_tau": sample_step,
        "reference_period_tau": reference_period,
        "single_training_rows": single_rows,
        "pooled_training_rows": rows,
        "pooled_baseline": baseline,
        "pooled_best_by_intersection_then_outside_closure": best,
    }
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    np.savez_compressed(
        ARTIFACT,
        training_initial_states=np.vstack(
            [record["initial_state"] for record in training_records]
        ),
        delay_samples=np.asarray(delays),
        delay_periods=np.asarray([row["delay_periods"] for row in rows]),
        pooled_intersection_total=np.asarray(
            [row["proper_intersections"]["total"] for row in rows]
        ),
        metadata_json=np.asarray(json.dumps(metrics, separators=(",", ":"))),
        **saved,
    )

    x_single = np.asarray([row["delay_periods"] for row in single_rows])
    x_pooled = np.asarray([row["delay_periods"] for row in rows])
    summary = [
        figure(
            "line",
            "Top-two intersections: one orbit versus five",
            "delay window / reference period",
            "proper intersections",
            [
                trace(
                    "one training orbit",
                    x_single,
                    np.asarray([row["proper_intersections"]["total"] for row in single_rows]),
                    color="#6b7280",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "five training orbits",
                    x_pooled,
                    np.asarray([row["proper_intersections"]["total"] for row in rows]),
                    color="#0891b2",
                    marker=True,
                    width=2.4,
                ),
            ],
        ),
        figure(
            "line",
            "Held-out cycle closure",
            "delay window / reference period",
            "median relative closure error",
            [
                trace(
                    "one orbit — inside",
                    x_single,
                    np.asarray([row["inside_median_cycle_closure"] for row in single_rows]),
                    color="#93c5fd",
                    marker=True,
                    width=1.8,
                ),
                trace(
                    "one orbit — outside",
                    x_single,
                    np.asarray([row["outside_median_cycle_closure"] for row in single_rows]),
                    color="#fca5a5",
                    marker=True,
                    width=1.8,
                ),
                trace(
                    "five orbits — inside",
                    x_pooled,
                    np.asarray([row["inside_median_cycle_closure"] for row in rows]),
                    color="#2563eb",
                    marker=True,
                    width=2.4,
                ),
                trace(
                    "five orbits — outside",
                    x_pooled,
                    np.asarray([row["outside_median_cycle_closure"] for row in rows]),
                    color="#dc2626",
                    marker=True,
                    width=2.4,
                ),
            ],
        ),
        figure(
            "line",
            "Training energy in the top two coordinates",
            "delay window / reference period",
            "energy fraction",
            [
                trace(
                    "one training orbit",
                    x_single,
                    np.asarray([row["training_first_two_energy"] for row in single_rows]),
                    color="#6b7280",
                    marker=True,
                    width=2.0,
                ),
                trace(
                    "five training orbits",
                    x_pooled,
                    np.asarray([row["training_first_two_energy"] for row in rows]),
                    color="#7c3aed",
                    marker=True,
                    width=2.4,
                ),
            ],
        ),
        phase_tile(
            f"Pooled best: {best['delay_periods']:.2f} periods, "
            f"{best['proper_intersections']['total']} intersections",
            cycles_by_delay[best["delay_samples"]],
        ),
    ]
    render_viewer(
        TRAINING,
        summary,
        name=SUMMARY_PLOT,
        label="Pooled-training injectivity comparison",
        columns=2,
        width=1900,
        height=1900,
    )

    selected = []
    for target in (baseline["delay_periods"], 2.0, 3.0, 4.0):
        candidate = min(rows, key=lambda item: abs(item["delay_periods"] - target))
        if candidate not in selected:
            selected.append(candidate)
    render_viewer(
        TRAINING,
        [
            phase_tile(
                f"Five-orbit training — {row['delay_periods']:.2f} periods — "
                f"{row['proper_intersections']['total']} intersections",
                cycles_by_delay[row["delay_samples"]],
            )
            for row in selected
        ],
        name=PHASE_PLOT,
        label="Pooled-training latent phase comparison",
        columns=2,
        width=1900,
        height=1900,
    )
    print(METRICS)
    print(ARTIFACT)


if __name__ == "__main__":
    main()
