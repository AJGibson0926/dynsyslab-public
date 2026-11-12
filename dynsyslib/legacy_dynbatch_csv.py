"""Read the legacy MATLAB ``data_node.save`` CSV hierarchy.

The historical MATLAB suite serializes a tree as tagged CSV records rather
than as one rectangular table.  This module decodes that format without
calling MATLAB's interactive ``data_node.load`` method.  It intentionally
keeps parsing separate from dynamical-system construction: imported numbers
are authoritative data, while mapping a legacy dynamics name to an executable
DynSysLab model is an explicit downstream decision.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from infrastructure.controlled_data import ControlledDataset, ControlledTrajectory
from time_series import TimeSeries


LEGACY_TAGS = frozenset(
    {"dynbatch", "dynsys", "dynclst", "inst", "time", "sol", "ctrl", "u"}
)


def _numeric_prefix(row: list[str], *, line: int, record: str) -> np.ndarray:
    values: list[float] = []
    encountered_empty = False
    for column, raw in enumerate(row[1:], start=2):
        token = raw.strip()
        if not token:
            encountered_empty = True
            continue
        if encountered_empty:
            raise ValueError(
                f"{record} record at line {line} has data after an empty column ({column})."
            )
        try:
            value = float(token)
        except ValueError as error:
            raise ValueError(
                f"{record} record at line {line}, column {column} is not numeric: {token!r}."
            ) from error
        if not np.isfinite(value):
            raise ValueError(f"{record} record at line {line} contains a non-finite value.")
        values.append(value)
    return np.asarray(values, dtype=float)


def _text_prefix(row: list[str]) -> tuple[str, ...]:
    values: list[str] = []
    for raw in row[1:]:
        token = raw.strip()
        if not token:
            break
        values.append(token)
    return tuple(values)


@dataclass(frozen=True)
class LegacyControlHistory:
    """Sample-aligned control values stored at one legacy ``dynsys`` node."""

    time_series: TimeSeries
    forcing: str = "off"
    subforcing: str = "off"

    @property
    def values(self) -> np.ndarray:
        return self.time_series.values

    @property
    def interval_values(self) -> np.ndarray:
        """Return the 200 interval commands after the legacy zero initializer.

        ``CTRL.hist`` was initialized with a zero sample before the 200 random
        commands were appended.  MATLAB's multistep constant-B regression uses
        ``hist.ts(2:end)``; this property gives that exact alignment.
        """

        if self.time_series.sample_count < 2:
            raise ValueError("A legacy control history needs an initializer and an interval sample.")
        return np.asarray(self.time_series.values[1:], dtype=float)


@dataclass(frozen=True)
class LegacyTrajectory:
    """One decoded legacy ``inst`` record and its state history."""

    system_index: int
    cluster_index: int
    trajectory_index: int
    time_series: TimeSeries
    source_lines: tuple[int, int]
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def initial_condition(self) -> np.ndarray:
        return np.asarray(self.time_series.values[0], dtype=float)


@dataclass(frozen=True)
class LegacySystem:
    """One decoded legacy ``dynsys`` record."""

    index: int
    dynamics_name: str
    delta_t_declared: float
    transition_count_declared: int
    trailing_fields: tuple[str, ...]
    time: np.ndarray
    trajectories: tuple[LegacyTrajectory, ...]
    control_history: LegacyControlHistory | None
    source_line: int

    def __post_init__(self):
        if not self.trajectories:
            raise ValueError(f"Legacy system {self.index} contains no trajectories.")
        if self.delta_t_declared <= 0.0 or self.transition_count_declared < 1:
            raise ValueError(f"Legacy system {self.index} has invalid declared time geometry.")
        if self.time.size and self.time.size != self.transition_count_declared + 1:
            raise ValueError(
                f"Legacy system {self.index} declares {self.transition_count_declared} transitions "
                f"but stores {self.time.size} time samples."
            )


@dataclass(frozen=True)
class LegacyDynBatchCSV:
    """Decoded, read-only representation of a legacy MATLAB CSV hierarchy."""

    source_path: Path
    source_sha256: str
    header: tuple[str, ...]
    declared_system_count: int | None
    systems: tuple[LegacySystem, ...]
    warnings: tuple[str, ...] = ()

    @property
    def trajectory_count(self) -> int:
        return sum(len(system.trajectories) for system in self.systems)

    def summary(self) -> dict[str, Any]:
        state_dimensions = sorted(
            {trajectory.time_series.values.shape[1] for system in self.systems for trajectory in system.trajectories}
        )
        dynamics = sorted({system.dynamics_name for system in self.systems})
        controlled = sum(system.control_history is not None for system in self.systems)
        return {
            "schema": "dynsyslab.legacy-matlab-dynbatch-csv/v1",
            "source_path": str(self.source_path),
            "source_sha256": self.source_sha256,
            "declared_system_count": self.declared_system_count,
            "decoded_system_count": len(self.systems),
            "trajectory_count": self.trajectory_count,
            "controlled_system_count": controlled,
            "dynamics_names": dynamics,
            "state_dimensions": state_dimensions,
            "warnings": list(self.warnings),
        }

    def pooled_time_series(self) -> TimeSeries:
        """Pool every trajectory when all decoded members share one mesh."""

        trajectories = [item for system in self.systems for item in system.trajectories]
        reference = trajectories[0].time_series.time
        shape = trajectories[0].time_series.values.shape
        for trajectory in trajectories[1:]:
            if trajectory.time_series.values.shape != shape or not np.array_equal(
                trajectory.time_series.time, reference
            ):
                raise ValueError("Legacy trajectories do not share one poolable time/state geometry.")
        values = np.stack([trajectory.time_series.values for trajectory in trajectories], axis=1)
        return TimeSeries(reference, values, pooled_members=len(trajectories))

    def controlled_dataset(self, *, input_scale: float = 1.0):
        """Convert controlled systems to the canonical regression container.

        Controls use MATLAB's multistep-regression alignment: the stored zero
        initializer is removed and samples 2:end become the interval inputs.
        ``input_scale`` is explicit because some external examples encode a
        noncanonical input gain in their equations of motion.
        """

        if not np.isfinite(input_scale):
            raise ValueError("input_scale must be finite.")
        converted = []
        for system in self.systems:
            history = system.control_history
            if history is None:
                raise ValueError(f"Legacy system {system.index} has no control history.")
            controls = history.interval_values * float(input_scale)
            for trajectory in system.trajectories:
                if controls.shape[0] != trajectory.time_series.sample_count - 1:
                    raise ValueError(
                        f"Legacy system {system.index} control/state histories are not interval aligned."
                    )
                converted.append(
                    ControlledTrajectory(
                        trajectory.time_series.time,
                        trajectory.time_series.values,
                        controls,
                        {
                            "source_path": str(self.source_path),
                            "source_sha256": self.source_sha256,
                            "legacy_system_index": system.index,
                            "legacy_cluster_index": trajectory.cluster_index,
                            "legacy_trajectory_index": trajectory.trajectory_index,
                            "legacy_forcing": history.forcing,
                            "legacy_subforcing": history.subforcing,
                            "control_alignment": "stored samples 2:end are interval inputs",
                            "input_scale": float(input_scale),
                        },
                    )
                )
        return ControlledDataset(
            converted,
            {
                **self.summary(),
                "control_alignment": "legacy CTRL.hist zero initializer removed",
                "input_scale": float(input_scale),
            },
        )


@dataclass
class _SystemBuilder:
    index: int
    dynamics_name: str
    delta_t: float
    transition_count: int
    trailing_fields: tuple[str, ...]
    source_line: int
    time: np.ndarray | None = None
    forcing: str = "off"
    subforcing: str = "off"
    control_values: np.ndarray | None = None
    control_time: np.ndarray | None = None
    cluster_index: int | None = None
    trajectory_index: int | None = None
    trajectory_time: np.ndarray | None = None
    trajectories: list[LegacyTrajectory] = field(default_factory=list)

    def finish(self) -> LegacySystem:
        time = np.asarray(self.time if self.time is not None else [], dtype=float)
        control = None
        if self.control_values is not None:
            control_time = self.control_time if self.control_time is not None else time
            if control_time is None or len(control_time) != len(self.control_values):
                raise ValueError(f"Legacy system {self.index} has mismatched control time and values.")
            values = np.asarray(self.control_values, dtype=float)
            if values.ndim == 1:
                values = values[:, None]
            control = LegacyControlHistory(
                TimeSeries(np.asarray(control_time, dtype=float), values), self.forcing, self.subforcing
            )
        return LegacySystem(
            self.index,
            self.dynamics_name,
            self.delta_t,
            self.transition_count,
            self.trailing_fields,
            time,
            tuple(self.trajectories),
            control,
            self.source_line,
        )


def read_legacy_dynbatch_csv(path: str | Path, *, strict: bool = True) -> LegacyDynBatchCSV:
    """Decode a legacy tagged CSV without invoking MATLAB or any file picker."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    digest = sha256(source.read_bytes()).hexdigest()
    with source.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream))
    if len(rows) < 2:
        raise ValueError("Legacy CSV is empty or contains only a header.")
    header = tuple(rows[0])
    declared_system_count: int | None = None
    warnings: list[str] = []
    systems: list[LegacySystem] = []
    current: _SystemBuilder | None = None
    index = 1
    while index < len(rows):
        row = rows[index]
        line = index + 1
        tag = row[0].strip() if row else ""
        if not tag:
            if strict:
                raise ValueError(f"Unexpected continuation row at line {line}.")
            warnings.append(f"Ignored unattached continuation row {line}.")
            index += 1
            continue
        if tag not in LEGACY_TAGS:
            message = f"Unknown legacy record tag {tag!r} at line {line}."
            if strict:
                raise ValueError(message)
            warnings.append(message)
            index += 1
            continue
        if tag == "dynbatch":
            fields = _text_prefix(row)
            if fields:
                declared_system_count = int(float(fields[0]))
        elif tag == "dynsys":
            if current is not None:
                systems.append(current.finish())
            fields = _text_prefix(row)
            if len(fields) < 4:
                raise ValueError(f"dynsys record at line {line} needs index, dynamics, dt, and transitions.")
            current = _SystemBuilder(
                int(float(fields[0])), fields[1], float(fields[2]), int(float(fields[3])), tuple(fields[4:]), line
            )
        else:
            if current is None:
                raise ValueError(f"{tag} record at line {line} precedes the first dynsys record.")
            if tag == "dynclst":
                fields = _text_prefix(row)
                if not fields:
                    raise ValueError(f"dynclst record at line {line} has no index.")
                current.cluster_index = int(float(fields[0]))
                current.trajectory_index = None
                current.trajectory_time = None
            elif tag == "inst":
                fields = _text_prefix(row)
                if current.cluster_index is None or not fields:
                    raise ValueError(f"inst record at line {line} has no active cluster/index.")
                current.trajectory_index = int(float(fields[0]))
                current.trajectory_time = None
            elif tag == "time":
                values = _numeric_prefix(row, line=line, record=tag)
                if values.size < 2 or np.any(np.diff(values) <= 0.0):
                    raise ValueError(f"time record at line {line} must be strictly increasing.")
                if current.trajectory_index is None:
                    current.time = values
                else:
                    current.trajectory_time = values
            elif tag == "ctrl":
                fields = _text_prefix(row)
                current.forcing = fields[0] if fields else "off"
                current.subforcing = fields[1] if len(fields) > 1 else "off"
            elif tag == "u":
                current.control_values = _numeric_prefix(row, line=line, record=tag)
                current.control_time = None if current.time is None else current.time.copy()
            elif tag == "sol":
                if current.cluster_index is None or current.trajectory_index is None:
                    raise ValueError(f"sol record at line {line} has no active cluster/trajectory.")
                component_rows = [_numeric_prefix(row, line=line, record=tag)]
                end_line = line
                lookahead = index + 1
                while lookahead < len(rows):
                    next_row = rows[lookahead]
                    next_tag = next_row[0].strip() if next_row else ""
                    if next_tag:
                        break
                    component_rows.append(
                        _numeric_prefix(next_row, line=lookahead + 1, record="sol-continuation")
                    )
                    end_line = lookahead + 1
                    lookahead += 1
                widths = {values.size for values in component_rows}
                if len(widths) != 1 or next(iter(widths)) < 2:
                    raise ValueError(f"sol block at lines {line}-{end_line} has inconsistent component lengths.")
                trajectory_time = current.trajectory_time
                if trajectory_time is None:
                    trajectory_time = current.time
                if trajectory_time is None or trajectory_time.size != next(iter(widths)):
                    raise ValueError(f"sol block at lines {line}-{end_line} has no matching time record.")
                states = np.vstack(component_rows).T
                current.trajectories.append(
                    LegacyTrajectory(
                        current.index,
                        current.cluster_index,
                        current.trajectory_index,
                        TimeSeries(trajectory_time, states),
                        (line, end_line),
                        {"dynamics_name": current.dynamics_name},
                    )
                )
                index = lookahead - 1
        index += 1
    if current is not None:
        systems.append(current.finish())
    if not systems:
        raise ValueError("No dynsys records were decoded from the legacy CSV.")
    if declared_system_count is not None and declared_system_count != len(systems):
        message = (
            f"dynbatch declares {declared_system_count} systems but {len(systems)} were decoded."
        )
        if strict:
            raise ValueError(message)
        warnings.append(message)
    return LegacyDynBatchCSV(
        source, digest, header, declared_system_count, tuple(systems), tuple(warnings)
    )


def export_legacy_dynbatch_csv(
    dataset: LegacyDynBatchCSV, output: str | Path, *, include_controlled_dataset: bool = False
) -> tuple[Path, Path]:
    """Write a lossless numerical NPZ plus a compact JSON provenance manifest."""

    destination = Path(output).expanduser().resolve()
    if destination.suffix.lower() != ".npz":
        destination = destination.with_suffix(".npz")
    destination.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    for system_position, system in enumerate(dataset.systems):
        prefix = f"system_{system_position:04d}"
        arrays[f"{prefix}_time"] = system.time
        if system.control_history is not None:
            arrays[f"{prefix}_control_time"] = system.control_history.time_series.time
            arrays[f"{prefix}_control_values"] = system.control_history.values
        for trajectory_position, trajectory in enumerate(system.trajectories):
            arrays[f"{prefix}_trajectory_{trajectory_position:04d}_time"] = trajectory.time_series.time
            arrays[f"{prefix}_trajectory_{trajectory_position:04d}_states"] = trajectory.time_series.values
    if include_controlled_dataset:
        controlled = dataset.controlled_dataset()
        arrays["controlled_states"] = np.stack([item.states for item in controlled.trajectories])
        arrays["controlled_interval_inputs"] = np.stack([item.inputs for item in controlled.trajectories])
    np.savez_compressed(destination, **arrays)
    manifest = destination.with_suffix(".json")
    payload = dataset.summary()
    payload["npz_path"] = str(destination)
    payload["npz_sha256"] = sha256(destination.read_bytes()).hexdigest()
    payload["array_keys"] = sorted(arrays)
    manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return destination, manifest


def _main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--allow-unknown-records", action="store_true")
    parser.add_argument("--include-controlled-dataset", action="store_true")
    arguments = parser.parse_args(list(argv) if argv is not None else None)
    dataset = read_legacy_dynbatch_csv(arguments.source, strict=not arguments.allow_unknown_records)
    payload = dataset.summary()
    if arguments.output is not None:
        archive, manifest = export_legacy_dynbatch_csv(
            dataset, arguments.output, include_controlled_dataset=arguments.include_controlled_dataset
        )
        payload.update({"archive": str(archive), "manifest": str(manifest)})
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())


__all__ = [
    "LEGACY_TAGS",
    "LegacyControlHistory",
    "LegacyDynBatchCSV",
    "LegacySystem",
    "LegacyTrajectory",
    "export_legacy_dynbatch_csv",
    "read_legacy_dynbatch_csv",
]
