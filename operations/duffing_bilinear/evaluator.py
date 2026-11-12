"""Per-attempt extraction, fitting, error fields, and control comparison."""
from __future__ import annotations

from pathlib import Path
from time import perf_counter
from typing import Any
import numpy as np

from rom import (KROMCandidate, KROMControlAttempt, KROMControlCampaignPolicy,
    control_approximation_error_fields, fit_bilinear_control_map,
    fit_bilinear_control_matrix, fit_conjugate_bilinear_control_map,
    fit_constant_control_matrix, load_krom_archive,
    resample_trajectory_by_arclength, reserve_control_representations_archive,
    reserve_krom_archive, save_control_approximation_error_fields,
    trajectory_tube_states)
from .common import ACTIVE_WINDOW, _write_json
from .control_cases import _solve_case_bounded
from .context import DuffingSearchContext
from .path_support import prior_successful_exact_path

class DuffingKROMSearch(DuffingSearchContext):
    """Evaluate one scheduled sub-KROM using every retained control-map form."""

    def __call__(self, attempt_number: int, candidate: KROMCandidate) -> KROMControlAttempt:
        started = perf_counter()
        attempt_root = self.output_root / "attempts" / f"attempt-{attempt_number:02d}"
        attempt_root.mkdir(parents=True, exist_ok=True)
        indices = list(candidate.retained_indices)
        label = "-".join(str(item) for item in indices)
        krom_archive = reserve_krom_archive(
            self.parent_archive,
            indices,
            label=f"duffing-bilinear-search-a{attempt_number:02d}-modes-{label}",
            selection_metadata={
                "campaign_attempt": attempt_number,
                "candidate": candidate.summary(),
            },
        )
        model = load_krom_archive(krom_archive)
        domain = model.training_domain
        if domain is None:
            raise ValueError("Campaign KROM lacks a free-training rectangle.")
        lower = np.asarray(domain["minimum"], dtype=float)
        upper = np.asarray(domain["maximum"], dtype=float)
        timings: dict[str, float] = {}

        stage = perf_counter()
        constant = fit_constant_control_matrix(
            model,
            self.controlled_dataset,
            method="multistep-output",
            ridge=1e-12,
        )
        timings["constant_regression_seconds"] = perf_counter() - stage
        stage = perf_counter()
        dynamics = fit_bilinear_control_matrix(
            model, self.controlled_dataset, ridge=1e-12
        )
        timings["bilinear_dynamics_regression_seconds"] = perf_counter() - stage
        stage = perf_counter()
        direct = fit_bilinear_control_map(
            model, self.controlled_states[::5], ridge=1e-12
        )
        timings["bilinear_map_regression_seconds"] = perf_counter() - stage
        stage = perf_counter()
        conjugate_direct = fit_conjugate_bilinear_control_map(
            model, self.controlled_states[::5], ridge=1e-12
        )
        timings["bilinear_conjugate_map_regression_seconds"] = perf_counter() - stage
        corridor_mask = (
            (self.controlled_states[:, 0] >= -0.6)
            & (self.controlled_states[:, 0] <= 0.6)
            & (np.abs(self.controlled_states[:, 1]) <= 0.5)
        )
        corridor_states = self.controlled_states[corridor_mask]
        if corridor_states.shape[0] < max(50, 10 * (model.lifted_dimension + 1)):
            raise ValueError("Controlled dataset has insufficient samples in the target corridor.")
        stage = perf_counter()
        local_direct = fit_bilinear_control_map(
            model, corridor_states, ridge=1e-12
        )
        timings["bilinear_local_map_regression_seconds"] = perf_counter() - stage
        path_conditioned = bool(candidate.scores.get("path_conditioned_candidate_stage", 0.0))
        path_states = tube_states = None
        path_source_archive = None
        path_direct = tube_direct = None
        arc3_states = arc7_states = arc9_states = None
        arc3_direct = arc7_conjugate = arc9_conjugate = None
        if path_conditioned:
            path_states, path_source_archive = prior_successful_exact_path(
                self.output_root, indices
            )
            tube_states = trajectory_tube_states(
                path_states, lower, upper, radii=(0.025, 0.05)
            )
            stage = perf_counter()
            path_direct = fit_bilinear_control_map(model, path_states, ridge=1e-12)
            timings["bilinear_path_map_regression_seconds"] = perf_counter() - stage
            stage = perf_counter()
            tube_direct = fit_bilinear_control_map(model, tube_states, ridge=1e-12)
            timings["bilinear_path_tube_map_regression_seconds"] = perf_counter() - stage
            if candidate.scores.get("path_arclength_collocation_stage", 0.0):
                arc3_states = resample_trajectory_by_arclength(path_states, 3)
                arc7_states = resample_trajectory_by_arclength(path_states, 7)
                arc9_states = resample_trajectory_by_arclength(path_states, 9)
                stage = perf_counter()
                arc3_direct = fit_bilinear_control_map(model, arc3_states, ridge=0.0)
                timings["bilinear_arc3_map_regression_seconds"] = perf_counter() - stage
                stage = perf_counter()
                arc7_conjugate = fit_conjugate_bilinear_control_map(
                    model, arc7_states, ridge=0.0
                )
                timings["bilinear_arc7_conjugate_map_regression_seconds"] = perf_counter() - stage
                stage = perf_counter()
                arc9_conjugate = fit_conjugate_bilinear_control_map(
                    model, arc9_states, ridge=0.0
                )
                timings["bilinear_arc9_conjugate_map_regression_seconds"] = perf_counter() - stage

        dynamics_model = constant.attach(model).with_control_representation(
            "bilinear", dynamics.representation
        )
        map_model = constant.attach(model).with_control_representation(
            "bilinear", direct.representation
        )
        local_map_model = constant.attach(model).with_control_representation(
            "bilinear", local_direct.representation
        )
        conjugate_map_model = constant.attach(model).with_control_representation(
            "bilinear", conjugate_direct.representation
        )
        path_map_model = (
            None
            if path_direct is None
            else constant.attach(model).with_control_representation(
                "bilinear", path_direct.representation
            )
        )
        tube_map_model = (
            None
            if tube_direct is None
            else constant.attach(model).with_control_representation(
                "bilinear", tube_direct.representation
            )
        )
        arc3_map_model = (
            None if arc3_direct is None else constant.attach(model).with_control_representation(
                "bilinear", arc3_direct.representation
            )
        )
        arc7_map_model = (
            None if arc7_conjugate is None else constant.attach(model).with_control_representation(
                "bilinear", arc7_conjugate.representation
            )
        )
        arc9_map_model = (
            None if arc9_conjugate is None else constant.attach(model).with_control_representation(
                "bilinear", arc9_conjugate.representation
            )
        )
        dynamics_archive = reserve_control_representations_archive(
            krom_archive,
            {"constant": constant.representation, "bilinear": dynamics.representation},
            label=f"duffing-bilinear-search-a{attempt_number:02d}-dynamics-maps",
        )
        map_archive = reserve_control_representations_archive(
            krom_archive,
            {"constant": constant.representation, "bilinear": direct.representation},
            label=f"duffing-bilinear-search-a{attempt_number:02d}-direct-map",
        )
        local_map_archive = reserve_control_representations_archive(
            krom_archive,
            {"constant": constant.representation, "bilinear": local_direct.representation},
            label=f"duffing-bilinear-search-a{attempt_number:02d}-target-corridor-map",
        )
        conjugate_map_archive = reserve_control_representations_archive(
            krom_archive,
            {"constant": constant.representation, "bilinear": conjugate_direct.representation},
            label=f"duffing-bilinear-search-a{attempt_number:02d}-conjugate-complete-map",
        )
        path_map_archive = tube_map_archive = None
        arc3_map_archive = arc7_map_archive = arc9_map_archive = None
        if path_direct is not None and tube_direct is not None:
            path_map_archive = reserve_control_representations_archive(
                krom_archive,
                {"constant": constant.representation, "bilinear": path_direct.representation},
                label=f"duffing-bilinear-search-a{attempt_number:02d}-successful-path-map",
            )
            tube_map_archive = reserve_control_representations_archive(
                krom_archive,
                {"constant": constant.representation, "bilinear": tube_direct.representation},
                label=f"duffing-bilinear-search-a{attempt_number:02d}-successful-path-tube-map",
            )
        if arc3_direct is not None and arc7_conjugate is not None and arc9_conjugate is not None:
            arc3_map_archive = reserve_control_representations_archive(
                krom_archive,
                {"constant": constant.representation, "bilinear": arc3_direct.representation},
                label=f"duffing-bilinear-search-a{attempt_number:02d}-arc3-map",
            )
            arc7_map_archive = reserve_control_representations_archive(
                krom_archive,
                {"constant": constant.representation, "bilinear": arc7_conjugate.representation},
                label=f"duffing-bilinear-search-a{attempt_number:02d}-arc7-conjugate-map",
            )
            arc9_map_archive = reserve_control_representations_archive(
                krom_archive,
                {"constant": constant.representation, "bilinear": arc9_conjugate.representation},
                label=f"duffing-bilinear-search-a{attempt_number:02d}-arc9-conjugate-map",
            )
        stage = perf_counter()
        dynamics_fields = control_approximation_error_fields(
            dynamics_model,
            (lower, upper),
            grid_shape=(self.grid_size, self.grid_size),
            representations=("constant", "bilinear"),
        )
        map_fields = control_approximation_error_fields(
            map_model,
            (lower, upper),
            grid_shape=(self.grid_size, self.grid_size),
            representations=("bilinear",),
        )
        local_map_fields = control_approximation_error_fields(
            local_map_model,
            (lower, upper),
            grid_shape=(self.grid_size, self.grid_size),
            representations=("bilinear",),
        )
        conjugate_map_fields = control_approximation_error_fields(
            conjugate_map_model,
            (lower, upper),
            grid_shape=(self.grid_size, self.grid_size),
            representations=("bilinear",),
        )
        fields = {
            "constant": dynamics_fields["constant"],
            "bilinear-dynamics": dynamics_fields["bilinear"],
            "bilinear-map": map_fields["bilinear"],
            "bilinear-local-map": local_map_fields["bilinear"],
            "bilinear-conjugate-map": conjugate_map_fields["bilinear"],
        }
        if path_map_model is not None and tube_map_model is not None:
            path_fields = control_approximation_error_fields(
                path_map_model,
                (lower, upper),
                grid_shape=(self.grid_size, self.grid_size),
                representations=("bilinear",),
            )
            tube_fields = control_approximation_error_fields(
                tube_map_model,
                (lower, upper),
                grid_shape=(self.grid_size, self.grid_size),
                representations=("bilinear",),
            )
            fields["bilinear-path-map"] = path_fields["bilinear"]
            fields["bilinear-path-tube-map"] = tube_fields["bilinear"]
        if arc3_map_model is not None and arc7_map_model is not None and arc9_map_model is not None:
            for field_name, field_model in (
                ("bilinear-arc3-map", arc3_map_model),
                ("bilinear-arc7-conjugate-map", arc7_map_model),
                ("bilinear-arc9-conjugate-map", arc9_map_model),
            ):
                arc_fields = control_approximation_error_fields(
                    field_model,
                    (lower, upper),
                    grid_shape=(self.grid_size, self.grid_size),
                    representations=("bilinear",),
                )
                fields[field_name] = arc_fields["bilinear"]
        timings["error_field_seconds"] = perf_counter() - stage
        error_archive = save_control_approximation_error_fields(
            attempt_root / "control-map-errors.npz",
            fields,
            metadata={
                "system": "duffing",
                "attempt_number": attempt_number,
                "retained_indices": indices,
            },
        )
        arrays_path = attempt_root / "control-representations.npz"
        representation_arrays = {
            "retained_indices": np.asarray(indices, dtype=int),
            "full_spectrum_indices": np.asarray(
                [self.parent._retained_indices[index] for index in indices], dtype=int
            ),
            "continuous_eigenvalues": np.asarray(
                [self.parent.eigenvalues[index].continuous_time for index in indices]
            ),
            "constant_B": constant.representation.matrix,
            "constant_B_discrete": constant.representation.discrete_matrix,
            "bilinear_dynamics_operators": dynamics.representation.operators,
            "bilinear_dynamics_offsets": dynamics.representation.offsets,
            "bilinear_map_operators": direct.representation.operators,
            "bilinear_map_offsets": direct.representation.offsets,
            "bilinear_local_map_operators": local_direct.representation.operators,
            "bilinear_local_map_offsets": local_direct.representation.offsets,
            "bilinear_local_map_training_states": corridor_states,
            "bilinear_conjugate_map_operators": conjugate_direct.representation.operators,
            "bilinear_conjugate_map_conjugate_operators": conjugate_direct.representation.conjugate_operators,
            "bilinear_conjugate_map_offsets": conjugate_direct.representation.offsets,
        }
        if path_direct is not None and tube_direct is not None:
            representation_arrays.update({
                "bilinear_path_map_operators": path_direct.representation.operators,
                "bilinear_path_map_offsets": path_direct.representation.offsets,
                "bilinear_path_map_training_states": path_states,
                "bilinear_path_tube_map_operators": tube_direct.representation.operators,
                "bilinear_path_tube_map_offsets": tube_direct.representation.offsets,
                "bilinear_path_tube_map_training_states": tube_states,
            })
        if arc3_direct is not None and arc7_conjugate is not None and arc9_conjugate is not None:
            representation_arrays.update({
                "bilinear_arc3_map_operators": arc3_direct.representation.operators,
                "bilinear_arc3_map_offsets": arc3_direct.representation.offsets,
                "bilinear_arc3_map_training_states": arc3_states,
                "bilinear_arc7_conjugate_map_operators": arc7_conjugate.representation.operators,
                "bilinear_arc7_conjugate_map_conjugate_operators": arc7_conjugate.representation.conjugate_operators,
                "bilinear_arc7_conjugate_map_offsets": arc7_conjugate.representation.offsets,
                "bilinear_arc7_conjugate_map_training_states": arc7_states,
                "bilinear_arc9_conjugate_map_operators": arc9_conjugate.representation.operators,
                "bilinear_arc9_conjugate_map_conjugate_operators": arc9_conjugate.representation.conjugate_operators,
                "bilinear_arc9_conjugate_map_offsets": arc9_conjugate.representation.offsets,
                "bilinear_arc9_conjugate_map_training_states": arc9_states,
            })
        np.savez_compressed(arrays_path, **representation_arrays)

        case_specs = [
            ("exact", "Exact analytic B_psi", krom_archive, "exact"),
            ("constant", "Controlled-data constant B", dynamics_archive, "constant"),
            (
                "bilinear-dynamics",
                "Controlled-data bilinear dynamics fit",
                dynamics_archive,
                "bilinear",
            ),
            (
                "bilinear-map",
                "Direct pointwise bilinear map fit",
                map_archive,
                "bilinear",
            ),
            (
                "bilinear-local-map",
                "Target-corridor pointwise bilinear map fit",
                local_map_archive,
                "bilinear",
            ),
            (
                "bilinear-conjugate-map",
                "Conjugate-complete real-split bilinear map fit",
                conjugate_map_archive,
                "bilinear",
            ),
        ]
        if path_map_archive is not None and tube_map_archive is not None:
            case_specs.extend((
                (
                    "bilinear-path-map",
                    "Successful-path pointwise bilinear map fit",
                    path_map_archive,
                    "bilinear",
                ),
                (
                    "bilinear-path-tube-map",
                    "Successful-path tube pointwise bilinear map fit",
                    tube_map_archive,
                    "bilinear",
                ),
            ))
        if arc3_map_archive is not None and arc7_map_archive is not None and arc9_map_archive is not None:
            case_specs.extend((
                (
                    "bilinear-arc3-map",
                    "Three-node arc-length collocation bilinear map",
                    arc3_map_archive,
                    "bilinear",
                ),
                (
                    "bilinear-arc7-conjugate-map",
                    "Seven-node arc-length conjugate-complete bilinear map",
                    arc7_map_archive,
                    "bilinear",
                ),
                (
                    "bilinear-arc9-conjugate-map",
                    "Nine-node arc-length conjugate-complete bilinear map",
                    arc9_map_archive,
                    "bilinear",
                ),
            ))
        records: dict[str, dict[str, Any]] = {}
        for key, case_label, case_archive, representation in case_specs:
            record = _solve_case_bounded(
                {
                    "attempt_number": attempt_number,
                    "key": key,
                    "label": case_label,
                    "krom_archive": str(Path(case_archive).resolve()),
                    "representation": representation,
                    "output_dt": self.output_dt,
                },
                attempt_root / f"{key}.worker.json",
                timeout_seconds=self.case_timeout_seconds,
            )
            records[key] = record
            _write_json(attempt_root / "cases.partial.json", records)

        artifacts = {
            "krom_archive": str(Path(krom_archive).resolve()),
            "dynamics_control_krom_archive": str(Path(dynamics_archive).resolve()),
            "map_control_krom_archive": str(Path(map_archive).resolve()),
            "local_map_control_krom_archive": str(Path(local_map_archive).resolve()),
            "conjugate_map_control_krom_archive": str(Path(conjugate_map_archive).resolve()),
            "control_representation_arrays": str(arrays_path.resolve()),
            "control_error_archive": str(Path(error_archive).resolve()),
        }
        if path_map_archive is not None and tube_map_archive is not None:
            artifacts.update({
                "path_map_control_krom_archive": str(Path(path_map_archive).resolve()),
                "path_tube_map_control_krom_archive": str(Path(tube_map_archive).resolve()),
                "path_reference_solution_archive": path_source_archive,
            })
        if arc3_map_archive is not None and arc7_map_archive is not None and arc9_map_archive is not None:
            artifacts.update({
                "arc3_map_control_krom_archive": str(Path(arc3_map_archive).resolve()),
                "arc7_conjugate_map_control_krom_archive": str(Path(arc7_map_archive).resolve()),
                "arc9_conjugate_map_control_krom_archive": str(Path(arc9_map_archive).resolve()),
            })
        regression_diagnostics = {
            "constant": constant.diagnostics,
            "bilinear_dynamics": dynamics.diagnostics,
            "bilinear_map": direct.diagnostics,
            "bilinear_local_map": {
                **local_direct.diagnostics,
                "corridor": {"x": [-0.6, 0.6], "absolute_velocity_maximum": 0.5},
                "corridor_sample_count": int(corridor_states.shape[0]),
            },
            "bilinear_conjugate_map": conjugate_direct.diagnostics,
        }
        if path_direct is not None and tube_direct is not None:
            regression_diagnostics.update({
                "bilinear_path_map": {
                    **path_direct.diagnostics,
                    "reference_solution_archive": path_source_archive,
                    "active_window": list(ACTIVE_WINDOW),
                    "path_sample_count": int(path_states.shape[0]),
                },
                "bilinear_path_tube_map": {
                    **tube_direct.diagnostics,
                    "reference_solution_archive": path_source_archive,
                    "active_window": list(ACTIVE_WINDOW),
                    "path_sample_count": int(path_states.shape[0]),
                    "tube_sample_count": int(tube_states.shape[0]),
                    "tube_radii": [0.025, 0.05],
                },
            })
        if arc3_direct is not None and arc7_conjugate is not None and arc9_conjugate is not None:
            regression_diagnostics.update({
                "bilinear_arc3_map": {
                    **arc3_direct.diagnostics,
                    "sampling": "uniform state-space arc length",
                    "collocation_states": arc3_states.tolist(),
                },
                "bilinear_arc7_conjugate_map": {
                    **arc7_conjugate.diagnostics,
                    "sampling": "uniform state-space arc length",
                    "collocation_states": arc7_states.tolist(),
                },
                "bilinear_arc9_conjugate_map": {
                    **arc9_conjugate.diagnostics,
                    "sampling": "uniform state-space arc length",
                    "collocation_states": arc9_states.tolist(),
                },
            })
        diagnostics = {
            "selected_full_spectrum_indices": [
                int(self.parent._retained_indices[index]) for index in indices
            ],
            "continuous_eigenvalues": [
                self.parent.eigenvalues[index].continuous_time for index in indices
            ],
            "regressions": regression_diagnostics,
            "error_fields": {name: field.summary() for name, field in fields.items()},
            "timings": {**timings, "attempt_total_seconds": perf_counter() - started},
        }
        attempt = KROMControlAttempt(
            attempt_number, candidate, records, artifacts, diagnostics
        )
        _write_json(attempt_root / "attempt.json", attempt.summary(KROMControlCampaignPolicy()))
        return attempt
