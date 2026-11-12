from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from dynsyslib import load_system
from postprocessing.koopman import (
    compare_hankel_observables,
    hankel_dmd,
    observable_spectrum_comparison_report,
)
from postprocessing.cli import build_parser, main
from provenance import serializable
from time_series import TimeSeries


def _rpe_training_data():
    time = np.arange(72, dtype=float) * 0.01
    radius = 1.0 + 0.08 * np.sin(2.0 * time) + 0.025 * np.cos(7.0 * time)
    velocity = 0.16 * np.cos(2.0 * time) - 0.175 * np.sin(7.0 * time)
    return time, np.column_stack((radius, velocity))


def _metadata(system):
    return {
        "uniform_time_steps": True,
        "system": system.name,
        "system_parameters": serializable(system),
        "integration_form": "dimensionless",
        "postprocessing_preprocessing": {"coordinate_system": "dimensionless"},
    }


class HDMDScalarObservableTests(unittest.TestCase):
    def test_named_archived_observable_matches_supplied_scalar_time_series(self):
        system = load_system("rpe")
        time, states = _rpe_training_data()
        settings = {
            "time": time,
            "matrix_size": 8,
            "sampling_window": 48,
            "rank": 5,
            "matrix_layout": "window-snapshots",
            "mode_sort": "terminal-amplitude-energy",
            "metadata": _metadata(system),
        }
        named = hankel_dmd(states, observable="far_field_acoustic_echo", **settings)
        signal = system.observable(
            "far_field_acoustic_echo", coordinate_system="dimensionless"
        ).evaluate(states)[:, 0]
        sampled = hankel_dmd(
            states,
            observable_series=TimeSeries(time, signal[:, None]),
            observable_name="recorded bubble echo",
            **settings,
        )
        named_contract = named.metadata["preprocessing"]["observable_contract"]
        sampled_contract = sampled.metadata["preprocessing"]["observable_contract"]
        self.assertEqual(named_contract["name"], "far_field_acoustic_echo")
        self.assertEqual(named_contract["descriptor"]["kind"], "rpe-observable")
        self.assertEqual(sampled_contract["encoding"], "sampled-scalar-observable/v1")
        self.assertEqual(sampled_contract["name"], "recorded bubble echo")
        self.assertFalse(sampled_contract["descriptor"]["state_callable"])
        np.testing.assert_allclose(
            np.sort_complex(named.eigenvalues),
            np.sort_complex(sampled.eigenvalues),
            rtol=1e-8,
            atol=1e-8,
        )
        self.assertIsNone(sampled.krom.observable_reconstruction)

    def test_rpe_radius_echo_and_energy_can_be_compared_on_one_state_dataset(self):
        system = load_system("rpe")
        time, states = _rpe_training_data()
        comparison = compare_hankel_observables(
            states,
            {
                "radius": "radius",
                "echo": "far_field_acoustic_echo",
                "energy": "inertial_kinetic_energy",
            },
            time=time,
            metadata=_metadata(system),
            matrix_size=8,
            sampling_window=48,
            rank=5,
        )
        self.assertEqual([item.name for item in comparison.entries],
                         ["radius", "echo", "energy"])
        self.assertTrue(all(item.discrete_eigenvalues.size == 5
                            for item in comparison.entries))
        distances = comparison.pairwise_distances("discrete")
        np.testing.assert_allclose(distances, distances.T)
        np.testing.assert_array_equal(np.diag(distances), np.zeros(3))
        report = observable_spectrum_comparison_report(comparison)
        self.assertEqual(report.tile_plan.name, "observable-spectrum-comparison")
        self.assertEqual(len(report.figures), 1)
        self.assertEqual(report.figures[0].kind, "koopman_spectrum")
        self.assertEqual([trace["label"] for trace in report.figures[0].traces],
                         ["unit circle", "radius", "echo", "energy"])

    def test_sampled_series_must_use_the_exact_source_time_grid(self):
        system = load_system("rpe")
        time, states = _rpe_training_data()
        shifted = TimeSeries(time + 1e-6, states[:, :1])
        with self.assertRaisesRegex(ValueError, "exact H-DMD source time grid"):
            hankel_dmd(
                states,
                time=time,
                observable_series=shifted,
                matrix_size=8,
                sampling_window=48,
                matrix_layout="window-snapshots",
                metadata=_metadata(system),
            )

    def test_nonscalar_and_complex_catalog_observables_are_rejected_explicitly(self):
        time, states = _rpe_training_data()
        rpe = load_system("rpe")
        with self.assertRaisesRegex(ValueError, "one scalar Observable"):
            hankel_dmd(
                states,
                time=time,
                observable="state",
                matrix_size=8,
                sampling_window=48,
                matrix_layout="window-snapshots",
                metadata=_metadata(rpe),
            )
        duffing = load_system("duffing")
        with self.assertRaisesRegex(ValueError, "real-valued scalar Observable"):
            hankel_dmd(
                states,
                time=time,
                observable="random_observation",
                matrix_size=8,
                sampling_window=48,
                matrix_layout="window-snapshots",
                metadata=_metadata(duffing),
            )

    def test_marmottant_surface_tension_is_an_hdmd_observable(self):
        system = load_system("marmottant")
        time = np.arange(72, dtype=float) * 0.002
        states = np.column_stack((
            1.0 + 0.04 * np.sin(5.0 * time),
            0.2 * np.cos(5.0 * time),
        ))
        result = hankel_dmd(
            states,
            time=time,
            observable="surface_tension",
            matrix_size=8,
            sampling_window=48,
            rank=4,
            matrix_layout="window-snapshots",
            metadata=_metadata(system),
        )
        contract = result.metadata["preprocessing"]["observable_contract"]
        self.assertEqual(contract["name"], "surface_tension")
        self.assertEqual(contract["descriptor"]["kind"], "marmottant-observable")
        self.assertEqual(contract["descriptor"]["units"], "N/m")

    def test_cli_accepts_named_observables_and_spectrum_comparison(self):
        parser = build_parser()
        planned = parser.parse_args([
            "plan", "solution.npz", "--observable", "far_field_acoustic_echo"
        ])
        self.assertEqual(planned.observable, "far_field_acoustic_echo")
        compared = parser.parse_args([
            "compare-observables", "solution.npz", "--observables",
            "0", "far_field_acoustic_echo", "inertial_kinetic_energy", "--execute"
        ])
        self.assertEqual(compared.observables,
                         [0, "far_field_acoustic_echo", "inertial_kinetic_energy"])
        listed = parser.parse_args(["observables", "solution.npz"])
        self.assertEqual(listed.command, "observables")

    def test_archive_cli_lists_and_compares_named_observables_in_integration_coordinates(self):
        system = load_system("rpe")
        tau, model_states = _rpe_training_data()
        physical_time = np.asarray(system.to_dimensional_time(tau))
        physical_states = np.asarray([
            system.to_dimensional_state(state) for state in model_states
        ])
        metadata = {
            "uniform_time_steps": True,
            "delta_t": float(physical_time[1] - physical_time[0]),
            "system": system.name,
            "system_parameters": serializable(system),
            "integration_form": "dimensionless",
        }
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "rpe.npz"
            np.savez(archive, time=physical_time, values=physical_states,
                     metadata_json=json.dumps(metadata))
            output = io.StringIO()
            with redirect_stdout(output):
                code = main(["observables", str(archive)])
            self.assertEqual(code, 0)
            catalog = json.loads(output.getvalue())
            by_key = {item["key"]: item for item in catalog["observables"]}
            self.assertTrue(by_key["far_field_acoustic_echo"]["hdmd_scalar_eligible"])
            self.assertFalse(by_key["state"]["hdmd_scalar_eligible"])

            output = io.StringIO()
            with redirect_stdout(output):
                code = main([
                    "compare-observables", str(archive), "--observables",
                    "radius", "far_field_acoustic_echo", "--matrix-size", "8",
                    "--sampling-window", "48", "--rank", "5", "--no-archive",
                    "--execute",
                ])
            self.assertEqual(code, 0)
            comparison = json.loads(output.getvalue())
            self.assertEqual(comparison["observable_names"],
                             ["radius", "far_field_acoustic_echo"])
            self.assertEqual(comparison["metadata"]["analysis_coordinates"],
                             "integration")


if __name__ == "__main__":
    unittest.main()
