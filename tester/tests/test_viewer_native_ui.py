from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from viewer.ui import build_catalog, build_payload, write_catalog
from viewer.ui_cli import parser


class NativeViewerUiTests(unittest.TestCase):
    def _archive(self, root: Path, relative: str, **arrays: object) -> Path:
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        np.savez(target, **arrays)
        return target

    def test_catalog_classifies_linked_scientific_objects_without_loading_payloads(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            time = np.linspace(0.0, 1.0, 5)
            training = self._archive(
                root, "archives/solutions/duffing/training.npz", time=time,
                values=np.column_stack((time, time**2)),
                metadata_json=json.dumps({"system": "duffing", "label": "free training"}),
            )
            self._archive(
                root, "archives/solutions/duffing/attempt.npz", time=time,
                values=np.column_stack((time, -time)), control_values=time[:, None],
                metadata_json=json.dumps({"system": "duffing", "label": "controlled attempt",
                                          "controller": {"type": "klqr"}}),
            )
            self._archive(
                root, "archives/derived/duffing/krom.npz",
                eigenvalues=np.asarray((-0.2 + 1.0j, -0.2 - 1.0j)),
                coefficients=np.ones((2, 3), dtype=complex),
                training_states=np.column_stack((time, time**2)),
                metadata_json=json.dumps({"system": "duffing", "label": "two-mode KROM",
                                          "source_archive": str(training.relative_to(root)),
                                          "analysis": "koopman-with-sindy",
                                          "analysis_metadata": {"krom": {
                                              "method": "koopman-with-sindy",
                                              "eigenfunction_count": 2,
                                              "retained_indices": [0, 1],
                                          }}}),
            )
            self._archive(
                root, "archives/derived/duffing/errors.npz", x_values=np.arange(3.0),
                y_values=np.arange(4.0), field_0_relative_error=np.ones((4, 3)),
                metadata_json=json.dumps({"format": "control-approximation-error/v1",
                                          "label": "map error"}),
            )
            self._archive(
                root, "archives/experiments/attempt-2/errors.npz", x_values=np.arange(3.0),
                y_values=np.arange(4.0), field_0_relative_error=np.ones((4, 3)),
                metadata_json=json.dumps({"format": "control-approximation-error/v1",
                                          "label": "map error"}),
            )
            self._archive(root, "archives/derived/duffing/pod.npz",
                          singular_values=np.arange(4.0), modes=np.eye(4))
            campaign = root / "archives/campaigns/duffing.json"
            campaign.parent.mkdir(parents=True)
            campaign.write_text(json.dumps({"name": "Duffing sweep", "system": "duffing",
                                            "campaign_type": "parameter_bifurcation", "jobs": [{}, {}]}),
                                encoding="utf-8")
            report = root / "viewer_outputs/result.pdf"
            report.parent.mkdir(parents=True)
            report.write_bytes(b"%PDF-1.4\n")

            catalog = build_catalog(root)
            by_label = {item["label"]: item for item in catalog["items"]}

            self.assertEqual(by_label["free training"]["kind"], "trajectory")
            self.assertEqual(by_label["free training"]["system"], "Duffing")
            self.assertEqual(by_label["controlled attempt"]["kind"], "control_attempt")
            self.assertIn("control", by_label["controlled attempt"]["views"])
            self.assertEqual(by_label["two-mode KROM"]["kind"], "krom")
            error_labels = [item["label"] for item in catalog["items"]
                            if item["kind"] == "control_error"]
            self.assertEqual(len(error_labels), 2)
            self.assertTrue(all(label.startswith("map error — ") for label in error_labels))
            self.assertEqual(by_label["pod"]["category"], "POD and reduced models")
            self.assertIn("bifurcation", by_label["Duffing sweep"]["views"])
            self.assertEqual(by_label["result"]["views"], ["open"])
            relation = by_label["two-mode KROM"]["relations"][0]
            self.assertEqual(relation["kind"], "derived from")
            self.assertEqual(relation["label"], "free training")

    def test_krom_catalog_collapses_lineages_and_separates_spectra(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            states = np.zeros((5, 2))
            base = self._archive(
                root, "archives/derived/koopman-with-sindy/base.npz",
                coefficients=np.ones((4, 3), dtype=complex), training_states=states,
                metadata_json=json.dumps({
                    "system": "duffing", "analysis": "koopman-with-sindy",
                    "analysis_metadata": {"krom": {
                        "method": "koopman-with-sindy", "eigenfunction_count": 4,
                        "retained_indices": [0, 1, 2, 3],
                    }},
                }),
            )
            extracted = self._archive(
                root, "archives/derived/koopman-krom-extract/minimal.npz",
                eigenvalues=np.asarray((-0.2 + 1j, -0.2 - 1j)),
                coefficients=np.ones((2, 3), dtype=complex), training_states=states,
                metadata_json=json.dumps({
                    "system": "duffing", "analysis": "koopman-krom-extract",
                    "analysis_metadata": {
                        "source_krom_archive": str(base),
                        "reservation_label": "successful-minimal-krom",
                        "krom": {"method": "koopman-with-sindy-krom-extract",
                                 "eigenfunction_count": 2, "retained_indices": [0, 1]},
                    },
                }),
            )
            self._archive(
                root, "archives/derived/koopman-krom-control-representations/maps.npz",
                eigenvalues=np.asarray((-0.2 + 1j, -0.2 - 1j)),
                coefficients=np.ones((2, 3), dtype=complex), training_states=states,
                constant_B=np.ones((2, 1)),
                metadata_json=json.dumps({
                    "system": "duffing", "analysis": "koopman-krom-control-representations",
                    "analysis_metadata": {
                        "source_krom_archive": str(extracted),
                        "reservation_label": "controlled-data-fit",
                        "krom": {"method": "koopman-with-sindy-krom-extract",
                                 "eigenfunction_count": 2, "retained_indices": [0, 1]},
                    },
                }),
            )
            self._archive(
                root, "archives/derived/koopman-hankel-dmd/spectrum.npz",
                eigenvalues=np.asarray((0.9 + 0.1j, 0.9 - 0.1j)), modes=np.eye(2),
                metadata_json=json.dumps({
                    "system": "duffing", "analysis": "koopman-hankel-dmd",
                    "analysis_metadata": {"krom": {"method": "hankel-dmd",
                                                       "eigenfunction_count": 2}},
                }),
            )
            self._archive(
                root, "archives/experiments/regression/control-map.npz",
                eigenvalues_continuous=np.asarray((-0.2 + 1j, -0.2 - 1j)),
                constant_B=np.ones((2, 1)),
            )

            catalog = build_catalog(root)
            by_label = {item["label"]: item for item in catalog["items"]}
            base_item = next(item for item in catalog["items"]
                             if item["source"].endswith("base.npz"))
            extracted_item = by_label["successful-minimal-krom"]
            maps_item = by_label["controlled-data-fit"]
            spectrum = next(item for item in catalog["items"]
                            if item["source"].endswith("spectrum.npz"))
            control_map = next(item for item in catalog["items"]
                               if item["source"].endswith("control-map.npz"))

            self.assertEqual(catalog["counts"]["KROMs"], 1)
            self.assertEqual(catalog["archive_counts"]["KROMs"], 3)
            self.assertTrue(base_item["primary"])
            self.assertEqual(base_item["variants"], 2)
            self.assertEqual(extracted_item["parent"], base_item["id"])
            self.assertEqual(maps_item["parent"], extracted_item["id"])
            self.assertEqual(maps_item["lineage"], base_item["id"])
            self.assertEqual(maps_item["role"], "control-map variant")
            self.assertEqual(spectrum["kind"], "koopman_spectrum")
            self.assertEqual(spectrum["category"], "Koopman spectra")
            self.assertEqual(control_map["category"], "Control maps")

    def test_training_catalog_groups_campaigns_by_physical_dynamics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            time = np.linspace(0.0, 1.0, 5)
            campaign_id = "campaign-pressure-sweep"
            campaign = root / "archives/campaigns/rpe-pressure.json"
            campaign.parent.mkdir(parents=True)
            campaign.write_text(json.dumps({
                "format": "dynsyslab.campaign/v1", "campaign_id": campaign_id,
                "label": "RPE pressure sweep", "jobs": ["one"],
                "metadata": {"campaign_type": "bifurcation", "system_name": "rpe",
                             "pressure_range_kpa": [90.0, 100.0],
                             "frequency_hz": 155000.0},
            }), encoding="utf-8")

            def parameters(viscosity: float) -> dict[str, object]:
                return {
                    "_dynsyslab_class": "RayleighPlessetSystem",
                    "ambient_pressure_pa": 100000.0,
                    "equilibrium_radius_m": 1e-5,
                    "preset": {"mu_pa_s": viscosity, "rho_kg_m3": 998.0,
                               "sigma_n_m": 0.0725, "kappa": 4 / 3,
                               "vapor_pressure_pa": 2330.0},
                }

            def member(name: str, viscosity: float, amplitude: float,
                       formulation: str) -> None:
                self._archive(
                    root, f"archives/solutions/rayleigh-plesset/{name}.npz",
                    time=time, values=np.column_stack((time, -time)),
                    metadata_json=json.dumps({
                        "system": "rayleigh-plesset", "label": name,
                        "system_parameters": parameters(viscosity),
                        "integration_form": formulation,
                        "forcing": {"sources": [{"source": "driver", "parameters": {
                            "amplitude": amplitude, "frequency_hz": 155000.0,
                            "phase_rad": 0.0}}]},
                        "job": {"tags": {"campaign_id": campaign_id}},
                    }),
                )

            member("pressure-90", 0.001, 90000.0, "dimensionless")
            member("pressure-91", 0.001, 91000.0, "dimensionless")
            member("pressure-92-dimensional", 0.001, 92000.0, "dimensional")
            member("different-viscosity", 0.002, 90000.0, "dimensionless")

            catalog = build_catalog(root)
            parameter_sets = [item for item in catalog["items"]
                              if item["kind"] == "physical_parameter_set"]
            campaigns = [item for item in catalog["items"]
                         if item["kind"] == "forcing_campaign"]
            forms = [item for item in catalog["items"]
                     if item["kind"] == "formulation_group"]

            self.assertEqual(len(parameter_sets), 2)
            self.assertEqual(len(campaigns), 2)
            self.assertEqual(len(forms), 2)
            self.assertEqual(catalog["counts"]["Systems"], 2)
            self.assertEqual(catalog["archive_counts"]["Training data"], 4)
            self.assertEqual(catalog["member_counts"]["Training data"], 4)
            three_member_campaign = next(item for item in campaigns
                                         if item["properties"]["dataset count"] == "3")
            self.assertEqual(three_member_campaign["properties"]["planned · pressure range"],
                             "90–100 kPa")
            self.assertEqual(three_member_campaign["properties"]["planned · frequency"],
                             "155 kHz")
            self.assertEqual(three_member_campaign["properties"]["forcing · amplitude"],
                             "90–92 kPa · 3 values")
            self.assertEqual({item["parent"] for item in forms}, {three_member_campaign["id"]})

    def test_compact_free_ensemble_members_are_cataloged_under_one_ensemble(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            grid = root / "archives/time_grids/shared.npz"
            grid.parent.mkdir(parents=True)
            np.savez(grid, time=np.linspace(0.0, 1.0, 5))
            parameters = {"_dynsyslab_class": "RayleighPlessetSystem",
                          "ambient_pressure_pa": 100000.0,
                          "equilibrium_radius_m": 1e-5,
                          "preset": {"mu_pa_s": 0.001, "rho_kg_m3": 998.0}}
            for index in range(3):
                self._archive(
                    root, f"archives/solutions/rayleigh-plesset/free-c000-t{index:03d}.npz",
                    values=np.column_stack((np.arange(5.0), np.full(5, index))),
                    raw_time=np.linspace(0.0, 1.0, 3), raw_values=np.zeros((3, 2)),
                    metadata_json=json.dumps({
                        "system": "rayleigh-plesset",
                        "label": f"rpe-free-ensemble-c000-t{index:03d}",
                        "system_parameters": parameters,
                        "integration_form": "dimensionless",
                        "time_grid_ref": "time_grids/shared.npz",
                        "initial_state_override": [1e-5 * (index + 1), 0.0],
                    }),
                )

            catalog = build_catalog(root)
            trajectories = [item for item in catalog["items"] if item["kind"] == "trajectory"]
            ensembles = [item for item in catalog["items"]
                         if item["kind"] == "trajectory_ensemble"]
            parameter = next(item for item in catalog["items"]
                             if item["kind"] == "physical_parameter_set")

            self.assertEqual(len(trajectories), 3)
            self.assertEqual(len(ensembles), 1)
            self.assertEqual(ensembles[0]["label"], "rpe-free-ensemble")
            data_stage = next(item for item in catalog["items"]
                              if item["parent"] == parameter["id"] and
                              item["label"] == "Data")
            self.assertEqual(ensembles[0]["parent"], data_stage["id"])
            self.assertEqual({item["parent"] for item in trajectories}, {ensembles[0]["id"]})
            self.assertEqual(ensembles[0]["properties"]["dataset count"], "3")
            self.assertEqual(ensembles[0]["properties"]["initial state x0"], "1e-05–3e-05")

    def test_system_workspace_links_data_krom_and_control_by_physical_realization(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            time = np.linspace(0.0, 1.0, 5)
            parameters = {
                "ambient_pressure_pa": 100000.0, "equilibrium_radius_m": 1e-5,
                "preset": {"description": "Air in water", "kappa": 4 / 3,
                           "mu_pa_s": 0.001, "rho_kg_m3": 998.0,
                           "sigma_n_m": 0.0725, "vapor_pressure_pa": 2330.0},
            }
            training = self._archive(
                root, "archives/solutions/rayleigh-plesset/free.npz", time=time,
                values=np.column_stack((time, -time)),
                metadata_json=json.dumps({
                    "system": "rayleigh-plesset", "label": "free data",
                    "system_parameters": parameters, "integration_form": "dimensionless",
                    "dimensionless_groups": {"viscosity": 0.12},
                }))
            krom = self._archive(
                root, "archives/derived/rpe/krom.npz",
                eigenvalues=np.asarray((-0.2 + 1j, -0.2 - 1j)),
                coefficients=np.ones((2, 3), dtype=complex),
                training_states=np.column_stack((time, -time)),
                metadata_json=json.dumps({
                    "system": "rayleigh-plesset", "label": "control KROM",
                    "source_archive": str(training.relative_to(root)),
                    "analysis": "koopman-with-sindy",
                    "analysis_metadata": {"krom": {"method": "koopman-with-sindy",
                                                      "eigenfunction_count": 2,
                                                      "retained_indices": [0, 1]}},
                }))
            self._archive(
                root, "archives/solutions/rayleigh-plesset/control.npz", time=time,
                values=np.column_stack((time, -time)), control_values=time[:, None],
                metadata_json=json.dumps({
                    "system": "rayleigh-plesset", "label": "KLQR test",
                    "system_parameters": parameters, "integration_form": "dimensionless",
                    "controller": {"controller": "koopman-lqr", "target": "multisine",
                                   "rom": {"source_archive": str(krom)}},
                }))

            catalog = build_catalog(root)
            parameter = next(item for item in catalog["items"]
                             if item["kind"] == "physical_parameter_set")
            self.assertEqual(parameter["category"], "Systems")
            self.assertTrue(parameter["label"].startswith("Air in water"))
            self.assertEqual(parameter["properties"]["dimensionless group · viscosity"],
                             "0.12")
            stages = {item["label"] for item in catalog["items"]
                      if item["parent"] == parameter["id"]}
            self.assertEqual(stages, {"Data", "KROMs and spectra", "Control attempts"})
            free = next(item for item in catalog["items"] if item["label"] == "free data")
            attempt = next(item for item in catalog["items"] if item["label"] == "KLQR test")
            self.assertIn("free", free["tags"])
            self.assertIn("KLQR", attempt["tags"])
            self.assertEqual(attempt["properties"]["control target"], "multisine")
            self.assertTrue(any(relation["kind"] == "controlled with"
                                for relation in attempt["relations"]))

    def test_generic_payload_is_bounded_and_separates_state_control_and_phase(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "controlled.npz"
            samples = 3000
            np.savez(source, time=np.arange(samples, dtype=float),
                     states=np.column_stack((np.arange(samples), -np.arange(samples))),
                     controls=np.ones((samples, 1)))
            payload = build_payload(source, "data_overview")

        self.assertEqual([figure["kind"] for figure in payload["figures"]],
                         ["time", "phase", "time"])
        self.assertLessEqual(len(payload["figures"][0]["traces"][0]["x"]), 2200)
        self.assertTrue(payload["figures"][1]["square"])
        self.assertEqual(payload["figures"][2]["title"], "Control excitation")

    def test_observable_spectrum_comparison_is_cataloged_and_rendered_natively(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self._archive(
                root,
                "archives/derived/koopman-observable-spectrum-comparison/echo-energy.npz",
                observable_names=np.asarray(("echo", "energy")),
                observable_contracts_json=np.asarray((
                    json.dumps({"name": "far_field_acoustic_echo"}),
                    json.dumps({"name": "inertial_kinetic_energy"}),
                )),
                mode_counts=np.asarray((2, 2)),
                discrete_eigenvalues=np.asarray(((0.9 + 0.1j, 0.9 - 0.1j),
                                                  (0.8 + 0.2j, 0.8 - 0.2j))),
                continuous_eigenvalues=np.asarray(((-0.1 + 1.0j, -0.1 - 1.0j),
                                                    (-0.2 + 2.0j, -0.2 - 2.0j))),
                metadata_json=json.dumps({
                    "format": "dynsyslab.derived/v1",
                    "analysis": "koopman-observable-spectrum-comparison",
                    "analysis_metadata": {"comparison": {"metadata": {
                        "same_state_training_data": True
                    }}},
                }),
            )
            catalog = build_catalog(root)
            item = next(entry for entry in catalog["items"]
                        if entry["kind"] == "observable_spectrum_comparison")
            self.assertEqual(item["views"], ["observable_spectrum_comparison"])
            payload = build_payload(source, "observable_spectrum_comparison")

        self.assertEqual(payload["figures"][0]["kind"], "koopman_spectrum")
        self.assertTrue(payload["figures"][0]["square"])
        self.assertEqual([trace["label"] for trace in payload["figures"][0]["traces"]],
                         ["unit circle", "echo", "energy"])

    def test_catalog_cache_is_not_rewritten_when_sources_are_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._archive(root, "archives/solutions/free.npz", time=np.arange(2.0),
                          values=np.zeros((2, 2)))
            output = root / "cache/catalog.json"
            write_catalog(root, output)
            first = output.stat().st_mtime_ns
            write_catalog(root, output)
            self.assertEqual(output.stat().st_mtime_ns, first)

    def test_cli_contract_exposes_catalog_payload_build_and_launch(self):
        command = parser()
        self.assertEqual(command.parse_args(["catalog", "--output", "x.json"]).command,
                         "catalog")
        self.assertEqual(command.parse_args(["payload", "--source", "x.npz", "--view",
                                             "phase", "--output", "y.json"]).view, "phase")
        self.assertTrue(command.parse_args(["launch", "--rebuild"]).rebuild)
        self.assertEqual(command.parse_args(["build"]).configuration, "Release")

    def test_native_viewer_modules_remain_small_and_cmake_declares_the_app(self):
        root = Path(__file__).resolve().parents[2]
        sources = list((root / "viewer/native/src").glob("*.cpp")) + \
                  list((root / "viewer/native/include").rglob("*.hpp")) + \
                  list((root / "viewer/ui").glob("*.py")) + [root / "viewer/ui_cli.py"]
        oversized = {path.relative_to(root).as_posix(): len(path.read_text(encoding="utf-8").splitlines())
                     for path in sources
                     if len(path.read_text(encoding="utf-8").splitlines()) >= 500}
        self.assertEqual(oversized, {})
        cmake = (root / "viewer/native/CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("add_executable(dynsyslab_viewer WIN32", cmake)
        self.assertIn("add_subdirectory(../../native_ui", cmake)
        self.assertIn("dynsyslab_ui_support", cmake)
        self.assertIn("src/PlotExport.cpp", cmake)
        self.assertIn("src/KromControlCampaignPayload.cpp", cmake)
        self.assertIn("src/HierarchyPayload.cpp", cmake)
        self.assertIn("krom_control_campaign_payload_tests.cpp", cmake)
        self.assertIn("comdlg32", cmake)
        self.assertIn("DynSysLabViewerIcon.rc", cmake)
        self.assertIn("/utf-8", cmake)
        workbench = (root / "viewer/native/src/Workbench.cpp").read_text(encoding="utf-8")
        self.assertIn("Help / About", workbench)
        self.assertIn("Copy plot", workbench)
        self.assertIn("Save PNG", workbench)
        self.assertIn("IDI_DYNSYSLAB_APP", workbench)
        self.assertIn("kViewSelectId", workbench)
        self.assertNotIn("viewButtons_.size()", workbench)
        self.assertIn("Search archive tree", workbench)
        self.assertNotIn("Filter objects", workbench)
        self.assertNotIn("linked archive objects", workbench)
        self.assertNotIn("lazy native plots", workbench)
        self.assertNotIn('L"Loaded "', workbench)
        self.assertNotIn("status_", workbench)
        self.assertEqual(workbench.count(
            'theme::scalePixels(window_, 10), L"DynSysLab Viewer", 16'), 1)
        self.assertIn("theme::metrics(window_)", workbench)
        plot_canvas = (root / "viewer/native/src/PlotFigureDrawing.cpp").read_text(
            encoding="utf-8")
        self.assertIn("automaticDetailInset", plot_canvas)
        self.assertIn("automatic detail", plot_canvas)
        campaign = (root / "viewer/native/src/KromControlCampaignPayload.cpp").read_text(
            encoding="utf-8")
        for title in ("Target error at control release", "Final error after release",
                      "Peak control demand", "Mean B_psi approximation error",
                      "Candidate pre-screening", "Extracted KROM size"):
            self.assertIn(title, campaign)


if __name__ == "__main__":
    unittest.main()
