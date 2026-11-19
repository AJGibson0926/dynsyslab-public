from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from compute_workbench.contracts import (JOB_SCHEMA, postprocess_template,
                                         solve_template, validate_job)
from compute_workbench.execution import plan_job, plan_reference_job
from compute_workbench.progress import PROGRESS_SCHEMA, ProgressReporter
from compute_workbench.initial_conditions import (generate_initial_conditions,
                                                  make_clusters)
from compute_workbench.registry import build_registry
from dynsyslib import load_system


ROOT = Path(__file__).resolve().parents[2]


class ComputeWorkbenchTests(unittest.TestCase):
    def test_templates_are_versioned_and_independently_validated(self):
        solve = validate_job(solve_template())
        self.assertEqual(solve["schema"], JOB_SCHEMA)
        self.assertEqual(solve["kind"], "solve")
        post = postprocess_template()
        post["sources"] = ["one.npz"]
        self.assertEqual(validate_job(post)["kind"], "postprocess")

    def test_registry_exposes_systems_presets_parameters_targets_methods_and_kroms(self):
        registry = build_registry(ROOT)
        self.assertEqual(registry["schema"], "dynsyslab.compute-registry/v1")
        names = {item["name"] for item in registry["systems"]}
        self.assertTrue({"rpe", "duffing", "marmottant", "reddy_szeri"}.issubset(names))
        rpe = next(item for item in registry["systems"] if item["name"] == "rpe")
        parameters = {item["name"] for item in rpe["parameters"]}
        self.assertTrue({"equilibrium_radius_m", "ambient_pressure_pa", "rho_kg_m3",
                         "mu_pa_s", "sigma_n_m", "kappa"}.issubset(parameters))
        self.assertIn("hdmd", registry["postprocessing_methods"])
        self.assertIn("havok", registry["postprocessing_methods"])
        self.assertTrue(registry["targets"])
        self.assertLess(len(json.dumps(registry)), 300_000)

    def test_initial_condition_generation_is_deterministic_and_clustered(self):
        system = load_system("duffing")
        settings = {"mode": "uniform", "count": 7, "clusters": 3, "seed": 19,
                    "ranges": [[-1.0, 1.0], [-2.0, 2.0]]}
        first = generate_initial_conditions(system, settings, controlled=False)
        second = generate_initial_conditions(system, settings, controlled=False)
        np.testing.assert_allclose(first, second)
        clusters = make_clusters(first, 3)
        self.assertEqual([len(item.trajectories) for item in clusters], [3, 2, 2])

    def test_rpe_controller_campaign_plans_batch_system_cluster_trajectory_counts(self):
        job = solve_template()
        job["system"] = {"name": "rpe", "preset": "air_water_20c", "parameters": {}}
        job["time"] = {"start": 0.0, "stop": 1.0, "output_step": 0.01,
                       "coordinate_system": "dimensionless"}
        job["initial_conditions"].update({"count": 4, "clusters": 2})
        job["controller"] = {"kind": "waveform", "settings": {
            "frequencies_hz": [100_000.0], "amplitudes": [50_000.0],
            "phases_rad": [0.0]}, "target": {"kind": "constant", "state": [1.0, 0.0]}}
        job["campaign"] = {"enabled": True, "mode": "cartesian", "axes": [
            {"path": "controller.settings.frequencies_hz.0", "values": [100_000.0, 200_000.0]},
            {"path": "controller.settings.amplitudes.0", "values": [25_000.0, 50_000.0, 75_000.0]},
        ]}
        plan = plan_job(job, ROOT)
        reference = plan_reference_job(job, ROOT)
        self.assertTrue(plan["campaign"])
        self.assertEqual(plan["system_count"], 6)
        self.assertEqual(plan["trajectory_count"], 24)
        self.assertEqual(plan["cluster_count"], 12)
        self.assertEqual(plan["clusters_per_system"], 2)
        for field in ("system_count", "trajectory_count", "cluster_count",
                      "clusters_per_system", "state_dimension",
                      "output_samples_per_trajectory"):
            self.assertEqual(plan[field], reference[field])

    def test_cli_executes_one_small_solve_and_independent_postprocess(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory); archive = work / "archives"
            job = solve_template(); job["label"] = "compute-ui-test"
            job["time"].update({"stop": 1.0, "output_step": 0.01})
            solve_path = work / "solve.json"; solve_path.write_text(json.dumps(job), encoding="utf-8")
            solve_output = work / "solve-result.json"; progress_output = work / "progress.json"
            environment = dict(os.environ); environment["DYNSYSLAB_ARCHIVE_DIR"] = str(archive)
            command = [sys.executable, "-m", "compute_workbench.cli", "--project-root", str(ROOT),
                       "run", "--job", str(solve_path), "--execute", "--output", str(solve_output),
                       "--progress-output", str(progress_output)]
            subprocess.run(command, cwd=ROOT, env=environment, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            solved = json.loads(solve_output.read_text(encoding="utf-8"))
            self.assertEqual(solved["status"], "complete")
            progress = json.loads(progress_output.read_text(encoding="utf-8"))
            self.assertEqual(progress["schema"], PROGRESS_SCHEMA)
            self.assertEqual(progress["state"], "complete")
            self.assertEqual(progress["percent"], 100.0)
            source = solved["solution_archives"][0]

            post = postprocess_template(); post["sources"] = [source]
            post["methods"] = [{"name": "dmd", "enabled": True, "settings": {}}]
            post_path = work / "post.json"; post_path.write_text(json.dumps(post), encoding="utf-8")
            post_output = work / "post-result.json"
            command = [sys.executable, "-m", "compute_workbench.cli", "--project-root", str(ROOT),
                       "run", "--job", str(post_path), "--execute", "--output", str(post_output)]
            subprocess.run(command, cwd=ROOT, env=environment, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            processed = json.loads(post_output.read_text(encoding="utf-8"))
            self.assertEqual(processed["status"], "complete")
            self.assertEqual(len(processed["derived_archives"]), 1)

    def test_cli_control_sweep_archives_one_batch_with_per_system_controllers(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory); archive = work / "archives"
            job = solve_template(); job["label"] = "duffing-waveform-sweep"
            job["time"].update({"stop": 0.02, "output_step": 0.01})
            job["controller"] = {"kind": "waveform", "settings": {
                "frequencies_hz": [1.0], "amplitudes": [0.1], "phases_rad": [0.0]},
                "target": {"kind": "constant", "state": [0.0, 0.0]}}
            job["campaign"] = {"enabled": True, "mode": "cartesian", "axes": [{
                "path": "controller.settings.amplitudes.0", "values": [0.1, 0.2]}]}
            job_path = work / "campaign.json"; job_path.write_text(json.dumps(job), encoding="utf-8")
            output = work / "result.json"; environment = dict(os.environ)
            environment["DYNSYSLAB_ARCHIVE_DIR"] = str(archive)
            completed = subprocess.run([sys.executable, "-m", "compute_workbench.cli", "--project-root", str(ROOT),
                                        "run", "--job", str(job_path), "--execute", "--output", str(output)],
                                       cwd=ROOT, env=environment, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr[-3000:])
            result = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(result["system_count"], 2)
            self.assertEqual(result["trajectory_count"], 2)
            self.assertTrue(Path(result["campaign_manifest"]).is_file())
            self.assertEqual(len(result["hierarchy_archives"]), 2)
            self.assertTrue(all(Path(path).is_file()
                                for path in result["hierarchy_archives"]))
            amplitudes = []
            for path in result["solution_archives"]:
                with np.load(path, allow_pickle=False) as data:
                    metadata = json.loads(str(data["metadata_json"]))
                self.assertEqual(metadata["campaign_id"], result["campaign_id"])
                amplitudes.append(metadata["controller_configuration"]["settings"]["amplitudes"][0])
            self.assertEqual(amplitudes, [0.1, 0.2])

    def test_native_app_reuses_shared_ui_library_and_all_new_modules_are_small(self):
        source = (ROOT / "compute_workbench/native/CMakeLists.txt").read_text(encoding="utf-8")
        support = (ROOT / "native_ui/CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("add_subdirectory(../../native_ui", source)
        self.assertIn("dynsyslab_ui_support", source)
        for component in ("Theme.cpp", "HierarchyCatalog.cpp", "RuntimePaths.cpp",
                          "Process.cpp", "NavigationPane.cpp", "InspectorPane.cpp",
                          "Model.cpp", "Help.cpp"):
            self.assertIn(component, support)
        queue = (ROOT / "compute_workbench/native/src/QueueDialog.cpp").read_text(
            encoding="utf-8")
        for label in ("Durable native job queue", "Queue active job", "Run next",
                      "Cancel selected", "Recover stale after"):
            self.assertIn(label, queue)
        self.assertLessEqual(len(queue.splitlines()), 500)
        workflow = (ROOT / "workflow/native/src/Queue.cpp").read_text(encoding="utf-8")
        self.assertIn("milliseconds(100)", workflow)
        self.assertIn("ExecutionCancelled", workflow)
        archive_tools = (ROOT / "compute_workbench/native/src/ArchiveToolsDialog.cpp").read_text(
            encoding="utf-8")
        for label in ("Canonical archive maintenance", "Run read-only audit",
                      "Preview deduplication", "Apply recoverable move"):
            self.assertIn(label, archive_tools)
        self.assertLessEqual(len(archive_tools.splitlines()), 500)
        self.assertIn("Help / About", (ROOT / "compute_workbench/native/src/Workbench.cpp").read_text(
            encoding="utf-8"))
        self.assertIn("Open Viewer", (ROOT / "compute_workbench/native/src/Workbench.cpp").read_text(
            encoding="utf-8"))
        backend = (ROOT / "compute_workbench/native/src/Backend.cpp").read_text(
            encoding="utf-8")
        self.assertIn("dynsyslab::workflow::executeJob", backend)
        self.assertIn("DynSysLabViewer.exe", backend)
        self.assertIn("dynsyslab.compute-launch/v1", (ROOT / "compute_workbench/cli.py").read_text(
            encoding="utf-8"))
        workbench = (ROOT / "compute_workbench/native/src/Workbench.cpp").read_text(
            encoding="utf-8")
        workbench_arguments = (ROOT / "compute_workbench/native/src/WorkbenchArguments.cpp").read_text(
            encoding="utf-8")
        self.assertIn("discoverRuntimePaths", workbench_arguments)
        self.assertNotIn("--python", workbench)
        self.assertIn("DynSysLabAppIcon.rc", source)
        self.assertIn("src/KromTrainingDialog.cpp", source)
        self.assertIn("src/LegacyImportDialog.cpp", source)
        self.assertIn("src/KromControlSearchDialog.cpp", source)
        self.assertIn("src/SolverExecutionPane.cpp", source)
        training_dialog = (ROOT / "compute_workbench/native/src/KromTrainingDialog.cpp").read_text(
            encoding="utf-8")
        for label in ("KROM training-domain expansion", "Preview exact plan",
                      "Apply to Solver"):
            self.assertIn(label, training_dialog)
        self.assertLessEqual(len(training_dialog.splitlines()), 500)
        regression_dialog = (ROOT / "compute_workbench/native/src/ControlRegressionDialog.cpp").read_text(
            encoding="utf-8")
        for label in ("Additional controlled trajectory archives", "Add files…",
                      "Korda–Mezić multistep output regression",
                      "Observable reconstruction C", "Exact-map sample stride"):
            self.assertIn(label, regression_dialog)
        self.assertLessEqual(len(regression_dialog.splitlines()), 500)
        import_dialog = (ROOT / "compute_workbench/native/src/LegacyImportDialog.cpp").read_text(
            encoding="utf-8")
        for label in ("Import legacy tagged trajectory CSV", "Preview / validate",
                      "Import and archive"):
            self.assertIn(label, import_dialog)
        self.assertLessEqual(len(import_dialog.splitlines()), 500)
        self.assertIn("Import data…", workbench)
        search_dialog = (ROOT / "compute_workbench/native/src/KromControlSearchDialog.cpp").read_text(
            encoding="utf-8")
        for label in ("Bounded KROM control-candidate campaign",
                      "Controlled regression ensemble",
                      "Maximum attempts (hard limit 25)", "Run campaign",
                      "Plan only", "Save job…", "Queue"):
            self.assertIn(label, search_dialog)
        self.assertLessEqual(len(search_dialog.splitlines()), 500)
        post_panel = (ROOT / "compute_workbench/native/src/PostprocessPanel.cpp").read_text(
            encoding="utf-8")
        self.assertIn("KROM control search…", post_panel)
        self.assertIn("Cancel active", workbench)
        self.assertIn("cancelRequested_.load()", workbench)
        solver_execution = (ROOT / "compute_workbench/native/src/SolverExecutionPane.cpp").read_text(
            encoding="utf-8")
        for label in ("Solver and execution", "Internal step", "Maximum step",
                      "Maximum steps", "Backend"):
            self.assertIn(label, solver_execution)
        self.assertLessEqual(len(solver_execution.splitlines()), 500)
        solver_panel = (ROOT / "compute_workbench/native/src/SolverPanel.cpp").read_text(
            encoding="utf-8")
        self.assertIn('caption(L"x*")', solver_panel)
        controller_panel = (ROOT /
            "compute_workbench/native/src/SolverPanelController.cpp").read_text(
                encoding="utf-8")
        self.assertIn("scalarOrArray", controller_panel)
        cli = (ROOT / "compute_workbench/native/src/Cli.cpp").read_text(encoding="utf-8")
        self.assertIn("import legacy-csv", cli)
        self.assertIn("krom_control_search", cli)
        self.assertTrue((ROOT / "native_ui/assets/dynsyslab-app.ico").is_file())
        self.assertTrue((ROOT / "native_ui/assets/dynsyslab-viewer.ico").is_file())
        resource = (ROOT / "native_ui/windows/DynSysLabAppIcon.rc").read_text(encoding="utf-8")
        self.assertIn("IDI_DYNSYSLAB_APP", resource)
        files = list((ROOT / "compute_workbench").rglob("*.py")) + \
                list((ROOT / "compute_workbench/native").rglob("*.cpp")) + \
                list((ROOT / "compute_workbench/native").rglob("*.hpp")) + \
                [ROOT / "native_ui/CMakeLists.txt"]
        oversized = {str(path.relative_to(ROOT)): len(path.read_text(encoding="utf-8").splitlines())
                     for path in files if len(path.read_text(encoding="utf-8").splitlines()) > 500}
        self.assertEqual(oversized, {})

    def test_progress_reporter_combines_trajectory_and_in_trajectory_progress(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "progress.json"
            reporter = ProgressReporter(path, kind="solve", planned_trajectories=2)
            reporter.event({"event": "heartbeat", "progress_fraction": 0.5,
                            "simulation_time": 0.25})
            halfway_first = json.loads(path.read_text(encoding="utf-8"))
            self.assertTrue(halfway_first["determinate"])
            self.assertEqual(halfway_first["percent"], 25.0)
            reporter.event({"event": "hierarchy_progress",
                            "completed_trajectories": 1,
                            "total_trajectories": 2})
            self.assertEqual(json.loads(path.read_text(encoding="utf-8"))["percent"], 50.0)
            reporter.complete()
            self.assertEqual(json.loads(path.read_text(encoding="utf-8"))["state"], "complete")

    def test_packaged_desktop_apps_have_headless_initialization_gate(self):
        smoke = (ROOT / "cmake/PackageSmoke.cmake").read_text(encoding="utf-8")
        compute = (ROOT / "compute_workbench/native/src/Workbench.cpp").read_text(
            encoding="utf-8")
        compute_arguments = (ROOT / "compute_workbench/native/src/WorkbenchArguments.cpp").read_text(
            encoding="utf-8")
        viewer = (ROOT / "viewer/native/src/WorkbenchArguments.cpp").read_text(encoding="utf-8")
        self.assertIn("foreach(_ui DynSysLab.exe DynSysLabViewer.exe)", smoke)
        self.assertIn("--smoke-test", compute_arguments)
        self.assertIn("--smoke-test", viewer)


if __name__ == "__main__":
    unittest.main()
