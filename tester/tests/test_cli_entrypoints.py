from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT=Path(__file__).resolve().parents[2]


class PackageCLIEntrypointTests(unittest.TestCase):
    def _run(self,*arguments: str,env: dict[str,str]|None=None):
        return subprocess.run(
            (sys.executable,*arguments),cwd=PROJECT_ROOT,
            capture_output=True,text=True,timeout=30,check=False,env=env,
        )

    def test_performance_module_routes_to_benchmark_help(self):
        completed=self._run("-m","performance","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("DynSysLab numerical benchmarks",completed.stdout)

    def test_postprocessing_module_routes_to_workflow_help(self):
        completed=self._run("-m","postprocessing","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("DynSysLab post-processing workflows",completed.stdout)
        self.assertIn("observables",completed.stdout)
        self.assertIn("compare-observables",completed.stdout)

    def test_postprocessing_cli_exposes_sindy_observable_coordinates(self):
        completed=self._run("-m","postprocessing","describe","--sindy-coordinates","c","k5",
                            "--coordinate-system","dimensionless")
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        sindy=payload["methods"][1]["settings"]
        self.assertEqual(sindy["coordinate_observable"],["c","k5"])
        self.assertEqual(sindy["coordinate_system"],"dimensionless")

    def test_archiver_module_routes_to_archive_help(self):
        completed=self._run("-m","archiver","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("audit",completed.stdout)
        self.assertIn("deduplicate",completed.stdout)

    def test_rom_module_exposes_observable_selection_help(self):
        completed=self._run("-m","rom","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("select-observable",completed.stdout)
        self.assertIn("fit-observation",completed.stdout)

    def test_estimation_module_exposes_observability_and_klqe_benchmark_help(self):
        completed=self._run("-m","estimation","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("observability",completed.stdout)
        self.assertIn("design",completed.stdout)
        self.assertIn("benchmark",completed.stdout)

    def test_legacy_gui_help_does_not_launch_a_window(self):
        completed=self._run("-m","gui","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("compatibility Tk interface",completed.stdout)
        self.assertIn("launch",completed.stdout)

    def test_primary_front_door_help_does_not_launch_an_application(self):
        completed=self._run("-m","dynsyslab","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("Solver/Post-processor",completed.stdout)
        self.assertIn("doctor",completed.stdout)
        self.assertIn("viewer",completed.stdout)

    def test_primary_front_door_doctor_is_machine_readable(self):
        completed=self._run("-m","dynsyslab","doctor")
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.workspace-diagnostics/v1")
        self.assertEqual(Path(payload["project_root"]),PROJECT_ROOT)
        self.assertTrue(payload["paths"]["start_document"]["exists"])

    def test_primary_system_catalog_delegates_to_compiled_registry(self):
        completed=self._run("-m","dynsyslab","systems")
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.system-catalog/v1")
        families={item["family"] for item in payload["systems"]}
        self.assertTrue({"ode","dde","lattice","pde"}.issubset(families))

    def test_primary_archive_deduplication_is_native_and_dry_by_default(self):
        from operations.front_door import _native_executable
        with tempfile.TemporaryDirectory() as directory:
            environment={**os.environ,"DYNSYSLAB_INSTALL_DIR":str(
                _native_executable("dynsyslab-cli.exe").parent.parent)}
            completed=self._run(
                "-m","dynsyslab","archive","deduplicate","--archive-root",
                directory,"--all-ages",env=environment,
            )
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.archive-deduplication/v1")
        self.assertTrue(payload["dry_run"])
        self.assertTrue(payload["recoverable"])

    def test_primary_front_door_routes_normal_workflows_without_execution(self):
        expectations={
            "solve":"Solver and Post-processor",
            "postprocess":"post-processing workflows",
            "archive":"usage:",
            "jobs":"campaigns",
            "estimate":"observability",
            "krom-training":"system-agnostic",
        }
        for command, marker in expectations.items():
            with self.subTest(command=command):
                completed=self._run("-m","dynsyslab",command,"--help")
                self.assertEqual(completed.returncode,0,completed.stderr)
                self.assertIn(marker.lower(),completed.stdout.lower())

    def test_primary_postprocess_front_door_plans_through_native_abi(self):
        from dynsyslab_client.workflow import find_library
        library=find_library(PROJECT_ROOT)
        job={
            "schema":"dynsyslab.compute-job/v1","kind":"postprocess",
            "label":"native-front-door-test","sources":["not-loaded-by-plan.npz"],
            "methods":[{"name":"hdmd","enabled":True,"settings":{}}],
            "archive":{"enabled":True},
        }
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"job.json"
            path.write_text(json.dumps(job),encoding="utf-8")
            completed=self._run(
                "-m","dynsyslab","postprocess","--native-library",str(library),
                "plan",str(path),
            )
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.workflow-plan/v1")
        self.assertEqual(payload["kind"],"postprocess")
        self.assertEqual(payload["methods"],["hdmd"])

    def test_primary_estimation_front_door_is_a_native_job_adapter(self):
        completed=self._run("-m","dynsyslab","estimate","template")
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["kind"],"postprocess")
        self.assertEqual(payload["methods"][0]["name"],"observability")
        self.assertEqual(payload["methods"][0]["settings"]["measurement"],"echo")

    def test_primary_jobs_front_door_delegates_to_native_queue(self):
        from operations.front_door import _native_executable
        environment={**os.environ,"DYNSYSLAB_INSTALL_DIR":str(
            _native_executable("dynsyslab-cli.exe").parent.parent)}
        with tempfile.TemporaryDirectory() as directory:
            completed=self._run(
                "-m","dynsyslab","jobs","list","--archive-root",directory,
                env=environment,
            )
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.job-queue/v1")
        self.assertEqual(payload["operation"],"list")

    def test_primary_archive_audit_delegates_to_compiled_cli(self):
        completed=self._run("-m","dynsyslab","archive","audit","--days","0.001")
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.archive-audit/v1")
        self.assertTrue(payload["read_only"])

    def test_gui_package_is_lazy_for_help_and_import(self):
        code=("import json,sys,gui; print(json.dumps({"
              "'desktop':'gui.desktop' in sys.modules,"
              "'has_launch':'launch' in dir(gui)}))")
        completed=self._run("-c",code)
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertEqual(json.loads(completed.stdout),
                         {"desktop":False,"has_launch":True})

    def test_krom_training_planner_help_is_system_agnostic(self):
        completed=self._run("-m","scripts.design_krom_training_campaign","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("system-agnostic",completed.stdout)
        self.assertIn("--points-per-axis",completed.stdout)

    def test_duffing_bilinear_search_help_exposes_bounded_campaign_controls(self):
        completed=self._run("-m","scripts.search_duffing_bilinear_krom","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("bounded search",completed.stdout.lower())
        self.assertIn("--max-attempts",completed.stdout)
        self.assertIn("--resume",completed.stdout)

    def test_duffing_bilinear_report_help_is_lightweight(self):
        completed=self._run("-m","scripts.report_duffing_bilinear_campaign","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("Viewer reports",completed.stdout)

    def test_duffing_bilinear_pdf_help_is_available_without_reporting_extra(self):
        completed=self._run("-m","scripts.export_duffing_bilinear_campaign_pdf","--help")
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertIn("formal static PDF",completed.stdout)

    def test_solver_module_emits_backend_diagnostics_json(self):
        completed=self._run("-m","solver")
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload["schema"],"dynsyslab.backend-diagnostics/v1")
        native=payload["solver"]["native"]
        if native["available"]:
            self.assertIn("rpe_dimensionless",native["fixed_step_models"])
        else:
            # Pybind is optional; native CLI/ABI integration is tested separately.
            self.assertEqual(native["fixed_step_models"],[])
            self.assertIn("not installed",native["reason"])
            self.assertTrue(payload["solver"]["python"]["available"])
        self.assertEqual(payload["postprocessing"]["linear_algebra"]["scipy"]["name"],"scipy-lapack")

    def test_solver_facade_does_not_eagerly_import_numerical_engines(self):
        code=(
            "import json,sys,solver; "
            "print(json.dumps({'scipy_integrate':'scipy.integrate' in sys.modules,"
            "'ode':'solver.ode' in sys.modules,'pde':'solver.pde' in sys.modules,"
            "'hierarchy':'solver.hierarchy' in sys.modules,'solve':callable(solver.solve)}))"
        )
        completed=self._run("-c",code)
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload,{"scipy_integrate":False,"ode":False,"pde":False,
                                  "hierarchy":False,"solve":True})

    def test_postprocessing_facade_keeps_method_families_lazy(self):
        code=(
            "import json,sys,postprocessing; "
            "print(json.dumps({'scipy_interpolate':'scipy.interpolate' in sys.modules,"
            "'archiver_core':'archiver.core' in sys.modules,"
            "'processor':'postprocessing.processor' in sys.modules,"
            "'pca':'postprocessing.pca_pod.algorithms' in sys.modules,"
            "'load_solution':callable(postprocessing.load_solution)}))"
        )
        completed=self._run("-c",code)
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload,{"scipy_interpolate":False,"archiver_core":False,"processor":False,
                                  "pca":False,"load_solution":True})

    def test_viewer_facade_does_not_eagerly_load_rendering_or_archives(self):
        code=(
            "import json,sys,viewer; "
            "print(json.dumps({'core':'viewer.core' in sys.modules,"
            "'archives':'archives' in sys.modules,"
            "'has_viewer':'Viewer' in dir(viewer),"
            "'has_routing':'classify_source' in dir(viewer)}))"
        )
        completed=self._run("-c",code)
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload,{"core":False,"archives":False,
                                  "has_viewer":True,"has_routing":True})

    def test_umbrella_facade_does_not_import_independent_subsystems(self):
        code=(
            "import json,sys,dynsyslab; "
            "print(json.dumps({'solver':'solver' in sys.modules,"
            "'control':'control' in sys.modules,'jobs':'job_manager' in sys.modules,"
            "'viewer':'viewer' in sys.modules,'has_solver':'solver' in dir(dynsyslab),"
            "'has_timeseries':'TimeSeries' in dir(dynsyslab)}))"
        )
        completed=self._run("-c",code)
        self.assertEqual(completed.returncode,0,completed.stderr)
        payload=json.loads(completed.stdout)
        self.assertEqual(payload,{"solver":False,"control":False,"jobs":False,
                                  "viewer":False,"has_solver":True,"has_timeseries":True})


if __name__=="__main__":
    unittest.main()
