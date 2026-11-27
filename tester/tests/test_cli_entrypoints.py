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
            "archive":"dynsyslab_cli archive recompress",
            "jobs":"campaigns",
            "estimate":"observability",
            "krom-training":"krom-training-plan",
        }
        # archive is native: use the CLI this checkout builds, not an older installed release.
        build=PROJECT_ROOT/"build"/"native-suite-fast"
        native={**os.environ,"DYNSYSLAB_INSTALL_DIR":str(build)} if (build/"bin"/"dynsyslab-cli.exe").is_file() else None
        for command, marker in expectations.items():
            with self.subTest(command=command):
                completed=self._run("-m","dynsyslab",command,"--help",
                                    **({"env":native} if command in {"archive","krom-training"} else {}))
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

if __name__=="__main__":
    unittest.main()
