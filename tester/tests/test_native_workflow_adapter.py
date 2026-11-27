from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest

from compute_workbench.contracts import solve_template
from compute_workbench.native_workflow import (
    NativeWorkflow, NativeWorkflowCancelled, NativeWorkflowError,
    NativeWorkflowResourceLimit, find_library,
)


ROOT = Path(__file__).resolve().parents[2]


class NativeWorkflowAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.library = find_library(ROOT)
        except FileNotFoundError as error:
            raise unittest.SkipTest(str(error)) from error

    def test_adapter_runs_native_solve_and_reports_progress(self):
        job = solve_template()
        job["label"] = "python-native-adapter"
        job["time"].update({"stop": 0.02, "output_step": 0.01})
        events = []
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = NativeWorkflow(self.library).run(
                job, root / "archives", progress_path=root / "progress.json",
                progress=events.append,
            )
            self.assertEqual(result["execution_backend"], "native-cpp")
            self.assertFalse(result["python_required"])
            self.assertEqual(result["status"], "complete")
            self.assertEqual(len(result["solution_archives"]), 1)
            self.assertTrue(Path(result["solution_archives"][0]).is_file())
            self.assertTrue(Path(result["hierarchy_archive"]).is_file())
            progress = json.loads((root / "progress.json").read_text(encoding="utf-8"))
            self.assertEqual(progress["state"], "complete")
            self.assertTrue(events)

    def test_adapter_plans_from_the_compiled_registry(self):
        job = solve_template()
        job["controller"]["kind"] = "pid"
        plan = NativeWorkflow(self.library).plan(job)
        self.assertEqual(plan["schema"], "dynsyslab.workflow-plan/v1")
        self.assertEqual(plan["controller"], "pid")
        self.assertEqual(plan["output_samples_per_trajectory"], 1001)

    def test_native_validation_error_keeps_status(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(NativeWorkflowError) as caught:
                NativeWorkflow(self.library).run(
                    {"schema": "wrong", "kind": "solve"}, directory)
        self.assertEqual(caught.exception.status, 1)

    def test_resource_limit_is_a_distinct_exception(self):
        job = solve_template()
        job["time"].update({"stop": 0.02, "output_step": 0.01})
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(NativeWorkflowResourceLimit) as caught:
                NativeWorkflow(self.library).run(
                    job, directory, maximum_memory_bytes=32)
        self.assertEqual(caught.exception.status, 6)

    def test_cancellation_callback_stops_native_execution(self):
        job = solve_template()
        job["time"].update({"stop": 10.0, "output_step": 0.001})
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(NativeWorkflowCancelled) as caught:
                NativeWorkflow(self.library).run(
                    job, directory, cancelled=lambda: True)
        self.assertEqual(caught.exception.status, 4)


if __name__ == "__main__":
    unittest.main()
