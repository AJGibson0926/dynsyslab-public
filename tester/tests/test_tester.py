from __future__ import annotations
import tempfile
import unittest
import json
from pathlib import Path
from tester import Tester

class TesterContractTests(unittest.TestCase):
    def test_reports_selected_group_and_writes_both_formats(self):
        with tempfile.TemporaryDirectory() as directory:
            report=Tester(report_root=directory).run("architecture")
            self.assertTrue(report.successful)
            self.assertGreaterEqual(report.tests_run,1)
            self.assertTrue(Path(report.report_json).is_file())
            self.assertTrue(Path(report.report_markdown).is_file())

    def test_all_is_an_explicit_available_selection(self):
        self.assertIn("native_workflow",Tester.available_groups())
        suite,selection,_=Tester()._suite(("all",),())
        self.assertEqual(selection,("all",))
        self.assertGreater(suite.countTestCases(),10)

    def test_unknown_group_is_rejected_with_actionable_error(self):
        with self.assertRaises(KeyError): Tester().run("not-a-test-group",write_report=False)

    def test_run_all_is_an_explicit_api(self):
        tester=Tester()
        self.assertEqual(tester.run_all.__name__,"run_all")

    def test_explicit_unittest_id_runs_one_requested_test(self):
        with tempfile.TemporaryDirectory() as directory:
            report=Tester(report_root=directory).run(test_ids=["tester.tests.test_tester.TesterContractTests.test_run_all_is_an_explicit_api"])
            self.assertTrue(report.successful); self.assertEqual(report.tests_run,1)

    def test_report_can_be_kept_in_memory_only(self):
        report=Tester().run("architecture",write_report=False)
        self.assertTrue(report.successful); self.assertIsNone(report.report_json); self.assertIsNone(report.report_markdown)

    def test_multiple_groups_are_combined(self):
        with tempfile.TemporaryDirectory() as directory:
            report=Tester(report_root=directory).run(["architecture","front_door"])
            self.assertTrue(report.successful); self.assertEqual(report.selection,("architecture","front_door")); self.assertGreaterEqual(report.tests_run,10)

    def test_report_json_contains_machine_readable_result_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            report=Tester(report_root=directory).run("architecture")
            payload=json.loads(Path(report.report_json).read_text(encoding="utf-8"))
            self.assertTrue(payload["successful"]); self.assertEqual(payload["selection"],["architecture"])

    def test_markdown_report_contains_runner_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            report=Tester(report_root=directory).run("architecture")
            content=Path(report.report_markdown).read_text(encoding="utf-8")
            self.assertIn("dynsyslab test report",content); self.assertIn("Tests:",content)
