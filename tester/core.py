"""Selection, execution, and durable reporting for automated tests."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from time import perf_counter
from typing import Iterable
import json
import platform
import sys
import unittest

TEST_GROUPS={
    "archive": ("test_archive_contracts.py","test_gateway_contracts.py"),
    "campaigner": ("test_campaign_contracts.py",),
    "control_estimation": ("test_control_estimation.py",),
    "job_manager": ("test_job_claims.py","test_pause_watchdog.py","test_history_policy.py"),
    "history_policy": ("test_history_policy.py",),
    "pause_watchdog": ("test_pause_watchdog.py",),
    "solver": ("test_solver_checkpoint_pause.py",),
    "viewer": ("test_viewer_contracts.py",),
}

@dataclass(frozen=True)
class TestReport:
    started_utc: str
    completed_utc: str
    selection: tuple[str,...]
    test_ids: tuple[str,...]
    tests_run: int
    failures: int
    errors: int
    skipped: int
    successful: bool
    elapsed_seconds: float
    python: str
    platform: str
    output: str
    report_json: str|None=None
    report_markdown: str|None=None

    def payload(self): return asdict(self)

class Tester:
    """Run registered groups, all tests, or explicit unittest test IDs.

    The runner does not alter production archives.  Its only writes are compact
    JSON/Markdown reports in a caller-selected report directory.
    """
    def __init__(self, *, project_root: str|Path|None=None, report_root: str|Path|None=None):
        self.project_root=Path(project_root or Path(__file__).resolve().parents[1]).resolve()
        self.test_root=self.project_root/"tester"/"tests"
        self.report_root=Path(report_root or self.project_root/"test_reports").resolve()
    @staticmethod
    def available_groups() -> dict[str,tuple[str,...]]:
        return {key:tuple(value) for key,value in TEST_GROUPS.items()}
    def _suite(self, selection: Iterable[str], test_ids: Iterable[str]) -> tuple[unittest.TestSuite,tuple[str,...],tuple[str,...]]:
        selected=tuple(selection)
        explicit=tuple(test_ids)
        if explicit:
            return unittest.defaultTestLoader.loadTestsFromNames(list(explicit)),selected,explicit
        if not selected or selected==("all",):
            return unittest.defaultTestLoader.discover(str(self.test_root),pattern="test_*.py",top_level_dir=str(self.project_root)),("all",),()
        unknown=[item for item in selected if item not in TEST_GROUPS]
        if unknown: raise KeyError(f"Unknown test group(s): {', '.join(unknown)}. Available: all, {', '.join(TEST_GROUPS)}.")
        suite=unittest.TestSuite()
        for group in selected:
            for filename in TEST_GROUPS[group]:
                suite.addTests(unittest.defaultTestLoader.discover(str(self.test_root),pattern=filename,top_level_dir=str(self.project_root)))
        return suite,selected,()
    def run(self, selection: str|Iterable[str]="all", *, test_ids: Iterable[str]=(), write_report: bool=True, verbosity: int=1) -> TestReport:
        groups=(selection,) if isinstance(selection,str) else tuple(selection)
        suite,groups,explicit=self._suite(groups,test_ids); started=datetime.now(timezone.utc); stream=StringIO(); timer=perf_counter()
        result=unittest.TextTestRunner(stream=stream,verbosity=verbosity).run(suite); completed=datetime.now(timezone.utc)
        report=TestReport(started.isoformat(),completed.isoformat(),groups,explicit,result.testsRun,len(result.failures),len(result.errors),len(result.skipped),result.wasSuccessful(),perf_counter()-timer,sys.version,platform.platform(),stream.getvalue())
        return self._write(report) if write_report else report
    def run_all(self, **kwargs) -> TestReport:
        """Explicit convenience entry point for the entire registered suite."""
        return self.run("all",**kwargs)
    def _write(self, report: TestReport) -> TestReport:
        self.report_root.mkdir(parents=True,exist_ok=True); stem=f"dynsyslab-tests-{report.completed_utc.replace(':','-').replace('+00:00','Z')}"
        json_path=self.report_root/f"{stem}.json"; markdown_path=self.report_root/f"{stem}.md"
        enriched=TestReport(**{**report.payload(),"report_json":str(json_path),"report_markdown":str(markdown_path)})
        json_path.write_text(json.dumps(enriched.payload(),indent=2,sort_keys=True),encoding="utf-8")
        summary="PASS" if enriched.successful else "FAIL"
        markdown_path.write_text(f"# dynsyslab test report — {summary}\n\n- Selection: `{', '.join(enriched.selection)}`\n- Tests: {enriched.tests_run}\n- Failures: {enriched.failures}; errors: {enriched.errors}; skipped: {enriched.skipped}\n- Elapsed: {enriched.elapsed_seconds:.3f} s\n- Started: {enriched.started_utc}\n- Completed: {enriched.completed_utc}\n\n## Runner output\n\n```text\n{enriched.output}\n```\n",encoding="utf-8")
        return enriched
