"""Command line interface for :mod:`tester`."""
from __future__ import annotations
import argparse
from .core import Tester

def main() -> int:
    parser=argparse.ArgumentParser(description="Run dynsyslab automated tests and write a report.")
    parser.add_argument("selection",nargs="*",default=["all"],help="Test groups, or 'all'.")
    parser.add_argument("--all",action="store_true",help="Explicitly run the complete test suite.")
    parser.add_argument("--test-id",action="append",default=[],help="Explicit unittest test ID; may be repeated.")
    parser.add_argument("--report-root")
    parser.add_argument("--verbosity",type=int,default=1)
    args=parser.parse_args(); report=Tester(report_root=args.report_root).run("all" if args.all else args.selection,test_ids=args.test_id,verbosity=args.verbosity)
    print(f"{'PASS' if report.successful else 'FAIL'}: {report.tests_run} tests in {report.elapsed_seconds:.3f}s")
    print(report.report_markdown)
    return 0 if report.successful else 1
