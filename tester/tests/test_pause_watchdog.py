from __future__ import annotations
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path
import tempfile
import unittest

from job_manager.pause import _pid_cpu_seconds, sample_worker_activity, watch_paused_jobs


class PauseWatchdogTests(unittest.TestCase):
    def _watch(self, event_time: datetime, *, grace: float=30.):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); (root/"worker.telemetry.jsonl").write_text(json.dumps({
                "event":"heartbeat", "job_id":"job-1", "time_utc":event_time.isoformat(),
                "simulation_time":2., "progress_fraction":.2,
            })+"\n",encoding="utf-8")
            rows={"job-1":{"status":"running","pause_after_current_solve":True,"pause_generation":1,"worker_pid":os.getpid()}}
            events=[]
            def status(job_id): return dict(rows.get(job_id,{}))
            def set_status(job_id, value, **details): rows[job_id].update({"status":value,**details})
            report=watch_paused_jobs(["job-1"],status=status,set_status=set_status,event=lambda *args,**kwargs: events.append((args,kwargs)),stamp=lambda:"now",telemetry_root=root,grace_seconds=grace)
            return report,rows,events

    def test_fresh_heartbeat_is_classified_as_finishing(self):
        report,rows,_=self._watch(datetime.now(timezone.utc))
        self.assertEqual(report["job-1"]["classification"],"finishing_current_trajectory")
        self.assertEqual(rows["job-1"]["status"],"running")

    def test_stale_live_worker_is_quarantined_not_released(self):
        report,rows,events=self._watch(datetime.now(timezone.utc)-timedelta(minutes=10))
        self.assertEqual(report["job-1"]["classification"],"rogue_suspected")
        self.assertEqual(rows["job-1"]["status"],"quarantined")
        self.assertTrue(any(args[1]=="pause_watchdog_quarantined" for args,_ in events))

    def test_activity_probe_reports_a_local_cpu_sample(self):
        rows={"job-1":{"worker_pid":os.getpid()}}
        report=sample_worker_activity(["job-1"],status=lambda job_id:rows[job_id],sample_seconds=0.)
        self.assertIsNotNone(_pid_cpu_seconds(os.getpid()))
        self.assertIn("actively_computing",report["job-1"])

    def test_watch_ignores_jobs_without_pause_request(self):
        report=watch_paused_jobs(["job-1"],status=lambda _: {"status":"running"},set_status=lambda *args,**kwargs: None,event=lambda *args,**kwargs: None,stamp=lambda:"now",telemetry_root=".")
        self.assertEqual(report,{})

    def test_missing_worker_identity_awaits_a_heartbeat(self):
        with tempfile.TemporaryDirectory() as directory:
            rows={"job-1":{"status":"running","pause_after_current_solve":True}}
            report=watch_paused_jobs(["job-1"],status=lambda key:dict(rows[key]),set_status=lambda key,value,**details:rows[key].update({"status":value,**details}),event=lambda *args,**kwargs:None,stamp=lambda:"now",telemetry_root=directory)
            self.assertEqual(report["job-1"]["classification"],"awaiting_worker_heartbeat")

    def test_paused_status_is_recognized_as_clean(self):
        with tempfile.TemporaryDirectory() as directory:
            rows={"job-1":{"status":"paused","pause_after_current_solve":True,"worker_pid":os.getpid()}}
            report=watch_paused_jobs(["job-1"],status=lambda key:dict(rows[key]),set_status=lambda key,value,**details:rows[key].update({"status":value,**details}),event=lambda *args,**kwargs:None,stamp=lambda:"now",telemetry_root=directory)
            self.assertEqual(report["job-1"]["classification"],"paused_cleanly")

    def test_activity_probe_handles_unknown_worker_without_error(self):
        report=sample_worker_activity(["missing"],status=lambda _: {},sample_seconds=0.)
        self.assertIsNone(report["missing"]["actively_computing"])

    def test_telemetry_evidence_is_preserved_in_watchdog_report(self):
        report,_,_=self._watch(datetime.now(timezone.utc))
        self.assertEqual(report["job-1"]["latest_telemetry"]["simulation_time"],2.)
