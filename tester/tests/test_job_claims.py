from __future__ import annotations
import tempfile
import unittest
from job_manager.manager import JobManager
from job_manager.manifests import JobManifest
from archives.schema import ArchiveSchemaError

class JobClaimTests(unittest.TestCase):
    def test_live_claim_excludes_a_second_job_and_releases_cleanly(self):
        spec={"system_name":"test","preset":None,"parameters":{},"driver":{"kind":"zero"},"time_span":[0.,1.],"initial_state":[0.],"solver":{},"output":{}}
        with tempfile.TemporaryDirectory() as directory:
            manager=JobManager(directory)
            manager._set_status("first","claiming")
            self.assertTrue(manager._claim_specs("first",(spec,))[0])
            manager._set_status("second","claiming")
            acquired,conflicts=manager._claim_specs("second",(spec,))
            self.assertFalse(acquired); self.assertTrue(conflicts)
            manager._set_status("first","completed"); manager._release_spec_claims("first")
            self.assertTrue(manager._claim_specs("second",(spec,))[0])

    def test_pause_preserves_an_immutable_resume_snapshot(self):
        from job_manager.manifests import JobManifest
        with tempfile.TemporaryDirectory() as directory:
            manager=JobManager(directory); job=JobManifest(label="pause",systems=({"system_name":"test"},)); manager.submit(job)
            manager.pause([job.job_id]); state=manager.status(job.job_id)
            self.assertEqual(state["status"],"paused")
            self.assertTrue(__import__("pathlib").Path(state["pause_snapshot_manifest"]).is_file())
            self.assertEqual(manager.resume([job.job_id]),[job.job_id])
            self.assertEqual(manager.status(job.job_id)["status"],"pending")

    def test_job_manifest_rejects_missing_format(self):
        import json
        from pathlib import Path
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"bad.json"; path.write_text(json.dumps({"label":"bad","systems":[]}),encoding="utf-8")
            with self.assertRaises(ArchiveSchemaError): JobManifest.read(path)

    def test_progress_estimate_uses_completed_job_durations(self):
        with tempfile.TemporaryDirectory() as directory:
            manager=JobManager(directory)
            manager._set_status("done","completed",duration_seconds=10.,planned_system_count=2)
            manager._set_status("next","pending",planned_system_count=4)
            summary=manager.progress()
            self.assertEqual(summary["mean_completed_system_seconds"],5.)
            self.assertEqual(summary["estimated_remaining_seconds"],20.)

    def test_pause_is_idempotent_for_an_already_paused_job(self):
        with tempfile.TemporaryDirectory() as directory:
            manager=JobManager(directory); job=JobManifest(label="idempotent",systems=({"system_name":"test"},)); manager.submit(job)
            self.assertEqual(manager.pause([job.job_id])["paused"],[job.job_id])
            self.assertEqual(manager.pause([job.job_id])["paused"],[])

    def test_progress_for_empty_manager_is_zero_and_unestimated(self):
        with tempfile.TemporaryDirectory() as directory:
            summary=JobManager(directory).progress()
            self.assertEqual(summary["total"],0); self.assertEqual(summary["fraction_complete"],0.)
            self.assertIsNone(summary["estimated_remaining_seconds"])

    def test_state_is_transactional_sqlite_and_legacy_json_is_imported_once(self):
        import json
        from pathlib import Path
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); legacy={"format":"legacy","jobs":{"old":{"status":"completed"}},"spec_claims":{}}
            (root/"state.json").write_text(json.dumps(legacy),encoding="utf-8")
            manager=JobManager(root)
            self.assertTrue((root/"state.sqlite").is_file())
            self.assertEqual(manager.status("old")["status"],"completed")
            manager._set_status("new","pending")
            self.assertEqual(JobManager(root).status("new")["status"],"pending")
