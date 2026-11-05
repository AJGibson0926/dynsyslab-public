from __future__ import annotations
import tempfile, unittest
from pathlib import Path
from job_manager import DeadlinePolicy, Dependency, JobArray, JobManager, JobManifest, ResourceRequest, SchedulerPolicy, StagingPlan
from job_manager.scheduling import admission_reason, dependency_satisfied
from job_manager.staging import stage_inputs
from job_manager.backends import pbs_script, slurm_script
from job_manager.backends.local import LocalProcessBackend
from job_manager.resource_limits import THREAD_VARIABLES, worker_thread_environment
from unittest.mock import patch

class SchedulerContractTests(unittest.TestCase):
    def test_resource_request_supports_accelerators_and_features(self):
        request=ResourceRequest(cores=4,memory_mb=8192,gpus=1,gpu_type="A100",features=("avx2",))
        self.assertEqual(request.payload()["gpus"],1)
    def test_invalid_resource_request_is_rejected(self):
        with self.assertRaises(ValueError): ResourceRequest(cores=0)
    def test_dependency_conditions_have_expected_terminal_semantics(self):
        self.assertTrue(dependency_satisfied("afterok","completed")); self.assertTrue(dependency_satisfied("afterany","failed")); self.assertTrue(dependency_satisfied("afternotok","interrupted"))
    def test_unknown_dependency_condition_is_rejected(self):
        with self.assertRaises(ValueError): Dependency("parent","later")
    def test_deadline_policy_validates_nonnegative_requeue_limit(self):
        with self.assertRaises(ValueError): DeadlinePolicy(max_requeues=-1)
    def test_staging_plan_defaults_to_primary_output(self):
        self.assertEqual(StagingPlan().output_location,"primary")
    def test_submission_records_resource_and_reproducibility_metadata(self):
        job=JobManifest(label="scheduler",systems=({"label":"spec"},),resources={"cores":2,"gpus":0},staging={"inputs":["input.dat"]},deadline={"max_requeues":1})
        with tempfile.TemporaryDirectory() as directory:
            manager=JobManager(directory); manager.submit(job); row=manager.status(job.job_id)
            self.assertEqual(row["resources"]["cores"],2); self.assertIn("python",row["provenance"]); self.assertEqual(row["staging"]["inputs"],["input.dat"])
            self.assertEqual(row["provenance"]["numerical_thread_policy"]["threads_per_process"],2)
    def test_worker_environment_caps_nested_numerical_threads(self):
        environment=worker_thread_environment(3,{"PATH":"kept","OMP_NUM_THREADS":"99"})
        self.assertEqual(environment["PATH"],"kept")
        self.assertEqual(environment["DYNSYSLAB_ALLOCATED_CORES"],"3")
        self.assertTrue(all(environment[name]=="3" for name in THREAD_VARIABLES))
    def test_local_process_receives_job_specific_thread_environment(self):
        completed=type("Completed",(),{"returncode":0,"stdout":"ok","stderr":""})()
        context={"manager_root":"manager","project_root":"project","resources_by_job":{"job":{"cores":2}}}
        with patch("job_manager.backends.local.subprocess.run",return_value=completed) as run:
            result=LocalProcessBackend().execute("job",execute_local=lambda _: {},context=context)
        self.assertEqual(run.call_args.kwargs["env"]["OPENBLAS_NUM_THREADS"],"2")
        self.assertEqual(result.details["numerical_thread_policy"]["allocated_cores"],2)
    def test_admission_detects_core_limit(self):
        reason=admission_reason(ResourceRequest(cores=2),running=[{"status":"running","resources":{"cores":2}}],policy=SchedulerPolicy(max_cores=3))
        self.assertEqual(reason,"max_cores")
    def test_array_submission_assigns_stable_indices(self):
        jobs=tuple(JobManifest(label=f"array-{index}",systems=({"label":f"spec-{index}"},)) for index in range(2))
        with tempfile.TemporaryDirectory() as directory:
            manager=JobManager(directory); path=manager.submit_array(JobArray("array",jobs,max_concurrent=2)); payload=__import__("json").loads(path.read_text())
            self.assertEqual(payload["max_concurrent"],2)
            self.assertEqual(manager.status(jobs[1].job_id)["tags"]["array_index"],1)
    def test_staging_copies_and_hash_verifies_declared_input(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/"input.txt"; source.write_text("scientific input",encoding="utf-8")
            job=JobManifest(label="stage",systems=({"label":"spec"},),staging={"inputs":[str(source)]})
            manifest=stage_inputs(job,Path(directory)/"scratch")
            self.assertTrue(Path(manifest["inputs"][0]["staged"]).is_file())
            self.assertEqual(len(manifest["inputs"][0]["sha256"]),64)
    def test_remote_scripts_translate_resources_and_requeue_policy(self):
        job=JobManifest(label="remote",systems=({"label":"spec"},),resources={"cores":4,"memory_mb":2048,"gpus":1,"wall_time_s":3600},deadline={"checkpoint_before_deadline_s":120,"requeue_on_timeout":True})
        self.assertIn("#SBATCH --gpus=1",slurm_script(job)); self.assertIn("#SBATCH --requeue",slurm_script(job)); self.assertIn("#PBS -l select=1:ncpus=4",pbs_script(job))
    def test_deadline_requeue_enables_checkpoint_capable_solver_metadata(self):
        job=JobManifest(label="deadline",systems=({"solver":{"options":{},"history":{}},},),resources={"wall_time_s":60},deadline={"checkpoint_before_deadline_s":10,"requeue_on_timeout":True})
        with tempfile.TemporaryDirectory() as directory:
            options,history,_=JobManager(directory)._options(job)
            self.assertEqual(options.adaptive_engine,"low_level_stepper"); self.assertEqual(history.mode,"checkpoint"); self.assertTrue(history.pause_at_checkpoint)
    def test_aftercorr_is_a_supported_dependency_contract(self):
        self.assertTrue(dependency_satisfied("aftercorr","completed"))
        self.assertFalse(dependency_satisfied("aftercorr","failed"))
    def test_scheduler_policy_accepts_qos_and_reservation_contracts(self):
        policy=SchedulerPolicy(allowed_qos=("normal","urgent"),reservation_windows={"night":("20:00","08:00")})
        self.assertIn("urgent",policy.allowed_qos); self.assertIn("night",policy.reservation_windows)
