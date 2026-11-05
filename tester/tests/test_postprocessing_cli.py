from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from archives.store import serializable_system
from dynsyslib import load_system
from postprocessing.cli import main
from postprocessing.processor import processor


class PostprocessingCLITests(unittest.TestCase):
    def _source(self, directory: str, samples: int=22000) -> Path:
        path=Path(directory)/"source.npz"; system=load_system("rpe")
        time=np.arange(samples,dtype=float)*.001
        values=np.column_stack((1.+.01*np.cos(time),-.01*np.sin(time)))
        metadata={"system":system.name,"system_parameters":serializable_system(system),"uniform_time_steps":True,"delta_t":.001}
        np.savez(path,time=time,values=values,metadata_json=json.dumps(metadata))
        return path

    def test_describe_is_json_and_does_not_require_a_source(self):
        output=StringIO()
        with redirect_stdout(output): status=main(["describe"])
        payload=json.loads(output.getvalue())
        self.assertEqual(status,0)
        self.assertEqual(payload["format"],"dynsyslab.postprocessing-workflow/v1")
        self.assertEqual(payload["workflow_policy"],"trajectory-hdmd-sindy")
        self.assertEqual(payload["methods"][0]["settings"]["max_memory_mb"],8192.)
        self.assertIsNone(payload["methods"][0]["settings"]["svd_threads"])
        self.assertIsNone(payload["methods"][0]["settings"]["eigendecomposition_threads"])
        self.assertEqual(payload["methods"][0]["settings"]["linear_algebra_backend"],"auto")
        self.assertEqual(payload["methods"][1]["settings"]["null_space_backend"],"gram-eigh")
        self.assertEqual(payload["methods"][1]["settings"]["archive_storage"],"coefficient-first")
        self.assertFalse(payload["methods"][1]["settings"]["compute_training_values"])

    def test_plan_reports_full_reference_geometry_without_computing(self):
        with tempfile.TemporaryDirectory() as directory:
            source=self._source(directory); output=StringIO()
            with redirect_stdout(output): status=main(["plan",str(source)])
        payload=json.loads(output.getvalue())
        self.assertEqual(status,0); self.assertTrue(payload["valid"])
        self.assertEqual(payload["hdmd"]["requested_shape"],[20000,2000])
        self.assertEqual(payload["hdmd"]["realized_shape"],[20000,2000])
        self.assertEqual(payload["hdmd"]["retained_rank_upper_bound"],2000)
        self.assertIsNone(payload["hdmd"]["svd_threads"])
        self.assertIsNone(payload["hdmd"]["eigendecomposition_threads"])
        self.assertEqual(payload["sindy"]["feature_count"],15)
        self.assertEqual(payload["sindy"]["archive_storage"],"coefficient-first")
        self.assertFalse(payload["sindy"]["algorithm_materializes_training_values"])
        self.assertFalse(payload["sindy"]["publication_materializes_training_values"])
        self.assertGreater(payload["sindy"]["estimated_training_values_mib"],600.)
        self.assertGreater(payload["sindy"]["estimated_full_values_mib"],600.)
        self.assertLess(payload["sindy"]["estimated_coefficient_first_mib"],2.)

    def test_run_requires_explicit_execute_confirmation(self):
        errors=StringIO()
        with self.assertRaises(SystemExit) as raised, redirect_stderr(errors): main(["run","missing.npz"])
        self.assertEqual(raised.exception.code,2)
        self.assertIn("requires --execute",errors.getvalue())

    def test_plan_reports_full_storage_materialization_only_when_requested(self):
        with tempfile.TemporaryDirectory() as directory:
            source=self._source(directory); output=StringIO()
            with redirect_stdout(output): status=main(["plan",str(source),"--sindy-storage","full"])
        payload=json.loads(output.getvalue())
        self.assertEqual(status,0)
        self.assertEqual(payload["sindy"]["archive_storage"],"full")
        self.assertFalse(payload["sindy"]["algorithm_materializes_training_values"])
        self.assertTrue(payload["sindy"]["publication_materializes_training_values"])

    def test_workflow_description_reflects_reference_and_override_backends(self):
        workflow=processor.preset("trajectory-hdmd-sindy",sources=[],archive=False,
                                  sindy_settings={"null_space_backend":"svd"},svd_driver="gesdd",svd_threads=2,
                                  eigendecomposition_threads=3,rank=8)
        payload=workflow.describe()
        self.assertFalse(payload["archive"])
        self.assertEqual(payload["methods"][0]["settings"]["rank"],8)
        self.assertEqual(payload["methods"][0]["settings"]["svd_driver"],"gesdd")
        self.assertEqual(payload["methods"][0]["settings"]["svd_threads"],2)
        self.assertEqual(payload["methods"][0]["settings"]["eigendecomposition_threads"],3)
        self.assertEqual(payload["methods"][1]["settings"]["null_space_backend"],"svd")

    def test_backends_reports_capabilities_without_a_source(self):
        output=StringIO()
        with redirect_stdout(output): status=main(["backends"])
        payload=json.loads(output.getvalue())
        self.assertEqual(status,0)
        self.assertEqual(payload["format"],"dynsyslab.postprocessing-backends/v1")
        self.assertTrue(payload["linear_algebra"]["scipy"]["available"])
        self.assertFalse(payload["linear_algebra"]["torch-cuda"]["auto_admitted"])


if __name__=="__main__":
    unittest.main()
