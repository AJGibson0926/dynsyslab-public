from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from dynsyslib import load_system
from postprocessing.core import Trajectory
from postprocessing.workflow import SourceResolver, WorkflowContext
from postprocessing.processor import PostprocessingPipeline, processor
from postprocessing.policies import TRAJECTORY_HDMD_SINDY_POLICY
from postprocessing.policies import RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY
from time_series import TimeSeries


class PostprocessingWorkflowTests(unittest.TestCase):
    def setUp(self):
        time=np.linspace(0.,1.,5)
        self.series=TimeSeries(time,np.column_stack((time,time**2)))
        self.source=Trajectory(self.series,{"uniform_time_steps":True},None)

    def test_source_resolver_flattens_in_memory_hierarchy_with_provenance(self):
        system=load_system("rpe")
        trajectory=SimpleNamespace(time_series=self.series)
        cluster=SimpleNamespace(trajectories=[trajectory])
        hierarchy=SimpleNamespace(dynamical_system=system,clusters=[cluster])
        resolved=SourceResolver().resolve([hierarchy])
        self.assertEqual(len(resolved),1)
        self.assertIs(resolved[0].time_series,self.series)
        self.assertEqual(resolved[0].metadata["system"],system.name)
        self.assertTrue(resolved[0].metadata["uniform_time_steps"])

    def test_conventional_pipeline_class_preserves_lowercase_compatibility_alias(self):
        self.assertIs(processor,PostprocessingPipeline)
        self.assertIsInstance(PostprocessingPipeline(archive=False),PostprocessingPipeline)

    def test_standard_policy_is_named_immutable_and_copied_per_pipeline(self):
        self.assertEqual(TRAJECTORY_HDMD_SINDY_POLICY.name,"trajectory-hdmd-sindy")
        with self.assertRaises(TypeError):
            TRAJECTORY_HDMD_SINDY_POLICY.hdmd["rank"]=12
        first=PostprocessingPipeline.preset(); second=PostprocessingPipeline.preset()
        first.methods[0].settings["rank"]=12
        self.assertIsNone(second.methods[0].settings["rank"])
        self.assertEqual(first.workflow_policy,"trajectory-hdmd-sindy")
        self.assertEqual(first.methods[1].settings["archive_storage"],"coefficient-first")

    def test_rpe_clean_isolation_policy_preserves_validated_sindy_conventions(self):
        policy=RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY
        self.assertEqual(policy.name,"rpe-clean-isolation-hdmd-sindy")
        self.assertEqual(policy.hdmd["mode_sort"],"terminal-amplitude-energy")
        self.assertIsNone(policy.hdmd["rank"])
        self.assertEqual(policy.sindy["polynomial_degree"],4)
        self.assertFalse(policy.sindy["normalize"])
        self.assertEqual(policy.sindy["threshold"],1e-6)
        self.assertEqual(policy.sindy["max_iterations"],1)
        self.assertFalse(policy.sindy["compute_training_values"])
        with self.assertRaises(TypeError): policy.sindy["threshold"]=.05

    def test_source_resolver_rejects_cycles_explicitly(self):
        cyclic=SimpleNamespace()
        cyclic.systems=[cyclic]
        with self.assertRaisesRegex(ValueError,"Cyclic post-processing source graph"):
            SourceResolver().resolve([cyclic])

    def test_preparation_controls_are_removed_without_mutating_attachment(self):
        converted=Trajectory(self.series,{"coordinate_system":"dimensionless"},None)
        settings={"analysis_coordinates":"integration","rank":3}
        context=WorkflowContext()
        with patch("postprocessing.workflow.integration_coordinates",return_value=converted) as transform:
            first=context.prepare(self.source,settings)
            second=context.prepare(self.source,settings)
        self.assertIs(first.source,converted)
        self.assertIs(second.source,converted)
        self.assertEqual(first.settings,{"rank":3})
        self.assertEqual(settings,{"analysis_coordinates":"integration","rank":3})
        transform.assert_called_once_with(self.source)

    def test_invalid_resampling_control_is_not_silently_ignored(self):
        with self.assertRaisesRegex(TypeError,"ResamplingOptions"):
            WorkflowContext().prepare(self.source,{"resampling":object()})

    def test_spectrum_and_archive_lineage_advance_atomically(self):
        context=WorkflowContext(); first=object(); second=object()
        context.observe(SimpleNamespace(spectrum=first),Path("first.npz"))
        self.assertIs(context.spectrum,first)
        self.assertEqual(context.spectrum_archive,Path("first.npz"))
        context.observe(SimpleNamespace(spectrum=second),None)
        self.assertIs(context.spectrum,second)
        self.assertIsNone(context.spectrum_archive)


if __name__=="__main__":
    unittest.main()
