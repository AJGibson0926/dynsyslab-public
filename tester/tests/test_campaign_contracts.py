from __future__ import annotations
import tempfile
import unittest
import json
from pathlib import Path
from campaigner import Campaigner, CampaignRequest
from archives.schema import ArchiveSchemaError

class CampaignContractTests(unittest.TestCase):
    def test_generic_bifurcation_design_uses_registered_system_and_preset(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="bifurcation",label="pendulum-test",system_name="nonlinear_pendulum",preset="holmes_chaotic",quality="low",parameters={"driver":{"kind":"sinusoid","frequency_hz":.33226853021918074,"amplitude":1.2},"sweep":{"target":"driver.amplitude","start":.9,"stop":1.5}},integration={"transient_cycles":2,"retained_cycles":4,"samples_per_cycle":16})
            campaign,collection=Campaigner(directory).design_request(request)
            systems=[spec for job in collection.jobs for spec in job.systems]
            self.assertEqual(campaign.metadata["system_name"],"nonlinear_pendulum")
            self.assertEqual(len(systems),31)
            self.assertEqual(systems[0]["driver"]["amplitude"],.9)
            self.assertEqual(systems[-1]["driver"]["amplitude"],1.5)

    def test_campaign_manifest_rejects_unversioned_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"bad.json"; path.write_text(json.dumps({"label":"bad"}),encoding="utf-8")
            with self.assertRaises(ArchiveSchemaError):
                __import__("campaigner").CampaignManifest.read(path)

    def test_parameter_sweep_can_target_a_system_parameter(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="parameter_sweep",label="parameter-test",system_name="nonlinear_pendulum",preset="holmes_chaotic",quality="low",parameters={"driver":{"kind":"sinusoid","frequency_hz":.3,"amplitude":1.},"sweep":{"target":"parameters.damping_n_m_s","values":[.1,.2,.3]}},integration={"time_span":[0.,1.],"delta_t":.1})
            _,collection=Campaigner(directory).design_request(request)
            values=[spec["parameters"]["damping_n_m_s"] for job in collection.jobs for spec in job.systems]
            self.assertEqual(values,[.1,.2,.3])

    def test_campaign_rejects_an_unaddressable_sweep_target(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="bifurcation",label="bad",system_name="nonlinear_pendulum",preset="holmes_chaotic",parameters={"driver":{"kind":"sinusoid","frequency_hz":.3},"sweep":{"target":"unknown.value","start":0.,"stop":1.}},integration={"time_span":[0.,1.],"delta_t":.1})
            with self.assertRaises(ValueError): Campaigner(directory).design_request(request)

    def test_generic_trajectory_ensemble_preserves_initial_conditions(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="trajectory_ensemble",label="ensemble",system_name="nonlinear_pendulum",preset="earth_small_angle",parameters={"initial_conditions":[[0.,0.],[.1,0.]]},integration={"time_span":[0.,1.],"delta_t":.1})
            _,collection=Campaigner(directory).design_request(request)
            initial=collection.jobs[0].systems[0]["hierarchy"]["clusters"][0]["initial_conditions"]
            self.assertEqual(initial,[[0.,0.],[.1,0.]])

    def test_campaign_carries_checkpoint_policy_as_job_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            policy={"enabled":True,"mode":"checkpoint","checkpoint_interval":10.}
            request=CampaignRequest(kind="bifurcation",label="policy",system_name="nonlinear_pendulum",preset="holmes_chaotic",parameters={"driver":{"kind":"sinusoid","frequency_hz":.3},"sweep":{"target":"driver.amplitude","start":.9,"stop":1.}},integration={"time_span":[0.,1.],"delta_t":.1,"checkpoint_policy":policy})
            campaign,collection=Campaigner(directory).design_request(request)
            resolved=campaign.metadata["checkpoint_policy"]
            self.assertEqual(resolved["checkpoint_interval"],10.)
            self.assertEqual(resolved["mode"],"checkpoint")
            self.assertEqual(collection.jobs[0].systems[0]["solver"]["history_policy"],resolved)

    def test_bifurcation_campaign_uses_fast_policy_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="bifurcation",label="default-policy",system_name="nonlinear_pendulum",preset="holmes_chaotic",parameters={"driver":{"kind":"sinusoid","frequency_hz":.3},"sweep":{"target":"driver.amplitude","start":.9,"stop":1.}},integration={"time_span":[0.,1.],"delta_t":.1})
            campaign,collection=Campaigner(directory).design_request(request)
            policy=campaign.metadata["checkpoint_policy"]
            self.assertEqual(policy,{})
            self.assertEqual(campaign.metadata["checkpoint_policy_profile"],"fast")
            self.assertEqual(collection.jobs[0].systems[0]["solver"]["history_policy"],policy)

    def test_bifurcation_can_select_recoverable_checkpoint_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="bifurcation",label="recoverable-policy",system_name="nonlinear_pendulum",preset="holmes_chaotic",parameters={"driver":{"kind":"sinusoid","frequency_hz":.3},"sweep":{"target":"driver.amplitude","start":.9,"stop":1.}},integration={"time_span":[0.,1.],"delta_t":.1,"checkpoint_policy_profile":"recoverable"})
            campaign,collection=Campaigner(directory).design_request(request)
            policy=campaign.metadata["checkpoint_policy"]
            self.assertTrue(policy["enabled"])
            self.assertEqual(policy["checkpoint_interval"],600.)
            self.assertEqual(campaign.metadata["checkpoint_policy_profile"],"recoverable")
            self.assertEqual(collection.jobs[0].systems[0]["solver"]["history_policy"],policy)

    def test_non_bifurcation_campaign_keeps_empty_default_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="trajectory_ensemble",label="generic-policy",system_name="nonlinear_pendulum",preset="earth_small_angle",parameters={"initial_conditions":[[0.,0.]]},integration={"time_span":[0.,1.],"delta_t":.1})
            campaign,collection=Campaigner(directory).design_request(request)
            self.assertEqual(campaign.metadata["checkpoint_policy"],{})
            self.assertEqual(collection.jobs[0].systems[0]["solver"]["history_policy"],{})

    def test_non_bifurcation_campaign_can_select_recoverable_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="trajectory_ensemble",label="recoverable-ensemble",system_name="nonlinear_pendulum",preset="earth_small_angle",parameters={"initial_conditions":[[0.,0.]]},integration={"time_span":[0.,1.],"delta_t":.1,"checkpoint_policy_profile":"recoverable"})
            campaign,collection=Campaigner(directory).design_request(request)
            self.assertEqual(campaign.metadata["checkpoint_policy_profile"],"recoverable")
            self.assertTrue(campaign.metadata["checkpoint_policy"]["enabled"])
            self.assertEqual(collection.jobs[0].systems[0]["solver"]["history_policy"],campaign.metadata["checkpoint_policy"])

    def test_checkpoint_policy_profiles_are_advertised_for_all_campaigns(self):
        profiles=Campaigner.checkpoint_policy_profiles()
        self.assertEqual(profiles["fast"],{})
        self.assertTrue(profiles["recoverable"]["enabled"])

    def test_campaign_rejects_unknown_dynsyslib_system(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="bifurcation",label="unknown",system_name="not-a-system",parameters={"driver":{"kind":"sinusoid","frequency_hz":1.},"sweep":{"target":"driver.amplitude","start":0.,"stop":1.}},integration={"time_span":[0.,1.],"delta_t":.1})
            with self.assertRaises(KeyError): Campaigner(directory).design_request(request)

    def test_campaign_rejects_unknown_preset(self):
        with tempfile.TemporaryDirectory() as directory:
            request=CampaignRequest(kind="bifurcation",label="preset",system_name="nonlinear_pendulum",preset="missing",parameters={"driver":{"kind":"sinusoid","frequency_hz":1.},"sweep":{"target":"driver.amplitude","start":0.,"stop":1.}},integration={"time_span":[0.,1.],"delta_t":.1})
            with self.assertRaises(KeyError): Campaigner(directory).design_request(request)
