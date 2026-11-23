"""Compiled packet map and nonlinear history-feature identification contracts."""
import math
from pathlib import Path
import unittest
from dynsyslab_client.periodic import identify_history_features, packet_features, fit_delay_model

ROOT=Path(__file__).resolve().parents[1]
CLI=ROOT/'build/native-suite-fast/bin/dynsyslab-cli.exe'

@unittest.skipUnless(CLI.is_file(),'native CLI unavailable')
class PeriodicFeatures(unittest.TestCase):
    def test_compiled_identification_and_causal_rejection(self):
        knots=[dict(parameter=p,channels=[dict(angular_frequency=p,coefficients=[[-1,p/2,0],[1,p/2,0]])]) for p in [1,1.25,1.5,1.75,2]]
        offsets=[-1.3,-.4,-.1,0]
        f=dict(channels=[0],time_offsets=offsets,outputs=4,coefficients=[float(i==j) for i in range(4) for j in range(4)])
        q=dict(knots=knots,features=f,histories=[1.42*math.cos(2.3+1.42*t) for t in offsets],
               noise_standard_deviations=[.1]*4,phase_seeds=64,forecast_times=[0,1])
        result=identify_history_features(q,executable=CLI)
        self.assertEqual(result['noise_rank'],4)
        self.assertAlmostEqual(result['estimates'][0]['parameters'][0],1.42,places=5)
        self.assertAlmostEqual(result['estimates'][0]['parameters'][1],2.3,places=5)
        self.assertEqual(len(result['estimates'][0]['forecast']),2)
        self.assertEqual(len(result['features_sha256']),64)
        q['features']['time_offsets'][-1]=.01
        with self.assertRaises(RuntimeError):identify_history_features(q,executable=CLI)

    def test_packet_contract(self):
        c=dict(angular_frequency=1.2,coefficients=[[-1,.5,0],[1,.5,0]])
        q=dict(training=[dict(weight=1,channels=[c])],channels=[0],time_offsets=[-.4,0,.3],
               phase_count=64,sample_interval=.05,method='mpedmd',assume_measure_preserving=True,include_model=True)
        fit=fit_delay_model(q,executable=CLI)
        packet=packet_features(dict(model=fit['model'],channel=0,packet_centers=[.04,.12],bandwidth=.1,order=2),executable=CLI)
        self.assertEqual(packet['model_identity'],fit['model_identity'])
        self.assertEqual(packet['features']['time_offsets'],[-.7,-.3,0])
        self.assertEqual(packet['features']['outputs'],4)
        self.assertEqual(len(packet['features']['coefficients']),12)
        self.assertEqual(packet['sample_interval'],.05)

if __name__=='__main__':unittest.main()
