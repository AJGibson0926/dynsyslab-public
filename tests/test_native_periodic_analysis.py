"""Public periodic/torus CLI and canonical snapshot archive contracts."""
import json
import math
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
from dynsyslab_client.periodic import analyze_torus, sample, snapshot_pairs, predict_orbit, check_coefficient_shift
from dynsyslab_client.workflow import NativeWorkflow
from dynsyslab_client.periodic import evaluate_quartic_orbit, predict_controlled

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("DYNSYSLAB_NATIVE_BIN", ROOT / "build/native-suite-fast/bin"))
CLI = BIN / ("dynsyslab-cli.exe" if os.name == "nt" else "dynsyslab-cli")
LIB = BIN / ("dynsyslab_workflow_c_api.dll" if os.name == "nt" else "libdynsyslab_workflow_c_api.so")


@unittest.skipUnless(CLI.is_file(), "Build the native CLI first")
class PeriodicContracts(unittest.TestCase):
    def test_controlled_chart_and_reference(self):
        knots=[]
        for a in np.linspace(1,2,33):
            knots.append(dict(parameter=float(a),channels=[
                dict(angular_frequency=1,coefficients=[[-1,a/2,0],[1,a/2,0]]),
                dict(angular_frequency=1,coefficients=[[-1,0,-a/2],[1,0,a/2]])]))
        request=dict(system=dict(name='harmonic_oscillator'),knots=knots,parameter=1.4,
            initial_phase=.3,final_time=3,output_step=.01,table_samples=512,
            input_maximum_harmonic=63,input=dict(bias=.02,state_gain=[0,-.02]),
            maximum_step=.02)
        result=predict_controlled(request,executable=CLI)
        self.assertEqual(result['method'],'controlled-periodic-chart')
        np.testing.assert_allclose(np.array(result['predictions']).reshape(-1,2),
                                   np.array(result['reference']).reshape(-1,2),atol=3e-6)
        self.assertLess(max(result['diagnostics']['drift_relative_l2']),1e-6)
        request['input_maximum_harmonic']=256
        with self.assertRaises(RuntimeError):predict_controlled(request,executable=CLI)

    def test_quartic_parameter_clock_is_explicit(self):
        r=evaluate_quartic_orbit(dict(coefficients=[1,0,-1,0,0],initial=1,velocity=0,
            elapsed=[0,.3]),executable=CLI)
        self.assertIn("not necessarily",r["clock"])
        np.testing.assert_allclose(r["values"],[[1,0],[math.cos(.3),-math.sin(.3)]],atol=1e-11)

    def test_coefficient_shift_identity(self):
        result=check_coefficient_shift(dict(system=dict(name="morse_oscillator"),
            initial_state=[.7,.3], samples=512, maximum_harmonic=31,
            observables=["position","momentum"], shifts=[.13,-.7], initial_time=1.3), executable=CLI)
        self.assertFalse(result["certified"])
        for row in result["checks"]:
            self.assertLess(max(row["coefficient_l2_error"]),1e-10)
            self.assertLess(row["initial_state_error"],1e-12)

    def test_frozen_atlas_prediction_clock(self):
        s=dict(angular_frequency=2,coefficients=[[-1,.5,0],[1,.5,0]])
        result=predict_orbit(dict(knots=[dict(parameter=1,channels=[s]),dict(parameter=2,channels=[s])],
            parameter=1.4,initial_phase=.3,initial_time=7,time=[7,7.2],table_samples=8192),executable=CLI)
        np.testing.assert_allclose(np.array(result['predictions'])[:,0],[math.cos(.3),math.cos(.7)],atol=1e-7)

    def test_torus_qualification(self):
        request = dict(system=dict(name="coupled_duffing_hamiltonian"),
                       initial_state=[.7,.4,0,0], shape=[8,9], maximum_harmonics=[3,3])
        result = analyze_torus(request, executable=CLI)
        self.assertEqual(result["torus_dimension"], 2)
        self.assertTrue(result["sampled_checks_pass"])
        self.assertFalse(result["kam_certified"])
        request["system"]["parameters"] = {"coupling": .1}
        self.assertEqual(analyze_torus(request, executable=CLI)["status"], "not_applicable")

    def test_elementary_periodic_observable(self):
        d = sample(dict(system=dict(name="isotonic_oscillator"), initial_state=[1.4,0],
                        samples=64, observables=["position_squared"], maximum_harmonic=5), executable=CLI)
        self.assertAlmostEqual(d["angular_frequency"], 2)
        c = {int(n):complex(re,im) for n,re,im in d["fourier_series"][0]["coefficients"]}
        self.assertLess(sum(abs(v) for n,v in c.items() if abs(n)>1), 1e-12)

    def test_pairs_npz_fit_and_immutability(self):
        s = dict(angular_frequency=1.3, coefficients=[[-1,.5,0],[1,.5,0]])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"pairs.npz"
            request = dict(orbits=[dict(weight=2,channels=[s])], phase_count=64,
                           time_offsets=[0,-.4],sample_interval=.15,output_archive=str(path))
            response=snapshot_pairs(request,executable=CLI)
            self.assertEqual(response["dimension"],2)
            with np.load(path) as archive:
                x,y,w=archive["x"],archive["y"],archive["weights"]
                metadata=json.loads(str(archive["metadata_json"].item()))
            self.assertEqual(metadata["schema"],"dynsyslab.koopman-snapshot-pairs/v1")
            np.testing.assert_allclose(x[0],[1,math.cos(-.52)],atol=1e-14)
            np.testing.assert_allclose(y[0],[math.cos(.195),math.cos(-.325)],atol=1e-14)
            self.assertAlmostEqual(float(w.sum()),2)
            with self.assertRaises(RuntimeError):snapshot_pairs(request,executable=CLI)
            if LIB.is_file():
                native=NativeWorkflow(library=LIB)
                result=native.koopman_spectral(dict(data=dict(archive=str(path)),settings=dict(
                    method="mpedmd",assume_measure_preserving=True,sample_interval=.15,
                    dictionary=dict(kind="linear",constant=False),
                    rigged_measures=[dict(name="projected",reconstruction=dict(name="coordinates",output=0),
                                          angles=[.195],bandwidth=.05,order=2)])))
                self.assertEqual(result["rank"],2)
                self.assertIn("explicitly selected",result["spectral"]["measures"][0]["assumptions"])


if __name__=="__main__": unittest.main()
