"""Native periodic spectral workflow, frozen archives, and Viewer contracts."""
import math
import json
import os
from pathlib import Path
import tempfile
import unittest
from dynsyslab_client.periodic import spectral_sweep,spectral_compare,quadrature_audit,fit_delay_model,audit_delay_forecast
from viewer.periodic_spectra import comparison_payload,write_comparison_dashboard,render_comparison_pages

ROOT=Path(__file__).resolve().parents[1]
BIN=Path(os.environ.get('DYNSYSLAB_NATIVE_BIN',ROOT/'build/native-suite-fast/bin'))
CLI=BIN/('dynsyslab-cli.exe' if os.name=='nt' else 'dynsyslab-cli')

def request():
    s=dict(angular_frequency=1,coefficients=[[-1,.5,0],[1,.5,0]])
    v=dict(angular_frequency=1,coefficients=[[-1,0,.5],[1,0,-.5]])
    return dict(training=[dict(channels=[s,v],weight=1)],validation=[dict(channels=[s,v],weight=1,parameter=1,family='circle')],
        groups=[dict(name='single',channels=[0]),dict(name='joint',channels=[0,1])],delay_counts=[4,6],phase_count=32,
        assume_measure_preserving=True,sample_interval=.2,bandwidth=.08,order=4,
        angles=[-math.pi+2*math.pi*j/64 for j in range(64)],packet_centers=[.2,.4],channel_names=['position','velocity'],audit_phase_counts=[16,32])

@unittest.skipUnless(CLI.is_file(),'Build native CLI first')
class PeriodicSpectralContracts(unittest.TestCase):
    def test_delay_fit_methods_and_immutable_provenance(self):
        import numpy as np
        from dynsyslab_client.workflow import NativeWorkflow
        api=NativeWorkflow(BIN/('dynsyslab_workflow_c_api.dll' if os.name=='nt' else 'libdynsyslab_workflow_c_api.so'))
        with tempfile.TemporaryDirectory() as directory:
            for method in ('edmd','resdmd','mpedmd'):
                path=Path(directory)/(method+'.npz')
                q=dict(training=request()['training'],channels=[0,1],time_offsets=[0],phase_count=32,
                    sample_interval=.2,method=method,storage='streaming_qr',assume_measure_preserving=True,
                    output_archive=str(path),source='independent rotation test')
                result=fit_delay_model(q,executable=CLI)
                self.assertEqual(result['fitting_method'],method);self.assertEqual(result['rank'],2)
                with np.load(path) as archive:
                    provenance=json.loads(archive['periodic_fit_provenance_json'].item())
                    self.assertEqual(provenance['source'],q['source'])
                    self.assertEqual(provenance['training_sha256'],result['provenance']['training_sha256'])
                prediction=api.koopman_spectral(dict(action='rollout',archive=str(path),initial=[1,0],steps=7))
                expected=np.column_stack((np.cos(.2*np.arange(8)),np.sin(.2*np.arange(8))))
                np.testing.assert_allclose(np.array(prediction['values']).reshape(-1,2),expected,atol=2e-12)
                audit=audit_delay_forecast(dict(archive=str(path),validation=q['training'],steps_per_orbit=[7]),executable=CLI)
                row=audit['orbits'][0]
                np.testing.assert_allclose(row['reference_mean_square'],[.5,.5],atol=1e-13)
                self.assertLess(max(row['error_mean_square']),1e-24)
                np.testing.assert_allclose(row['prediction_mean_square'],.5,atol=2e-12)
                np.testing.assert_allclose(np.asarray(row['persistence_error_mean_square']).reshape(-1,2),
                                          np.repeat((1-np.cos(.2*np.arange(8)))[:,None],2,axis=1),atol=1e-13)
                sparse=audit_delay_forecast(dict(archive=str(path),validation=q['training'],evaluation_steps=[[0,1,7,100]]),executable=CLI)
                self.assertEqual(sparse['orbits'][0]['step_indices'],[0,1,7,100])
                self.assertLess(max(sparse['orbits'][0]['error_mean_square']),1e-23)
                with self.assertRaises(RuntimeError):
                    audit_delay_forecast(dict(archive=str(path),validation=q['training'],evaluation_steps=[[0,7]],steps_per_orbit=[7]),executable=CLI)
                with self.assertRaises(RuntimeError):
                    audit_delay_forecast(dict(archive=str(path),validation=q['training'],steps_per_orbit=[]),executable=CLI)
                with self.assertRaises(RuntimeError):fit_delay_model(q,executable=CLI)
            q.pop('output_archive');q['method']='tlsdmd'
            with self.assertRaises(RuntimeError):fit_delay_model(q,executable=CLI)
            q['method']='mpedmd';q['assume_measure_preserving']=False
            with self.assertRaises(RuntimeError):fit_delay_model(q,executable=CLI)

    def test_streaming_storage_provenance(self):
        q=request();q['storage']='streaming_qr'
        result=spectral_sweep(q,executable=CLI)
        for model in result['models']:
            contract=json.loads(model['coordinate_contract'])
            self.assertEqual(contract['storage'],'streaming_qr')
            self.assertEqual(contract['original_samples'],32)
            self.assertLess(contract['fit_rows'],32)
        for row in result['comparisons']:
            self.assertLess(max(row['packet_complex_relative_l2']),1e-8)
        q['storage']='bad'
        with self.assertRaises(RuntimeError):spectral_sweep(q,executable=CLI)

    def test_registered_analytic_sampling_adapter(self):
        from postprocessing.koopman.periodic_spectral import sample_bank
        bank=sample_bank(dict(name='isotonic_oscillator'),[[1.4,0]],['position_squared'],
            weights=[2],samples=64,maximum_harmonic=5,families=['radial'],parameters=[1.4],executable=CLI)
        self.assertAlmostEqual(bank[0]['channels'][0]['angular_frequency'],2)
        self.assertEqual(bank[0]['family'],'radial')
        self.assertEqual(bank[0]['parameter'],1.4)
        self.assertIn('equations_fingerprint',bank[0]['source'])
        with self.assertRaises(ValueError):sample_bank({},[[0]],['x'],weights=[],samples=64,maximum_harmonic=5,executable=CLI)

    def test_frozen_sweep_and_separate_viewer(self):
        with tempfile.TemporaryDirectory() as directory:
            q=request();q['output_directory']=str(Path(directory)/'models')
            result=spectral_sweep(q,executable=CLI)
            self.assertEqual(len(result['models']),4)
            for row in result['comparisons']:
                self.assertLess(row['density_relative_l2'],1e-9)
                self.assertLess(max(row['packet_complex_relative_l2']),1e-8)
            for model in result['models']:self.assertTrue(Path(model['archive']).is_file())
            frozen=spectral_compare(dict(q,models=result['models']),executable=CLI)
            self.assertEqual(result['comparisons'],frozen['comparisons'])
            with self.assertRaises(RuntimeError):spectral_sweep(q,executable=CLI)
            payload=comparison_payload(result)
            self.assertEqual(len(payload['page_labels']),4)
            self.assertEqual(payload['figures'][0]['value_scale'],'symlog')
            self.assertEqual(payload['figures'][1]['value_scale'],'linear')
            path=write_comparison_dashboard(result,Path(directory)/'viewer')
            with self.assertRaises(FileExistsError):write_comparison_dashboard(result,path.parent)
            renderer=BIN/('DynSysLabViewerCli.exe' if os.name=='nt' else 'DynSysLabViewerCli')
            if renderer.is_file():
                outputs=render_comparison_pages(path,executable=renderer,pages=[0],width=1100,height=550)
                self.assertGreater(outputs[0].stat().st_size,1000)

    def test_resolution_and_clock_rejections(self):
        q=request();fast=dict(angular_frequency=1,coefficients=[[-19,.5,0],[19,.5,0]])
        q['training'][0]['channels']=[fast,fast]
        with self.assertRaises(RuntimeError):spectral_sweep(q,executable=CLI)
        audit=quadrature_audit(dict(training=q['training'],audit_phase_counts=[19,40]),executable=CLI)['audit']
        self.assertEqual(audit['alias_free_phase_count'],39)
        self.assertAlmostEqual(audit['sampled_mass'][0],1)
        self.assertAlmostEqual(audit['sampled_mass'][2],.5)
        q=request();q['include_models']=True;result=spectral_sweep(q,executable=CLI)
        q['models']=result['models'];q['models'][0]['time_offsets'][0]-=.1
        with self.assertRaises(RuntimeError):spectral_compare(q,executable=CLI)
        q['models'][0]['time_offsets'][0]+=.1;q['sample_interval']=.3
        with self.assertRaises(RuntimeError):spectral_compare(q,executable=CLI)

    def test_validation_weights_and_empty_rejection(self):
        q=request();q['validation'][0]['weight']=0
        with self.assertRaises(RuntimeError):spectral_sweep(q,executable=CLI)
        q=request();q['groups']=[]
        with self.assertRaises(RuntimeError):spectral_sweep(q,executable=CLI)

if __name__=='__main__':unittest.main()
