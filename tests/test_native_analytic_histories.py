"""Analytic history sampling, complex channels, native clocks and archive contracts."""
import json,math,os,tempfile,unittest
from pathlib import Path
import numpy as np
from dynsyslab_client.periodic import analytic_histories,evaluate_trajectory
ROOT=Path(__file__).resolve().parents[1]
CLI=Path(os.environ.get('DYNSYSLAB_NATIVE_BIN',ROOT/'build/native-suite-fast/bin'))/('dynsyslab-cli.exe' if os.name=='nt' else 'dynsyslab-cli')

def request():
    return dict(system=dict(name='nonlinear_pendulum',parameters={'gravity_m_s2':1,'damping_n_m_s':0}),
        initial_states=[[.4,.5],[0,3],[0,-3]],time_offsets=[-.4,0,.8],observables=['phase_sech','phase_sech_components'])

@unittest.skipUnless(CLI.is_file(),'Build native CLI first')
class AnalyticHistoryContracts(unittest.TestCase):
    def test_complex_channels_and_signed_clock(self):
        q=request();result=analytic_histories(q,executable=CLI)
        a=np.array(result['values']);z=(a[:,0]+1j*a[:,1]).reshape(3,3,3)
        self.assertEqual(result['channels'],3)
        for i,state in enumerate(q['initial_states']):
            truth=evaluate_trajectory(dict(system=q['system'],initial_state=state,time=q['time_offsets']),executable=CLI)
            x=np.array(truth['values']).reshape(3,2);expected=np.exp(1j*x[:,0])/np.cosh(x[:,1])
            np.testing.assert_allclose(z[i,:,0],expected,rtol=1e-13,atol=1e-13)
            np.testing.assert_allclose(z[i,:,1]+1j*z[i,:,2],expected,rtol=1e-13,atol=1e-13)
            self.assertAlmostEqual(abs(z[i,1,0]-complex(math.cos(state[0]),math.sin(state[0]))/math.cosh(state[1])),0,places=13)
        q['initial_time']=12.5
        shifted=analytic_histories(q,executable=CLI)
        np.testing.assert_allclose(shifted['values'],result['values'],atol=1e-13,rtol=1e-13)

    def test_archive_and_rejections(self):
        with tempfile.TemporaryDirectory() as d:
            q=request();q['output_archive']=str(Path(d)/'histories.npz')
            result=analytic_histories(q,executable=CLI)
            with np.load(q['output_archive']) as z:
                self.assertEqual(z['values'].shape,(3,3,3));self.assertEqual(z['values'].dtype,np.dtype('complex128'))
                metadata=json.loads(z['metadata_json'].item());self.assertEqual(metadata['schema'],'dynsyslab.analytic-histories/v1')
                self.assertEqual(len(metadata['solution_families']),3)
            with self.assertRaises(RuntimeError):analytic_histories(q,executable=CLI)
        q=request();q['initial_states']=[[1]]
        with self.assertRaises(RuntimeError):analytic_histories(q,executable=CLI)
        q=request();q['system']['parameters']['damping_n_m_s']=.1
        with self.assertRaises(RuntimeError):analytic_histories(q,executable=CLI)
        q=request();q['observables']=[]
        with self.assertRaises(RuntimeError):analytic_histories(q,executable=CLI)

if __name__=='__main__':unittest.main()
