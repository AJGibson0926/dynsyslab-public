"""End-to-end checks for independent publication, native tools, deploy/replay and comparisons."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
side=Path(__file__).resolve().parent
sys.path.insert(0,str(side.parent))
from miscellaneous_computes.modal import ModalModel
from miscellaneous_computes.modal_study import run_manifest
from miscellaneous_computes.witness_checks import match_visits,difference_weights
from miscellaneous_computes.comparison import comparison_pages
exe,dll=[Path(sys.argv.pop(1)).resolve() for _ in range(2)]
closure_cli=Path(sys.argv.pop(1)) if len(sys.argv)>1 else None
class CLI(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.d=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def call(self,*args):
        subprocess.run([sys.executable,side/'toolkit.py',*map(str,args)],check=True,capture_output=True)
    def records(self):
        t=np.arange(0,10,.025)
        for i,tag in enumerate(['train','validate','test']):
            u=t+.31*i;pv=[np.sin(u)+.22*np.sin(2*u)+.15*np.cos(3*u),
                np.cos(u)+.44*np.cos(2*u)-.45*np.sin(3*u),-np.sin(u)-.88*np.sin(2*u)-1.35*np.cos(3*u)]
            np.savez(self.d/(tag+'.npz'),time=t,values=np.stack(pv[:2],axis=1),rates=np.stack(pv[1:],axis=1),
                     state=np.column_stack([np.sin(u),np.cos(u)]))
        return dict(training={'A':'train.npz'},validation={'B':'validate.npz'},test={'C':'test.npz'})
    def test_study_and_replay(self):
        cfg=dict(schema='dynsyslab.miscellaneous-computes/modal-study/v1',**self.records(),
                 geometries=[dict(id='three-taps',lags=[0,.57,1.103])],interpolation='quintic_pair',
                 candidates=[dict(family='polynomial',degree=1,time_weighted=True,ridge=1e-10,threshold=1e-10)],
                 rate_channel=1,solver=dict(method='exact-affine'),peak_prominence=.1,fit_stride=2)
        cfg['observation_provenance']={'A':{'rates':'analytic harmonic time derivatives'},'C':{'rates':'held-out'}}
        manifest=self.d/'spec.json';manifest.write_text(json.dumps(cfg))
        self.call('study',manifest,self.d/'first')
        folder=self.d/'first'/'three-taps';model=ModalModel.load(folder/'model.npz')
        result=json.loads((self.d/'first/results.json').read_text())['geometries']['three-taps']
        self.assertEqual(model.metadata['origin'],'fresh');self.assertEqual(model.metadata['training'],['A'])
        self.assertEqual(model.metadata['observation_provenance'],{'A':cfg['observation_provenance']['A']})
        self.assertEqual(set(model.metadata['source_files']),{'A','B'})
        self.assertLess(max(result['results']['C']['output']['errors']['relative']),1e-6)
        with np.load(folder/'curves/C.npz') as z:
            values=z['predicted'];t=z['time'];initial=z['xi'][0]
        np.savez(self.d/'query.npz',time=t,initial=initial)
        self.call('replay',folder/'model.npz',self.d/'query.npz',self.d/'free.npz','--method','rk4','--library',dll,'--step','.005')
        with np.load(self.d/'free.npz') as z:
            np.testing.assert_allclose(z['values'],values,atol=1e-7)
            np.testing.assert_array_equal(z['values'],z['windows'][:,0])
        self.call('equations',folder/'model.npz',self.d/'eq.json');self.assertEqual(len(json.loads((self.d/'eq.json').read_text())['A']),6)
        self.call('score',folder/'curves/C.npz',self.d/'free.npz',self.d/'scores.json',
                  '--reference-key','reference','--reference-time','reference_time',
                  '--reference-columns','0','--prediction-columns','0','--prominence','.1')
        self.assertTrue(json.loads((self.d/'scores.json').read_text())['complete'])
        with self.assertRaises(FileExistsError):run_manifest(manifest,self.d/'first')
        cfg['test']={'C':'train.npz'};manifest.write_text(json.dumps(cfg))
        with self.assertRaises(ValueError):run_manifest(manifest,self.d/'failed')
        self.assertFalse((self.d/'failed').exists())
    def test_structured_study_explicit_origin(self):
        cfg=dict(schema='dynsyslab.miscellaneous-computes/modal-study/v1',**self.records(),
                 geometries=[dict(id='endpoints',lags=[0,1.1])],interpolation='quintic_pair',experiment='structured',
                 candidates=[dict(family='second-order',degree=1,ridge=.001)],rate_channel=1)
        manifest=self.d/'spec.json';manifest.write_text(json.dumps(cfg));run_manifest(manifest,self.d/'end')
        cfg['geometries']=[dict(id='interior',lags=[0,.55,1.1])]
        cfg['candidates']=[dict(family='damped-interior',gamma=4.,
            inherited_model='end/endpoints/model.npz',drive_options=dict(degree=1),feedback_options=dict(degree=1))]
        manifest.write_text(json.dumps(cfg))
        result=run_manifest(manifest,self.d/'structured')
        model=ModalModel.load(self.d/'structured/interior/model.npz')
        self.assertEqual(model.metadata['origin'],'inherited');self.assertTrue(result['interior']['results']['C']['complete'])
        cfg['experiment']='direct';manifest.write_text(json.dumps(cfg))
        with self.assertRaises(ValueError):run_manifest(manifest,self.d/'illegal')
    def test_native_commands_and_crossing_match(self):
        t=np.linspace(0,4,500);x=np.column_stack([np.sin(t),np.cos(t)]);db=x[:,[1,0]]*[1,-1]
        src=self.d/'coords.npz';model=self.d/'model.npz';np.savez(src,time=t,b=x,db=db)
        subprocess.run([exe,'fit',src,'1','0.0001',model],check=True)
        q=self.d/'query.npz';np.savez(q,time=t[:100],initial=x[0]);p=self.d/'out.npz'
        subprocess.run([exe,'rk4',model,q,'.002',p],check=True)
        with np.load(p) as z:np.testing.assert_allclose(z['b'],x[:100],atol=1e-8)
        t=np.linspace(-.4,3.5,300);xy=np.column_stack([np.sin(t),np.sin(2*t)])
        plane=self.d/'plane.npz';np.savez(plane,time=t,values=xy)
        subprocess.run([exe,'crossings',plane,self.d/'cross.npz'],check=True)
        scalar=self.d/'scalar.npz';np.savez(scalar,time=t,values=t)
        subprocess.run([exe,'integral',scalar,'.245',self.d/'integrated.npz'],check=True)
        with np.load(self.d/'integrated.npz') as z:
            np.testing.assert_allclose(z['values'][:,1],.245*(z['time']-.245/2),atol=1e-10)
        with np.load(self.d/'cross.npz') as z:self.assertEqual(len(z['times']),1)
        self.call('crossings',plane,self.d/'hits.json','--library',dll,'--refine-with',plane)
        h=json.loads((self.d/'hits.json').read_text())['candidates'];self.assertTrue(h[0]['refinement']['matched'])
        paired=match_visits(h,[dict(times=[0,np.pi]),dict(times=[.2,3.3])],.001)
        self.assertEqual(paired['unmatched_second'],[1]);self.assertEqual(match_visits([],h,1.)['unmatched_second'],[0])
        np.testing.assert_allclose(difference_weights([-1,0,1],2),[1,-2,1])
        self.assertEqual(difference_weights([-1,0,1],2,exact=True),[1,-2,1])
    def test_receiver_cli(self):
        splits=self.records()
        cfg=dict(schema='dynsyslab.miscellaneous-computes/reconstruction/v1',training=splits['training'],test=splits['test'],
                 designs={'two delays':[dict(channel=0,delay=0),dict(channel=0,delay=.41)]},options=dict(train_stride=2,steps_ahead=True))
        p=self.d/'r.json';p.write_text(json.dumps(cfg));self.call('reconstruct',p,self.d/'recon')
        with np.load(self.d/'recon/two delays/C.npz') as z:
            self.assertEqual(z['actual'].shape,z['neighbor'].shape);self.assertEqual(z['forecast'].shape[1],2)
    @unittest.skipUnless(closure_cli,'Native suite CLI unavailable')
    def test_exact_closure(self):
        query=side.parents[2]/'examples/native/duffing-position-global-jet.json'
        p=self.d/'proof.json'
        self.call('closure',query,p,'--cli',closure_cli)
        proof=json.loads(p.read_text())
        self.assertEqual(proof['certified_order'],1)
        self.assertTrue(proof['minimality_proved'])
if __name__=='__main__':unittest.main()
