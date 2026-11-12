"""Independent oscillator, delayed multichannel covariance, frozen-basis and free-forecast regressions."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from miscellaneous_computes.history import fit_basis,HistoryBasis,separation,scalar_delay_pairs
from miscellaneous_computes.tapped_history import sample_delay_jet,full_tap_svd,delay_consistency
from miscellaneous_computes.models import fit_field,free_rollout,SparseField
from miscellaneous_computes.sparse_stlsq import sparse_stlsq
from miscellaneous_computes.local_field import LocalField
from miscellaneous_computes.paired import PairedField
from miscellaneous_computes.readouts import NeighborReadout,TangentReadout
from miscellaneous_computes.prehistory import backward_history
from miscellaneous_computes.audit import replay_controls
from miscellaneous_computes.study import rank_sweep
exe=Path(sys.argv.pop(1)).resolve();dll=Path(sys.argv.pop(1)).resolve()

def oscillator(amplitude):
    t=np.arange(0,8,.015);y=amplitude*np.cos(2*t+.23);v=-2*amplitude*np.sin(2*t+.23)
    return t,y,v
class Regression(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.d=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def test_sparse_modal_derivative_fit(self):
        from scipy.sparse import csr_matrix
        rng=np.random.default_rng(126)
        theta=rng.normal(size=(500,16));theta[:,-1]=0.
        target=np.column_stack([3*theta[:,4],-2*theta[:,1],np.zeros(500)])
        coef=sparse_stlsq(csr_matrix(theta.T@theta),theta.T@target,np.linalg.norm(target,axis=0),
                          ridge=.0001,threshold=.01)
        self.assertEqual(np.count_nonzero(coef),2)
        np.testing.assert_allclose(coef[[4,1],[0,1]],[3,-2],atol=.0005)
    def test_fractional_delay_window(self):
        t=np.arange(150)*.01
        jet=lambda t:np.stack([t**5-2*t,5*t**4-2,20*t**3],axis=-1)
        lags=np.array([0.,.1234,.251])
        x=sample_delay_jet(t,jet(t),t[40:],lags)
        np.testing.assert_allclose(x,jet(t[40:,None]-lags),atol=3e-9)
        with self.assertRaises(ValueError):sample_delay_jet(t,jet(t),[.403],[.001])
        rng=np.random.default_rng(52);h=rng.normal(size=(230,6,2))
        s=full_tap_svd([h[:110],h[110:]],[.7,9.])
        v=h.reshape(len(h),-1)@s['encoder']
        np.testing.assert_allclose(v.T@v,np.eye(12),atol=1e-12)
        np.testing.assert_allclose(v@s['decoder'],h.reshape(len(h),-1),atol=1e-12)
        clock=t[40:]
        pack=lambda a:np.stack([a**3,3*a*a],axis=-1)
        pairs=pack(clock[:,None]-lags).reshape(len(clock),-1)
        np.testing.assert_allclose(delay_consistency(clock,pairs,pairs,lags),0,atol=1e-12)
        altered=pairs.copy();altered[:,2]+=.1
        self.assertGreater(delay_consistency(clock,altered,pairs,lags)[0][0],.01)
    def test_causal_scalar_pairs(self):
        t=np.arange(100)*.02;signal=t*t+3*t
        ends=np.array([50,80,98]);lags=np.array([0,12]);half=3
        raw,x,a=scalar_delay_pairs(signal,ends,lags,.02,half)
        newest=ends[:,None]-lags
        midpoint=t[newest-half]
        np.testing.assert_array_equal(raw[:,1::2],signal[newest])
        np.testing.assert_array_equal(raw[:,::2],signal[newest-2*half])
        np.testing.assert_allclose(x[:,::2],midpoint**2+3*midpoint+.06**2)
        np.testing.assert_allclose(x[:,1::2],2*midpoint+3)
        np.testing.assert_allclose(a,2.)
        with self.assertRaises(ValueError):scalar_delay_pairs(signal,[17],lags,.02,half)
    def test_covariance_projection_and_locked_modes(self):
        rng=np.random.default_rng(4)
        series=[rng.normal(size=(n,2)) for n in (61,39)]
        records=[(np.arange(len(x))*.1,x) for x in series]
        gap=3;L=4
        H=np.vstack([np.array([x[i:i+gap*L:gap].ravel() for i in range(len(x)-(L-1)*gap)]) for x in series])
        scales=[1.,1.]
        basis=fit_basis(records,L,8,dll,gap=gap,scales=scales)
        np.testing.assert_allclose(basis.sigma,np.linalg.svd(H,compute_uv=False),rtol=5e-12)
        fixed=(basis.U[:,:3],basis.sigma[:3])
        extended=fit_basis(records,L,6,dll,gap=gap,scales=scales,locked=fixed)
        np.testing.assert_array_equal(extended.U[:,:3],fixed[0])
        self.assertLess(extended.gram_error,1e-10)
        t,b,_=basis.project(*records[0]);N=len(b)
        np.testing.assert_allclose(t,records[0][0][gap*(L-1):])
        np.testing.assert_allclose(b,H[:N]@basis.U*np.sqrt(len(H))/basis.sigma,rtol=1e-12,atol=2e-13)
        np.testing.assert_allclose(b@basis.endpoint(8,1),series[0][-N:,1],atol=1e-12)
        with self.assertRaises(ValueError):fit_basis(records,L,6,dll,gap=gap,scales=[1,2],locked=fixed)
        p=self.d/'basis.npz';extended.save(p);np.testing.assert_array_equal(HistoryBasis.load(p).U,extended.U)
        with self.assertRaises(ValueError):basis.project(t[:2],[1,2])
    def test_free_forecast_and_control(self):
        data=[oscillator(a) for a in (1.,1.3,1.6)]
        basis=fit_basis([(data[0][0],data[0][1])],24,2,dll)
        rec=[]
        for t,y,v in data:
            tm,b,db=basis.project(t,y,v);rec.append(dict(time=tm,b=b,db=db,truth=np.column_stack([y[-len(tm):],v[-len(tm):]])))
        entry,field,_=next(rank_sweep(basis,[2],{'train':rec[0]},{'valid':rec[1]},{'test':rec[2]},
                              [dict(family='polynomial',degree=1,ridge=1e-10,threshold=.0001)],
                              solver=dict(max_step=.1)))
        self.assertLess(max(entry['results']['test']['test']['free']),1e-6)
        checks=replay_controls(field,rec[2],2,basis.endpoint(2),restart_time=1.2,compare_from=2.)
        self.assertLess(max(checks['restart_on_common']),1e-6)
        model_path=self.d/'weak.json'
        weak=fit_field([rec[0]],family='polynomial',degree=1,weak_window=8);weak.save(model_path)
        self.assertLess(np.linalg.norm(SparseField.load(model_path)(rec[2]['b'])-rec[2]['db'])/np.linalg.norm(rec[2]['db']),.001)
        for family in ['bounded','localized','harmonic','gated','paired-gated']:
            fitted=fit_field([rec[0]],family=family,centers=6)
            self.assertTrue(np.isfinite(fitted(rec[2]['b'][:10])).all())
        paired=PairedField.fit(rec,harmonics=[-1],weak_windows=[8],instant_weight=.8)
        self.assertTrue(np.isfinite(paired(rec[2]['b'])).all())
        d,dist=separation(rec[0]['time'],rec[0]['b'],1.5,1.5)
        np.testing.assert_array_equal(dist,0.)
        with self.assertRaises(ValueError):list(rank_sweep(basis,[2],{'x':rec[0]},{'x':rec[1]},{'z':rec[2]},[{}]))
    def test_local_and_readouts(self):
        t=np.linspace(0,7,501);records=[]
        for a in (1.,1.3,1.6):
            u=a*np.column_stack([np.sin(t),np.cos(t)]);v=u[:,[1,0]]*[1,-1]
            records.append(dict(time=t,b=u,db=v))
        data=records[1];b=data['b'];v=data['db']
        local=LocalField.fit(records,stride=9,samples=35,neighbors=6,library=dll)
        p=self.d/'local.npz';local.save(p)
        np.testing.assert_allclose(local(b),LocalField.load(p)(b),rtol=1e-9,atol=1e-10)
        np.savez(self.d/'query.npz',b=b)
        subprocess.run([exe,'local-eval',p,self.d/'query.npz',self.d/'rates.npz'],check=True)
        with np.load(self.d/'rates.npz') as z:np.testing.assert_allclose(z['db'],local(b),rtol=1e-10,atol=1e-10)
        self.assertLess(np.linalg.norm(local(b)-v)/np.linalg.norm(v),.04)
        readout=NeighborReadout(b,b);np.testing.assert_allclose(readout.predict(b,1),b)
        f=TangentReadout(b,v,b[:,0],v[:,0]);q=b[20:31]+.002;dq=v[20:31];eps=1e-7
        numerical=(f.phase(q+eps*dq,dq,8)[:,0]-f.phase(q-eps*dq,dq,8)[:,0])/(2*eps)
        np.testing.assert_allclose(f.phase(q,dq,8)[:,1],numerical,rtol=2e-5,atol=1e-7)
        ts=np.linspace(-.75,-.01,70)
        obs,err=backward_history(lambda _,x:[x[1],-x[0]],lambda _,x:x[0],ts,0.,[1.,0.])
        np.testing.assert_allclose(obs,np.cos(ts),atol=1e-9);self.assertLess(err,1e-9)
    def test_native_cli_and_manifest(self):
        t,y,v=oscillator(1.)
        b=np.column_stack([y,v]);db=np.column_stack([v,-4*y])
        src=self.d/'train.npz';np.savez(src,time=t,b=b,db=db,values=y[:,None],rates=v[:,None])
        dest=self.d/'gram.npz';subprocess.run([exe,'gram',src,'5',dest],check=True)
        h=np.lib.stride_tricks.sliding_window_view(y,5)
        with np.load(dest) as z:np.testing.assert_allclose(z['gram'],h.T@h,atol=1e-9)
        model=self.d/'native.npz';subprocess.run([exe,'fit',src,'1','.00001',model],check=True)
        query=self.d/'query.npz';np.savez(query,time=t,initial=b[0])
        out=self.d/'forecast.npz';subprocess.run([exe,'rollout',model,query,out],check=True)
        with np.load(out) as z:np.testing.assert_allclose(z['b'],b,atol=3e-7)
        for a,k in [(1.3,'valid'),(1.6,'test')]:
            t,y,v=oscillator(a);np.savez(self.d/(k+'.npz'),time=t,values=y,rates=v)
        plan=dict(schema='dynsyslab.miscellaneous-computes/study/v1',basis=dict(taps=12),ranks=[2],
                  training={'a':'train.npz'},validation={'b':'valid.npz'},test={'c':'test.npz'},
                  candidates=[dict(family='polynomial',degree=1,ridge=1e-10,threshold=0.)])
        manifest=self.d/'study.json';manifest.write_text(json.dumps(plan))
        plan['controls']=dict(neighbor_count=8,tangent_neighbors=[8,16],restart_time={'c':2.},visits={'a':[1.,4.]})
        manifest.write_text(json.dumps(plan))
        subprocess.run([sys.executable,Path(__file__).with_name('run_study.py'),manifest,self.d/'results','--library',dll],check=True)
        result=json.loads((self.d/'results/rank-2/scores.json').read_text())
        self.assertLess(max(result['results']['test']['c']['free']),1e-6)
        self.assertIn('tangent_free',result['controls']['readouts']['c'])
        self.assertEqual(len(result['controls']['visits']['a']['distance']),2)
if __name__=='__main__':unittest.main()
