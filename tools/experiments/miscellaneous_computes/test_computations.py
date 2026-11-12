"""Numerical regressions for reusable kernels, derivative fits, charts and diagnostics."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.linalg import expm
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from miscellaneous_computes.integration import exact_affine,polynomial_rk4,affine_equations
from miscellaneous_computes.models import SparseField,fit_field
from miscellaneous_computes.modal import ModalChart,ModalModel
from miscellaneous_computes.compact_field import CompactField
from miscellaneous_computes.structured import PairNetwork,SecondOrderField,InteriorField,EndpointAcceleration
from miscellaneous_computes.calibration import calibrate
from miscellaneous_computes.sampling import weak_rows,delay_samples,trapezoid_weights,sample_indices
from miscellaneous_computes.waveform_fit import filtered_library,fit_damped_drive
from miscellaneous_computes.reconstruction import compare_reconstructions
from miscellaneous_computes.crossings import candidates,refine
from miscellaneous_computes.diagnostics import forecast,kinematics
from miscellaneous_computes.delay_dmd import DiscreteDelayModel
from miscellaneous_computes.replay import replay
from miscellaneous_computes.comparison import comparison_pages
from miscellaneous_computes.equations import modal_equations
DLL=Path(sys.argv.pop(1)).resolve()
def linear(A,bias):
    n=np.shape(A)[1]
    return SparseField(dict(n=n,family='polynomial',degree=1,powers=[[i] for i in range(n)]),
                       np.column_stack([bias,A]))
class Tests(unittest.TestCase):
    def test_nested_physical_window_projection(self):
        from miscellaneous_computes.history import fit_basis
        rng = np.random.default_rng(13)
        t = np.arange(50) * .1
        x = rng.normal(size=(len(t), 2))
        basis = fit_basis([(t, x)], 5, 6, DLL, scales=[.5, 3.])
        windows = np.stack([x[j:j+5] for j in range(len(x)-4)])
        _, expected, _ = basis.project(t, x, rank=4)
        np.testing.assert_allclose(basis.project_windows(windows, 4), expected, atol=1e-12)
        np.testing.assert_allclose(basis.project_windows(windows[0])[:4], expected[0], atol=1e-12)
        with self.assertRaises(ValueError):
            basis.project_windows(windows[:, :-1], 4)
        with self.assertRaises(ValueError):
            basis.window_encoder(7)

    def test_waveform_objective_and_resolved_support(self):
        from miscellaneous_computes.waveform_objective import waveform_weights,resolved_directions
        t=np.array([0.,.1,.5,1.2,2.]);y=np.column_stack([np.ones(5)*4,np.ones(5)*2])
        w=waveform_weights(t,y,[.6,.2],[.2,0.],.4,.1)
        self.assertAlmostEqual(np.sum((w*y)**2),1.)
        x=np.array([[1.,2.,0.],[2.,4.,0.],[3.,6.,0.]])
        basis=np.eye(3)
        dirs,response=resolved_directions(x,basis)
        self.assertEqual(dirs.shape,(3,1))
        np.testing.assert_allclose(x@dirs,response)
        self.assertEqual(dirs[2,0],0.)
    def test_autonomous_delay_rollout(self):
        from miscellaneous_computes.delay_rollout import delay_rollout
        # x'=x(t-L), x(t)=1+t before zero. The method-of-steps solution is cubic on [L,2L].
        step=.004;lag=.13;past=np.arange(-34,1)*step
        run=delay_rollout(DLL,lambda window:window[1],(1+past)[:,None],[0.,lag],step,66)
        t=run['time'];extra=np.maximum(t-lag,0.)
        expected=1.+(1.-lag)*t+t*t/2.-lag*extra*extra/2.+extra**3/6.
        self.assertTrue(run['complete'])
        np.testing.assert_allclose(run['values'][:,0],expected,atol=3e-7,rtol=0.)
        with self.assertRaises(ValueError):
            delay_rollout(DLL,lambda window:window[1],(1+past[-5:])[:,None],[0.,lag],step,66)
        limited=delay_rollout(DLL,lambda window:np.ones(1),(1+past)[:,None],[0.,lag],step,66,limits=[1.05])
        self.assertFalse(limited['complete'])
        self.assertGreater(len(limited['time']),1)
        self.assertLess(limited['values'][-1,0],1.05)
    def test_signed_equations_and_partial_reference(self):
        from miscellaneous_computes.field_codec import decode
        from miscellaneous_computes.modal_study import evaluate_candidate
        chart=ModalChart(np.eye(2),np.eye(2),[0],[1,1],3.)
        quadratic=SparseField(dict(n=2,family='polynomial',powers=[[0],[1],[0,1]]),
                              [[1.,2.,0.,-.5],[0.,0.,-3.,.75]])
        m=ModalModel(chart,quadratic)
        point=np.array([.4,-.7])
        exported=modal_equations(m)
        value=[sum(term['coefficient']*np.prod(point**term['exponents']) for term in row)
               for row in exported['rows']]
        np.testing.assert_allclose(value,m.rhs(point))
        struct=ModalModel(chart,SecondOrderField(linear([[-2.,.1]],[.3])),'decoded')
        eq=modal_equations(struct)
        field=decode(eq['field'],{k:np.asarray(v) for k,v in eq['arrays'].items()})
        np.testing.assert_allclose(field(point@eq['decoder'])@eq['encoder'],struct.rhs(point))
        t=np.linspace(0,2,21);source=np.column_stack([t,np.ones_like(t)])
        ds=np.column_stack([np.ones_like(t),np.zeros_like(t)])
        history=source[:,None,:];rec=chart.record(t,history,ds[:,None,:])
        rec.update(truth=source,truth_rate=ds,taps=history)
        moving=ModalModel(chart,linear(np.zeros((2,2)),[chart.gain,0.]))
        score,curves=evaluate_candidate(moving,rec,0,1,dict(method='rk4',library=DLL,max_step=.05,max_steps=8))
        self.assertFalse(score['complete'])
        self.assertLess(len(curves['time']),len(curves['reference_time']))
        np.testing.assert_array_equal(curves['reference_time'],t)
        np.testing.assert_array_equal(curves['reference'],source)
        np.testing.assert_array_equal(curves['actual'],source[:len(curves['time'])])
    def test_nonuniform_weak_rows_and_weighted_fit(self):
        t=np.array([0,.07,.3,.42,.8,1.,1.4]);duration=.36
        avg,delta=weak_rows(t,duration,[0,2,4])
        start=t[[0,2,4]]
        np.testing.assert_allclose(avg@np.ones(len(t)),1.)
        np.testing.assert_allclose(avg@t,start+duration/2)
        np.testing.assert_allclose(delta@t,1.)
        self.assertAlmostEqual(trapezoid_weights(t).sum(),1.4)
        np.testing.assert_array_equal(sample_indices(np.array([[0],[2],[5],[1],[0]]),4,[(3,1)]),[0,2,4])
        rec=dict(time=t,b=t[:,None],db=np.ones((len(t),1)))
        f=fit_field([rec],degree=1,weak_durations=[duration],time_weighted=True,threshold=.02,ridge=1e-8)
        np.testing.assert_allclose(f([[.1],[.8]]),1.,atol=1e-6)
        with self.assertRaises(ValueError):weak_rows(t,1.1,[3])
        from miscellaneous_computes.regression import weighted_moments
        design=lambda x:np.column_stack([np.ones(len(x)),x,x*x])
        rows=np.flatnonzero(t+duration<=t[-1])
        W,V=weak_rows(t,duration,rows);aug=np.column_stack([W@design(t[:,None]),V@t])
        np.testing.assert_allclose(weighted_moments([rec],design,weak_durations=[duration],block_size=2),
                                   aug.T@aug,atol=1e-12)
        with self.assertRaises(ValueError):delay_samples(t,t,[.4],[0.])
    def test_native_rollout_and_complete_prefix(self):
        f=linear(np.array([[0.,1.],[-3.,-.2]]),[.1,.2]);A,c=affine_equations(f)
        t=np.r_[0.,.013,.21,.8,1.1];x0=[.4,.6]
        ck=exact_affine(A,c,t,x0)
        r=polynomial_rk4(f,t,x0,DLL,max_step=.005)
        self.assertTrue(r['complete']);np.testing.assert_allclose(ck,r['b'],atol=1e-9)
        regular=np.linspace(0,2,200)
        np.testing.assert_allclose(exact_affine(A,c,regular,x0)[-1],exact_affine(A,c,[0,2],x0)[-1],atol=1e-12)
        r=polynomial_rk4(f,t,x0,DLL,max_step=.01,max_steps=3)
        self.assertEqual(r['message'],'step budget');np.testing.assert_array_equal(r['time'],t[:2])
        r=polynomial_rk4(f,t,x0,DLL,max_step=.01,limits=[.41,2])
        self.assertFalse(r['complete']);self.assertLess(len(r['b']),len(t))
    def test_bundle_fields_and_kinematics(self):
        chart=ModalChart(np.eye(2),np.eye(2),[0],[2.,2.],7.)
        f=linear([[0,1],[-1,0]],[0,0]);m=ModalModel(chart,f)
        window=[[1.,0.]];v=m.chart.project(window)
        np.testing.assert_allclose(v,[2,0])
        q=replay(m,np.linspace(0,2.,69),window=window,method='exact-affine')
        np.testing.assert_allclose(q['values'][:,0,0],np.cos(q['time']),atol=1e-13)
        self.assertLess(kinematics(m,q['xi'],value_channel=0,rate_channel=1)['relative_mismatch'],1e-12)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'model.npz';m.save(p);m1=ModalModel.load(p)
            np.testing.assert_allclose(m1.rhs(v),m.rhs(v))
            with self.assertRaises(FileExistsError):m.save(p)
        with self.assertRaises(ValueError):ModalChart(np.eye(2),2*np.eye(2),[0],[1,1])
    def test_compact_field_multioutput_and_readout(self):
        t=np.linspace(0,5,340);x=np.column_stack([np.sin(t),np.cos(t)])
        dy=np.column_stack([x[:,1],-x[:,0]])
        f=CompactField.fit([dict(time=t,b=x,db=dy)],stride=7,neighbors=7,affine_tail=True,
                           time_weighted=True,speed_scale=4.,asinh_scale=1.,ridge=.0005,threshold=.00001)
        self.assertLess(np.linalg.norm(f(x)-dy)/np.linalg.norm(dy),.015)
        g=CompactField.from_labels([(t,x,dy[:,:1])],stride=7,neighbors=5)
        self.assertEqual(g(x).shape,(len(x),1))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'c.npz';f.save(p);np.testing.assert_allclose(CompactField.load(p)(x),f(x))
    def test_structured_maps_and_filtered_forcing(self):
        t=np.linspace(0,10,1501);lags=np.array([0.,.5,1.])
        ts=t[:,None]-lags;b=np.stack([np.sin(ts),np.cos(ts)],axis=-1).reshape(len(t),6)
        v=np.stack([np.cos(ts),-np.sin(ts)],axis=-1).reshape(len(t),6);r=[dict(time=t,b=b,db=v)]
        weights=(np.ones((3,3))-np.eye(3))/2
        pair=PairNetwork.fit(r,lags,weights,degree=1,threshold=0.,ridge=1e-8)
        self.assertEqual(len(pair.fields),2)
        interior=InteriorField.fit(r,3.,drive_options=dict(degree=1,threshold=0.,ridge=1e-8))
        simple=SecondOrderField.fit(r,degree=1,ridge=1e-8,threshold=0.)
        with tempfile.TemporaryDirectory() as tmp:
            for i,field in enumerate([pair,interior,simple]):
                y=field(b[::31]);np.testing.assert_allclose(y,v[::31],atol=2e-6)
                model=ModalModel(ModalChart(np.eye(6),np.eye(6),lags,[1.,1.]),field,'decoded')
                p=Path(tmp)/f'{i}.npz';model.save(p);np.testing.assert_allclose(ModalModel.load(p).rhs(b[::31]),y)
            base=ModalModel(ModalChart(np.eye(4),np.eye(4),[0,1],[1,1]),
                SecondOrderField(linear([[-1,0,0,0],[0,0,-1,0]],[0,0])),'decoded')
            hybrid=InteriorField.fit(r,3.,outer=EndpointAcceleration(base),feedback_options=dict(degree=1,ridge=1e-8))
            z=ModalModel(model.chart,hybrid,'decoded');p=Path(tmp)/'inherited.npz';z.save(p)
            np.testing.assert_allclose(ModalModel.load(p).rhs(b[::31]),z.rhs(b[::31]))
        gamma=2.;P,V=filtered_library(t,np.ones((len(t),1)),gamma)
        np.testing.assert_allclose(P[:,0],1-(1+gamma*t)*np.exp(-gamma*t),atol=1e-10)
        np.testing.assert_allclose(V[:,0],gamma**2*t*np.exp(-gamma*t),atol=1e-10)
        force=fit_damped_drive([(t,b[:,:2],P,V)],gamma,family='polynomial',degree=1,ridge=1e-9,threshold=.01)
        np.testing.assert_allclose(force(b[:10,:2]),gamma**2,atol=1e-6)
    def test_crossing_witness_and_partial_peaks(self):
        plane=lambda t:np.array([np.sin(t),np.sin(2*t)])
        rate=lambda t:np.array([np.cos(t),2*np.cos(2*t)])
        t=np.linspace(-.4,3.5,300)
        h=candidates(t,plane(t).T,DLL,min_separation=1.)
        self.assertEqual(len(h),1)
        ref=refine(h[0],plane,velocity=rate,lift=np.cos)
        self.assertTrue(ref['matched']);np.testing.assert_allclose(ref['times'],[0,np.pi],atol=1e-9)
        self.assertGreater(ref['acute_angle_degrees'],40)
        self.assertAlmostEqual(ref['lift_difference'][0],-2.)
        t=np.arange(50)*.1;x=np.sin(4*t)[:,None];p=forecast(t,x,t[:32],x[:32],prominence=.4)
        self.assertFalse(p['complete']);self.assertEqual(p['errors']['relative'],[0.])
        self.assertTrue(any(not m['full_interval'] for m in p['peaks']['matches']))
        empty=forecast(t,x,[],np.empty((0,1)));self.assertIsNone(empty['errors'])
    def test_discrete_advance_and_receiver_delays(self):
        t=np.arange(500)*.02
        def record(a,shift):
            p=a*np.cos(t+shift);v=-a*np.sin(t+shift);return t,p[:,None],np.column_stack([p,v])
        train={'a':record(1.,0.),'b':record(1.4,.4)};test={'test':record(1.2,.1)}
        designs={'two samples':[dict(channel=0,delay=0),dict(channel=0,delay=.4)]}
        r=compare_reconstructions(train,test,designs,neighbors=6,train_stride=2,steps_ahead=True)['two samples']
        self.assertAlmostEqual(r['spectral_radius'],1.)
        self.assertLess(max(r['cases']['test']['scores']['linear']['relative']),1e-10)
        self.assertLess(max(r['cases']['test']['scores']['autonomous']['errors']['relative']),.002)
        with self.assertRaises(ValueError):DiscreteDelayModel.fit([(t,record(1,0)[2]),(2*t,record(1,.5)[2])])
    def test_calibrate_and_plot_contracts(self):
        t=np.linspace(0,3,40);y=np.exp(-t)[:,None];r=[dict(time=t,b=y,db=-y)]
        f=linear([[-.8]],[0]);g,audit=calibrate(f,r,duration=1.,max_evaluations=10)
        self.assertTrue(audit['accepted']);self.assertLess(abs(g.coef[0,1]+1),.002)
        self.assertEqual(f.coef[0,1],-.8)
        x=np.column_stack([y,-y]);pages=comparison_pages({'test':{'reference':(t,x),'old':(t,x),'fresh':(t[:20],x[:20])}},
                     ['fresh'],baseline='old')
        curves=pages['test'][0]['curves'];self.assertEqual([c[0] for c in curves[:2]],['old','reference'])
        self.assertIn('to',curves[2][0]);self.assertEqual(len(pages),2)

if __name__=='__main__':unittest.main()
