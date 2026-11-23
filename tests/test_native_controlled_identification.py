"""Native controlled-window identification: physical and chart comparators."""
import os
from pathlib import Path
import unittest
import numpy as np
from dynsyslab_client.periodic import identify_controlled, identify_controlled_sequential, replay_controlled, simulate_feedback
from dynsyslab_client.observables import evaluate_batch

ROOT=Path(__file__).resolve().parents[1]
CLI=Path(os.environ.get('DYNSYSLAB_NATIVE_BIN',ROOT/'build/native-suite-fast/bin'))/('dynsyslab-cli.exe' if os.name=='nt' else 'dynsyslab-cli')


@unittest.skipUnless(CLI.is_file(),'Build native CLI first')
class ControlledIdentification(unittest.TestCase):
    def test_native_batch_energy_readout(self):
        state=[[1.,0.],[2.,.3],[.4,-3.]]
        result=evaluate_batch(dict(system=dict(name='rpe',parameters=dict(mu_pa_s=0)),
            states=np.ravel(state).tolist(),observables=['radius','mechanical_energy','energy_dissipation','control_power_coefficient']),executable=CLI)
        values=np.asarray(result['values']).reshape(-1,4)
        self.assertEqual(result['samples'],3)
        np.testing.assert_array_equal(values[:,0],np.asarray(state)[:,0])
        self.assertEqual(values[0,1],0)
        self.assertTrue(np.all(values[1:,1]>0))
        np.testing.assert_array_equal(values[:,2],0)
        self.assertLess(values[1,3],0)
        self.assertGreater(values[2,3],0)
        with self.assertRaises(RuntimeError):evaluate_batch(dict(system=dict(name='rpe'),states=[1.,2.,3.],observables=['energy']),executable=CLI)

    def test_online_feedback_and_recorded_replay_agree(self):
        knots=[dict(parameter=float(a),channels=[
            dict(angular_frequency=1,coefficients=[[-1,a/2,0],[1,a/2,0]]),
            dict(angular_frequency=1,coefficients=[[-1,0,-a/2],[1,0,a/2]])]) for a in np.linspace(.7,2,33)]
        request=dict(system=dict(name='harmonic_oscillator'),knots=knots,observables=['position','velocity'],
            time=[0.,1.2],sample_time=.03,first_fit_index=20,fit_every=10,prediction_substeps=4,
            initial_parameters=[1.4,.7],initial_candidates=[[1.2,.3],[1.6,.9]],measurement_channels=[0],
            measurement_variances=[.001],measurement_feedthrough=[.3],table_samples=512,input_maximum_harmonic=63,
            relative_tolerance=1e-10,absolute_tolerance=1e-12,controller=dict(state_gain=[0,-.25],limit=.4))
        for model in ('physical','atlas'):
            result=simulate_feedback({**request,'model':model},executable=CLI)
            self.assertTrue(result['completed'],result['termination_reason'])
            self.assertEqual(result['estimate_ready'],[False]*20+[True]*21)
            self.assertEqual(result['information_time'][:20],[None]*20)
            truth=np.asarray(result['sample_truth']).reshape(-1,2)
            estimated=np.asarray(result['estimates']).reshape(-1,2)
            np.testing.assert_allclose(estimated[20:],truth[20:],atol=3e-6)
            np.testing.assert_allclose(result['measurements'],truth[:,0]+.3*np.array(result['previous_applied_input']),atol=1e-12)
            np.testing.assert_allclose(result['commands'][20:],np.clip(-.25*estimated[20:,1],-.4,.4),atol=1e-12)
            np.testing.assert_array_equal(result['commands'][:20],0)
            np.testing.assert_array_equal(result['previous_applied_input'][1:],result['commands'][:-1])
            self.assertTrue(np.all(np.asarray(result['dense_information_time'])<=result['dense_time']))
            common={k:v for k,v in request.items() if k not in ('sample_time','first_fit_index','fit_every','initial_parameters','controller')}
            replay=identify_controlled_sequential({**common,'model':model,'time':result['sample_time'],
                'measurements':result['measurements'],'update_indices':[20,30,40],
                'input':dict(time=result['sample_time'],values=result['commands'],before_first=0)},executable=CLI)
            np.testing.assert_allclose(replay['time'],result['dense_time'],atol=1e-12)
            np.testing.assert_allclose(replay['predictions'],result['dense_estimates'],atol=3e-6)
        noisy={**request,'measurement_noise':dict(kind='white',relative=False,level=.001,seed=31)}
        measured=simulate_feedback({**noisy,'model':'physical'},executable=CLI)
        baseline=simulate_feedback({**noisy,'model':'physical','feedback_source':'truth'},executable=CLI)
        np.testing.assert_allclose(measured['noise'],baseline['noise'],atol=1e-15)
        np.testing.assert_array_equal(measured['measurements'][:21],baseline['measurements'][:21])
        initialized=simulate_feedback({**request,'model':'physical','coarse_initialization':
            dict(parameters=np.linspace(.7,2,33).tolist(),phase_count=64)},executable=CLI)
        self.assertTrue(initialized['completed'],initialized['termination_reason'])
        np.testing.assert_allclose(initialized['updates'][0]['fit']['parameters'][0],1.4,atol=2e-6)
        with self.assertRaises(RuntimeError):simulate_feedback({**request,'input':dict(bias=.1)},executable=CLI)
        escaped=simulate_feedback({**request,'model':'atlas','time':[0.,3.],
            'controller':dict(state_gain=[-40.,0.],limit=5.)},executable=CLI)
        self.assertFalse(escaped['completed'])
        self.assertTrue(escaped['termination_reason'])
        self.assertGreater(len(escaped['time']),1)

    def test_recorded_switches_and_measurement_timing(self):
        knots=[dict(parameter=float(a),channels=[
            dict(angular_frequency=1,coefficients=[[-1,a/2,0],[1,a/2,0]]),
            dict(angular_frequency=1,coefficients=[[-1,0,-a/2],[1,0,a/2]])]) for a in np.linspace(.7,2,33)]
        switch=np.array([0.,.173,.1731,.5,1.]);inputs=np.array([.1,200.,-.3,.4,7.])
        radius,phase=1.4,.7;x0,v0=radius*np.cos(phase),-radius*np.sin(phase)
        def exact(time):
            out=[]
            for t in time:
                x,v=x0,v0
                for k in range(len(switch)-1):
                    dt=max(0,min(t,switch[k+1])-switch[k]);b=inputs[k]
                    x,v=b+(x-b)*np.cos(dt)+v*np.sin(dt),-(x-b)*np.sin(dt)+v*np.cos(dt)
                i=np.searchsorted(switch,t,side='left')-1
                u=inputs[i] if i>=0 else -.2
                out.append([x,v,x+.3*u])
            return np.array(out)
        t=np.linspace(0,.4,41);data=exact(t)
        request=dict(system=dict(name='harmonic_oscillator'),knots=knots,observables=['position','velocity'],
            time=t.tolist(),measurements=data[:,2].tolist(),measurement_channels=[0],measurement_variances=[.001],
            measurement_feedthrough=[.3],initial_candidates=[[1.2,.3],[1.6,.9]],
            input=dict(time=switch.tolist(),values=inputs.tolist(),before_first=-.2),
            table_samples=512,input_maximum_harmonic=63,relative_tolerance=1e-10,absolute_tolerance=1e-12,
            forecast_final_time=1.,output_step=.01)
        for model in ['physical','atlas']:
            r=identify_controlled({**request,'model':model},executable=CLI)
            np.testing.assert_allclose(r['fit']['parameters'],[radius,phase],atol=3e-6)
            expected=exact(r['time'])
            np.testing.assert_allclose(np.array(r['predictions']).reshape(-1,2),expected[:,:2],atol=3e-6)
            np.testing.assert_allclose(r['predicted_measurements'],expected[:,2],atol=3e-6)
        changed={**request['input'],'values':[.1,200.,-.3,-.7,7.]}
        future=identify_controlled({**request,'model':'atlas','input':changed},executable=CLI)
        self.assertEqual(r['fit'],future['fit'])
        self.assertGreater(np.max(np.abs(np.array(r['predictions'])-future['predictions'])),.1)
        replay={k:v for k,v in request.items() if k not in ('measurements','initial_candidates','forecast_final_time')}
        replay.update(time=[0,1],initial_parameters=[radius,phase],model='physical')
        forward=replay_controlled(replay,executable=CLI)
        np.testing.assert_allclose(np.array(forward['predictions']).reshape(-1,2),exact(forward['time'])[:,:2],atol=2e-8)
        self.assertEqual(forward['sensor_input'][0],-.2)
        self.assertEqual(forward['applied_input'][0],.1)
        with self.assertRaises(RuntimeError):replay_controlled({**replay,'time':[0,1.1]},executable=CLI)
        with self.assertRaises(RuntimeError):replay_controlled({**replay,'input':{**request['input'],'bias':.2}},executable=CLI)

    def test_exact_forced_oscillator_and_future_independence(self):
        knots=[]
        for a in np.linspace(.7,2,33):
            knots.append(dict(parameter=float(a),channels=[
                dict(angular_frequency=1,coefficients=[[-1,a/2,0],[1,a/2,0]]),
                dict(angular_frequency=1,coefficients=[[-1,0,-a/2],[1,0,a/2]])]))
        time=np.linspace(2,5,101);elapsed=time-2;radius,phase=1.4,.7;bias=.04
        x0,v0=radius*np.cos(phase),-radius*np.sin(phase)
        y=bias+(x0-bias)*np.cos(elapsed)+v0*np.sin(elapsed)+.3*bias
        request=dict(system=dict(name='harmonic_oscillator'),knots=knots,observables=['position','velocity'],
            time=time.tolist(),measurements=y.tolist(),measurement_channels=[0],measurement_variances=[.001],
            measurement_feedthrough=[.3],initial_candidates=[[1.2,.3],[1.6,.9]],input=dict(bias=bias),
            table_samples=512,input_maximum_harmonic=63,maximum_step=.04,relative_tolerance=1e-10,absolute_tolerance=1e-12,
            forecast_final_time=8,output_step=.02)
        results={}
        for model in ['physical','atlas']:
            r=identify_controlled({**request,'model':model},executable=CLI)
            np.testing.assert_allclose(r['fit']['parameters'],[radius,phase],atol=2e-6)
            t=np.array(r['time'])-2
            expected=np.column_stack((bias+(x0-bias)*np.cos(t)+v0*np.sin(t),-(x0-bias)*np.sin(t)+v0*np.cos(t)))
            np.testing.assert_allclose(np.array(r['predictions']).reshape(-1,2),expected,atol=3e-6)
            np.testing.assert_allclose(r['predicted_measurements'],expected[:,0]+.3*bias,atol=3e-6)
            self.assertEqual(r['measurement_cutoff'],5)
            results[model]=r
        shorter=identify_controlled({**request,'model':'atlas','forecast_final_time':5},executable=CLI)
        self.assertEqual(shorter['fit'],results['atlas']['fit'])
        with self.assertRaises(RuntimeError):
            identify_controlled({**request,'input':dict(state_gain=[0,-.1])},executable=CLI)
        sequential={k:v for k,v in request.items() if k not in ('forecast_final_time','output_step')}
        sequential.update(model='atlas',update_indices=[50,80],prediction_substeps=2)
        tracked=identify_controlled_sequential(sequential,executable=CLI)
        t=np.array(tracked['time'])-2
        expected=np.column_stack((bias+(x0-bias)*np.cos(t)+v0*np.sin(t),-(x0-bias)*np.sin(t)+v0*np.cos(t)))
        np.testing.assert_allclose(np.array(tracked['predictions']).reshape(-1,2),expected,atol=3e-6)
        self.assertTrue(np.all(np.array(tracked['information_time'])<=tracked['time']))
        changed=y.copy();changed[80:]+=.04
        perturbed=identify_controlled_sequential({**sequential,'measurements':changed.tolist()},executable=CLI)
        before=np.array(tracked['time'])<time[80]
        np.testing.assert_array_equal(np.array(tracked['predictions']).reshape(-1,2)[before],
                                      np.array(perturbed['predictions']).reshape(-1,2)[before])
        self.assertGreater(np.max(np.abs(np.array(tracked['predictions'])-perturbed['predictions'])),.001)


if __name__=='__main__':unittest.main()
