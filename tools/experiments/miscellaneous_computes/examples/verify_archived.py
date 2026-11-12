"""Optional local-archive parity check. Runtime tools have no dependency on the experiment folders."""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import numpy as np
side=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(side.parent))
from miscellaneous_computes.modal import ModalModel
from miscellaneous_computes.models import fit_field
from miscellaneous_computes.replay import replay
from miscellaneous_computes.diagnostics import forecast
from miscellaneous_computes.toolkit import save_json

def historical_evaluator(path):
    """Load an original experiment evaluator by path as an independent numerical reference."""
    spec=importlib.util.spec_from_file_location('archived_reference_only',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def verify(root,kernel):
    """Rebuild four selected direct fits and compare direct rollouts and structured fields to the RPE archives."""
    direct=root/'rpe-direct-modal-sindy-20261107';structured=root/'rpe-full-hsvd-delay-orders-20261107'
    summary=json.loads((direct/'summary.json').read_text())
    original=historical_evaluator(structured/'identified.py')
    figure=historical_evaluator(root/'rpe-two-echo-hsvd-sindy-20261107'/'identified.py')
    checks={'direct':{},'structured':{},'four_mode_baseline':{}}
    for rank in (6,8,10,12):
        m=ModalModel.load(direct/f'final-model-{rank}.npz')
        with np.load(direct/f'training-{rank}.npz') as packed:
            B,V,t,offset=[packed[k] for k in ('b','db','time','offsets')]
        records=[dict(time=t[a:b],b=B[a:b],db=V[a:b]) for a,b in zip(offset[:-1],offset[1:])]
        cfg=m.metadata['selection']
        fresh=fit_field(records,family='polynomial',degree=cfg['degree'],ridge=cfg['ridge'],threshold=cfg['threshold'],
                        time_weighted=True,speed_scale=24. if cfg['balanced'] else None)
        np.testing.assert_array_equal(fresh.coef!=0,m.field.coef!=0)
        coeff_error=float(np.linalg.norm(fresh.coef-m.field.coef)/np.linalg.norm(m.field.coef))
        assert coeff_error<1e-9,(rank,coeff_error)
        gate={'refit_relative_coefficients':coeff_error,'nonzero':int(np.count_nonzero(fresh.coef))}
        for ic in (0,4):
            with np.load(structured/f'record-{rank}-{ic}.npz') as z:t=z['time']
            with np.load(direct/f'final-forecast-{rank}-{ic}.npz') as z:
                xi=z['free_xi'];tm=z['modal_time'];saved=z['predicted'];initial=z['actual_xi'][0]
            q=replay(m,tm,initial=initial,method='exact-affine')
            fixed=replay(m,tm,initial=initial,method='rk4',library=kernel,max_step=.00005)
            ix=np.searchsorted(t,tm);scale=max(np.linalg.norm(saved[ix]),1e-15)
            exact=float(np.linalg.norm(q['values'][:,0]-saved[ix])/scale)
            compiled=float(np.linalg.norm(fixed['values'][:,0]-saved[ix])/scale)
            assert exact<1e-6 and compiled<1e-5 and fixed['complete'],(rank,ic,exact,compiled)
            with np.load(structured/f'forecast-4-{ic}.npz') as z:
                reference=z['actual']
                if rank==6:
                    baseline=forecast(z['time'],reference,z['prediction_time'],z['predicted'],prominence=1.)
                    checks['four_mode_baseline'][str(ic)]=dict(relative=baseline['errors']['relative'],
                        coverage=baseline['coverage'],peaks=baseline['peaks']['forecast_count'])
            e=forecast(t,reference,t,saved,prominence=1.)
            ref=summary['fits'][str(rank)][str(ic)]
            np.testing.assert_allclose(e['errors']['relative'],ref['error'],rtol=1e-10)
            assert e['peaks']['forecast_count']==ref['pulses']
            gate[str(ic)]=dict(exact_trace_error=exact,native_trace_error=compiled,pulses=ref['pulses'],
                               physical_relative_error=e['errors']['relative'])
        checks['direct'][str(rank)]=gate
    for rank in (4,6,8,10,12):
        path=structured/f'portable-{rank}.npz';m=ModalModel.load(path);reference=original.DelayModel(path);gate={}
        for ic in (0,4):
            with np.load(structured/f'forecast-{rank}-{ic}.npz') as z:xi=z['predicted_V'][::89]
            got=m.rhs(xi);expected=reference.rhs(xi)
            diff=float(np.linalg.norm(got-expected)/max(np.linalg.norm(expected),1e-20))
            np.testing.assert_allclose(m.chart.decode(xi).reshape(len(xi),-1),reference.decode(xi),rtol=1e-10,atol=1e-9)
            assert diff<1e-8,(rank,ic,diff)
            gate[str(ic)]=dict(relative_rhs_difference=diff,sampled_states=len(xi))
            if rank==4:
                old_rates=figure.fv(xi)
                figure_error=float(np.linalg.norm(got-old_rates)/np.linalg.norm(old_rates))
                assert figure_error<1e-8,(ic,figure_error)
                np.testing.assert_allclose(m.chart.decode(xi)[:,0],figure.decode_echo(xi),rtol=1e-10,atol=1e-9)
                gate[str(ic)]['original_figure_rhs_difference']=figure_error
        checks['structured'][str(rank)]=gate
    return checks

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('viewer_outputs',type=Path);ap.add_argument('kernel',type=Path);ap.add_argument('output',type=Path)
    a=ap.parse_args();save_json(a.output,verify(a.viewer_outputs.resolve(),a.kernel.resolve()))
