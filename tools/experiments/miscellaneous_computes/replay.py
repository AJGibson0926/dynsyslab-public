"""Replay a deployed modal field from one window using its declared coordinate convention."""
import numpy as np
from .models import free_rollout,SparseField
from .integration import affine_equations,exact_affine,polynomial_rk4
from .sampling import clock

def replay(model,time,*,initial=None,window=None,method='adaptive',library=None,**options):
    """Advance one signed state or one measurement window and return the accepted clock, signed states and decoded taps."""
    if (initial is None)==(window is None):raise ValueError('Supply one initial modal vector or one initial window')
    xi=np.asarray(model.chart.project(window) if window is not None else initial,float)
    t=clock(time);gain=model.chart.gain
    if xi.shape!=(model.chart.encoder.shape[1],) or not np.isfinite(xi).all():raise ValueError('Invalid initial state')
    if method=='adaptive':
        # Integrate unit-training-RMS amplitudes for consistent absolute tolerances across record lengths.
        result=free_rollout(lambda b:gain*model.rhs(b/gain),t,xi*gain,**options)
    elif method in ('rk4','exact-affine'):
        if model.space!='direct' or model.feature_asinh or not isinstance(model.field,SparseField):
            raise ValueError(f'{method} expects a direct polynomial field')
        if method=='rk4':
            if library is None:raise ValueError('Provide the compiled kernel library for RK4')
            result=polynomial_rk4(model.field,t,xi*gain,library,**options)
        else:
            if options:raise ValueError('Exact affine propagation takes no solver options')
            A,c=affine_equations(model.field);b=exact_affine(A,c,t,xi*gain)
            if not np.isfinite(b).all():raise ArithmeticError('Nonfinite affine trajectory')
            result=dict(time=t,b=b,complete=True,evaluations=0,message='matrix exponential')
    else:raise ValueError(f'Unknown integration method {method}')
    b=result.pop('b');v=b.reshape(-1,len(xi))/gain
    result.update(xi=v,values=model.chart.decode(v))
    return result
