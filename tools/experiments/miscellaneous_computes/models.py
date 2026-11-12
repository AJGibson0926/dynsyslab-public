"""Fit, evaluate and serialize autonomous sparse fields for arbitrary modal records."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from .history import matrix
from .features import prepare_library,library,regression_moment,thresholded_least_squares

class SparseField:
    """Evaluate a frozen feature dictionary with one sparse coefficient row per output coordinate."""
    def __init__(self,spec,coefficients):
        self.spec=spec;self.coef=matrix(coefficients)
    def __call__(self,x):
        result=library(np.atleast_2d(x),self.spec)@self.coef.T
        return result[0] if np.ndim(x)==1 else result
    def save(self,path):
        with open(path,'x',encoding='utf-8') as f:
            json.dump(dict(schema='dynsyslab.miscellaneous-computes/sparse-field/v1',
                           library=self.spec,coef=self.coef.tolist()),f,allow_nan=False)
    @classmethod
    def load(cls,path):
        with open(path,encoding='utf-8') as f:data=json.load(f)
        if data.get('schema')!='dynsyslab.miscellaneous-computes/sparse-field/v1':
            raise ValueError('Unsupported sparse field schema')
        return cls(data['library'],data['coef'])

def records_checked(records):
    """Validate independent timestamped state and derivative blocks and return normalized array shapes."""
    if not records:raise ValueError('At least one independent trajectory is required')
    out=[]
    for rec in records:
        b=matrix(rec['b']);db=matrix(rec['db']);t=np.asarray(rec['time'],dtype=float)
        if db.shape!=b.shape or len(t)<2 or t.shape!=(len(b),) or not np.isfinite(t).all() or np.any(np.diff(t)<=0):
            raise ValueError('Expected aligned, finite coordinates, rates and increasing time')
        if out and b.shape[1]!=out[0]['b'].shape[1]:raise ValueError('State dimensions disagree')
        out.append(dict(time=t,b=b,db=db))
    return out

def fit_field(records,*,family='polynomial',degree=2,centers=24,harmonics=5,ridge=1e-6,threshold=.003,weak_window=0,
              time_weighted=False,speed_scale=None,weak_durations=(),window_stride=1,rational_scales=(.25,1.,4.)):
    """records are independent {time,b,db} mappings; no windows span two records."""
    records=records_checked(records)
    spec=prepare_library(np.vstack([r['b'] for r in records]),family,degree,centers,harmonics,rational_scales)
    if time_weighted or speed_scale is not None or weak_durations or window_stride!=1:
        from .regression import weighted_moments
        if weak_window:raise ValueError('Choose sample-count or physical-duration windows')
        gram=weighted_moments(records,lambda x:library(x,spec),time_weighted=time_weighted,
                              speed_scale=speed_scale,weak_durations=weak_durations,window_stride=window_stride)
        return SparseField(spec,thresholded_least_squares(gram,spec['n'],ridge,threshold))
    gram=0.
    for rec in records:
        if weak_window:
            gram+=regression_moment(library(rec['b'],spec),rec['b'],rec['db'],rec['time'],weak_window)
        else:
            for i in range(0,len(rec['time']),2048):
                sl=slice(i,i+2048);x=library(rec['b'][sl],spec)
                gram+=regression_moment(x,rec['b'][sl],rec['db'][sl],rec['time'][sl],0)
    return SparseField(spec,thresholded_least_squares(gram,spec['n'],ridge,threshold))

def free_rollout(field,time,initial,*,rtol=1e-8,atol=1e-10,max_step=np.inf,max_abs=np.inf):
    """Integrate once from the supplied initial state. Preserve partial coverage and solver status."""
    t=np.asarray(time,dtype=float);x=np.asarray(initial,dtype=float)
    if t.ndim!=1 or len(t)<2 or not np.isfinite(t).all() or np.any(np.diff(t)<=0):
        raise ValueError('Expected at least two strictly increasing finite output times')
    if x.ndim!=1 or not x.size or not np.isfinite(x).all() or np.max(abs(x))>=max_abs:
        raise ValueError('Invalid initial state or guard')
    def guard(_,b):return max_abs-np.max(abs(b))
    guard.terminal=True;guard.direction=-1
    sol=solve_ivp(lambda _,b:field(b),(t[0],t[-1]),x,t_eval=t,method='DOP853',rtol=rtol,atol=atol,
                  max_step=max_step,events=guard if np.isfinite(max_abs) else None)
    return dict(time=np.asarray(sol.t),b=np.asarray(sol.y).T,complete=len(sol.t)==len(t) and sol.success,
                message=sol.message,evaluations=sol.nfev)
