"""Identify critically damped forced responses from waveforms using a filtered SINDy library."""
import numpy as np
from scipy.signal import lfilter
from .features import prepare_library,library,thresholded_least_squares
from .models import SparseField
from .history import uniform_clock,matrix

def filtered_library(time,features,gamma):
    """Return response columns of x''+2γx'+γ²x=γ²φ from zero initial conditions.

    Each input interval holds the arithmetic midpoint of its endpoint features.
    The filter update is exact for that held input. Caller supplies one record.
    """
    F=matrix(features);t,h=uniform_clock(time,len(F))
    if not np.isfinite(gamma) or gamma<=0:raise ValueError('A positive finite damping rate is required')
    x=gamma*h;q=np.exp(-x)
    h0=(-np.expm1(-x)-x*q) if x>1e-4 else x*x*(.5-x/3+x*x/8)
    h1=gamma*x*q
    A00=q*(1+x);A01=h*q;A10=-gamma*x*q;A11=(1-x)*q
    mid=.5*(F[:-1]+F[1:]);den=[1.,-2*q,q*q]
    def filter_pair(num):
        return np.vstack([np.zeros(F.shape[1]),lfilter(num,den,mid,axis=0)])
    return filter_pair([h0,A01*h1-A11*h0]),filter_pair([h1,A10*h0-A00*h1])

def fit_damped_drive(records,gamma,*,family='localized',degree=1,centers=32,ridge=1e-5,threshold=.002):
    """Fit physical forcing for InteriorField from (time, input, position, velocity).
    All records share input and output dimensions; each has its own initial conditions.
    """
    data=[(np.asarray(t,float),matrix(x),matrix(p),matrix(v)) for t,x,p,v in records]
    if not data or any(len(t)!=len(x) or len(t)!=len(p) or p.shape!=v.shape or p.shape[1]!=data[0][2].shape[1] for t,x,p,v in data):
        raise ValueError('Expected aligned inputs and position/velocity labels')
    spec=prepare_library(np.concatenate([x for _,x,_,_ in data]),family,degree,centers)
    gram=0.
    for t,x,p,v in data:
        pos,vel=filtered_library(t,library(x,spec),gamma)
        dt=(t-t[0])[:,None];e=np.exp(-gamma*dt)
        transient_p=e*((1+gamma*dt)*p[0]+dt*v[0])
        transient_v=e*((1-gamma*dt)*v[0]-gamma**2*dt*p[0])
        for design,target,observed in [(pos,p-transient_p,p),(vel,v-transient_v,v)]:
            scale=np.sqrt(np.mean(observed**2))
            rows=np.column_stack([design,target])/(scale if scale else 1.)
            gram+=rows.T@rows
    return SparseField(spec,gamma**2*thresholded_least_squares(gram,data[0][2].shape[1],ridge,threshold))
