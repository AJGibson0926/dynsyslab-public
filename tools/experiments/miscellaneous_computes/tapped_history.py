"""Sample fractional causal delays and fit an SVD to the selected observable/rate windows."""
import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicHermiteSpline
from .history import uniform_clock

def sample_delay_jet(time,jet,endpoints,lags):
    """Interpolate (y,y',y'') at endpoint-lag using one C² quintic Hermite curve.

    The result is endpoint by tap by derivative order. Extrapolation is rejected.
    At a fractional delay the two bracketing source samples must both be available
    at the endpoint. Zero lag on a native tick returns that exact sample.
    """
    jets=np.asarray(jet,dtype=float)
    if jets.ndim!=2 or jets.shape[1]!=3 or not np.isfinite(jets).all():
        raise ValueError('Expected aligned finite value, rate and acceleration columns')
    time,dt=uniform_clock(time,len(jets))
    t=np.asarray(endpoints,dtype=float);lags=np.asarray(lags,dtype=float)
    if (t.ndim!=1 or not len(t) or lags.ndim!=1 or not len(lags) or np.any(lags<0) or
        not np.isfinite(t).all() or not np.isfinite(lags).all()):
        raise ValueError('Invalid endpoints or delays')
    ticks=(t[:,None]-lags-time[0])/dt
    nearest=np.rint(ticks)
    exact=abs(ticks-nearest)<=2e-8
    ticks=np.where(exact,nearest,ticks)
    if np.min(ticks)<0 or np.max(ticks)>len(time)-1:
        raise ValueError('Insufficient source history')
    if np.any(time[np.ceil(ticks).astype(int)]>t[:,None]+dt*2e-8):
        raise ValueError('Fractional interpolation would use an unobserved future sample')
    i=np.minimum(ticks.astype(np.int64),len(time)-2)
    u=ticks-i
    p,v,a=np.moveaxis(jets[i],-1,0);q,w,b=np.moveaxis(jets[i+1],-1,0)
    C1=dt*v;C2=.5*dt**2*a
    d=(q-p)-C1-C2;s=dt*w-C1-2*C2;c=dt**2*b-2*C2
    c3=10*d-4*s+.5*c;c4=-15*d+7*s-c;c5=6*d-3*s+.5*c
    y=p+u*(C1+u*(C2+u*(c3+u*(c4+u*c5))))
    dy=(C1+u*(2*C2+u*(3*c3+u*(4*c4+5*u*c5))))/dt
    ddy=(2*C2+u*(6*c3+u*(12*c4+20*u*c5)))/dt**2
    out=np.stack([y,dy,ddy],axis=-1)
    if np.any(exact):out[exact]=jets[nearest[exact].astype(np.int64)]
    return out

def full_tap_svd(blocks,scales,rank=None):
    """Fit all signed right-singular coordinates on disjoint row blocks of measured histories.
    Blocks have shape (windows, taps, channels). Return row-vector maps, with
    xi = raw @ encoder and raw = xi @ decoder. xi is signed V on training rows.
    """
    sc=np.asarray(scales,dtype=float)
    if sc.ndim!=1 or not sc.size or not np.isfinite(sc).all() or np.any(sc<=0):
        raise ValueError('Positive finite channel scales required')
    G=None;N=0
    for a in blocks:
        a=np.asarray(a,dtype=float)
        if a.ndim!=3 or a.shape[2]!=len(sc) or not a.shape[1] or not np.isfinite(a).all():
            raise ValueError('Invalid training windows')
        width=a.shape[1]*len(sc);H=(a/sc).reshape(len(a),width)
        if G is None:G=np.zeros((width,width))
        if width!=len(G):raise ValueError('Incompatible tap counts')
        G+=H.T@H;N+=len(H)
    if not N:raise ValueError('No training windows')
    ev,U=eigh(G);ev=ev[::-1];U=U[:,::-1]
    width=len(ev)
    rank=width if rank is None else rank
    if not isinstance(rank,(int,np.integer)) or not 1<=rank<=width or ev[rank-1]<=np.finfo(float).eps*width*ev[0]:
        raise ValueError('History is numerically rank deficient at the requested rank')
    ev=ev[:rank];U=U[:,:rank]
    U*=np.sign(U[np.argmax(abs(U),axis=0),np.arange(len(ev))])
    sigma=np.sqrt(ev);units=np.tile(sc,width//len(sc))
    return dict(U=U,sigma=sigma,scales=sc,windows=N,encoder=(U/sigma)/units[:,None],
                decoder=(sigma[:,None]*U.T)*units[None,:],gram_error=np.max(abs((U.T@G@U)/sigma[:,None]/sigma-np.eye(len(ev)))))


def delay_consistency(time, prediction, reference, lags, stride=1):
    """Relative RMS discrepancy between each free tap and a shifted copy of the free current pair.

    Columns are (value, rate) per tap, starting with zero delay. Reference values
    set each error scale. The comparison excludes unavailable forecast history.
    Undefined relative errors are returned as None.
    """
    t=np.asarray(time,dtype=float);x=np.asarray(prediction,dtype=float);r=np.asarray(reference,dtype=float)
    d=np.asarray(lags,dtype=float)
    if (t.ndim!=1 or len(t)<2 or not np.isfinite(t).all() or np.any(np.diff(t)<=0)
       or d.ndim!=1 or not len(d) or d[0]!=0. or np.any(d<0) or not np.isfinite(d).all()
       or x.shape!=(len(t),2*len(d)) or r.shape!=x.shape or not np.isfinite(x).all()
       or not np.isfinite(r).all() or not isinstance(stride,(int,np.integer)) or stride<1):
        raise ValueError('Expected an increasing clock, finite aligned (value, rate) pairs and nonnegative lags')
    interp=CubicHermiteSpline(t,x[:,0],x[:,1])
    errors=[]
    for j,lag in enumerate(d[1:],1):
        ids=np.flatnonzero(t-lag>=t[0])[::stride]
        if len(ids)<2:
            errors.append(None)
            continue
        clock=t[ids]-lag;delayed=np.column_stack([interp(clock),interp(clock,1)])
        sl=slice(2*j,2*j+2)
        numerator=np.linalg.norm(delayed-x[ids,sl],axis=0)
        denom=np.linalg.norm(r[ids,sl],axis=0)
        errors.append([float(a/b) if b>0 else None for a,b in zip(numerator,denom)])
    return errors
