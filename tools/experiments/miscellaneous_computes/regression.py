"""Blockwise derivative and weak-form regression moments on separate trajectory clocks."""
import numpy as np
from .sampling import trapezoid_weights

def _primitive(time,states,design,at,block_size):
    """Piecewise-linear library integral at sorted query times with bounded feature workspace."""
    from scipy.integrate import cumulative_trapezoid
    at=np.asarray(at);out=None;carry=None
    for i in range(0,len(time)-1,block_size):
        end=min(i+block_size+1,len(time));t=time[i:end];f=design(states[i:end])
        integral=cumulative_trapezoid(f,t,axis=0,initial=0)
        if carry is not None:integral+=carry
        if out is None:out=np.empty((len(at),f.shape[1]))
        lo=at.searchsorted(t[0]);hi=at.searchsorted(t[-1],side='right')
        q=at[lo:hi]
        if len(q):
            j=np.clip(np.searchsorted(t,q,side='right')-1,0,len(t)-2)
            h=(q-t[j])[:,None];dt=(t[j+1]-t[j])[:,None]
            out[lo:hi]=integral[j]+h*(f[j]+.5*h/dt*(f[j+1]-f[j]))
        carry=integral[-1].copy()
    return out

def weighted_moments(records, design, *, time_weighted=False, speed_scale=None,
                     weak_durations=(), window_stride=1, block_size=1024):
    """Accumulate an augmented Gram matrix from instantaneous derivative rows or physical-duration state increments."""
    if speed_scale is not None and (not np.isfinite(speed_scale) or speed_scale<=0):
        raise ValueError("Invalid speed scale")
    if not isinstance(window_stride,int) or window_stride<1 or block_size<1:
        raise ValueError("Invalid regression stride or block size")
    durations=tuple(weak_durations)
    if any(not np.isfinite(v) or v<=0 for v in durations):
        raise ValueError("Weak durations must be positive and finite")
    gram=None
    def add(features, targets, mass):
        nonlocal gram
        weight=np.sqrt(mass)
        if speed_scale is not None:
            weight=weight/np.sqrt(1+np.linalg.norm(targets,axis=1)/speed_scale)
        rows=np.column_stack([features,targets])*weight[:,None]
        value=rows.T@rows
        if gram is None:gram=value
        else:gram+=value
    for rec in records:
        t,x,v=rec['time'],rec['b'],rec['db']
        if not durations:
            mass=trapezoid_weights(t) if time_weighted else np.ones(len(t))
            for start in range(0,len(x),block_size):
                sl=slice(start,start+block_size)
                add(design(x[sl]),v[sl],mass[sl])
            continue
        for duration in durations:
            starts=np.flatnonzero(t+duration<=np.nextafter(t[-1],np.inf))[::window_stride]
            mass=trapezoid_weights(t[starts]) if time_weighted and len(starts)>1 else np.ones(len(starts))
            for i in range(0,len(starts),block_size):
                selected=starts[i:i+block_size]
                lo=selected[0];hi=min(len(t),np.searchsorted(t,t[selected[-1]]+duration)+1)
                begin=t[selected];end=np.minimum(begin+duration,t[-1])
                both=np.r_[begin,end];order=np.argsort(both)
                primitive=_primitive(t[lo:hi],x[lo:hi],design,both[order],block_size)
                integrals=np.empty_like(primitive);integrals[order]=primitive
                right=np.clip(t.searchsorted(end),1,len(t)-1);left=right-1
                frac=((end-t[left])/(t[right]-t[left]))[:,None]
                final=x[left]+frac*(x[right]-x[left])
                add((integrals[len(begin):]-integrals[:len(begin)])/duration,(final-x[selected])/duration,
                    mass[i:i+len(selected)]/len(durations))
    if gram is None:raise ValueError("No complete regression windows")
    return gram
