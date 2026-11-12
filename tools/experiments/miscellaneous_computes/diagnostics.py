"""Free forecast auditing on the reference clock, including extrema, coverage and kinematics."""
import numpy as np
from scipy.signal import find_peaks
from .history import matrix
from .sampling import clock,trapezoid_weights

def errors(reference,predicted,weights=None):
    """Return channelwise RMS and relative errors and a pooled relative norm using the supplied sample weights."""
    x=matrix(reference);y=matrix(predicted)
    if x.shape!=y.shape:raise ValueError('Reference and prediction shapes disagree')
    w=np.ones(len(x)) if weights is None else np.asarray(weights,float)
    if w.shape!=(len(x),) or np.any(w<0) or not np.isfinite(w).all() or w.sum()<=0:
        raise ValueError('Invalid sample weights')
    scale=np.sqrt(w@np.square(x)/w.sum());rmse=np.sqrt(w@np.square(y-x)/w.sum())
    norm=np.linalg.norm(scale)
    return dict(rmse=rmse.tolist(),reference_rms=scale.tolist(),
                relative=[float(a/b) if b else None for a,b in zip(rmse,scale)],
                total=float(np.linalg.norm(rmse)/norm) if norm else None)

def forecast(time,actual,prediction_time,predicted,*,time_weighted=False,prominence=None,channel=0):
    """Interpolate predictions only onto reference times within the available forecast.
    Events use their original native clocks. Unfinished peak intervals are flagged.
    """
    t=clock(time);x=matrix(actual);pt=np.asarray(prediction_time,float);p=np.asarray(predicted,float)
    if p.ndim==1:p=p[:,None]
    if (len(x)!=len(t) or p.shape!=(len(pt),x.shape[1]) or not np.isfinite(p).all()
        or not 0<=channel<x.shape[1]):raise ValueError('Invalid forecast arrays')
    if len(pt):pt=clock(pt,1)
    mask=np.zeros(len(t),dtype=bool) if not len(pt) else (t>=pt[0])&(t<=pt[-1])
    shared=t[mask]
    out=dict(requested=[float(t[0]),float(t[-1])],sample_coverage=float(mask.mean()),complete=bool(mask.all()),
             reference_count=len(t),compared_count=len(shared),prediction_count=len(pt))
    if not len(shared):return dict(out,coverage=0.,compared=None,errors=None)
    aligned=np.column_stack([np.interp(shared,pt,col) for col in p.T])
    covered=max(0.,min(pt[-1],t[-1])-max(pt[0],t[0]))
    out.update(coverage=float(covered/(t[-1]-t[0])),compared=[shared[0],shared[-1]],
               errors=errors(x[mask],aligned,trapezoid_weights(shared) if time_weighted and len(shared)>1 else None))
    if prominence is not None:
        if not np.isfinite(prominence) or prominence<0:raise ValueError('Invalid peak prominence')
        true=find_peaks(x[:,channel],prominence=prominence)[0]
        peaks=find_peaks(p[:,channel],prominence=prominence)[0]
        edges=np.r_[t[0],(t[true[:-1]]+t[true[1:]])/2,t[-1]]
        matches=[]
        for j,ref in enumerate(true):
            left,right=edges[j:j+2];ids=peaks[(pt[peaks]>=left)&(pt[peaks]<right)]
            item=dict(reference_time=t[ref],reference_value=x[ref,channel],count=len(ids),
                      full_interval=bool(pt[0]<=left and pt[-1]>=right))
            if len(ids):
                i=ids[np.argmax(p[ids,channel])];height=x[ref,channel]
                item.update(time=pt[i],value=p[i,channel],time_error=pt[i]-t[ref],
                            relative_height=float(p[i,channel]/height-1) if height else None)
            matches.append(item)
        out['peaks']=dict(reference_count=len(true),forecast_count=int(np.sum((pt[peaks]>=t[0])&(pt[peaks]<=t[-1]))),matches=matches)
    return out

def kinematics(model,xi,*,value_channel,rate_channel,tap=0,true_rate=None):
    """Decode the modeled rate and independently differentiate the modeled value."""
    v=matrix(xi);pairs=model.chart.decode(v)
    derivative=np.vstack([model.rhs(v[i:i+1024]) for i in range(0,len(v),1024)])
    chain=model.chart.decode(derivative)[:,tap,value_channel]
    rate=pairs[:,tap,rate_channel];gap=chain-rate
    reference=rate if true_rate is None else np.asarray(true_rate,float)
    if reference.shape!=rate.shape or not np.isfinite(reference).all():raise ValueError('Invalid reference rate')
    scale=np.linalg.norm(reference)
    return dict(decoded_rate=rate,chain_rate=chain,
                relative_mismatch=float(np.linalg.norm(gap)/scale) if scale else None)

def refine_rollout(integrate, time, initial, coarse_options, fine_options, reference):
    """Compare two independent calls on their common prefix, normalized by reference RMS."""
    first=integrate(time,initial,**coarse_options);second=integrate(time,initial,**fine_options)
    n=min(len(first['time']),len(second['time']))
    if not np.array_equal(first['time'][:n],second['time'][:n]):raise ValueError('Integrator output clocks differ')
    target=matrix(reference)
    if len(target)!=len(time):raise ValueError('Reference clock length mismatch')
    discrepancy=None if not n else errors(target[:n],target[:n]+second['b'][:n]-first['b'][:n])
    return dict(first_complete=bool(first['complete']),second_complete=bool(second['complete']),shared_samples=n,
                difference=discrepancy)
