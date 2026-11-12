"""Record-local sampling weights, causal delays and physical-duration integral rows."""
import numpy as np
from scipy import sparse

def clock(time, minimum=2):
    """Validate the minimum length, finiteness and strict ordering of a one-dimensional time grid."""
    time = np.asarray(time, dtype=float)
    if time.ndim != 1 or len(time) < minimum or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError("Expected a finite strictly increasing clock")
    return time

def trapezoid_weights(time):
    """Return nodal integration weights for the piecewise-linear interpolant on an increasing clock."""
    t = clock(time)
    dt = np.diff(t)
    return .5 * (np.r_[0., dt] + np.r_[dt, 0.])

def sample_indices(velocity, stride=1, refinements=()):
    """Union of regular strides at user-specified speed thresholds, including both endpoints."""
    v = np.asarray(velocity, dtype=float)
    if v.ndim != 2 or not len(v) or not np.isfinite(v).all():
        raise ValueError("Expected finite velocity rows")
    index = np.arange(len(v))
    speed = np.linalg.norm(v, axis=1)
    keep = np.zeros(len(v), dtype=bool)
    for threshold, step in [(0., stride), *refinements]:
        if not np.isfinite(threshold) or threshold < 0 or not isinstance(step,(int,np.integer)) or step < 1:
            raise ValueError("Invalid speed threshold or stride")
        keep |= (speed >= threshold) & (index % step == 0)
    keep[[0,-1]] = True
    return index[keep]

def delay_samples(time, values, endpoints, lags):
    """Piecewise-linear sampled observations on each completed causal window, in supplied lag order."""
    t = clock(time)
    y = np.asarray(values, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    ends, delays = map(lambda a: np.asarray(a,dtype=float), (endpoints,lags))
    if (y.ndim != 2 or len(y) != len(t) or not np.isfinite(y).all() or ends.ndim != 1
        or delays.ndim != 1 or not len(delays) or not np.isfinite(ends).all()
        or not np.isfinite(delays).all() or np.any(delays < 0)):
        raise ValueError("Invalid delayed measurement inputs")
    q = ends[:, None] - delays
    # Snap only roundoff-sized differences to native ticks (including the start boundary).
    i = np.searchsorted(t,q).clip(0,len(t)-1)
    left = (i-1).clip(0)
    near = np.where(abs(t[i]-q)<abs(t[left]-q),i,left)
    eps = np.minimum(8*np.spacing(np.maximum(1.,abs(q))), np.min(np.diff(t))*1e-6)
    q = np.where(abs(t[near]-q)<=eps,t[near],q)
    right = np.searchsorted(t,q).clip(0,len(t)-1)
    if np.any(q<t[0]) or np.any(q>t[-1]) or np.any(t[right]>ends[:,None]+eps):
        raise ValueError("Incomplete causal history at these endpoints")
    left = (right-1).clip(0)
    alpha = np.divide(q-t[left],t[right]-t[left],out=np.zeros_like(q),where=right!=left)
    return y[left]+alpha[:,:,None]*(y[right]-y[left])

def weak_rows(time, duration, starts):
    """Sparse exact averages of a piecewise-linear library, with physical-time windows.

    Starts are observed row indices. Return the averaging operator and a second
    operator for (linearly interpolated final state - initial state)/duration.
    """
    t=clock(time)
    idx=np.asarray(starts)
    if (not np.isfinite(duration) or duration<=0 or idx.ndim!=1 or not np.issubdtype(idx.dtype,np.integer)
        or np.any(idx<0) or np.any(idx>=len(t))):
        raise ValueError("Invalid weak window or start indices")
    end=t[idx]+duration
    if np.any(end>np.nextafter(t[-1],np.inf)):
        raise ValueError("Weak window crosses the record boundary")
    end=np.minimum(end,t[-1])
    stop=np.searchsorted(t,end,side='left')
    if np.any(stop<=idx):
        raise ValueError("Window is below clock resolution")
    counts=stop-idx+1
    indptr=np.r_[0,np.cumsum(counts)]
    ids=np.arange(indptr[-1])+np.repeat(idx-indptr[:-1],counts)
    rows=np.repeat(np.arange(len(idx)),counts)
    w=np.zeros(len(ids))
    valid=ids<stop[rows]
    j=np.flatnonzero(valid)
    length=np.minimum(t[ids[j]+1],end[rows[j]])-t[ids[j]]
    right=.5*length**2/(t[ids[j]+1]-t[ids[j]])
    w[j]+=(length-right)/duration
    w[j+1]+=right/duration
    avg=sparse.csr_matrix((w,ids,indptr),shape=(len(idx),len(t)))
    f=(end-t[stop-1])/(t[stop]-t[stop-1])
    delta=sparse.csr_matrix((np.column_stack([-np.ones(len(idx)),1-f,f]).ravel()/duration,
                np.column_stack([idx,stop-1,stop]).ravel(),3*np.arange(len(idx)+1)),shape=avg.shape)
    return avg,delta
