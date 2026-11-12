"""Training-only physical-output weights and parameter directions for fixed-support waveform fitting."""
import numpy as np
from .sampling import trapezoid_weights


def waveform_weights(time, observations, global_weights, local_weights, duration, rms_floors, *, global_floors=None):
    """Return sample-by-channel multipliers for a squared waveform residual.

    The sum of squared weighted differences combines global-RMS and fixed-interval local-RMS
    normalized errors. Time quadrature works on irregular clocks. Weights are loss weights,
    and floors have the units of their respective observation channels.
    """
    t=np.asarray(time,float)
    y=np.asarray(observations,float)
    if y.ndim!=2 or len(y)!=len(t) or not np.isfinite(y).all() or not np.isfinite(duration) or duration<=0.:
        raise ValueError("Supply aligned finite waveform channels and a positive interval")
    global_w,local_w,floor=[np.broadcast_to(np.asarray(v,float),(y.shape[1],)) for v in
                            (global_weights,local_weights,rms_floors)]
    global_floor=floor if global_floors is None else np.broadcast_to(np.asarray(global_floors,float),floor.shape)
    if (np.any(global_w<0.) or np.any(local_w<0.) or np.any(floor<=0.) or
        np.any(global_floor<=0.) or not all(np.isfinite(v).all() for v in (global_w,local_w,floor,global_floor))):
        raise ValueError("Waveform weights must be nonnegative and RMS floors strictly positive")
    mass=trapezoid_weights(t);mass/=mass.sum()
    total=np.maximum(mass@(y*y),global_floor**2)
    local=np.empty_like(y)
    blocks=np.floor((t-t[0])/duration)
    for block in np.unique(blocks):
        take=blocks==block
        local[take]=np.maximum(np.average(y[take]**2,axis=0,weights=mass[take]),floor**2)
    return np.sqrt(mass[:,None]*(global_w/total+local_w/local))


def resolved_directions(response, coefficient_directions, *, relative_cutoff=.001, singular_floor=.015, scale=.035):
    """Rotate proposed directions using a training response Jacobian, keeping only resolved combinations.

    Columns of response and coefficient_directions refer to the same perturbations. Rows of the latter
    refer to coefficients. Zero rows remain zero, so a fixed sparse support is preserved.
    Returns (new coefficient directions, corresponding linearized responses).
    """
    x=np.asarray(response,float);v=np.asarray(coefficient_directions,float)
    if (x.ndim!=2 or v.ndim!=2 or x.shape[1]!=v.shape[1] or not x.size or not v.size or
        not np.isfinite(x).all() or not np.isfinite(v).all() or not 0<=relative_cutoff<1 or
        not np.isfinite([singular_floor,scale]).all() or min(singular_floor,scale)<=0.):
        raise ValueError("Invalid training response or coefficient subspace")
    norms=np.linalg.norm(x,axis=0);active=norms>norms.max()*1e-10
    if not active.any():raise ValueError("Training data resolve no supplied coefficient direction")
    _,s,vt=np.linalg.svd(x[:,active]/norms[active],full_matrices=False)
    keep=s>s[0]*relative_cutoff
    w=(vt[keep].T/np.maximum(s[keep],singular_floor))*(scale/norms[active,None])
    return v[:,active]@w, x[:,active]@w
