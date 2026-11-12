"""Optional simulation-only recovery of a missing causal history."""
import numpy as np
from scipy.integrate import solve_ivp

def backward_history(flow,observe,time_before_initial,initial_time,initial_state,*,rtol=1e-11,atol=1e-13,max_step=np.inf,roundtrip_limit=1e-8):
    """Generate observations before the supplied initial state and verify that the flow returns to it."""
    times=np.array(time_before_initial,dtype=float);x=np.asarray(initial_state,dtype=float)
    if times.ndim!=1 or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times)<=0) or times[-1]>=initial_time:
        raise ValueError('Prehistory times must be finite, increasing and strictly before initialization')
    opts=dict(method='DOP853',rtol=rtol,atol=atol,max_step=max_step)
    past=solve_ivp(flow,(initial_time,times[0]),x,t_eval=times[::-1],**opts)
    if not past.success:raise RuntimeError(past.message)
    forward=solve_ivp(flow,(times[0],initial_time),past.y[:,-1],**opts)
    if not forward.success:raise RuntimeError(forward.message)
    error=float(np.max(abs(forward.y[:,-1]-x)))
    if error>roundtrip_limit:raise ValueError(f'Backward/forward state discrepancy {error:g}')
    return np.array([observe(t,state) for t,state in zip(times,past.y.T[::-1])]),error
