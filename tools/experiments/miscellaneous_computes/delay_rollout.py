"""Method-of-steps integration of an observation law on its own causal delay history."""
import ctypes as ct
from pathlib import Path
import numpy as np

DOUBLE=ct.POINTER(ct.c_double)
def _pointer(array):
    return array.ctypes.data_as(DOUBLE)

def delay_rollout(library, rhs, history, lags, step, samples, *, context=None, limits=None):
    """Advance from a uniformly sampled history ending at time zero.

    rhs may be a C function pointer (context, flattened window, output) or a Python
    callable taking a taps-by-channels array. The returned relative clock starts at zero.
    Delay interpolation is linear on accepted states, including the supplied prehistory.
    """
    hist=np.ascontiguousarray(history,dtype=float);lag=np.ascontiguousarray(lags,dtype=float)
    if hist.ndim!=2 or not hist.shape[1] or lag.ndim!=1 or not len(lag) or not isinstance(samples,int) or samples<1:
        raise ValueError("Expected aligned history, lag vector, and positive output count")
    dim=hist.shape[1];native=ct.CDLL(str(Path(library).resolve()))
    callback=ct.CFUNCTYPE(ct.c_int,ct.c_void_p,DOUBLE,DOUBLE)
    errors=[]
    if isinstance(rhs,ct._CFuncPtr):
        handle=ct.cast(rhs,ct.c_void_p)
    else:
        @callback
        def python_rhs(_,x,y):
            try:
                window=np.ctypeslib.as_array(x,shape=(len(lag)*dim,)).reshape(len(lag),dim)
                value=np.asarray(rhs(window),dtype=float)
                if value.shape!=(dim,):raise ValueError("Delay RHS returned the wrong state dimension")
                np.ctypeslib.as_array(y,shape=(dim,))[:]=value
                return 0
            except Exception as exc:
                errors.append(exc);return 1
        handle=ct.cast(python_rhs,ct.c_void_p)
    fn=native.mc_delay_rk4_v1
    fn.argtypes=[ct.c_void_p,ct.c_void_p,ct.c_size_t,DOUBLE,ct.c_size_t,DOUBLE,ct.c_size_t,ct.c_double,
                 ct.c_size_t,DOUBLE,DOUBLE,ct.POINTER(ct.c_size_t)]
    fn.restype=ct.c_int
    bound=None if limits is None else np.ascontiguousarray(np.broadcast_to(limits,(dim,)),dtype=float)
    out=np.empty((samples,dim));written=ct.c_size_t()
    status=fn(handle,context,dim,_pointer(lag),len(lag),_pointer(hist),len(hist),step,samples,
              None if bound is None else _pointer(bound),_pointer(out),ct.byref(written))
    if errors:raise errors[0]
    if status==3:raise ValueError("Invalid delay grid, history, limits or workspace")
    return dict(time=step*np.arange(written.value),values=out[:written.value],complete=status==0,
                message=("complete","state limit","nonfinite field")[status])
