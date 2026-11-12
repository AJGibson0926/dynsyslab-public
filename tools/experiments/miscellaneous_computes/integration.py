"""Compiled fixed-step polynomial rollouts and independent affine propagation."""
import ctypes as ct
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from .sampling import clock

_D=ct.POINTER(ct.c_double)
def _ptr(a):
    return a.ctypes.data_as(_D)

def polynomial_arrays(field):
    """Translate the public polynomial dictionary into portable integer exponents."""
    spec=field.spec
    if spec['family']!='polynomial':raise ValueError('Expected a polynomial library')
    n=spec['n'];powers=np.zeros((1+len(spec['powers']),n),dtype=np.int32)
    for i,term in enumerate(spec['powers'],1):
        for j in term:
            if not 0<=j<n:raise ValueError('Invalid feature variable')
            powers[i,j]+=1
    coef=np.ascontiguousarray(field.coef,dtype=float)
    if coef.shape!=(n,len(powers)):raise ValueError('Expected one complete coefficient row per state')
    return powers,coef

def affine_equations(field):
    """Return A, c in column convention x'=A x+c, accepting unused higher-order terms."""
    powers,coef=polynomial_arrays(field);n=coef.shape[0]
    A=np.zeros((n,n));bias=np.zeros(n)
    for exponents,row in zip(powers,coef.T):
        degree=sum(exponents)
        if degree==0:bias+=row
        elif degree==1:A[:,np.argmax(exponents)]+=row
        elif np.any(row):raise ValueError('Active nonlinear terms in this field')
    return A,bias

def exact_affine(A,bias,time,initial):
    """Propagate with an augmented matrix exponential, including nonuniform output times."""
    t=clock(time,1);A=np.asarray(A,float);b=np.asarray(bias,float);x=np.asarray(initial,float)
    n=len(b)
    if A.shape!=(n,n) or x.shape!=(n,) or b.ndim!=1 or not n or not all(np.isfinite(a).all() for a in (A,b,x)):
        raise ValueError('Invalid affine field or initial state')
    M=np.zeros((n+1,n+1));M[:n,:n]=A;M[:n,n]=b
    init=np.r_[x,1.]
    if len(t)==1:return x[None].copy()
    dt=(t[-1]-t[0])/(len(t)-1)
    if np.allclose(t-t[0],np.arange(len(t))*dt,rtol=0,atol=1e-9*dt):
        # Binary doubling produces all grid states with O(log N) matrix products.
        state=np.empty((len(t),n+1));state[0]=init;step=expm(dt*M);k=1
        while k<len(t):
            count=min(k,len(t)-k);state[k:k+count]=state[:count]@step.T
            k*=2
            if k<len(t):step=step@step
        return state[:,:n]
    return np.stack([(expm((s-t[0])*M)@init)[:n] for s in t])

def polynomial_rk4(field,time,initial,library,*,max_step,max_steps=10000000,limits=None):
    """Call the independent side-car RK4 kernel. Return only accepted output rows."""
    t=np.ascontiguousarray(clock(time,1));x=np.ascontiguousarray(initial,dtype=float)
    powers,coef=polynomial_arrays(field);n=coef.shape[0]
    if x.shape!=(n,) or not np.isfinite(x).all() or not isinstance(max_steps,int) or max_steps<1:
        raise ValueError('Invalid initial state or step budget')
    bound=None if limits is None else np.ascontiguousarray(np.broadcast_to(limits,(n,)),dtype=float)
    lib=ct.CDLL(str(Path(library).resolve()))
    lib.mc_poly_new_v1.argtypes=[ct.c_size_t]*3+[ct.POINTER(ct.c_int),_D];lib.mc_poly_new_v1.restype=ct.c_void_p
    lib.mc_poly_drop_v1.argtypes=[ct.c_void_p];lib.mc_poly_drop_v1.restype=None
    lib.mc_rk4_v1.argtypes=[ct.c_void_p,ct.c_void_p,ct.c_size_t,_D,_D,ct.c_size_t,ct.c_double,
                          ct.c_size_t,_D,_D,ct.POINTER(ct.c_size_t),ct.POINTER(ct.c_size_t)]
    lib.mc_rk4_v1.restype=ct.c_int
    handle=lib.mc_poly_new_v1(n,n,len(powers),powers.ctypes.data_as(ct.POINTER(ct.c_int)),_ptr(coef))
    if not handle:raise ValueError('Invalid polynomial or allocation failure')
    try:
        out=np.empty((len(t),n));written=ct.c_size_t();calls=ct.c_size_t()
        status=lib.mc_rk4_v1(ct.cast(lib.mc_poly_eval_v1,ct.c_void_p),handle,n,_ptr(x),_ptr(t),len(t),
                             max_step,max_steps,None if bound is None else _ptr(bound),_ptr(out),ct.byref(written),ct.byref(calls))
    finally:lib.mc_poly_drop_v1(handle)
    if status>=4:raise ValueError('Invalid fixed-step inputs or native allocation failure')
    return dict(time=t[:written.value],b=out[:written.value],complete=status==0,
                evaluations=calls.value,message=['complete','state guard','nonfinite field','step budget'][status])
