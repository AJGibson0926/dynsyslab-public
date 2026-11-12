"""Empirical planar crossing witnesses and independent branch refinement."""
import ctypes as ct
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .sampling import clock

class _Hit(ct.Structure):
    _fields_=[('i',ct.c_size_t),('j',ct.c_size_t),('u',ct.c_double),('v',ct.c_double)]

def candidates(time,plane,library,*,min_separation=0.,limit=100000):
    """Return visit brackets/times and segment intersection points on the supplied path."""
    t=clock(time);x=np.ascontiguousarray(plane,dtype=float)
    if (x.shape!=(len(t),2) or not np.isfinite(x).all() or min_separation<0 or
        not isinstance(limit,int) or limit<1):raise ValueError('Invalid planar path or witness limit')
    lib=ct.CDLL(str(Path(library).resolve()));fn=lib.mc_crossings_v1
    ptr=ct.POINTER(ct.c_double);fn.argtypes=[ptr,ct.c_size_t,ct.POINTER(_Hit),ct.c_size_t,ct.POINTER(ct.c_size_t)]
    fn.restype=ct.c_int;found=ct.c_size_t()
    if fn(x.ctypes.data_as(ptr),len(t),None,0,ct.byref(found)):raise ValueError('Native crossing search failed')
    if found.value>limit:raise ValueError(f'{found.value} candidates exceed limit {limit}')
    hits=(_Hit*found.value)()
    if fn(x.ctypes.data_as(ptr),len(t),hits,len(hits),ct.byref(found)):raise RuntimeError('Crossing enumeration failed')
    results=[]
    for h in hits:
        i,j=h.i,h.j
        times=[t[i]+h.u*(t[i+1]-t[i]),t[j]+h.v*(t[j+1]-t[j])]
        if times[1]-times[0]<min_separation:continue
        results.append(dict(segments=[i,j],times=times,brackets=[t[i:i+2].tolist(),t[j:j+2].tolist()],
                            point=(x[i]+h.u*(x[i+1]-x[i])).tolist()))
    return sorted(results,key=lambda h:h['times'])

def refine(witness,observe,*,velocity=None,lift=None,scales=None,tolerance=1e-9):
    """Refine two distinct visits on an independently supplied continuous trajectory.
    Optional velocity and lift are evaluated at both visits, in their original units.
    """
    box=np.asarray(witness['brackets'],float);guess=np.asarray(witness['times'],float)
    s=np.ones(2) if scales is None else np.asarray(scales,float)
    if (box.shape!=(2,2) or not np.isfinite(box).all() or np.any(box[:,0]>=box[:,1]) or
        box[0,1]>box[1,0] or s.shape!=(2,) or not np.isfinite(s).all() or np.any(s<=0) or tolerance<=0):
        raise ValueError('Two ordered, distinct visit brackets and positive coordinate scales required')
    def pair(t):return np.asarray([observe(float(v)) for v in t],float)
    def residual(t):
        obs=pair(t)
        if obs.shape!=(2,2):raise ValueError('Expected a continuous two-coordinate observation')
        return (obs[0]-obs[1])/s
    def jac(t):return np.column_stack([velocity(float(t[0])),-np.asarray(velocity(float(t[1])) )])/s[:,None]
    fit=least_squares(residual,np.clip(guess,box[:,0],box[:,1]),bounds=(box[:,0],box[:,1]),
                      jac=jac if velocity is not None else '2-point',gtol=1e-14,ftol=1e-14,xtol=1e-14)
    t=fit.x;gap=float(np.linalg.norm(residual(t)))
    result=dict(times=t.tolist(),values=pair(t).tolist(),scaled_gap=gap,
                matched=bool(fit.success and gap<=tolerance and t[0]<t[1]))
    if velocity is not None:
        v=np.array([velocity(float(a)) for a in t]);norms=np.linalg.norm(v/s,axis=1)
        result['velocities']=v.tolist();result['velocity_difference']=(v[1]-v[0]).tolist()
        result['acute_angle_degrees']=float(np.degrees(np.arccos(np.clip(abs((v[0]/s)@(v[1]/s))/np.prod(norms),0,1)))) if np.all(norms>0) else None
    if lift is not None:
        d=np.atleast_1d(lift(float(t[1])))-np.atleast_1d(lift(float(t[0])))
        result['lift_difference']=d.tolist();result['cumulative_distance']=np.sqrt(np.cumsum(d*d)).tolist()
    return result
