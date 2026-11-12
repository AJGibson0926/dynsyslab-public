"""Sparse polynomial vector fields in two-dimensional local tangent charts of a modal state."""
import ctypes as ct
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .models import records_checked

def quadratic(q):
    """Return the six constant, linear and quadratic monomials on a two-dimensional tangent chart."""
    a,b=q[...,0],q[...,1]
    return np.stack([np.ones_like(a),a,b,a*a,a*b,b*b],axis=-1)

class LocalField:
    """Blend local sparse quadratic derivatives in tangent charts with optional transverse damping."""
    def __init__(self,arrays,decay=20.,neighbors=14,library=None):
        self.arrays={k:np.ascontiguousarray(arrays[k],dtype=float) for k in ('center','mean','P','Q','coef','scale')}
        self.decay,self.neighbors=float(decay),int(neighbors)
        center=self.arrays['center'];count,n=center.shape
        expected={'mean':(count,n),'P':(count,n,2),'Q':(count,n,2),'coef':(count,6,n),'scale':(n,)}
        if any(self.arrays[k].shape!=shape for k,shape in expected.items()) or any(not np.isfinite(v).all() for v in self.arrays.values()):
            raise ValueError('Invalid local-field arrays')
        if n<2 or not np.all(self.arrays['scale']>0) or not 1<=self.neighbors<count or not np.isfinite(self.decay) or self.decay<0:
            raise ValueError('Invalid local-field dimension, scales, support size or transverse damping')
        self.tree=cKDTree(center)
        self.handle=None
        if library is not None and n<=64 and self.neighbors<=64:
            lib=ct.CDLL(str(Path(library).resolve()));d=ct.POINTER(ct.c_double)
            lib.mc_local_field_v1.argtypes=[ct.c_int]*3+[ct.c_double]+[d]*6
            lib.mc_local_field_v1.restype=ct.c_void_p
            lib.mc_local_drop_v1.argtypes=[ct.c_void_p];lib.mc_local_drop_v1.restype=None
            lib.mc_local_batch_v1.argtypes=[ct.c_void_p,ct.c_int,d,d];lib.mc_local_batch_v1.restype=None
            self.ptr=lambda a:a.ctypes.data_as(d)
            self.lib=lib
            self.handle=lib.mc_local_field_v1(n,count,self.neighbors,self.decay,*[self.ptr(a) for a in self.arrays.values()])
            if not self.handle:raise ValueError('Native local-field initialization failed')
    def __del__(self):
        if getattr(self,'handle',None):self.lib.mc_local_drop_v1(self.handle)
    def __call__(self,b):
        b=np.asarray(b,dtype=float)
        if b.ndim not in (1,2) or b.shape[-1]!=len(self.arrays['scale']) or not np.isfinite(b).all():
            raise ValueError('Invalid local-field query')
        if self.handle:
            x=np.ascontiguousarray(b);y=np.empty_like(x)
            self.lib.mc_local_batch_v1(self.handle,len(np.atleast_2d(x)),self.ptr(x),self.ptr(y))
            return y
        if b.ndim==2 and len(b)>512:
            return np.vstack([self(b[i:i+512]) for i in range(0,len(b),512)])
        a=self.arrays;x=b/a['scale']
        d,ids=self.tree.query(x,k=self.neighbors+1);idx=ids[...,:-1]
        ratio=np.minimum(d[...,:-1]/np.maximum(d[...,-1:],1e-12),1.)
        w=np.where(np.all(ratio>=1.-1e-6,axis=-1,keepdims=True),1.,(1-ratio)**4*(1+4*ratio))
        w/=w.sum(axis=-1,keepdims=True)
        delta=np.expand_dims(x,-2)-a['mean'][idx]
        coords=np.einsum('...i,...ij->...j',delta,a['P'][idx])
        velocity=np.einsum('...i,...ij->...j',quadratic(coords),a['coef'][idx])
        q=a['Q'][idx]
        tangents=np.einsum('...i,...ij->...j',delta,q)
        velocity-=self.decay*(delta-np.einsum('...j,...ij->...i',tangents,q))
        return np.einsum('...i,...ij->...j',w,velocity)*a['scale']

    def save(self,path):
        with open(path,'xb') as f:np.savez_compressed(f,schema='dynsyslab.miscellaneous-computes/local-field/v1',
            decay=self.decay,neighbors=self.neighbors,**self.arrays)
    @classmethod
    def load(cls,path,library=None):
        with np.load(path) as z:
            if str(z['schema'])!='dynsyslab.miscellaneous-computes/local-field/v1':raise ValueError('Unsupported field schema')
            return cls(z,float(z['decay']),int(z['neighbors']),library)
    @classmethod
    def fit(cls,records,*,stride=36,samples=160,neighbors=14,decay=20.,ridge=.0002,threshold=.006,library=None):
        recs=records_checked(records)
        b=np.vstack([r['b'] for r in recs]);v=np.vstack([r['db'] for r in recs])
        if b.shape[1]<2 or stride<1 or not 6<=samples<=len(b) or ridge<=0 or threshold<0:
            raise ValueError('Invalid local polynomial fit options or insufficient samples')
        scale=np.sqrt(np.mean(b*b,axis=0))
        if np.any(scale<1e-14):raise ValueError('Cannot normalize a zero-energy coordinate')
        x,y=b/scale,v/scale;centers=np.concatenate([r['b'][::stride]/scale for r in recs])
        distances,indices=cKDTree(x).query(centers,k=samples)
        means=[];maps=[];units=[];coefficients=[]
        for d,idx in zip(distances,indices):
            w=np.exp(-2.5*(d/max(d[-1],1e-12))**2);w/=w.sum()
            mean=w@x[idx];delta=x[idx]-mean
            ev,vec=np.linalg.eigh((w[:,None]*delta).T@delta)
            unit=vec[:,-2:];P=unit/np.sqrt(np.maximum(ev[-2:],1e-7))
            f=quadratic(delta@P);G=(w[:,None]*f).T@f
            G.flat[::7]+=ridge*np.maximum(np.diag(G),1e-10)
            rhs=f.T@(w[:,None]*y[idx]);coef=np.linalg.solve(G,rhs)
            keep=abs(coef)*np.sqrt(np.diag(G))[:,None]>threshold*np.sqrt(w@(y[idx]**2))
            keep[0]=True
            for active in np.unique(keep.T,axis=0):
                rows=np.flatnonzero(active);cols=np.flatnonzero(np.all(keep.T==active,axis=1))
                coef[:,cols]=0.
                coef[np.ix_(rows,cols)]=np.linalg.solve(G[np.ix_(rows,rows)],rhs[np.ix_(rows,cols)])
            means.append(mean);maps.append(P);units.append(unit);coefficients.append(coef)
        return cls(dict(center=centers,mean=means,P=maps,Q=units,coef=coefficients,scale=scale),decay,neighbors,library)
