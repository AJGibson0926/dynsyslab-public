"""Compact support radial dictionaries and global sparse regression for arbitrary input/output dimensions."""
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from .history import matrix
from .models import records_checked
from .sampling import trapezoid_weights
from .sparse_stlsq import sparse_stlsq

class CompactFeatures:
    """Build sparse neighbor-indexed radial features from fixed training centers in the declared coordinates."""
    def __init__(self, centers, neighbors, asinh_scale=0.):
        self.centers=matrix(centers)
        self.k=int(neighbors)
        self.scale=float(asinh_scale)
        if self.k!=neighbors or not 1<=self.k<len(self.centers) or not np.isfinite(self.scale) or self.scale<0:
            raise ValueError("Invalid compact support or input scale")
        self.tree=cKDTree(self.centers)

    def __call__(self, values):
        x=matrix(values)
        z=np.arcsinh(x/self.scale) if self.scale else x
        distances,indices=self.tree.query(z,self.k+1)
        r=distances[:,:-1]/np.maximum(distances[:,-1:],1e-12)
        w=(1-r)**4*(1+4*r)
        w[np.all(r>=1.-1e-6,axis=1)]=1.
        w/=w.sum(axis=1,keepdims=True)
        return sparse.csr_matrix((w.ravel(),indices[:,:-1].ravel(),
            self.k*np.arange(len(x)+1)),shape=(len(x),len(self.centers)))

class CompactField:
    """Fit or replay a vector-valued compact radial map, optionally adding an affine term."""
    def __init__(self, centers, coefficients, linear, neighbors, asinh_scale=0.):
        self.features=CompactFeatures(centers,neighbors,asinh_scale)
        self.coefficients=matrix(coefficients)
        self.linear=np.asarray(linear,dtype=float)
        if (len(self.coefficients)!=len(self.features.centers) or
            self.linear.shape!=(self.features.centers.shape[1],self.coefficients.shape[1]) or
            not np.isfinite(self.linear).all()):
            raise ValueError("Incompatible compact-field arrays")

    def __call__(self, values):
        x=np.atleast_2d(values)
        y=self.features(x)@self.coefficients+x@self.linear
        return y[0] if np.ndim(values)==1 else y

    @classmethod
    def fit(cls, records, *, stride=40, neighbors=10, asinh_scale=0., affine_tail=False,
            ridge=.01, threshold=1e-5, time_weighted=False, speed_scale=None, block_size=2048):
        return cls.from_labels([(r['time'],r['b'],r['db']) for r in records_checked(records)],
            stride=stride,neighbors=neighbors,asinh_scale=asinh_scale,affine_tail=affine_tail,
            ridge=ridge,threshold=threshold,time_weighted=time_weighted,speed_scale=speed_scale,block_size=block_size)

    @classmethod
    def from_labels(cls, records, *, stride, neighbors, asinh_scale=0., affine_tail=False,
                    ridge=.01, threshold=1e-5, time_weighted=False, speed_scale=None, block_size=2048):
        from .sampling import clock
        if not isinstance(stride,int) or stride<1 or block_size<1 or not records:
            raise ValueError("Invalid center or block stride")
        if speed_scale is not None and (not np.isfinite(speed_scale) or speed_scale<=0):
            raise ValueError("Invalid target speed scale")
        data=[(clock(t),matrix(x),matrix(y)) for t,x,y in records]
        if any(len(t)!=len(x) or len(x)!=len(y) or x.shape[1]!=data[0][1].shape[1]
               or y.shape[1]!=data[0][2].shape[1] for t,x,y in data):
            raise ValueError("Label shapes or clocks disagree")
        raw=np.unique(np.concatenate([x[::stride] for _,x,_ in data]),axis=0)
        centers=np.arcsinh(raw/asinh_scale) if asinh_scale else raw
        features=CompactFeatures(centers,neighbors,asinh_scale)
        gram=None;rhs=None;energy=np.zeros(data[0][2].shape[1])
        for t,x,y in data:
            weight=np.sqrt(trapezoid_weights(t)) if time_weighted else np.ones(len(t))
            if speed_scale is not None:weight/=np.hypot(1.,np.linalg.norm(y,axis=1)/speed_scale)
            for i in range(0,len(t),block_size):
                sl=slice(i,i+block_size);a=features(x[sl])
                if affine_tail:a=sparse.hstack([a,sparse.csr_matrix(x[sl]-a@raw)],format='csr')
                a=a.multiply(weight[sl,None]).tocsr();target=y[sl]*weight[sl,None]
                g=a.T@a;r=a.T@target;energy+=(target**2).sum(axis=0)
                gram=g if gram is None else gram+g
                rhs=r if rhs is None else rhs+r
        coef=sparse_stlsq(gram,rhs,np.sqrt(energy),ridge,threshold)
        linear=coef[len(raw):] if affine_tail else np.zeros((raw.shape[1],len(energy)))
        return cls(centers,coef[:len(raw)]-raw@linear,linear,neighbors,asinh_scale)

    def save(self,path):
        with open(path,'xb') as stream:
            np.savez_compressed(stream,schema='dynsyslab.miscellaneous-computes/compact-field/v1',
                centers=self.features.centers,coefficients=self.coefficients,linear=self.linear,
                neighbors=self.features.k,asinh_scale=self.features.scale)

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:
            if str(z['schema'])!='dynsyslab.miscellaneous-computes/compact-field/v1':
                raise ValueError("Unsupported compact field")
            return cls(z['centers'],z['coefficients'],z['linear'],int(z['neighbors']),float(z['asinh_scale']))
