"""Frozen delay bases, exact covariance extension, and causal projection on independent records."""
from dataclasses import dataclass
import ctypes
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.signal import oaconvolve
from scipy.sparse.linalg import LinearOperator,eigsh

def matrix(a):
    """Convert finite scalar or multichannel samples into a nonempty samples-by-channels matrix."""
    a=np.asarray(a,dtype=float)
    if a.ndim==1:
        a=a[:,None]
    if a.ndim!=2 or not a.size or not np.isfinite(a).all():
        raise ValueError("Expected a finite, nonempty sample-by-channel matrix")
    return a

def scalar_delay_pairs(values, ends, lags, dt, half_width=1):
    """Build two scalar taps per lag from data available at each end index.

    Raw tap columns are (older,newer) for each lag, in the supplied order.
    Return the taps, interleaved midpoint averages and secant slopes, and
    optional training labels for curvature from centered second differences.
    At zero lag the estimates refer to end-half_width, with half_width*dt latency.
    """
    y=matrix(values)
    rows,offsets=np.asarray(ends),np.asarray(lags)
    if y.shape[1]!=1 or not np.isfinite(dt) or dt<=0 or not isinstance(half_width,(int,np.integer)) or half_width<1:
        raise ValueError('Expected a scalar signal and positive sampling and half-width')
    for indices in (rows,offsets):
        if indices.ndim!=1 or not len(indices) or not np.issubdtype(indices.dtype,np.integer):
            raise ValueError('End indices and lags must be nonempty integer vectors')
    if offsets.min()<0 or rows.max()>=len(y) or int(rows.min())<int(offsets.max())+2*half_width:
        raise ValueError('Each pair must lie inside its supplied trajectory')
    latest=rows.astype(np.intp)[:,None]-offsets.astype(np.intp)
    p=y[:,0];older=p[latest-2*half_width];newer=p[latest]
    delta=dt*half_width
    raw=np.stack([older,newer],axis=-1).reshape(len(rows),-1)
    state=np.empty_like(raw)
    state[:,::2]=(newer+older)/2
    state[:,1::2]=(newer-older)/(2*delta)
    curvature=(newer-2*p[latest-half_width]+older)/delta**2
    return raw,state,curvature

def uniform_clock(t,length):
    """Validate record alignment and uniform spacing and return the clock with its step."""
    t=np.asarray(t,dtype=float)
    if t.shape!=(length,) or length<2 or not np.isfinite(t).all():
        raise ValueError('Expected a finite clock aligned with the record')
    dt=(t[-1]-t[0])/(length-1)
    if dt<=0 or not np.allclose(np.diff(t),dt,atol=abs(dt)*1e-8,rtol=1e-7):
        raise ValueError('A uniform, strictly increasing clock is required')
    return t,dt

@dataclass
class HistoryBasis:
    """Frozen signed Hankel modes, singular values, scales and sampling convention for independent trajectories."""
    U: np.ndarray
    sigma: np.ndarray
    scales: np.ndarray
    energy: np.ndarray
    windows: int
    tap_interval: float
    gram_error: float=0.

    @property
    def taps(self):
        return len(self.U)//len(self.scales)

    def window_encoder(self, rank=None):
        """Return oldest-first tap × channel × mode weights in physical input units for b=√N ξ."""
        rank = len(self.sigma) if rank is None else rank
        if not isinstance(rank, (int, np.integer)) or not 1 <= rank <= len(self.sigma):
            raise ValueError('Invalid retained rank')
        vectors = self.U[:, :rank].reshape(self.taps, len(self.scales), rank)
        return vectors * (np.sqrt(self.windows) / self.sigma[:rank]) / self.scales[None, :, None]

    def project_windows(self, windows, rank=None):
        """Project one or more complete physical windows; last axes are oldest-first taps and channels."""
        x = np.asarray(windows, dtype=float)
        encoder = self.window_encoder(rank)
        if x.shape[-2:] != encoder.shape[:2] or not np.isfinite(x).all():
            raise ValueError('Expected finite complete windows with matching tap and channel axes')
        return x.reshape(*x.shape[:-2], len(self.U)) @ encoder.reshape(len(self.U), -1)

    def project(self,t,observations,rates=None,rank=None):
        """Unit training RMS coordinates b at valid endpoints. Rates follow the same linear map."""
        y=matrix(observations);t,dt=uniform_clock(t,len(y))
        rank=len(self.sigma) if rank is None else rank
        gap=round(self.tap_interval/dt)
        if not 1<=rank<=len(self.sigma) or y.shape[1]!=len(self.scales):
            raise ValueError('Invalid rank or channel count')
        if gap<1 or not np.isclose(gap*dt,self.tap_interval,rtol=1e-7,atol=self.tap_interval*1e-9):
            raise ValueError('The tap interval must be an integer multiple of the sample interval')
        first=(self.taps-1)*gap
        if len(y)<=first:
            raise ValueError('Record has no complete window')
        kernel=self.window_encoder(rank)
        def apply(samples):
            samples=matrix(samples)
            if samples.shape!=y.shape:
                raise ValueError('Observation values and rates must align')
            out=np.zeros((len(y)-first,rank))
            for phase in range(min(gap,len(out))):
                seq=samples[phase::gap]
                for j in range(rank):
                    for c in range(seq.shape[1]):
                        out[phase::gap,j]+=oaconvolve(seq[:,c],kernel[::-1,c,j],mode='valid')
            return out
        return t[first:],apply(y),None if rates is None else apply(rates)

    def endpoint(self,rank,channel=0):
        """Linear decoder in observation units. Apply to b and bdot for a consistent phase pair."""
        if not 1<=rank<=len(self.sigma) or not 0<=channel<len(self.scales):
            raise ValueError('Invalid rank or output channel')
        return self.U[-len(self.scales)+channel,:rank]*self.sigma[:rank]*self.scales[channel]/np.sqrt(self.windows)

    def raw_scalar_xi(self,b):
        if len(self.scales)!=1:
            raise ValueError('The original raw scalar convention requires one channel')
        return np.asarray(b)*self.scales[0]/np.sqrt(self.windows)

    def save(self,path):
        with Path(path).open('xb') as stream:
            np.savez(stream,schema='dynsyslab.miscellaneous-computes/history-basis/v1',**vars(self))

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:
            if str(z['schema'])!='dynsyslab.miscellaneous-computes/history-basis/v1':
                raise ValueError('Unsupported basis schema')
            return cls(z['U'],z['sigma'],z['scales'],z['energy'],int(z['windows']),float(z['tap_interval']),float(z['gram_error']))

def fit_basis(records,taps,rank,library,*,gap=1,scales=None,locked=None):
    """Pool complete windows per (clock, observations). A locked (U,sigma) prefix stays bit-identical."""
    if not records or min(taps,rank,gap)<1:
        raise ValueError('Positive taps, rank, gap and at least one training record are required')
    series=[matrix(values) for _,values in records]
    channels=series[0].shape[1]
    dts=[uniform_clock(t,len(x))[1] for (t,_),x in zip(records,series)]
    if any(x.shape[1]!=channels for x in series) or not np.allclose(dts,dts[0],atol=1e-10*dts[0]):
        raise ValueError('All basis-training records must share channels and sample interval')
    if scales is None:
        scales=np.sqrt(sum(np.einsum('ij,ij->j',x,x) for x in series)/sum(map(len,series)))
    scales=np.atleast_1d(np.asarray(scales,dtype=float))
    if scales.shape!=(channels,) or not np.isfinite(scales).all() or np.any(scales<=0):
        raise ValueError('Positive finite per-channel scales are required')
    dim=taps*channels
    if rank>dim:
        raise ValueError('Requested rank exceeds the number of delay coordinates')
    gram=np.zeros((dim,dim));windows=0
    native=ctypes.CDLL(str(Path(library).resolve())).mc_hankel_gram_v1
    ptr=np.ctypeslib.ndpointer(dtype='float64',flags='C_CONTIGUOUS')
    native.argtypes=[ptr,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_size_t,ptr];native.restype=ctypes.c_int
    for x in series:
        if len(x)<1+(taps-1)*gap:
            raise ValueError('Each record must supply a complete window')
        for phase in range(min(gap,len(x))):
            v=np.ascontiguousarray(x[phase::gap]/scales)
            if len(v)<taps:continue
            if native(v,len(v),channels,taps,gram):
                raise ValueError('Native covariance kernel rejected its inputs')
            windows+=len(v)-taps+1
    U=np.empty((dim,0));s=np.empty(0)
    if locked is not None:
        U,s=map(np.asarray,locked)
        if U.shape!=(dim,len(s)) or len(s)>rank or np.any(s<=0) or not np.allclose(U.T@U,np.eye(len(s)),atol=1e-8):
            raise ValueError('Invalid frozen prefix')
        if np.max(np.linalg.norm(gram@U-U*s**2,axis=0)/s**2,initial=0)>1e-6:
            raise ValueError('Frozen eigenpairs do not match these training records, scales and taps')
    extra=rank-len(s)
    if extra:
        def project(v):
            return v-U@(U.T@v)
        if dim<160 or extra>=dim-2:
            # The complement projection preserves locked modes even if eigenvalues are repeated.
            p=project(np.eye(dim))
            ev,vec=eigh(p@gram@p,subset_by_index=[dim-extra,dim-1])
        else:
            op=LinearOperator((dim,dim),matvec=lambda v:project(gram@project(v)),dtype=float)
            ev,vec=eigsh(op,k=extra,which='LA',tol=1e-10,v0=np.random.default_rng(73).normal(size=dim))
        ev=ev[::-1];vec=vec[:,::-1]
        if ev[-1]<=max(1e-14,np.finfo(float).eps*dim)*max(ev[0],s[0]**2 if len(s) else 0):
            raise ValueError('Requested rank exceeds numerical rank')
        vec*=np.sign(vec[np.argmax(abs(vec),axis=0),np.arange(extra)])
        U=np.column_stack([U,vec]);s=np.r_[s,np.sqrt(ev)]
    product=(U.T@gram@U)/s[:,None]/s
    return HistoryBasis(U,s,scales,np.cumsum(s*s)/np.trace(gram),windows,dts[0]*gap,float(np.max(abs(product-np.eye(rank)))))

def separation(time,coordinates,first,second):
    """Per-mode signed differences and cumulative Euclidean distance at two specified times."""
    from scipy.interpolate import CubicSpline
    x=matrix(coordinates);t=np.asarray(time,dtype=float)
    if len(t)!=len(x) or not np.isfinite(t).all() or np.any(np.diff(t)<=0) or not t[0]<=min(first,second)<=max(first,second)<=t[-1]:
        raise ValueError('Visits must lie on the supplied increasing clock')
    a,b=CubicSpline(t,x)([first,second])
    return a-b,np.sqrt(np.cumsum((a-b)**2))
