"""Small, verified PCA/POD family based solely on NumPy linear algebra."""
from __future__ import annotations
import numpy as np
from collections.abc import Callable
from time_series import TimeSeries
from .result import HAVOKResult, PCAResult, SVDResult
from .modes import pod_modes_from_hankel_svd, pod_modes_from_svd

def _matrix(values):
    if isinstance(values,TimeSeries): values=values.values
    x=np.asarray(values,dtype=float)
    return x[:,None] if x.ndim==1 else x.reshape(x.shape[0],-1)

def pca(values, rank: int|None=None, center: bool=True) -> PCAResult:
    x=_matrix(values); mean=x.mean(axis=0) if center else np.zeros(x.shape[1]); centered=x-mean
    U,s,Vh=np.linalg.svd(centered,full_matrices=False); r=min(rank or len(s),len(s)); total=np.sum(s*s)
    return PCAResult("pca",mean,Vh[:r],s[:r],U[:,:r]*s[:r],(s[:r]**2/total if total else np.zeros(r)),{"rank":r,"center":center})

def snapshot_pod(values, rank: int|None=None, center: bool=True) -> PCAResult:
    """Snapshot POD: spatial modes are components and temporal coefficients are scores."""
    result=pca(values,rank=rank,center=center); result.method="snapshot-pod"; return result

def _soft_threshold(x,t): return np.sign(x)*np.maximum(np.abs(x)-t,0.)
def _svt(x,t):
    u,s,vh=np.linalg.svd(x,full_matrices=False); keep=s>t
    return (u[:,keep]*(s[keep]-t))@vh[keep]

def robust_pca(values, lam: float|None=None, tolerance: float=1e-7, max_iter: int=1000):
    """Principal-component pursuit via inexact augmented-Lagrange iterations.

    Returns low-rank and sparse components; convergence is recorded in metadata.
    """
    x=_matrix(values); lam=lam or 1/np.sqrt(max(x.shape)); norm=np.linalg.norm(x,"fro") or 1.
    low,sparse,dual=np.zeros_like(x),np.zeros_like(x),np.zeros_like(x); mu=1.25/(np.linalg.norm(x,2) or 1.)
    for iteration in range(1,max_iter+1):
        low=_svt(x-sparse+dual/mu,1/mu); sparse=_soft_threshold(x-low+dual/mu,lam/mu)
        residual=x-low-sparse; error=np.linalg.norm(residual,"fro")/norm; dual+=mu*residual
        if error<tolerance: break
    return {"method":"robust-pca","low_rank":low,"sparse":sparse,"iterations":iteration,"relative_residual":error,"lambda":lam}

def bagged_pca(values, rank: int|None=None, bags: int=32, seed: int|None=None) -> PCAResult:
    """Bootstrap aggregation of PCA bases, with component sign alignment."""
    x=_matrix(values); rng=np.random.default_rng(seed); reference=pca(x,rank=rank); components=[]; ratios=[]
    for _ in range(bags):
        sample=x[rng.integers(x.shape[0],size=x.shape[0])]; trial=pca(sample,rank=reference.components.shape[0])
        aligned=trial.components.copy()
        for index in range(aligned.shape[0]):
            if np.dot(aligned[index],reference.components[index])<0: aligned[index]*=-1
        components.append(aligned); ratios.append(trial.explained_variance_ratio)
    mean_components=np.mean(components,axis=0); mean_components/=np.linalg.norm(mean_components,axis=1,keepdims=True)
    scores=(x-reference.mean)@mean_components.T
    return PCAResult("bagged-pca",reference.mean,mean_components,reference.singular_values,scores,np.mean(ratios,axis=0),{"rank":mean_components.shape[0],"bags":bags,"seed":seed})

def svd(values: np.ndarray|TimeSeries, *, center: bool=False, full_matrices: bool=False,
        source_metadata: dict|None=None, source_archive: str|None=None) -> SVDResult:
    """Compute direct SVD; unlike PCA, input is uncentered by default."""
    if isinstance(values,TimeSeries) and values.is_pooled: raise ValueError("Regular SVD requires one individual TimeSeries; choose a pooling convention explicitly.")
    x=_matrix(values)
    if not x.size or not np.all(np.isfinite(x)): raise ValueError("SVD requires a finite, nonempty samples-by-features array.")
    mean=x.mean(axis=0) if center else np.zeros(x.shape[1],dtype=float)
    left,singular,right=np.linalg.svd(x-mean,full_matrices=full_matrices)
    metadata={"method":"svd","centered":bool(center),"full_matrices":bool(full_matrices),"samples":int(x.shape[0]),"features":int(x.shape[1])}
    modes=pod_modes_from_svd(left,singular,right,x,mean=mean,source_metadata=source_metadata,source_archive=source_archive,metadata=metadata)
    metadata["pod_modes"]=modes.summary()
    return SVDResult(left,singular,right,mean,metadata,modes)

def hankel_pod(values: np.ndarray|TimeSeries, *, matrix_size: int=10, rank: int|None=None,
               delay_stride: int=1, center: bool=False, source_metadata: dict|None=None,
               source_archive: str|None=None) -> SVDResult:
    """POD of a Hankel time-delay observable matrix.

    The samples-by-delayed-observables layout matches the delay-coordinate
    construction used by Arbabi and Mezić for Hankel-DMD.  A requested rank
    uses an economical truncated SVD, which is sufficient when only leading
    POD modes are required.
    """
    if isinstance(values,TimeSeries):
        if values.is_pooled: raise ValueError("Hankel POD requires an individual trajectory; pooled data lose delay boundaries.")
        x=values.values
    else: x=np.asarray(values,dtype=float)
    x=_matrix(x)
    if int(matrix_size)!=matrix_size or matrix_size<2: raise ValueError("matrix_size must be an integer of at least 2.")
    if int(delay_stride)!=delay_stride or delay_stride<1: raise ValueError("delay_stride must be a positive integer.")
    delays=int(matrix_size); stride=int(delay_stride); last=(delays-1)*stride
    if x.shape[0]<=last+2: raise ValueError("Hankel POD delay settings leave fewer than three embedded samples.")
    hankel=np.hstack([x[index*stride:x.shape[0]-last+index*stride] for index in range(delays)])
    mean=hankel.mean(axis=0) if center else np.zeros(hankel.shape[1],dtype=float); centered=hankel-mean
    maximum=min(centered.shape)
    if rank is not None and (int(rank)!=rank or rank<1): raise ValueError("rank must be a positive integer when specified.")
    requested=maximum if rank is None else min(int(rank),maximum)
    if requested<maximum:
        from scipy.sparse.linalg import svds
        left,singular,right=np.linalg.svd(centered,full_matrices=False) if requested==maximum else svds(centered,k=requested,which="LM")
        order=np.argsort(singular)[::-1]; left,singular,right=left[:,order],singular[order],right[order,:]
    else: left,singular,right=np.linalg.svd(centered,full_matrices=False)
    metadata={"method":"hankel-pod","centered":bool(center),"rank":int(singular.size),"delays":delays,
              "delay_stride_samples":stride,"state_dimension":int(x.shape[1]),"embedded_feature_dimension":int(hankel.shape[1]),
              "embedded_samples":int(hankel.shape[0]),"reference":"Arbabi & Mezić (2017), doi:10.1137/17M1125236; SVD/POD of a Hankel time-delay matrix."}
    modes=pod_modes_from_svd(left,singular,right,x[:hankel.shape[0]],mean=mean,source_metadata=source_metadata,source_archive=source_archive,metadata=metadata)
    metadata["pod_modes"]=modes.summary()
    return SVDResult(left,singular,right,mean,metadata,modes)

def havok(values: np.ndarray|TimeSeries, *, time: np.ndarray|None=None, delays: int=100, rank: int=15,
          observable: int|Callable[[np.ndarray],np.ndarray]=0, derivative_order: int=2,
          source_metadata: dict|None=None, source_archive: str|None=None) -> HAVOKResult:
    """Fit the HAVOK forced-linear delay-coordinate representation."""
    if isinstance(values,TimeSeries):
        if time is not None: raise ValueError("Pass either a TimeSeries or separate data and time, not both.")
        if values.is_pooled: raise ValueError("HAVOK requires an individual trajectory; pooled data lose delay boundaries.")
        x,timestamps=values.values,values.time
    else:
        x=np.asarray(values,dtype=float)
        if time is None: raise ValueError("Array HAVOK input requires an explicit time array.")
        timestamps=np.asarray(time,dtype=float)
    if x.ndim==1: x=x[:,None]
    timestamps=np.asarray(timestamps,dtype=float).reshape(-1)
    if x.ndim!=2 or x.shape[0]!=timestamps.size or x.shape[0]<4: raise ValueError("HAVOK requires at least four timestamped state samples.")
    steps=np.diff(timestamps)
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(timestamps)) or np.any(steps<=0): raise ValueError("HAVOK requires finite values and strictly increasing timestamps.")
    # Explicit uniform grids such as np.arange(...)*dt retain harmless binary
    # roundoff at long horizons; accept it while still rejecting true adaptive
    # solver grids.
    if not np.allclose(steps,steps[0],rtol=1e-10,atol=1e-13): raise ValueError("HAVOK requires uniform time steps. Resample explicitly before processing nonuniform data.")
    signal=x[:,observable] if isinstance(observable,int) else np.asarray(observable(x),dtype=float).reshape(-1)
    if signal.size!=timestamps.size or not np.all(np.isfinite(signal)): raise ValueError("The HAVOK observable must return one finite scalar per time sample.")
    if delays<2 or delays>=signal.size-1: raise ValueError("HAVOK delays must be at least 2 and leave at least two Hankel columns.")
    hankel=np.lib.stride_tricks.sliding_window_view(signal,int(delays)).T.copy(); left,singular,right=np.linalg.svd(hankel,full_matrices=False)
    usable=min(int(rank),singular.size)
    if usable<2: raise ValueError("HAVOK needs rank at least 2 after Hankel truncation.")
    coordinates=right[:usable].T; coordinate_time=timestamps[:coordinates.shape[0]]
    derivative=np.gradient(coordinates,coordinate_time,axis=0,edge_order=derivative_order)
    operator=np.linalg.lstsq(np.column_stack((coordinates[:-1,:-1],coordinates[:-1,-1])),derivative[:-1,:-1],rcond=None)[0]
    observable_metadata={"kind":"state_index","index":int(observable)} if isinstance(observable,int) else {"kind":"callable"}
    metadata={"method":"havok","delays":int(delays),"rank":usable,"delta_t":float(steps[0]),"observable":observable_metadata,"derivative_order":int(derivative_order),"state_dimension":int(x.shape[1]),"hankel_shape":list(hankel.shape),"reference":"Brunton et al., Chaos as an intermittently forced linear system (2017), doi:10.1038/s41467-017-00030-8"}
    prom=pod_modes_from_hankel_svd(left[:,:usable],singular[:usable],right[:usable],x[:coordinates.shape[0]],mean=np.zeros(left.shape[0]),source_metadata=source_metadata,source_archive=source_archive,metadata=metadata)
    metadata["prom"]=prom.summary()
    return HAVOKResult(operator[:-1,:].T,operator[-1:,:].T,coordinates[:,-1],coordinates,coordinate_time,left[:,:usable],singular[:usable],right[:usable],signal,metadata,prom)
