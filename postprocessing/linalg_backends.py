"""Dense linear-algebra backends for complete DMD factorization kernels."""
from __future__ import annotations

from contextlib import nullcontext
from dataclasses import asdict, dataclass
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from functools import lru_cache
from time import perf_counter
from typing import Any, Literal

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits


LinearAlgebraBackendPreference=Literal["auto","scipy","torch-cuda"]


@dataclass(frozen=True)
class LinearAlgebraBackend:
    name: str
    available: bool
    compiled: bool
    accelerator: str
    full_spectrum: bool
    svd_drivers: tuple[str,...]
    auto_admitted: bool
    version: str|None=None
    libraries: tuple[str,...]=()
    reason: str|None=None

    def provenance(self,requested: LinearAlgebraBackendPreference) -> dict[str,Any]:
        return {"requested":requested,**asdict(self)}


@dataclass(frozen=True)
class DMDKernelOutput:
    eigenvalues: np.ndarray
    modes: np.ndarray
    coordinates: np.ndarray|None
    operator: np.ndarray
    retained_rank: int
    timing: dict[str,float|str]
    backend: dict[str,Any]


def _blas_libraries() -> tuple[str,...]:
    names=[]
    for item in threadpool_info():
        if item.get("user_api")!="blas": continue
        identity="/".join(str(item.get(key) or "unknown") for key in ("internal_api","prefix","version"))
        if identity not in names: names.append(identity)
    return tuple(names)


@lru_cache(maxsize=1)
def scipy_backend() -> LinearAlgebraBackend:
    import scipy
    return LinearAlgebraBackend("scipy-lapack",True,True,"cpu",True,("gesvd","gesdd"),True,
                                scipy.__version__,_blas_libraries())


@lru_cache(maxsize=1)
def torch_cuda_backend() -> LinearAlgebraBackend:
    try:
        torch=import_module("torch")
        installed=version("torch")
    except (ImportError,PackageNotFoundError) as error:
        return LinearAlgebraBackend("torch-cuda",False,True,"cuda",True,("gesvd",),False,
                                    reason=f"optional PyTorch is not installed: {error}")
    try:
        available=bool(torch.cuda.is_available())
        cuda_version=None if torch.version.cuda is None else str(torch.version.cuda)
        device=None if not available else str(torch.cuda.get_device_name(torch.cuda.current_device()))
    except Exception as error:
        return LinearAlgebraBackend("torch-cuda",False,True,"cuda",True,("gesvd",),False,
                                    installed,reason=f"CUDA capability inspection failed: {error}")
    reason=None if available else "PyTorch is installed but no CUDA device/runtime is available."
    libraries=tuple(item for item in (f"cuda/{cuda_version}" if cuda_version else None,device) if item)
    return LinearAlgebraBackend("torch-cuda",available,True,"cuda",True,("gesvd",),False,
                                installed,libraries,reason)


def available_linear_algebra_backends() -> dict[str,LinearAlgebraBackend]:
    return {"scipy":scipy_backend(),"torch-cuda":torch_cuda_backend()}


def select_linear_algebra_backend(preference: LinearAlgebraBackendPreference) -> LinearAlgebraBackend:
    if preference not in {"auto","scipy","torch-cuda"}:
        raise ValueError("linear_algebra_backend must be 'auto', 'scipy', or 'torch-cuda'.")
    if preference in {"auto","scipy"}: return scipy_backend()
    candidate=torch_cuda_backend()
    if not candidate.available: raise RuntimeError(candidate.reason or "The requested CUDA linear-algebra backend is unavailable.")
    return candidate


def _retained_rank(singular_values: np.ndarray,rank: int|None,svd_tolerance: float|None) -> int:
    if rank is not None: return min(len(singular_values),int(rank))
    if svd_tolerance is not None: return int(np.sum(singular_values>float(svd_tolerance)))
    return int(np.sum(singular_values>np.finfo(float).eps*singular_values[0]*max(len(singular_values),1)))


def _scipy_kernel(X: np.ndarray,Y: np.ndarray,rank: int|None,exact: bool,svd_tolerance: float|None,
                  svd_driver: str,compute_coordinates: bool,svd_threads: int|None,
                  eigendecomposition_threads: int|None,backend: LinearAlgebraBackend,
                  requested_backend: LinearAlgebraBackendPreference) -> DMDKernelOutput:
    started=perf_counter(); requested=None if rank is None else int(rank)
    limit=nullcontext() if svd_threads is None else threadpool_limits(limits=int(svd_threads),user_api="blas")
    with limit:
        if requested is not None and 0<requested<min(X.shape):
            from scipy.sparse.linalg import svds
            U,s,Vh=svds(X,k=requested,which="LM")
            order=np.argsort(s)[::-1]; U,s,Vh=U[:,order],s[order],Vh[order,:]
        else:
            from scipy.linalg import svd
            U,s,Vh=svd(X,full_matrices=False,lapack_driver=svd_driver,check_finite=False)
    decomposed=perf_counter(); retained=_retained_rank(s,rank,svd_tolerance)
    if retained<1: raise ValueError("DMD retained no singular values; lower svd_tolerance or inspect the input data.")
    U,s,Vh=U[:,:retained],s[:retained],Vh[:retained]
    operator=((U.conj().T@Y)@Vh.conj().T)/s[None,:]
    projected_at=perf_counter()
    eig_limit=nullcontext() if eigendecomposition_threads is None else threadpool_limits(
        limits=int(eigendecomposition_threads),user_api="blas")
    with eig_limit: eigenvalues,W=np.linalg.eig(operator)
    eigendecomposed=perf_counter()
    modes=(((Y@Vh.conj().T)/s[None,:])@W) if exact else U@W
    modes_built=perf_counter()
    coordinates=(np.linalg.pinv(modes)@X).T if compute_coordinates else None
    completed=perf_counter()
    timing={"svd_seconds":decomposed-started,"projected_operator_seconds":projected_at-decomposed,
            "eigendecomposition_seconds":eigendecomposed-projected_at,"mode_construction_seconds":modes_built-eigendecomposed,
            "coordinate_recovery_seconds":completed-modes_built,"factorization_seconds":completed-started,
            "clock":"time.perf_counter"}
    return DMDKernelOutput(eigenvalues,modes,coordinates,operator,retained,timing,
                           backend.provenance(requested_backend))


def _torch_cuda_kernel(X: np.ndarray,Y: np.ndarray,rank: int|None,exact: bool,svd_tolerance: float|None,
                       svd_driver: str,compute_coordinates: bool,backend: LinearAlgebraBackend,
                       requested_backend: LinearAlgebraBackendPreference) -> DMDKernelOutput:
    if svd_driver!="gesvd":
        raise ValueError("torch-cuda currently supports only reference-compatible svd_driver='gesvd'.")
    torch=import_module("torch"); device=torch.device("cuda")
    rows,columns=X.shape; width=min(rows,columns)
    estimated_bytes=int(2.*(2*rows*columns*8+rows*width*8+width*columns*8+
                            width*width*8+width*width*16+rows*width*16))
    free_bytes,total_bytes=torch.cuda.mem_get_info()
    if estimated_bytes>int(.8*free_bytes):
        raise MemoryError(f"torch-cuda DMD requires an estimated {estimated_bytes/1024**3:.2f} GiB "
                          f"but only {free_bytes/1024**3:.2f} GiB is currently free.")
    synchronize=torch.cuda.synchronize; started=perf_counter()
    xt=torch.as_tensor(X,dtype=torch.float64,device=device)
    yt=torch.as_tensor(Y,dtype=torch.float64,device=device)
    U,s,Vh=torch.linalg.svd(xt,full_matrices=False,driver="gesvd"); synchronize(); decomposed=perf_counter()
    singular_values=s.detach().cpu().numpy(); retained=_retained_rank(singular_values,rank,svd_tolerance)
    if retained<1: raise ValueError("DMD retained no singular values; lower svd_tolerance or inspect the input data.")
    U,s,Vh=U[:,:retained],s[:retained],Vh[:retained]
    operator=((U.mH@yt)@Vh.mH)/s.unsqueeze(0); synchronize(); projected_at=perf_counter()
    eigenvalues,W=torch.linalg.eig(operator); synchronize(); eigendecomposed=perf_counter()
    complex_dtype=W.dtype
    modes=((yt@Vh.mH)/s.unsqueeze(0)).to(complex_dtype)@W if exact else U.to(complex_dtype)@W
    synchronize(); modes_built=perf_counter()
    coordinates=(torch.linalg.pinv(modes)@xt.to(complex_dtype)).mT if compute_coordinates else None
    if coordinates is not None: synchronize()
    operator_np=operator.detach().cpu().numpy(); eigenvalues_np=eigenvalues.detach().cpu().numpy()
    modes_np=modes.detach().cpu().numpy(); coordinates_np=None if coordinates is None else coordinates.detach().cpu().numpy()
    synchronize(); completed=perf_counter()
    timing={"svd_seconds":decomposed-started,"projected_operator_seconds":projected_at-decomposed,
            "eigendecomposition_seconds":eigendecomposed-projected_at,"mode_construction_seconds":modes_built-eigendecomposed,
            "coordinate_recovery_seconds":completed-modes_built,"factorization_seconds":completed-started,
            "clock":"time.perf_counter"}
    provenance=backend.provenance(requested_backend)
    provenance["device_memory"]={"estimated_working_set_bytes":estimated_bytes,
                                  "free_bytes_before":int(free_bytes),"total_bytes":int(total_bytes)}
    return DMDKernelOutput(eigenvalues_np,modes_np,coordinates_np,operator_np,retained,timing,provenance)


def run_dmd_kernel(X: np.ndarray,Y: np.ndarray,*,rank: int|None,exact: bool,
                   svd_tolerance: float|None,svd_driver: str,compute_coordinates: bool,
                   svd_threads: int|None,eigendecomposition_threads: int|None,
                   linear_algebra_backend: LinearAlgebraBackendPreference="auto") -> DMDKernelOutput:
    backend=select_linear_algebra_backend(linear_algebra_backend)
    if svd_driver not in backend.svd_drivers:
        raise ValueError(f"Backend {backend.name!r} does not support SVD driver {svd_driver!r}.")
    if backend.name=="scipy-lapack":
        return _scipy_kernel(X,Y,rank,exact,svd_tolerance,svd_driver,compute_coordinates,
                             svd_threads,eigendecomposition_threads,backend,linear_algebra_backend)
    if svd_threads is not None or eigendecomposition_threads is not None:
        raise ValueError("BLAS thread limits apply only to the SciPy CPU backend, not torch-cuda.")
    return _torch_cuda_kernel(X,Y,rank,exact,svd_tolerance,svd_driver,compute_coordinates,
                              backend,linear_algebra_backend)


__all__=["DMDKernelOutput","LinearAlgebraBackend","LinearAlgebraBackendPreference",
         "available_linear_algebra_backends","run_dmd_kernel","scipy_backend",
         "select_linear_algebra_backend","torch_cuda_backend"]
