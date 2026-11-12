"""Sparse generator nullspace kernels for Koopman-with-SINDy."""
from __future__ import annotations

import numpy as np

from .spectrum import KoopmanEigenvalue


def _matlab_single_absolute_null_vector(matrix: np.ndarray, threshold: float, *, diagnostics: dict | None=None) -> np.ndarray:
    """Reproduce the legacy MATLAB SINDy eigenfunction coefficient step.

    MATLAB computes the final right singular vector of the unscaled generator
    matrix once and then zeros entries whose absolute magnitude is below the
    fixed tolerance.  This opt-in policy intentionally performs no iterative
    refit and leaves the accepted default policy unchanged.
    """
    if threshold<0: raise ValueError("Koopman-with-SINDy threshold must be nonnegative.")
    _,_,vh=np.linalg.svd(np.asarray(matrix,dtype=complex),full_matrices=False)
    if diagnostics is not None: diagnostics["svd_subproblems"]=diagnostics.get("svd_subproblems",0)+1
    vector=vh.conj().T[:,-1]
    vector[np.abs(vector)<threshold]=0.
    return vector

def _sparse_null_vector(matrix: np.ndarray, threshold: float, max_iterations: int, *, diagnostics: dict | None=None) -> np.ndarray:
    """Sequentially threshold the smallest-singular-vector null-space estimate."""
    if threshold<0: raise ValueError("Koopman-with-SINDy threshold must be nonnegative.")
    scale=np.linalg.norm(matrix,axis=0); scale=np.where(scale>np.finfo(float).eps,scale,1.)
    scaled=matrix/scale; active=np.ones(matrix.shape[1],dtype=bool)
    for _ in range(max_iterations):
        _,_,vh=np.linalg.svd(scaled[:,active],full_matrices=False)
        if diagnostics is not None: diagnostics["svd_subproblems"]=diagnostics.get("svd_subproblems",0)+1
        vector=vh.conj().T[:,-1]
        cutoff=threshold*np.max(np.abs(vector))
        retained=np.abs(vector)>=cutoff
        if not np.any(retained): retained[np.argmax(np.abs(vector))]=True
        next_active=active.copy(); next_active[active]=retained
        if np.array_equal(next_active,active): break
        active=next_active
    _,_,vh=np.linalg.svd(scaled[:,active],full_matrices=False)
    if diagnostics is not None: diagnostics["svd_subproblems"]=diagnostics.get("svd_subproblems",0)+1
    compact=vh.conj().T[:,-1]; result=np.zeros(matrix.shape[1],dtype=complex); result[active]=compact/scale[active]
    return result/np.linalg.norm(result)

def _gram_smallest_vector(gram: np.ndarray, active: np.ndarray) -> tuple[np.ndarray|None,dict]:
    """Return a trustworthy smallest eigenvector or request exact-SVD fallback."""
    indices=np.flatnonzero(active); block=np.asarray(gram[np.ix_(indices,indices)],dtype=complex)
    block=(block+block.conj().T)/2
    if block.shape[0]==1: return np.ones(1,dtype=complex),{"relative_gap":float("inf"),"error_indicator":0.,"negative_floor":False}
    eigenvalues,eigenvectors=np.linalg.eigh(block); vector=eigenvectors[:,0]
    spectral_scale=max(float(np.max(np.abs(eigenvalues))),np.finfo(float).eps)
    relative_gap=float(max(0.,eigenvalues[1]-eigenvalues[0])/spectral_scale)
    residual=float(np.linalg.norm(block@vector-eigenvalues[0]*vector)/spectral_scale)
    error_indicator=float(residual/max(relative_gap,np.finfo(float).eps))
    negative_floor=bool(eigenvalues[0] < -1e-10*spectral_scale)
    trustworthy=not negative_floor and relative_gap>1e-10 and error_indicator<=1e-7
    return (vector if trustworthy else None),{"relative_gap":relative_gap,"error_indicator":error_indicator,"negative_floor":negative_floor}

def _sparse_null_vector_from_gram(gram: np.ndarray, threshold: float, max_iterations: int, *, matrix_factory, diagnostics: dict) -> np.ndarray:
    """Sequential thresholding from a small Gramian with guarded SVD fallback."""
    if threshold<0: raise ValueError("Koopman-with-SINDy threshold must be nonnegative.")
    gram=np.asarray(gram,dtype=complex); diagonal=np.maximum(np.real(np.diag(gram)),0.)
    scale=np.sqrt(diagonal); scale=np.where(scale>np.finfo(float).eps,scale,1.)
    normalized=gram/(scale[:,None]*scale[None,:]); active=np.ones(gram.shape[0],dtype=bool); scaled_matrix=None
    def smallest(mask):
        nonlocal scaled_matrix
        vector,quality=_gram_smallest_vector(normalized,mask)
        diagnostics["minimum_relative_gap"]=min(diagnostics.get("minimum_relative_gap",float("inf")),quality["relative_gap"])
        diagnostics["maximum_error_indicator"]=max(diagnostics.get("maximum_error_indicator",0.),quality["error_indicator"])
        if vector is not None:
            diagnostics["gram_eigh_subproblems"]=diagnostics.get("gram_eigh_subproblems",0)+1
            return vector
        diagnostics["gram_eigh_fallbacks"]=diagnostics.get("gram_eigh_fallbacks",0)+1
        if scaled_matrix is None: scaled_matrix=np.asarray(matrix_factory(),dtype=complex)/scale
        _,_,vh=np.linalg.svd(scaled_matrix[:,mask],full_matrices=False)
        diagnostics["svd_subproblems"]=diagnostics.get("svd_subproblems",0)+1
        return vh.conj().T[:,-1]
    for _ in range(max_iterations):
        vector=smallest(active); cutoff=threshold*np.max(np.abs(vector)); retained=np.abs(vector)>=cutoff
        if not np.any(retained): retained[np.argmax(np.abs(vector))]=True
        next_active=active.copy(); next_active[active]=retained
        if np.array_equal(next_active,active): break
        active=next_active
    compact=smallest(active); result=np.zeros(gram.shape[0],dtype=complex); result[active]=compact/scale[active]
    return result/np.linalg.norm(result)

def _generator_gram_terms(theta: np.ndarray, gamma: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Precompute the three small matrices defining every generator Gramian."""
    return gamma.conj().T@gamma,gamma.conj().T@theta,theta.conj().T@theta

def _generator_gram(terms: tuple[np.ndarray,np.ndarray,np.ndarray], eigenvalue: complex) -> np.ndarray:
    gamma_gamma,gamma_theta,theta_theta=terms
    return gamma_gamma-eigenvalue*gamma_theta-np.conjugate(eigenvalue)*gamma_theta.conj().T+abs(eigenvalue)**2*theta_theta

def _conjugate_fit_plan(values: list[KoopmanEigenvalue], theta: np.ndarray, gamma: np.ndarray, *, enabled: bool) -> list[tuple[int,bool]]:
    """Map each eigenvalue to a directly fitted or conjugate-reused coefficient row.

    For real-valued dictionaries and generator data, ``G(conj(lambda))`` is
    exactly ``conj(G(lambda))``. Its right singular vector and every
    magnitude-based sparsification decision are therefore conjugates too.
    """
    plan=[(index,False) for index in range(len(values))]
    if not enabled or not (np.allclose(theta.imag,0.,atol=1e-14,rtol=0.) and np.allclose(gamma.imag,0.,atol=1e-14,rtol=0.)):
        return plan
    assigned=np.zeros(len(values),dtype=bool)
    for index,value in enumerate(values):
        if assigned[index]: continue
        assigned[index]=True
        if value.is_real: continue
        target=np.conjugate(value.continuous_time)
        partner=next((candidate for candidate in range(index+1,len(values)) if not assigned[candidate]
                      and np.isclose(values[candidate].continuous_time,target,atol=1e-9,rtol=1e-7)),None)
        if partner is not None:
            plan[partner]=(index,True); assigned[partner]=True
    return plan
