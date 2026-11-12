"""Memory-bounded STLSQ for large sparse libraries with per-component active sets."""
import numpy as np
from scipy import sparse

def sparse_stlsq(gram, cross, target_norm, ridge, threshold, *, tolerance=2e-8,
                 max_cg=500, max_threshold=12):
    """Solve normalized ridge problems and threshold independently for each output.

    gram = X.T @ X, cross = X.T @ Y, and target_norm contains column norms of Y,
    including any chosen row weighting. Returns feature-by-output coefficients.
    Batched conjugate gradients reuse one sparse operator across different
    supports and need no sparse direct factorization or dense Gram conversion.
    """
    G=sparse.csr_matrix(gram,dtype=float);R=np.asarray(cross,dtype=float)
    norm=np.asarray(target_norm,dtype=float)
    if (G.shape[0]!=G.shape[1] or R.ndim!=2 or R.shape[0]!=G.shape[0]
        or norm.shape!=(R.shape[1],) or not np.isfinite(norm).all() or np.any(norm<0)
        or not np.isfinite(G.data).all() or not np.isfinite(R).all() or ridge<=0 or threshold<0
        or tolerance<=0 or max_cg<1 or max_threshold<1):
        raise ValueError('Invalid sparse regression moments, dimensions or solver options')
    diag=G.diagonal()
    scale=np.sqrt(np.maximum(diag,1e-30))
    scaled=G.multiply((1/scale)[:,None]).multiply((1/scale)[None,:]).tocsr()
    scaled.setdiag(scaled.diagonal()+ridge)
    rhs=R/scale[:,None];active=np.broadcast_to((diag>0)[:,None],R.shape).copy()
    weights=np.zeros_like(R)
    for pruning in range(max_threshold+1):
        rhs_mask=rhs*active;weights*=active
        residual=(rhs_mask-scaled@weights)*active
        p=residual.copy();squared=np.sum(residual**2,axis=0)
        limit=np.maximum(np.sum(rhs_mask**2,axis=0)*tolerance**2,1e-60)
        for _ in range(max_cg):
            working=squared>limit
            if not np.any(working):break
            p[:,~working]=0.
            product=(scaled@p)*active
            denom=np.sum(p*product,axis=0)
            if np.any(denom[working]<=0):raise ArithmeticError('Sparse Gram is not positive definite')
            step=np.divide(squared,denom,out=np.zeros_like(squared),where=working)
            weights+=p*step;residual-=product*step
            updated=np.sum(residual**2,axis=0)
            beta=np.divide(updated,squared,out=np.zeros_like(updated),where=working)
            p=residual+p*beta;squared=updated
        if np.any(squared>limit):raise ArithmeticError('Sparse regression did not converge')
        retained=active&(abs(weights)>=threshold*norm)
        if np.array_equal(retained,active):break
        if pruning==max_threshold:raise ArithmeticError('Sparse threshold support did not converge')
        active=retained
    return weights/scale[:,None]
