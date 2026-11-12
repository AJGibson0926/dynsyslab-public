"""Sampling-convergence controls for empirical closure witnesses."""
from fractions import Fraction
from math import factorial
import numpy as np
from scipy.optimize import linear_sum_assignment

def difference_weights(offsets, order, *, exact=False):
    """Weights at zero on distinct offsets. For integer-grid offsets divide by step**order."""
    values=list(offsets);n=len(values)
    if not isinstance(order,int) or not 0<=order<n or len(set(values))!=n or not np.isfinite(np.asarray(values,float)).all():
        raise ValueError('Invalid derivative order or stencil nodes')
    if not exact:
        x=np.asarray(values,float);scale=max(np.max(abs(x)),1.)
        A=np.vander(x/scale,N=n,increasing=True).T;b=np.zeros(n);b[order]=factorial(order)
        return np.linalg.solve(A,b)/scale**order
    x=[Fraction(v) for v in values]
    aug=[[v**i for v in x]+[Fraction(factorial(i) if i==order else 0)] for i in range(n)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if aug[i][j])
        aug[j],aug[pivot]=aug[pivot],aug[j]
        divisor=aug[j][j];aug[j]=[v/divisor for v in aug[j]]
        for i in range(n):
            if i!=j:
                factor=aug[i][j];aug[i]=[a-factor*b for a,b in zip(aug[i],aug[j])]
    return [row[-1] for row in aug]

def match_visits(first,second,time_tolerance):
    """One-to-one maximum-cardinality matching of ordered visit times under a per-visit tolerance."""
    tol=np.broadcast_to(np.asarray(time_tolerance,float),(2,))
    if np.any(tol<=0) or not np.isfinite(tol).all():raise ValueError('Positive finite visit tolerances required')
    def times(items):
        v=np.asarray([r['times'] for r in items],float).reshape(-1,2)
        if not np.isfinite(v).all() or np.any(v[:,0]>=v[:,1]):raise ValueError('Expected ordered visit pairs')
        return v
    a,b=times(first),times(second);m,n=len(a),len(b)
    # A dummy assignment costs more than all valid assignment distances combined.
    cost=np.max(abs(a[:,None]-b[None,:])/tol,axis=-1)
    penalty=float(m+n+1)
    augmented=np.column_stack([np.where(cost<=1,cost,2*penalty),np.full((m,m),penalty)])
    rows,cols=linear_sum_assignment(augmented);pairs=[(int(i),int(j)) for i,j in zip(rows,cols) if j<n]
    return dict(matches=[dict(first=i,second=j,time_shift=(b[j]-a[i]).tolist()) for i,j in pairs],
                unmatched_first=sorted(set(range(m))-{i for i,_ in pairs}),
                unmatched_second=sorted(set(range(n))-{j for _,j in pairs}))
