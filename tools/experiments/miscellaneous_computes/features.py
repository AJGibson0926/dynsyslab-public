"""State-only library constructors for derivative and finite-window sparse regression."""
from itertools import combinations_with_replacement
import numpy as np
from scipy.spatial.distance import cdist

def prepare_library(x,family,degree=2,centers=24,harmonics=5,rational_scales=(.25,1.,4.)):
    """Freeze all basis parameters on training rows before validation and forecasting."""
    if degree<1 or harmonics<1 or centers<1:
        raise ValueError('Positive degree, harmonic count and center count are required')
    n=x.shape[1]
    spec=dict(family=family,n=n,degree=degree,harmonics=harmonics)
    if family in ('polynomial','bounded','localized','rational'):
        spec['powers']=[list(indices) for d in range(1,degree+1)
                        for indices in combinations_with_replacement(range(n),d)]
    if family=='rational':
        widths=np.asarray(rational_scales,float)
        if widths.ndim!=1 or not len(widths) or not np.isfinite(widths).all() or np.any(widths<=0):
            raise ValueError('Regularized reciprocal widths must be positive and finite')
        spec['rational_scales']=widths.tolist()
    if family=='localized':
        raw=x[::max(1,len(x)//20000)]
        data=raw/np.hypot(1,raw/3)
        indices=[np.argmin(np.linalg.norm(data,axis=1))]
        dist=np.full(len(data),np.inf)
        while len(indices)<min(centers,len(data)):
            dist=np.minimum(dist,np.sum((data-data[indices[-1]])**2,axis=1))
            if max(dist)<=1e-14:break
            indices.append(dist.argmax())
        c=data[indices]
        k=min(3,len(c)-1)
        spec.update(centers=c.tolist(),width2=np.maximum(.05,np.partition(cdist(c,c,'sqeuclidean'),k,axis=1)[:,k]).tolist())
    return spec

def library(b,spec):
    """Evaluate the frozen feature specification on row-oriented states in the training convention."""
    x=np.atleast_2d(b);family=spec['family']
    if family in ('gated','paired-gated'):
        energy=np.mean(x*x,axis=1)
        gates=[np.ones(len(x))]+[energy/(energy+s) for s in (.04,.4,4.,16.)]
        if family=='paired-gated':
            for j in range(0,x.shape[1],2):
                p=np.mean(x[:,j:j+2]**2,axis=1);gates.append(p/(1+p))
        return (np.column_stack(gates)[:,:,None]*x[:,None,:]).reshape(len(x),-1)
    if family=='harmonic':
        if x.shape[1]<2:raise ValueError('Harmonics require a modal pair')
        p=x[:,0]**2+x[:,1]**2;energy=np.mean(x*x,axis=1)
        u=(x[:,0]+1j*x[:,1])/np.sqrt(p+.015)
        gates=np.column_stack([np.ones_like(p),energy/(energy+.1),energy/(energy+1),energy/(energy+6),p/(p+1)])
        amplitude=gates*np.sqrt(energy+.015)[:,None]
        result=[x]
        for k in range(1,spec['harmonics']+1):
            h=u**k;result.extend([amplitude*h.real[:,None],amplitude*h.imag[:,None]])
        return np.column_stack(result)
    if family not in ('polynomial','bounded','localized','rational'):
        raise ValueError(f'Unknown SINDy library {family!r}')
    z=x if family=='polynomial' else x/np.hypot(1,x/3)
    columns=[np.ones(len(x))]
    if family!='polynomial':columns.extend(x.T)
    columns.extend(np.prod(z[:,term],axis=1) for term in spec['powers'])
    if family=='rational':
        for width in spec['rational_scales']:
            h=np.hypot(width,x)
            even=(width/h)**2
            odd=(x/h)*(width/h)
            columns.extend([even,odd])
    if family=='localized':
        distance=cdist(z,np.array(spec['centers']),'sqeuclidean')
        weights=(1+distance/np.asarray(spec['width2']))**-2
        weights/=weights.sum(axis=1,keepdims=True)
        columns.append((weights[:,:,None]*np.column_stack([np.ones(len(x)),z])[:,None,:]).reshape(len(x),-1))
    return np.column_stack(columns)

def thresholded_least_squares(gram,outputs,ridge,threshold):
    """Fit one normalized ridge regression per output and iteratively remove terms below the coefficient threshold."""
    if ridge<=0 or threshold<0:raise ValueError('Use a positive ridge and a nonnegative threshold')
    nf=len(gram)-outputs
    norm=np.sqrt(np.maximum(np.diag(gram),1e-60))
    g=gram/norm[:,None]/norm
    system=g[:nf,:nf].copy();system.flat[::nf+1]+=ridge
    result=np.zeros((outputs,nf))
    for j in range(outputs):
        active=np.flatnonzero(norm[:nf]>1e-10*max(norm[:nf]))
        while len(active):
            weights=np.linalg.solve(system[np.ix_(active,active)],g[active,nf+j])
            keep=abs(weights)>=threshold
            if np.all(keep):
                result[j,active]=weights*norm[nf+j]/norm[active]
                break
            active=active[keep]
    return result

def regression_moment(features,values,rates,clock,weak_window):
    """Window integrals are computed inside one trajectory, with the supplied clock."""
    if weak_window:
        from scipy.integrate import cumulative_trapezoid
        length=int(weak_window)
        if length!=weak_window or length<1 or length>=len(values):
            raise ValueError('Invalid weak integration window in sample intervals')
        integral=cumulative_trapezoid(features,clock,axis=0,initial=0)
        duration=(clock[length:]-clock[:-length])[:,None]
        features=(integral[length:]-integral[:-length])/duration
        rates=(values[length:]-values[:-length])/duration
    x=np.column_stack([features,rates])
    return x.T@x
