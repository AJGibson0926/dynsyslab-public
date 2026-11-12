"""Trajectory calibration on a sparse field's selected coefficients, using training shots only."""
import numpy as np
from scipy.optimize import least_squares
from .models import SparseField,free_rollout,records_checked

def calibrate(field,training,*,duration,starts_per_record=2,active=None,max_evaluations=40,
              coefficient_penalty=1e-4,solver=None):
    """Refit a chosen subset of the STLSQ support against short free trajectories.

    Return a new field and provenance. The uncalibrated field stays usable.
    Failed shots carry a finite penalty; an incomplete or worse final candidate is rejected.
    """
    records=records_checked(training);C=field.coef;mask=C!=0 if active is None else np.asarray(active,bool)
    if (mask.shape!=C.shape or np.any(mask&(C==0)) or not mask.any() or not np.isfinite(duration)
        or duration<=0 or not isinstance(starts_per_record,int) or starts_per_record<1 or max_evaluations<1
        or coefficient_penalty<0):raise ValueError('Invalid shooting settings or selected support')
    shots=[]
    for r in records:
        t=r['time'];possible=np.flatnonzero(t+duration<=np.nextafter(t[-1],np.inf))
        if not len(possible):raise ValueError('Training record shorter than requested shot')
        for i in np.unique(np.linspace(0,possible[-1],starts_per_record).astype(int)):
            stop=t.searchsorted(t[i]+duration,side='right')
            if stop<=i+1:raise ValueError('Shot has fewer than two output samples')
            shots.append((t[i:stop],r['b'][i:stop]))
    scale=np.sqrt(sum(np.sum(b*b,axis=0) for _,b in shots)/sum(len(b) for _,b in shots))
    scale=np.where(scale>0,scale,1.)
    c0=C[mask];unit=np.maximum(abs(c0),1e-10);opts={} if solver is None else solver
    def evaluate(change):
        cc=C.copy();cc[mask]=c0+unit*change;model=SparseField(field.spec,cc)
        residual=[];complete=True
        for t,b in shots:
            part=np.full_like(b,1000.)
            try:
                pred=free_rollout(model,t,b[0],**opts);k=len(pred['time'])
                if k:part[:k]=(pred['b']-b[:k])/scale
                complete &= pred['complete']
            except (FloatingPointError,ArithmeticError,ValueError):complete=False
            residual.append(part.ravel()/np.sqrt(b.size))
        residual.append(np.sqrt(coefficient_penalty)*change)
        return np.concatenate(residual),model,complete
    zero=np.zeros_like(c0);before,_,_=evaluate(zero)
    sol=least_squares(lambda d:evaluate(d)[0],zero,max_nfev=max_evaluations)
    final,candidate,complete=evaluate(sol.x)
    accepted=bool(complete and final@final<=before@before)
    return (candidate if accepted else field),dict(recipe='training-shooting',shots=len(shots),
        fitted_terms=len(c0),initial_loss=float(before@before),candidate_loss=float(final@final),
        accepted=accepted,evaluations=sol.nfev)
