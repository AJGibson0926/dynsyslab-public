"""Train-only library fitting, independent validation, untouched test rollouts and rank diagnostics."""
import numpy as np
from .models import fit_field,free_rollout
from .local_field import LocalField
from .paired import PairedField
from .compact_field import CompactField

def relative_rms(actual,prediction):
    """Compute relative state error per component over corresponding samples of two trajectories."""
    if np.shape(actual)!=np.shape(prediction):raise ValueError('Metric shapes disagree')
    error=np.linalg.norm(prediction-actual,axis=0);scale=np.linalg.norm(actual,axis=0)
    return [float(e/s) if s>0 else None for e,s in zip(np.atleast_1d(error),np.atleast_1d(scale))]

def evaluate(field,record,rank,endpoint,**solver):
    """Only the first measured window initializes the ODE. Separate field, state and decoder errors."""
    t,b,db,truth=record['time'],record['b'][:,:rank],record['db'][:,:rank],record['truth']
    out=free_rollout(field,t,b[0],**solver);pred=out['b'];n=len(out['time'])
    metrics=dict(complete=out['complete'],message=out['message'],coverage=n/len(t))
    if not n:return metrics,dict(time=out['time'],b=pred)
    if not np.isfinite(pred).all():raise ValueError('Nonfinite rollout')
    def rates(x):
        return np.vstack([field(x[i:i+512]) for i in range(0,len(x),512)])
    learned=rates(b)
    measured=np.column_stack([b@endpoint,db@endpoint])
    free=np.column_stack([pred@endpoint,rates(pred)@endpoint])
    metrics.update(field=relative_rms(db,learned),field_total=relative_rms(db.ravel(),learned.ravel())[0],
        modal=relative_rms(b[:n],pred),modal_total=relative_rms(b[:n].ravel(),pred.ravel())[0],
        measured=relative_rms(truth,measured),free=relative_rms(truth[:n],free),
        evaluations=out['evaluations'])
    return metrics,dict(time=t[:n],actual=truth[:n],predicted=free,projected=measured[:n],b=pred)

def rank_sweep(basis,ranks,training,validation,test,candidates,*,channel=0,solver=None,validation_weights=None):
    """Each split maps distinct record IDs to {time,b,db,truth}. truth holds (observable, rate).
    Selection averages validation scores only. Every returned fitted field can be serialized.
    """
    if not training or not validation or not test or not candidates:
        raise ValueError('Nonempty training, validation and test splits and candidates are required')
    if set(training)&set(validation) or set(test)&(set(training)|set(validation)):
        raise ValueError('Training, validation and test identities must be disjoint')
    solver={} if solver is None else solver
    weights=np.ones(len(validation)) if validation_weights is None else np.asarray(validation_weights,dtype=float)
    if weights.shape!=(len(validation),) or not np.isfinite(weights).all() or np.any(weights<0) or weights.sum()<=0:
        raise ValueError('Invalid validation weights')
    weights=weights/weights.sum()
    for rank in ranks:
        ep=basis.endpoint(rank,channel)
        records=[dict(time=r['time'],b=r['b'][:,:rank],db=r['db'][:,:rank]) for r in training.values()]
        trials=[];best=None;best_score=np.inf
        for config in candidates:
            opts=dict(config);family=opts.pop('family','polynomial')
            try:
                if family=='local':field=LocalField.fit(records,**opts)
                elif family=='compact':field=CompactField.fit(records,**opts)
                elif family=='oscillator':field=PairedField.fit(records,**opts)
                else:field=fit_field(records,family=family,**opts)
                stats=[evaluate(field,r,rank,ep,**solver)[0] for r in validation.values()]
                good=all(s['complete'] and None not in [s['modal_total'],*s['free']] for s in stats)
                score=float(sum(w*(.5*s['modal_total']+.25*sum(s['free'])) for w,s in zip(weights,stats))) if good else None
                trials.append(dict(config=config,validation=stats,score=score))
                if score is not None and score<best_score:best_score=score;best=field
            except (ArithmeticError,ValueError,RuntimeError,np.linalg.LinAlgError) as err:
                trials.append(dict(config=config,score=None,error=str(err)))
        result=dict(rank=rank,energy=float(basis.energy[rank-1]),trials=trials)
        if best is None:
            yield result,None,{}
            continue
        results={};curves={}
        for role,split in [('training',training),('test',test)]:
            results[role]={}
            for name,r in split.items():
                results[role][name],curves[name]=evaluate(best,r,rank,ep,**solver)
        result.update(validation_score=best_score,results=results)
        yield result,best,curves
