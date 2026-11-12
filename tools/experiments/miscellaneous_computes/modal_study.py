"""Explicit delay-grid and model-family selection, with frozen training charts and test-only scoring."""
import json
from pathlib import Path
import numpy as np
from .history import matrix
from .sampling import clock,delay_samples,sample_indices
from .tapped_history import sample_delay_jet,delay_consistency
from .modal import ModalChart,ModalModel
from .models import fit_field
from .local_field import LocalField
from .compact_field import CompactField
from .paired import PairedField
from .structured import PairNetwork,SecondOrderField,InteriorField,EndpointAcceleration
from .diagnostics import errors,forecast,kinematics
from .replay import replay
from .manifest import split_paths,publish_directory,file_id

def load_record(path):
    """Load a trajectory clock, measurement matrix and aligned rate labels from its NPZ file."""
    with np.load(path,allow_pickle=False) as z:
        t=clock(z['time']);y=matrix(z['values']);dy=matrix(z['rates'])
    if y.shape!=dy.shape or len(y)!=len(t):raise ValueError(f'Unaligned observations in {path}')
    return t,y,dy

def window_blocks(record,lags,indices,interpolation='linear',block_size=8192):
    """Yield endpoint times, causal observation windows and their rate labels in bounded row blocks."""
    t,y,dy=record
    if interpolation not in ('linear','quintic_pair') or not isinstance(block_size,int) or block_size<1:
        raise ValueError('Invalid interpolation or block size')
    for i in range(0,len(indices),block_size):
        ids=indices[i:i+block_size];q=t[ids];lo=max(0,t.searchsorted(q[0]-max(lags))-1);hi=ids[-1]+1
        if hi-lo<2:lo=max(0,lo-1);hi=max(hi,2)
        ts=t[lo:hi]
        if interpolation=='quintic_pair':
            if y.shape[1]!=2 or not np.allclose(dy[lo:hi,0],y[lo:hi,1],rtol=1e-7,atol=1e-10):
                raise ValueError('Quintic pair requires [value, rate] and [rate, acceleration]')
            jet=np.column_stack([y[lo:hi],dy[lo:hi,1]])
            v=sample_delay_jet(ts,jet,q,lags)
            yield q,v[:,:,:2],v[:,:,1:]
        else:
            v=delay_samples(ts,np.column_stack([y[lo:hi],dy[lo:hi]]),q,lags)
            yield q,v[:,:,:y.shape[1]],v[:,:,y.shape[1]:]

def fit_candidate(cfg,chart,recs):
    """Fit one direct or decoded candidate on the supplied training records in its declared chart."""
    opts=dict(cfg);family=opts.pop('family','polynomial')
    calibration=opts.pop('calibration',None)
    if family in ('second-order','pair-network','damped-interior'):
        if calibration:raise ValueError('Shooting calibration requires a direct sparse dictionary')
        if chart.encoder.shape[0]!=chart.encoder.shape[1] or len(chart.scales)!=2:
            raise ValueError('Structured value/rate fields need full rank and two channels')
        raw=[dict(time=r['time'],b=r['b']@chart.decoder/chart.gain,db=r['db']@chart.decoder/chart.gain) for r in recs]
        if family=='pair-network':field=PairNetwork.fit(raw,chart.lags,**opts)
        elif family=='damped-interior':
            inherited=opts.pop('inherited_model',None);tol=opts.pop('lag_tolerance',1.e-12)
            if inherited is not None:
                ref=ModalModel.load(inherited)
                if tol<0 or not np.isfinite(tol) or not np.allclose(ref.chart.lags,chart.lags[[0,-1]],atol=tol,rtol=0):
                    raise ValueError('Inherited model endpoints differ from requested lags')
                opts['outer']=EndpointAcceleration(ref)
            field=InteriorField.fit(raw,**opts)
        else:field=SecondOrderField.fit(raw,**opts)
        return ModalModel(chart,field,'decoded')
    if family=='local':field=LocalField.fit(recs,**opts)
    elif family=='compact':field=CompactField.fit(recs,**opts)
    elif family=='oscillator':field=PairedField.fit(recs,**opts)
    else:field=fit_field(recs,family=family,**opts)
    info={}
    if calibration is not None:
        from .calibration import calibrate
        from .models import SparseField
        if not isinstance(field,SparseField):raise ValueError('Calibration requires a global sparse field')
        field,info=calibrate(field,recs,**calibration)
    return ModalModel(chart,field,metadata={'calibration':info} if info else {})

def evaluate_candidate(model,r,value,rate,solver,prominence=None):
    """Compare the field, measured decode and single-initial-state forecast against the full reference interval."""
    t=r['time'];xi=r['b']/model.chart.gain
    run=replay(model,t,initial=xi[0],**solver);n=len(run['time'])
    def rhs(q):
        return np.vstack([model.rhs(q[i:i+512]) for i in range(0,len(q),512)])
    truth=r['truth'];source_rate=r['truth_rate']
    ref=np.column_stack([truth[:,value],source_rate[:,value] if rate is None else truth[:,rate]])
    first=model.chart.decode(xi)[:,0,:]
    true_dx=r['db']/model.chart.gain
    ref_projection=np.column_stack([first[:,value],model.chart.decode(true_dx)[:,0,value] if rate is None else first[:,rate]])
    metrics=dict(complete=bool(run['complete']),message=run['message'],field=errors(true_dx,rhs(xi)),
                 measured_decode=errors(ref,ref_projection))
    q=run['xi'];current=run['values'][:,0,:]
    pred=np.empty((0,2));chain=np.empty(0)
    if n:
        metrics['modal']=errors(xi[:n],q)
        chain=model.chart.decode(rhs(q))[:,0,value]
        pred=np.column_stack([current[:,value],chain if rate is None else current[:,rate]])
        if rate is not None:
            denom=np.linalg.norm(ref[:n,1])
            metrics['rate_mismatch']=float(np.linalg.norm(chain-pred[:,1])/denom) if denom else None
            if len(model.chart.scales)==2 and value==0 and rate==1 and n>1:
                metrics['delay_consistency']=delay_consistency(t[:n],run['values'].reshape(n,-1),r['taps'][:n].reshape(n,-1),model.chart.lags)
    metrics['output']=forecast(t,ref,run['time'],pred,prominence=prominence)
    return metrics,dict(time=run['time'],xi=q,predicted=pred,actual=ref[:n],chain_rate=chain,
                        reference_time=t,reference=ref)

def run_manifest(manifest,destination,library_path=None):
    """Fit each requested geometry, select its law on validation records and publish models, scores and curves."""
    path=Path(manifest).resolve();cfg=json.loads(path.read_text(encoding='utf-8'))
    if cfg.get('schema')!='dynsyslab.miscellaneous-computes/modal-study/v1':raise ValueError('Unsupported modal study')
    splits=split_paths(cfg,path.parent);data={k:load_record(v) for role in splits.values() for k,v in role.items()}
    observations=cfg.get('observation_provenance',{})
    if not isinstance(observations,dict) or not set(observations)<=data.keys():
        raise ValueError('Observation provenance must be keyed by declared trajectory IDs')
    basis_ids=cfg.get('basis_training',list(splits['training']))
    if not basis_ids or not set(basis_ids)<=set(splits['training']) or len(set(basis_ids))!=len(basis_ids):
        raise ValueError('Basis sources must be unique training IDs')
    grids=cfg['geometries'];candidates=cfg['candidates']
    if not grids or not candidates:raise ValueError('Empty geometry or candidate list')
    space=cfg.get('experiment','direct')
    structured={'pair-network','second-order','damped-interior'}
    if space not in ('direct','structured','mixed') or any(
        (space=='direct' and c.get('family') in structured) or
        (space=='structured' and c.get('family') not in structured) for c in candidates):
        raise ValueError('Candidate field space disagrees with the stated experiment')
    names=[file_id(g['id']) for g in grids]
    if len(names)!=len(set(names)):raise ValueError('Duplicate geometry ID')
    all_lags=[np.asarray(g['lags'],float) for g in grids]
    if any(l.ndim!=1 or not len(l) or l[0]!=0 or np.any(np.diff(l)<=0) or not np.isfinite(l).all() for l in all_lags):
        raise ValueError('Supply increasing physical lags beginning at zero')
    span=max(l[-1] for l in all_lags)
    method=cfg.get('interpolation','linear');solver=dict(cfg.get('solver',{}))
    if solver.get('method')=='rk4':solver['library']=library_path
    value=cfg.get('value_channel',0);rate=cfg.get('rate_channel')
    nc=data[basis_ids[0]][1].shape[1]
    if any(y.shape[1]!=nc for _,y,_ in data.values()) or not 0<=value<nc or (rate is not None and not 0<=rate<nc):
        raise ValueError('Invalid observable channels')
    scales=cfg.get('basis_scales')
    if scales is None:
        scales=np.sqrt(sum(np.einsum('ij,ij->j',data[k][1],data[k][1]) for k in basis_ids)/sum(len(data[k][0]) for k in basis_ids))
    weights=np.array([cfg.get('validation_weights',{}).get(k,1.) for k in splits['validation']],float)
    if not np.isfinite(weights).all() or np.any(weights<0) or weights.sum()<=0:raise ValueError('Invalid validation weights')
    weights/=weights.sum()
    modal_weight=float(cfg.get('modal_weight',.5))
    if not 0<=modal_weight<=1:raise ValueError('Invalid modal weight')
    def indices(k,lag,role):
        t=data[k][0];lo,hi=cfg.get('windows',{}).get(k,[t[0],t[-1]])
        if not np.isfinite([lo,hi]).all() or hi<lo:raise ValueError('Invalid evaluation interval')
        return np.flatnonzero((t>=max(lo,t[0]+lag))&(t<=hi))
    summary={}
    with publish_directory(destination) as folder:
        for geometry,lags in zip(grids,all_lags):
            def blocks():
                for k in basis_ids:
                    indices0=indices(k,lags[-1],'basis')
                    stride=cfg.get('basis_stride',1)
                    if not isinstance(stride,int) or stride<1 or not len(indices0):raise ValueError('No usable basis windows or invalid stride')
                    for _,h,_ in window_blocks(data[k],lags,indices0[::stride],method):yield h
            chart=ModalChart.fit(blocks(),lags,scales,geometry.get('rank'))
            records={}
            for k in data:
                ids=indices(k,span,'evaluation')
                if len(ids)<2:raise ValueError(f'{k} has insufficient common-window history')
                h=[];v=[]
                for _,x,dx in window_blocks(data[k],lags,ids,method):h.append(x);v.append(dx)
                h=np.concatenate(h);r=chart.record(data[k][0][ids],h,np.concatenate(v))
                r.update(truth=data[k][1][ids],truth_rate=data[k][2][ids],taps=h);records[k]=r
            train=[]
            for k in splits['training']:
                r=records[k];ix=sample_indices(r['db'],cfg.get('fit_stride',1),cfg.get('speed_refinements',[]))
                train.append({key:r[key][ix] for key in ('time','b','db')})
            trials=[];best=None;best_score=np.inf
            for spec in candidates:
                try:
                    options=dict(spec)
                    if 'inherited_model' in options:
                        if options.get('family')!='damped-interior':raise ValueError('Unsupported inherited model recipe')
                        options['inherited_model']=(path.parent/options['inherited_model']).resolve()
                    if options.get('family')=='local' and library_path:options['library']=library_path
                    model=fit_candidate(options,chart,train)
                    metrics=[evaluate_candidate(model,records[k],value,rate,solver)[0] for k in splits['validation']]
                    scores=[]
                    for m in metrics:
                        out=m['output']['errors']
                        if not m['complete'] or out is None or m['modal']['total'] is None or None in out['relative']:break
                        scores.append(modal_weight*m['modal']['total']+(1-modal_weight)*np.mean(out['relative']))
                    score=float(weights@scores) if len(scores)==len(weights) else None
                    trials.append(dict(config=spec,validation=metrics,score=score))
                    if score is not None and score<best_score:
                        best_score=score;best=model;best.metadata.update(recipe=spec,geometry=geometry,training=list(splits['training']),
                            validation=list(splits['validation']),basis_training=basis_ids,origin='inherited' if 'inherited_model' in spec else 'fresh',
                            clock_unit=cfg.get('time_unit','model time'),channel_units=cfg.get('channel_units',[]),
                            basis_channel_rms=np.asarray(scales).tolist(),
                            source_files={k:str(p) for role in ('training','validation') for k,p in splits[role].items()},
                            observation_provenance={k:observations[k] for k in (*splits['training'],*splits['validation']) if k in observations})
                except (ValueError,ArithmeticError,RuntimeError,np.linalg.LinAlgError) as err:
                    trials.append(dict(config=spec,score=None,error=str(err)))
            dest=folder/geometry['id'];dest.mkdir();entry=dict(trials=trials,selected=None,results={})
            if best is not None:
                best.save(dest/'model.npz');entry['selected']=best.metadata;entry['score']=best_score
                (dest/'curves').mkdir()
                for role in ('training','test'):
                    for k in splits[role]:
                        metrics,curves=evaluate_candidate(best,records[k],value,rate,solver,cfg.get('peak_prominence'))
                        entry['results'][k]=dict(role=role,**metrics)
                        np.savez_compressed(dest/'curves'/(k+'.npz'),**curves)
            summary[geometry['id']]=entry
        (folder/'results.json').write_text(json.dumps(dict(schema='dynsyslab.miscellaneous-computes/modal-results/v1',
            manifest=cfg,source=str(path),geometries=summary),indent=2,allow_nan=False),encoding='utf-8')
    return summary
