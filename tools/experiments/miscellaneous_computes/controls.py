"""Optional rank-study diagnostics using the same fixed training and validation split."""
import numpy as np
from .audit import replay_controls
from .history import separation
from .readouts import NeighborReadout,TangentReadout
from .study import relative_rms

def run_controls(config,rank,sets,prepared,field,curves,endpoint,solver):
    """Evaluate restart, tolerance, readout and visit-separation controls for one fitted modal rank."""
    train,valid,test=(list(sets[role]) for role in ('training','validation','test'))
    report={};by_rank=lambda a:a[:,:rank]
    if config.get('neighbor_count') or config.get('tangent_neighbors'):
        records=[prepared[k] for k in train]
        b=np.vstack([by_rank(r['b']) for r in records]);db=np.vstack([by_rank(r['db']) for r in records])
        truth=np.vstack([r['truth'] for r in records])
        k=config.get('neighbor_count')
        lookup=NeighborReadout(b,truth) if k else None
        sizes=config.get('tangent_neighbors',[])
        if sizes:
            smooth=TangentReadout(b,db,truth[:,0],truth[:,1],config.get('amplitude_floor',1.))
            losses={}
            for size in sizes:
                per_record=[]
                for name in valid:
                    r=prepared[name]
                    phase=smooth.phase(by_rank(r['b']),by_rank(r['db']),size)
                    errors=relative_rms(r['truth'],phase)
                    if None in errors:raise ValueError('Undefined validation readout error on a zero-energy channel')
                    per_record.append(np.mean(errors))
                losses[size]=float(np.mean(per_record))
            size=min(losses,key=losses.get)
            report['tangent_selection']=dict(neighbors=size,validation_measured_phase=losses)
        report['readouts']={}
        for name in test:
            r=prepared[name];curve=curves[name];n=len(curve['time']);x=by_rank(r['b'][:n]);y=r['truth'][:n]
            if not n:continue
            prediction=curve['b'];scores={'coverage':n/len(r['time'])}
            if lookup:
                for label,state in [('measured',x),('free',prediction)]:
                    phase=lookup.predict(state,k)
                    scores['neighbor_'+label]=relative_rms(y,phase);curve['neighbor_'+label]=phase
            if sizes:
                dp=np.vstack([field(prediction[i:i+512]) for i in range(0,n,512)])
                for label,state,derivative in [('measured',x,by_rank(r['db'][:n])),('free',prediction,dp)]:
                    phase=smooth.phase(state,derivative,size)
                    scores['tangent_'+label]=relative_rms(y,phase);curve['tangent_'+label]=phase
            report['readouts'][name]=scores
    for name,time in config.get('restart_time',{}).items():
        report.setdefault('restart',{})[name]=replay_controls(field,prepared[name],rank,endpoint,
            restart_time=time,compare_from=config.get('compare_from',{}).get(name),**solver)
    for name,pair in config.get('visits',{}).items():
        r=prepared[name];difference,distance=separation(r['time'],by_rank(r['b']),*pair)
        report.setdefault('visits',{})[name]=dict(times=pair,difference=difference.tolist(),distance=distance.tolist())
    return report
