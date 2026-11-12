"""Track a rebound jump down to its training neighbors and compare independent solvers."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from study import errors, load, measurements, rms


def read_run(folder):
    stats=json.loads((folder/'results.json').read_text())
    with np.load(folder/'diagnostics.npz') as z:
        return stats,{key:z[key] for key in z.files}


def explain(exe, folder, replicate):
    stats,d=read_run(folder);other,s=read_run(replicate)
    train=load(Path(stats['training']));held=load(Path(stats['held_out']))
    for key in ('receiver_config','parameters'):
        if train['meta'][key]!=load_metadata(Path(other['training']))[key]:
            raise ValueError('Solver control uses different geometry or physics')
    if not np.array_equal(d['time'],s['time']) or not np.allclose(d['truth'][0],s['truth'][0]):
        raise ValueError('Solver control must score the same clock and trajectory')
    t=d['time']; x=d['truth']; y=d['monopole_two']; original=d['monopole_observations']
    bounds=d['bounds'];complete=bounds[(bounds>=t[0])&(bounds<=t[-1])]
    if len(bounds)<2:raise ValueError('A bounded ring is required to diagnose rebound')
    between=(t>=bounds[0])&(t<bounds[1])
    collapse=np.flatnonzero(between)[np.argmin(x[between,0])]
    # Skip the fastest physical turnaround. Identify the largest increment in
    # the velocity estimation error in the remaining part of the first rebound.
    start=t[collapse]+.03*(bounds[1]-bounds[0])
    steps=np.flatnonzero((t[:-1]>=start)&(t[1:]<bounds[1]))
    if not len(steps):raise ValueError('Too few samples after the collapse')
    i=int(steps[np.argmax(np.abs(np.diff(y[:,1]-x[:,1])[steps]))])
    train_bounds=train['ring_bounds'];train_t=d['train_time'];train_x=d['train_state']
    event=[]
    for j in (i,i+1):
        ids=d['neighbor_rows'][j]; weights=d['weights'][j]
        ring=np.searchsorted(train_bounds,train_t[ids],side='right')
        if not np.allclose(weights@train_x[ids],y[j],atol=1e-12):raise ValueError('Cannot replay lookup')
        event.append(dict(tau=float(t[j]),state=x[j].tolist(),estimate=y[j].tolist(),
                          pressure_Pa=(original[j]/held['ranges_m']).tolist(),
                          training_neighbors=[dict(row=int(row),tau=float(train_t[row]),ring=int(idx),
                                                  v=float(train_x[row,1]),weight=float(ww))
                                              for row,idx,ww in zip(ids,ring,weights)],
                          weight_by_ring={str(int(idx)):float(weights[ring==idx].sum()) for idx in np.unique(ring)},
                          replicate_estimate=s['monopole_two'][j].tolist()))
    differences={}
    for key in ('truth','monopole_two','near_field_two','near_field_aligned','monopole_observations'):
        if key in d and key in s:
            delta=d[key]-s[key]
            differences[key]={'rms':rms(delta).tolist(),'max_abs':np.max(abs(delta),axis=0).tolist()}
    ft,_=measurements(train,(0,));scale=rms(ft[:,0])
    # Change training density only. The held-out pressure samples stay bit-identical.
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        np.savez(root/'train.npz',features=ft[::2],state=train_x[::2],scales=[scale,scale])
        np.savez(root/'query.npz',features=original,time=t)
        subprocess.run([exe,'estimate',root/'train.npz',root/'query.npz',str(stats['k']),root/'est.npz'],check=True)
        with np.load(root/'est.npz') as z:sparse=z['prediction']
    coef=held['ranges_m']; error_rms=rms(y-x)
    result=dict(example_index=i,collapse_tau=float(t[collapse]),rebound_event=event,
                true_velocity_change=float(x[i+1,1]-x[i,1]),estimated_velocity_change=float(y[i+1,1]-y[i,1]),
                solver_differences=differences,lookup_error_rms=error_rms.tolist(),
                simulation_repeat_to_lookup_error_ratio=(rms(d['monopole_two']-s['monopole_two'])/error_rms).tolist(),
                half_training_density_same_input_rms_change=rms(sparse-y).tolist(),
                half_training_density_event_estimates=sparse[i:i+2].tolist(),
                pressure_solver_difference_rms_pa=(rms(d['monopole_observations']-s['monopole_observations'])/coef).tolist(),
                complete_rings=len(complete)-1, ring_intervals=len(bounds)-1,
                first_eight_complete_tau=None if len(complete)<9 else [complete[0],complete[8]],
                group_scores={})
    if len(complete)>=9:
        for group,edges in [('first_eight_complete',complete[[0,8]]),('last_eight_complete',complete[[-9,-1]])]:
            mask=(t>=edges[0])&(t<edges[1])
            result['group_scores'][group]={key:errors(x,d[key],mask) for key in
                ['monopole_one','monopole_two','near_field_two','near_field_aligned','monopole_linear'] if key in d}
    offsets=train_t[d['neighbor_rows'][:,0]]-t
    ring_stats=[]
    for idx,(a,b) in enumerate(zip(bounds[:-1],bounds[1:])):
        m=(t>=a)&(t<b)
        if a<t[0] or b>t[-1] or not m.any():continue
        median=float(np.median(offsets[m]))
        ring_stats.append(dict(ring=idx+1, start=a, end=b,
            velocity_error_pct=errors(x,y,m)['percent'][1],median_training_time_offset=median,
            fraction_outside_quarter_ring_of_median=float(np.mean(abs(offsets[m]-median)>.25*(b-a)))))
    pairs=[(a,b) for a,b in zip(ring_stats[:-1],ring_stats[1:]) if a['velocity_error_pct'] is not None and b['velocity_error_pct'] is not None]
    if pairs:result['largest_ring_error_drop']=max(pairs,key=lambda pair:pair[0]['velocity_error_pct']-pair[1]['velocity_error_pct'])
    return result


def load_metadata(path):
    with np.load(path) as z:return json.loads(z['metadata_json'].item())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('exe','study','replicate','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    result=explain(args.exe.resolve(),args.study,args.replicate)
    args.output.write_text(json.dumps(result,indent=2))
