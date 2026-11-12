"""Run a controlled receiver ablation with the standalone native sidecars."""
import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import expm

TERMS = ('monopole_pressure', 'spherical_near_field_pressure',
         'dipole_induction_pressure', 'dipole_radiation_pressure')
CASES = {'monopole': (0,), 'near_field': (0, 1), 'induction': (0, 2),
         'radiation': (0, 3), 'near_plus_induction': (0, 1, 2), 'all': (0, 1, 2, 3)}


def load(path):
    with np.load(path, allow_pickle=False) as z:
        data = {key: z[key] for key in z.files if key != 'metadata_json'}
        data['meta'] = json.loads(z['metadata_json'].item())
    if data['meta']['schema'] != 'dynsyslab.receiver-signals/v1':
        raise ValueError('Expected native receiver signals')
    return data


def measurements(d, indices):
    # Only the sum of the selected physical pressures enters the inverse.
    raw = sum(d[TERMS[i]] for i in indices)
    past = sum(d['near_at_'+TERMS[i]][:, 1] for i in indices)
    scaled = raw*d['ranges_m']
    contrast = scaled[:, 1]-past*d['ranges_m'][0]
    return scaled, contrast


def rms(x, axis=0):
    return np.sqrt(np.mean(np.square(x), axis=axis))


def errors(truth, prediction, mask):
    delta = rms((prediction-truth)[mask])
    ref = rms((truth-[1., 0.])[mask])
    percent = [float(100*a/b) if b>0 else None for a,b in zip(delta,ref)]
    return {'rmse': delta.tolist(), 'percent': percent}


def linear_inverse(d):
    # Linearized RPE at r=1, v=0. The spherical correction is quadratic
    # in velocity; its first derivative vanishes here. The dipoles vanish too.
    g = d['meta']['dimensionless_groups']; p = d['meta']['parameters']
    k = 3*p['kappa']*g['gas_pressure_ratio']-g['capillarity']
    mu = g['viscosity']; S = g['pressure_scale_pa']*p['equilibrium_radius_m']
    A = np.array([[0., 1.], [-k, -mu]])
    H = -S*np.array([k, mu])
    O = np.stack([H@expm(-delay*A) for delay in d['lag_tau']])
    if np.linalg.matrix_rank(O) != 2:
        raise ValueError('The two-receiver linear map is rank deficient')
    return np.linalg.inv(O)


def digest(*arrays):
    h = sha256()
    for a in arrays:
        h.update(str(a.shape).encode())
        h.update(np.ascontiguousarray(a).view(np.uint8))
    return h.hexdigest()


def _compute(exe, train_path, held_path, out, k):
    train, held = load(train_path), load(held_path)
    if len(train['ranges_m']) != 2 or train['ranges_m'][0] >= train['ranges_m'][1]:
        raise ValueError('This two-receiver comparison expects the nearer receiver in column zero')
    for key in ('receiver_config', 'parameters', 'time_scale_s'):
        if train['meta'][key] != held['meta'][key]:
            raise ValueError('Mismatched training and held-out setup: '+key)
    if train['meta']['source'] == held['meta']['source']:
        raise ValueError('Training and held-out trajectories must be distinct')
    t = held['time']; x = held['state']; bounds = held['ring_bounds']
    windows = {'full': np.ones(len(t), bool), 'early': t <= 5, 'late': t >= 76}
    if len(bounds)>=9:
        windows['rings_1_to_8']=(t>=bounds[0]) & (t<bounds[8])
    complete=bounds[(bounds>=t[0]) & (bounds<=t[-1])]
    if len(complete)>=9:
        windows['first_8_complete']=(t>=complete[0]) & (t<complete[8])
        windows['last_8_complete']=(t>=complete[-9]) & (t<complete[-1])
    windows = {name: mask for name,mask in windows.items() if mask.any()}
    summary = {'schema':'dynsyslab.receiver-comparison/v1', 'k': k, 'units':'Pa',
               'training': str(train_path), 'held_out':str(held_path),
               'n_train':len(train['time']),'n_test':len(t),'ring_intervals':max(0,len(bounds)-1),
               'complete_observed_rings':max(0,len(complete)-1),
               'source_clock_limits':[float(t[0]),float(t[-1])], 'cases':{}, 'term_levels':{}}
    for role,d in [('training',train),('held_out',held)]:
        summary['term_levels'][role]={}
        for key in TERMS:
            a=d[key]
            summary['term_levels'][role][key]={'rms_pa':rms(a).tolist(),'peak_abs_pa':abs(a).max(axis=0).tolist(),
                                               'exactly_zero':not bool(np.any(a))}

    traces={'time':t, 'truth':x, 'bounds':bounds, 'train_time':train['time'], 'train_state':train['state']}
    cache={}; ring_rows=[]; linear=linear_inverse(train)
    with tempfile.TemporaryDirectory(dir=out) as work:
        work=Path(work)
        reference=work/'reference.npz'
        np.savez(reference,time=t,state=x,ring_bounds=bounds,error_center=[1.,0.])

        def score(predict_file):
            result=work/'rings.npz'
            subprocess.run([exe,'rings',reference,predict_file,result],check=True)
            with np.load(result) as z:
                stats={key:z[key] for key in ('count','partial','relative_rms_percent','error_sum_squares','reference_sum_squares')}
            result.unlink()
            return stats

        for case,indices in CASES.items():
            ft,ct=measurements(train,indices); fq,cq=measurements(held,indices)
            shared=float(rms(ft[:,0])); scale_c=float(rms(ct))
            if shared<=0: raise ValueError('No training signal')
            present=scale_c>64*np.finfo(float).eps*shared
            rec={'terms':[TERMS[i] for i in indices], 'contrast_resolved':bool(present),
                 'estimates':{},'contrast_rms_Pa_m':float(rms(cq)),
                 'equivalent_per_channel_noise_rms_Pa':float(rms(cq)/np.linalg.norm(held['ranges_m']))}
            if present:
                ti=t-held['lag_tau'][1];ok=ti>=t[0]
                sampled=fq[ok,1]-np.interp(ti[ok],t,fq[:,0])
                rec['sampled_delay_interpolation_error_to_contrast_rms']=float(rms(sampled-cq[ok])/rms(cq[ok]))
            all_features=[('one',ft[:,:1],fq[:,:1],[shared]),('two',ft,fq,[shared,shared])]
            if present:
                all_features.append(('aligned',np.column_stack((ft,ct)),np.column_stack((fq,cq)),[shared,shared,scale_c]))
            all_features.append(('linear',ft,fq,[]))
            for mode,F,Q,scales in all_features:
                key=(mode,digest(F,Q,np.asarray(scales)))
                if key not in cache:
                    dest=work/'est.npz'
                    if mode=='linear':
                        pred=Q@linear.T+[1.,0.]
                        np.savez(dest,time=t,prediction=pred)
                    else:
                        np.savez(work/'train.npz',features=F,state=train['state'],scales=scales)
                        np.savez(work/'query.npz',features=Q,time=t)
                        subprocess.run([exe,'estimate',work/'train.npz',work/'query.npz',str(k),dest],check=True)
                        with np.load(dest) as z:
                            pred=z['prediction']
                            if case=='monopole' and mode=='two':
                                traces['neighbor_rows']=z['neighbor_rows'];traces['weights']=z['weights']
                                traces['nearest_distance']=z['nearest_distance']
                    scores={name:errors(x,pred,mask) for name,mask in windows.items()}
                    per_ring=score(dest) if len(bounds)>=2 else None
                    dest.unlink()
                    cache[key]=(case,scores,per_ring)
                    traces[case+'_'+mode]=pred
                canonical,scores,per_ring=cache[key]
                rec['estimates'][mode]={'computed_as':canonical,'observation_digest':key[1],'scores':scores}
                for i,(lo,hi) in enumerate(zip(bounds[:-1],bounds[1:])):
                    ring_rows.append((case,mode,i+1,lo,hi,int(per_ring['count'][i]),int(per_ring['partial'][i]),
                                      *per_ring['relative_rms_percent'][i]))
            if not present:
                rec['estimates']['aligned']={'computed_as':'two','reason':'common-source contrast below roundoff threshold'}
            if case in ('monopole','near_field'):
                traces[case+'_observations']=fq;traces[case+'_contrast']=cq
            summary['cases'][case]=rec
    with (out/'by-ring.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['case','estimator','ring','start_tau','end_tau','samples','partial','r_minus_1_error_pct','velocity_error_pct'])
        w.writerows(ring_rows)
    np.savez_compressed(out/'diagnostics.npz',**traces)
    (out/'results.json').write_text(json.dumps(summary,indent=2))


def run(exe,train_path,held_path,out,k):
    if out.exists():
        raise FileExistsError(f'Study destination already exists: {out}')
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as staging:
        ready=Path(staging)/'run';ready.mkdir()
        _compute(exe,train_path,held_path,ready,k)
        ready.rename(out)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('exe','train','held','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--neighbors',type=int,default=8)
    a=parser.parse_args();run(a.exe.resolve(),a.train.resolve(),a.held.resolve(),a.out.resolve(),a.neighbors)
