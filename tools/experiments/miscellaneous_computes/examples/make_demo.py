"""Generate three synthetic multiharmonic sensor records and runnable toolkit manifests."""
from pathlib import Path
import argparse,json
import numpy as np

def make(folder):
    """Create independent harmonic signal records and compatible study, reconstruction and Viewer manifests."""
    folder.mkdir(parents=True,exist_ok=False)
    t=np.arange(0,12,.02)
    for tag,phase,scale in [('train',0.,1.),('validate',.37,.88),('test',.61,1.15)]:
        u=t+phase
        y=scale*(np.sin(u)+.2*np.sin(2*u)+.08*np.cos(3*u))
        v=scale*(np.cos(u)+.4*np.cos(2*u)-.24*np.sin(3*u))
        a=scale*(-np.sin(u)-.8*np.sin(2*u)-.72*np.cos(3*u))
        np.savez_compressed(folder/(tag+'.npz'),time=t,values=np.column_stack([y,v]),
                            rates=np.column_stack([v,a]),state=np.column_stack([np.sin(u),np.cos(u)]))
    splits=dict(training={'A':'train.npz'},validation={'B':'validate.npz'},test={'C':'test.npz'})
    study=dict(schema='dynsyslab.miscellaneous-computes/modal-study/v1',**splits,experiment='direct',
        geometries=[dict(id='two-taps',lags=[0.,1.1]),dict(id='three-taps',lags=[0.,.55,1.1])],
        candidates=[dict(family='polynomial',degree=1,ridge=1e-10,threshold=1e-10,time_weighted=True)],
        interpolation='quintic_pair',value_channel=0,rate_channel=1,time_unit='s',channel_units=['value','value/s'],
        observation_provenance={tag:dict(values='Analytic harmonic signal and its rate',
                                         rates='Analytic first and second time derivatives') for tag in ('A','B','C')},
        peak_prominence=.15,solver=dict(method='exact-affine'))
    recon=dict(schema='dynsyslab.miscellaneous-computes/reconstruction/v1',training=splits['training'],test=splits['test'],
        designs={'present':[dict(channel=0,delay=0.)],'history':[dict(channel=0,delay=0.),dict(channel=0,delay=.55),dict(channel=0,delay=1.1)]},
        options=dict(neighbors=8,steps_ahead=True,train_stride=2))
    cases={}
    for role,tag in [('Training','A'),('Held out','C')]:
        curves={name:dict(file=f'results/{name}/curves/{tag}.npz',values='predicted') for name in ['two-taps','three-taps']}
        curves['reference']=dict(file=f'results/three-taps/curves/{tag}.npz',time='reference_time',values='reference')
        cases[role]=curves
    plot=dict(schema='dynsyslab.miscellaneous-computes/comparison/v1',models=['three-taps'],cases=cases,
              options=dict(baseline='two-taps',labels=['Observable','Rate (per second)'],expanded=False))
    for name,obj in [('study',study),('sensors',recon),('plot',plot)]:
        (folder/(name+'.json')).write_text(json.dumps(obj,indent=2),encoding='utf-8')
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('folder',type=Path);make(ap.parse_args().folder)
