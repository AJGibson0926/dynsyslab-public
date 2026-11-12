"""Compare sensor arrays and delayed copies with aligned targets on one common interval."""
import numpy as np
from .history import matrix
from .sampling import clock,delay_samples
from .diagnostics import errors,forecast
from .delay_dmd import LinearReadout,DiscreteDelayModel
from .readouts import NeighborReadout

def sensor_features(time,observations,ends,design):
    """Each feature declares a channel, nonnegative physical delay and optional unit conversion."""
    y=matrix(observations)
    if not design:raise ValueError('Empty sensor design')
    delays=[term['delay'] for term in design]
    # Gather one channel at a time to avoid allocating channels × all taps.
    columns=[]
    for term,d in zip(design,delays):
        c=term['channel'];scale=term.get('scale',1.)
        if not isinstance(c,int) or not 0<=c<y.shape[1] or not np.isfinite(scale):raise ValueError('Invalid receiver channel or scale')
        columns.append(delay_samples(time,y[:,c],ends,[d])[:,0,0]*scale)
    return np.column_stack(columns)

def compare_reconstructions(training,test,designs,*,neighbors=8,train_stride=1,steps_ahead=False):
    """Inputs are ID -> (time, measurements, state labels); free DMD uses a single query window."""
    if not training or not test or not designs or set(training)&set(test) or not isinstance(train_stride,int) or train_stride<1:
        raise ValueError('Expected disjoint nonempty training and test sets and a positive training stride')
    span=max(term['delay'] for design in designs.values() for term in design)
    out={}
    for name,design in designs.items():
        fitted=[]
        def prepare(record):
            t=clock(record[0]);x=matrix(record[1]);y=matrix(record[2])
            if len(y)!=len(t) or len(x)!=len(t):raise ValueError('Misaligned receiver or state record')
            i=np.flatnonzero(t>=t[0]+span)
            if len(i)<2:raise ValueError('Insufficient common history')
            return t[i],sensor_features(t,x,t[i],design),y[i]
        for r in training.values():fitted.append(prepare(r))
        pairs=[(x[::train_stride],y[::train_stride]) for _,x,y in fitted]
        linear=LinearReadout.fit(pairs)
        x=np.vstack([a for a,_ in pairs]);y=np.vstack([a for _,a in pairs])
        rms=np.sqrt(np.mean(x*x,axis=0));lookup=NeighborReadout(x,y,np.where(rms>0,rms,1.))
        step=DiscreteDelayModel.fit([(t,x) for t,x,_ in fitted],stride=train_stride) if steps_ahead else None
        cases={}
        for tag,r in test.items():
            t,q,ref=prepare(r);a=linear(q)
            nearest=np.vstack([lookup.predict(q[i:i+4096],neighbors) for i in range(0,len(q),4096)])
            curves=dict(time=t,actual=ref,linear=a,neighbor=nearest)
            scores=dict(linear=errors(ref,a),neighbor=errors(ref,nearest))
            if step is not None:
                count=int(np.floor((t[-1]-t[0])/step.dt+1e-10))
                free=step.rollout(q[0],count,start=t[0]);pred=linear(free['b'])
                curves.update(forecast_time=free['time'],forecast=pred)
                scores['autonomous']=forecast(t,ref,free['time'],pred)
            cases[tag]=dict(scores=scores,curves=curves)
        out[name]=dict(training_scales=lookup.scale.tolist(),spectral_radius=step.spectral_radius if step else None,cases=cases)
    return out
