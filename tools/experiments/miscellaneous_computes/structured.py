"""Explicit second-order SINDy constructions in decoded delay coordinates."""
import numpy as np
from .history import matrix
from .features import prepare_library, library, thresholded_least_squares
from .models import SparseField, records_checked
from .sampling import trapezoid_weights

def fit_map(samples, *, family='polynomial', degree=2, centers=24, harmonics=5,
            ridge=1e-6, threshold=.003, time_weighted=False, baseline=None,rational_scales=(.25,1.,4.)):
    """Fit Y=Theta(X)C; an explicit baseline produces a residual model with that baseline included."""
    samples=list(samples)
    if not samples:raise ValueError('Empty map training set')
    from .sampling import clock
    data=[(clock(t),matrix(x),matrix(y)) for t,x,y in samples]
    if any(len(t)!=len(x) or len(x)!=len(y) or x.shape[1]!=data[0][1].shape[1] or y.shape[1]!=data[0][2].shape[1]
           for t,x,y in data):raise ValueError('Map training dimensions disagree')
    spec=prepare_library(np.concatenate([x for _,x,_ in data]),family,degree,centers,harmonics,rational_scales)
    gram=0.
    for t,x,y in data:
        w=np.sqrt(trapezoid_weights(t)) if time_weighted else np.ones(len(t))
        for i in range(0,len(t),1024):
            sl=slice(i,i+1024);target=y[sl]
            if baseline is not None:
                ref=np.asarray(baseline(x[sl]))
                if ref.shape!=target.shape:raise ValueError('Baseline output shape mismatch')
                target=target-ref
            rows=np.column_stack([library(x[sl],spec),target])*w[sl,None]
            gram+=rows.T@rows
    model=SparseField(spec,thresholded_least_squares(gram,data[0][2].shape[1],ridge,threshold))
    return SumMap([baseline,model]) if baseline is not None else model

class SumMap:
    """Evaluate the sum of a baseline and fitted residual maps, checking their output dimensions."""
    def __init__(self, fields):
        self.fields=list(fields)
        if not self.fields:raise ValueError('Empty sum of fields')
    def __call__(self,x):
        base=np.asarray(self.fields[0](x),float).copy()
        for field in self.fields[1:]:
            value=np.asarray(field(x))
            if value.shape!=base.shape:raise ValueError('Incompatible residual map')
            base+=value
        return base

class EndpointAcceleration:
    """Use a specifically supplied full two-tap deployment as an endpoint acceleration map."""
    def __init__(self,model):
        if model.chart.encoder.shape!=(4,4) or len(model.chart.lags)!=2 or len(model.chart.scales)!=2:
            raise ValueError('Expected an invertible two-tap value/rate reference model')
        self.model=model
    def __call__(self,x):
        x=np.asarray(x,float)
        if x.shape[-1]!=4:raise ValueError('Expected four endpoint coordinates')
        xi=self.model.chart.project(x.reshape(*x.shape[:-1],2,2))
        return self.model.chart.decode(self.model.rhs(xi))[...,1]

class SecondOrderField:
    """One pressure/velocity, position/velocity, or other value/rate pair per tap."""
    def __init__(self,acceleration):
        self.acceleration=acceleration
    def __call__(self,x):
        b=np.atleast_2d(np.asarray(x,float));v=np.asarray(self.acceleration(b))
        if b.shape[1]%2 or v.shape!=b[:,1::2].shape:raise ValueError('Acceleration shape mismatch')
        result=np.empty_like(b);result[:,::2]=b[:,1::2];result[:,1::2]=v
        return result[0] if np.ndim(x)==1 else result

    @classmethod
    def fit(cls,records,**options):
        data=records_checked(records)
        if data[0]['b'].shape[1]%2:raise ValueError('Expected value/rate pairs')
        return cls(fit_map([(r['time'],r['b'],r['db'][:,1::2]) for r in data],**options))

class PairNetwork:
    """Time-translation sharing by actual lag separation, with explicit directed mixing weights."""
    def __init__(self,lags,weights,pairs,fields):
        self.lags=np.asarray(lags,float);self.weights=np.asarray(weights,float);n=len(self.lags)
        if (self.lags.ndim!=1 or n<2 or self.lags[0]!=0 or not np.isfinite(self.lags).all()
            or np.any(np.diff(self.lags)<=0) or self.weights.shape!=(n,n) or not np.isfinite(self.weights).all()
            or np.any(self.weights<0) or np.any(np.diag(self.weights)!=0) or not np.allclose(self.weights.sum(axis=1),1)):
            raise ValueError('Expected increasing delays from zero and normalized off-diagonal weights')
        self.fields=list(fields);self.pairs=[tuple(edge) for edge in pairs]
        if (len({(i,j) for i,j,_ in self.pairs})!=len(self.pairs) or
            any(not (0<=i<j<n and 0<=k<len(self.fields)) for i,j,k in self.pairs) or
            {(i,j) for i,j,_ in self.pairs}!={(i,j) for i in range(n) for j in range(i+1,n) if self.weights[i,j] or self.weights[j,i]}):
            raise ValueError('Missing, repeated or invalid active pair')
    def __call__(self,x):
        a=np.atleast_2d(np.asarray(x,float))
        if a.shape[1]!=2*len(self.lags):raise ValueError('Delay network state dimension mismatch')
        out=np.zeros_like(a);out[:,::2]=a[:,1::2]
        for i,j,k in self.pairs:
            v=self.fields[k](a[:,[2*i,2*i+1,2*j,2*j+1]])
            if v.shape!=(len(a),2):raise ValueError('Expected two learned accelerations')
            out[:,2*i+1]+=self.weights[i,j]*v[:,0]
            out[:,2*j+1]+=self.weights[j,i]*v[:,1]
        return out[0] if np.ndim(x)==1 else out
    @classmethod
    def fit(cls,records,lags,weights,*,tolerance=1e-12,**options):
        data=records_checked(records);lags=np.asarray(lags,float);weights=np.asarray(weights,float);n=len(lags)
        if weights.shape!=(n,n) or data[0]['b'].shape[1]!=2*n or not np.isfinite(tolerance) or tolerance<0:
            raise ValueError('Inconsistent pair network design')
        gaps=[];samples=[];edges=[]
        for i in range(n):
            for j in range(i+1,n):
                if not (weights[i,j] or weights[j,i]):continue
                gap=lags[j]-lags[i];k=next((k for k,d in enumerate(gaps) if abs(gap-d)<=tolerance),len(gaps))
                if k==len(gaps):gaps.append(gap);samples.append([])
                edges.append((i,j,k))
                for r in data:
                    samples[k].append((r['time'],r['b'][:,[2*i,2*i+1,2*j,2*j+1]],r['db'][:,[2*i+1,2*j+1]]))
        return cls(lags,weights,edges,[fit_map(group,**options) for group in samples])

class InteriorField:
    """Damped interiors driven by a fixed function of the evolving endpoints; optional all-tap feedback."""
    def __init__(self,outer,drive,gamma,feedback=None):
        self.outer=outer;self.drive=drive;self.gamma=np.asarray(gamma,float);self.feedback=feedback
        if not np.isfinite(self.gamma).all() or np.any(self.gamma<=0):raise ValueError('Damping must be positive')
    def __call__(self,x):
        z=np.atleast_2d(np.asarray(x,float));k=z.shape[1]//2-2
        if k<1 or z.shape[1]%2:raise ValueError('Expected two outer taps and at least one inner tap')
        ends=z[:,[0,1,-2,-1]];v=np.asarray(self.outer(ends)).copy();f=self.drive(ends)
        if v.shape!=(len(z),2) or f.shape!=(len(z),k):raise ValueError('Interior or endpoint map dimension mismatch')
        if self.feedback is not None:
            correction=self.feedback(z)
            if correction.shape!=v.shape:raise ValueError('Endpoint feedback dimension mismatch')
            v+=correction
        out=np.empty_like(z);out[:,::2]=z[:,1::2];out[:,[1,-1]]=v
        out[:,3:-1:2]=f-2*self.gamma*z[:,3:-1:2]-self.gamma**2*z[:,2:-2:2]
        return out[0] if np.ndim(x)==1 else out
    @classmethod
    def fit(cls,records,gamma,*,outer=None,drive_options=None,feedback_options=None):
        data=records_checked(records);ends=lambda x:x[:,[0,1,-2,-1]]
        gamma=np.asarray(gamma,float)
        options=dict(drive_options or {});waveform=options.pop('waveform_fit',False)
        if data[0]['b'].shape[1]<6 or data[0]['b'].shape[1]%2:raise ValueError('Insufficient delay pairs')
        outer_data=[(r['time'],ends(r['b']),r['db'][:,[1,-1]]) for r in data]
        if outer is None:outer=fit_map(outer_data,**options)
        feedback=None
        if feedback_options is not None:
            feedback=fit_map([(r['time'],r['b'],r['db'][:,[1,-1]]-outer(ends(r['b']))) for r in data],**feedback_options)
        if waveform:
            from .waveform_fit import fit_damped_drive
            if gamma.ndim:raise ValueError('Waveform fitting uses a common scalar damping rate')
            drive=fit_damped_drive([(r['time'],ends(r['b']),r['b'][:,2:-2:2],r['b'][:,3:-1:2])
                                   for r in data],float(gamma),**options)
        else:
            drive=fit_map([(r['time'],ends(r['b']),r['db'][:,3:-1:2]+2*gamma*r['b'][:,3:-1:2]+gamma**2*r['b'][:,2:-2:2])
                           for r in data],**options)
        return cls(outer,drive,gamma,feedback)
