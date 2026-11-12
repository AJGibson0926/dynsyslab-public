"""Discrete delay-coordinate forecast and supervised state readout on a frozen basis."""
import numpy as np
from .sampling import clock
from .history import matrix

class LinearReadout:
    """Least-squares map from a row of measured or modal features to one or more labeled outputs."""
    def __init__(self,weights,bias):
        self.weights=matrix(weights);self.bias=np.asarray(bias,float)
        if self.bias.shape!=(self.weights.shape[1],) or not np.isfinite(self.bias).all():raise ValueError('Invalid readout bias')
    def __call__(self,z):return np.asarray(z)@self.weights+self.bias

    @classmethod
    def fit(cls,records,*,affine=True,ridge=0.,rcond=None):
        pairs=[(matrix(x),matrix(y)) for x,y in records]
        if not pairs or any(len(x)!=len(y) for x,y in pairs):raise ValueError('Expected aligned state labels')
        x=np.vstack([a for a,_ in pairs]);y=np.vstack([a for _,a in pairs])
        xm=x.mean(axis=0) if affine else np.zeros(x.shape[1])
        ym=y.mean(axis=0) if affine else np.zeros(y.shape[1])
        if not np.isfinite(ridge) or ridge<0:raise ValueError('Invalid ridge')
        X=x-xm;Y=y-ym
        if ridge:
            X=np.vstack([X,np.sqrt(ridge)*np.eye(x.shape[1])]);Y=np.vstack([Y,np.zeros((x.shape[1],y.shape[1]))])
        B=np.linalg.lstsq(X,Y,rcond=rcond)[0]
        return cls(B,ym-xm@B)

class DiscreteDelayModel:
    """Frozen discrete delay-state transition with a physical sampling interval and a free rollout."""
    def __init__(self,advance,dt):
        if advance.weights.shape[0]!=advance.weights.shape[1] or not np.isfinite(dt) or dt<=0:raise ValueError('Invalid one-step map or clock')
        self.advance=advance;self.dt=float(dt)
    @classmethod
    def fit(cls,records,*,stride=1,affine=False,ridge=0.,rcond=None):
        """Use only within-record snapshot pairs. Input rows are frozen delay coordinates."""
        if not isinstance(stride,int) or stride<1:raise ValueError('Invalid advance stride')
        pairs=[];steps=[]
        for time,values in records:
            t=clock(time);z=matrix(values)
            if len(z)!=len(t) or len(t)<=stride:raise ValueError('Not enough aligned snapshots')
            dt=t[stride:]-t[:-stride];steps.append(dt[0])
            if not np.allclose(dt,steps[0],rtol=1e-8,atol=steps[0]*1e-10):raise ValueError('Incompatible shift clocks')
            pairs.append((z[:-stride],z[stride:]))
        if not steps:raise ValueError('No training trajectories')
        return cls(LinearReadout.fit(pairs,affine=affine,ridge=ridge,rcond=rcond),steps[0])
    @property
    def spectral_radius(self):
        return float(np.max(abs(np.linalg.eigvals(self.advance.weights))))
    def rollout(self,initial,steps,*,start=0.,limit=np.inf):
        x=np.asarray(initial,float)
        if (x.shape!=(len(self.advance.weights),) or not np.isfinite(x).all() or not isinstance(steps,int)
            or steps<0 or not np.isfinite(start) or not limit>0):raise ValueError('Invalid DMD rollout request')
        rows=[]
        for i in range(steps+1):
            if not np.isfinite(x).all() or np.max(abs(x))>=limit:break
            rows.append(x)
            if i<steps:x=self.advance(x)
        return dict(time=start+self.dt*np.arange(len(rows)),b=np.asarray(rows).reshape(-1,len(x)),complete=len(rows)==steps+1)
