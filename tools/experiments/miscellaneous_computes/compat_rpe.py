"""Compatibility for the sampled pressure dictionaries in rpe-full-delay-sindy/v1 archives."""
from fractions import Fraction
import json
import numpy as np
from scipy.spatial.distance import cdist
from .compact_field import CompactFeatures
from .structured import PairNetwork

def channels(x):
    """Convert archived source-equivalent pressure and rate pairs to the fitted RPE dictionary variables."""
    pressure=x[:,::2]
    p=.5*(np.hypot(1.,pressure)+pressure)
    u=x[:,1::2]/(6*(1+p*p))
    z=np.empty_like(x);z[:,::2]=np.arcsinh(pressure);z[:,1::2]=u/np.hypot(1.,u/3.)
    return z,40*(1+p**3)

def rational(p,u):
    """Evaluate the historical regularized pressure and rate interaction dictionary."""
    r=[np.ones_like(p),p/np.hypot(1.,p)]
    for v in (.25,.5,1,2,4,8,16,32):
        den=v*v+p*p;r.extend([v*v/den,v*p/den])
    b=np.column_stack(r)
    return np.stack([b-np.tile([1.,0.],9),b*u[:,None],b*(u*u)[:,None]],axis=2).reshape(len(p),-1)

def localized(z,c,width,zero):
    """Evaluate the origin-referenced local rational functions in the archived pressure basis."""
    weights=(1+cdist(z,c,'sqeuclidean')*width)**-2
    weights/=weights.sum(axis=1,keepdims=True)
    return np.concatenate([(weights-zero)[:,:,None],weights[:,:,None]*z[:,None,:]],axis=2).reshape(len(z),-1)

class PressurePair:
    """Evaluate both fitted acceleration rows of the archived two-tap pressure field."""
    def __init__(self,spec):
        self.mon=[[i for i in powers if i>=0] for powers in spec['mon']]
        self.C=np.asarray(spec['coef']);self.c=np.asarray(spec['centers'])
        self.width=np.asarray(spec['inv_width2']);self.zero=np.asarray(spec['b0'])
    def __call__(self,x):
        rows=[]
        for j in (0,1):
            a=np.roll(x,-2*j,axis=1);z,g=channels(a)
            cols=[np.prod(z[:,p],axis=1) for p in self.mon]+[rational(a[:,0],z[:,1])]
            if len(self.c):cols.append(localized(z,self.c,self.width,self.zero))
            rows.append(g[:,0]*(np.column_stack(cols)@self.C[j]))
        return np.column_stack(rows)

class HistoricalPressureField:
    """Replay an archived multi-tap law using its original pressure normalization and field family."""
    def __init__(self,law,lags):
        self.law=law;self.lags=np.asarray(lags,float);self.family=law['family']
        if self.family=='pairwise':
            K=len(lags);W=np.asarray(law['config']['weights']);edges=[];groups={};fields=[]
            for i in range(K):
                for j in range(i+1,K):
                    if not (W[i,j] or W[j,i]):continue
                    gap=Fraction(j-i,K-1);key=f'{gap.numerator}-{gap.denominator}'
                    if key not in groups:groups[key]=len(fields);fields.append(PressurePair(law['pairs'][key]))
                    edges.append((i,j,groups[key]))
            self.network=PairNetwork(lags,W,edges,fields)
        elif self.family in ('joint','local forcing','waveform'):
            self.outer=PressurePair(law['endpoint_field']) if law.get('endpoint_field') is not None else None
            if self.family=='local forcing':
                self.bank=CompactFeatures(law['centers'],int(law['neighbors']));self.force=np.asarray(law['force'])
        else:raise ValueError(f'Unknown historical pressure field {self.family}')
    def __call__(self,x):
        if self.family=='pairwise':return self.network(x)
        a=np.atleast_2d(np.asarray(x,float));s=self.law;z,g=channels(a);ends=a[:,[0,1,-2,-1]]
        out=np.empty_like(a);out[:,::2]=a[:,1::2];acc=out[:,1::2];i,j=np.triu_indices(a.shape[1])
        if self.family=='joint':
            b=s['bank'];parts=[z,z[:,i]*z[:,j]]
            parts.extend(rational(a[:,k],z[:,k+1]) for k in range(0,a.shape[1],2))
            parts.append(localized(z,b['c'],b['iw'],b['zero']))
            acc[:]=(np.column_stack(parts)@np.asarray(s['coef']).T)*g
            if self.outer is not None:acc[:,[0,-1]]+=self.outer(ends)
        else:
            bounded=z.copy();bounded[:,::2]=a[:,::2]/np.hypot(1.,a[:,::2]);bounded[:,1::2]/=3
            phi=np.column_stack([bounded,bounded[:,i]*bounded[:,j]])
            acc[:,[0,-1]]=self.outer(ends)+g[:,[0,-1]]*(phi@np.asarray(s['feedback']).T)
            key=z[:,[0,1,-2,-1]].copy();gamma=s['config']['gamma']
            if self.family=='local forcing':F=self.bank(key)@self.force
            else:
                key[:,::2]=np.tanh(key[:,::2]/3);key[:,1::2]/=3;i,j=np.triu_indices(4)
                phi=np.column_stack([np.ones(len(key)),key,key[:,i]*key[:,j],np.exp(-cdist(key,s['centers'],'sqeuclidean')*s['width'])])
                F=gamma**2*(phi@np.asarray(s['coef']).T)
            acc[:,1:-1]=F-2*gamma*a[:,3:-1:2]-gamma**2*a[:,2:-2:2]
        return out[0] if np.ndim(x)==1 else out
