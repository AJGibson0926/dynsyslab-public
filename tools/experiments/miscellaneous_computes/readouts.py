"""Supervised readouts from modal states for independent-label and chain-rule controls."""
import numpy as np
from scipy.spatial import cKDTree
from .history import matrix
class NeighborReadout:
    """Estimate labels from independent queries by distance-weighting training neighbors in scaled asinh coordinates."""
    def __init__(self,states,labels,scales=None):
        b=matrix(states);self.labels=matrix(labels)
        if len(b)!=len(self.labels):raise ValueError('Readout labels must align with states')
        self.scale=np.asarray(scales) if scales is not None else np.sqrt(np.mean(b*b,axis=0))
        if self.scale.shape!=(b.shape[1],) or np.any(self.scale<=0) or not np.isfinite(self.scale).all():
            raise ValueError('Positive finite per-channel scales are required')
        self.x=np.arcsinh(b/self.scale);self.tree=cKDTree(self.x)
    def predict(self,states,k=8):
        b=matrix(states)
        if not 1<=k<=len(self.x):raise ValueError('Neighbor count exceeds training data')
        d,idx=self.tree.query(np.arcsinh(b/self.scale),list(range(1,k+1)))
        w=np.maximum(d,1e-12)**-2;w/=w.sum(axis=1,keepdims=True)
        return np.einsum('ij,ijk->ik',w,self.labels[idx])

class TangentReadout(NeighborReadout):
    """Blend locally bounded first-order pressure approximations; return (g, Dg v)."""
    def __init__(self,states,rates,outputs,output_rates,amplitude_floor=1.):
        super().__init__(states,outputs)
        b,v,dy=matrix(states),matrix(rates),matrix(output_rates)
        if b.shape!=v.shape or dy.shape!=self.labels.shape or self.labels.shape[1]!=1 or amplitude_floor<=0:
            raise ValueError('Expected one scalar output and aligned velocities')
        dx=v/np.hypot(self.scale,b)
        self.slopes=dy*dx/np.maximum((dx*dx).sum(axis=1,keepdims=True),1e-14)
        self.amplitude=np.maximum(amplitude_floor,abs(self.labels[:,0]))
    def phase(self,states,rates,k=16):
        b,v=matrix(states),matrix(rates)
        if b.shape!=v.shape or not 1<=k<len(self.x):raise ValueError('Invalid tangent readout query')
        result=np.empty((len(b),2))
        for start in range(0,len(b),400):
            sl=slice(start,start+400)
            x=np.arcsinh(b[sl]/self.scale);dx=v[sl]/np.hypot(self.scale,b[sl])
            distances,ids=self.tree.query(x,k+1)
            off=x[:,None,:]-self.x[ids]
            dvel=np.einsum('ikn,in->ik',off,dx)/np.maximum(distances,1e-14)
            radius=np.maximum(distances[:,-1:],1e-13)
            ratio=distances[:,:-1]/radius
            d_ratio=(dvel[:,:-1]-ratio*dvel[:,-1:])/radius
            ids=ids[:,:-1];off=off[:,:-1]
            w=(1-ratio)**4*(1+4*ratio)
            dw=-20*ratio*(1-ratio)**3*d_ratio
            total=w.sum(axis=1)
            ties=total<1e-25
            w[ties]=1.;dw[ties]=0.;total[ties]=k
            slope=self.slopes[ids];amp=self.amplitude[ids]
            tangent=np.tanh(np.einsum('ikj,ikj->ik',slope,off)/amp)
            approximant=self.labels[ids,0]+amp*tangent
            value=(w*approximant).sum(axis=1)/total
            local_rates=(1-tangent*tangent)*np.einsum('ikj,ij->ik',slope,dx)
            velocity=(dw*(approximant-value[:,None])+w*local_rates).sum(axis=1)/total
            result[sl]=np.column_stack([value,velocity])
        return result
