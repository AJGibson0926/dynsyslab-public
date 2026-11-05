"""Neutral feature dictionaries shared by EDMD, SINDy, and executable ROMs."""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations_with_replacement
import numpy as np

class ObservableDictionary:
    def transform(self, states: np.ndarray) -> np.ndarray: raise NotImplementedError
    def gradient(self, states: np.ndarray) -> np.ndarray:
        """Dictionary Jacobian, shaped (samples, features, state_dimension)."""
        x=np.asarray(states,dtype=float); h=1e-5*np.maximum(1.,np.abs(x))
        columns=[]
        for axis in range(x.shape[1]):
            delta=np.zeros_like(x); delta[:,axis]=h[:,axis]
            # Five-point central stencil: robust for deterministic callable dictionaries.
            columns.append((-self.transform(x+2*delta)+8*self.transform(x+delta)-8*self.transform(x-delta)+self.transform(x-2*delta))/(12*h[:,axis,None]))
        return np.stack(columns,axis=2)

@dataclass(frozen=True)
class PolynomialDictionary(ObservableDictionary):
    degree: int=2
    include_constant: bool=True
    def transform(self, states):
        x=np.asarray(states); columns=[]
        if self.include_constant: columns.append(np.ones(x.shape[0]))
        for order in range(1,self.degree+1):
            for indices in combinations_with_replacement(range(x.shape[1]),order): columns.append(np.prod(x[:,indices],axis=1))
        return np.column_stack(columns)
    def gradient(self, states):
        x=np.asarray(states,dtype=float); pieces=[]
        if self.include_constant: pieces.append(np.zeros((x.shape[0],x.shape[1])))
        for order in range(1,self.degree+1):
            for indices in combinations_with_replacement(range(x.shape[1]),order):
                derivative=np.zeros((x.shape[0],x.shape[1]))
                for axis in set(indices):
                    product=np.ones(x.shape[0])
                    consumed=False
                    for item in indices:
                        if item==axis and not consumed: consumed=True
                        else: product*=x[:,item]
                    derivative[:,axis]=indices.count(axis)*product
                pieces.append(derivative)
        return np.stack(pieces,axis=1)

@dataclass(frozen=True)
class RBFDictionary(ObservableDictionary):
    centers: np.ndarray
    width: float=1.
    def transform(self, states):
        diff=np.asarray(states)[:,None,:]-np.asarray(self.centers)[None,:,:]
        return np.exp(-np.sum(diff*diff,axis=2)/(2*self.width*self.width))
    def gradient(self, states):
        x=np.asarray(states,dtype=float); diff=x[:,None,:]-np.asarray(self.centers)[None,:,:]
        values=self.transform(x)
        return values[:,:,None]*(-diff/(self.width*self.width))


__all__=["ObservableDictionary","PolynomialDictionary","RBFDictionary"]
