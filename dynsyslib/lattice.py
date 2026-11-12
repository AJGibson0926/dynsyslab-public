"""Statistical-mechanics backends, intentionally separate from ODE chaos tools."""
from __future__ import annotations
from dataclasses import dataclass, field as dataclass_field
import numpy as np
from .interfaces import LatticeSystem

@dataclass
class IsingModel(LatticeSystem):
    """2-D square-lattice Ising model with periodic boundary and Metropolis dynamics."""
    size: int=32; coupling: float=1.; field: float=0.; temperature: float=2.269; seed: int|None=None; name: str="ising"
    spins: np.ndarray=dataclass_field(init=False); rng: np.random.Generator=dataclass_field(init=False)
    def __post_init__(self):
        self.rng=np.random.default_rng(self.seed); self.spins=self.rng.choice((-1,1), size=(self.size,self.size))
    def energy(self):
        s=self.spins
        return float(-self.coupling*np.sum(s*(np.roll(s,-1,0)+np.roll(s,-1,1)))-self.field*np.sum(s))
    def magnetization(self): return float(self.spins.mean())
    def sweep(self):
        for _ in range(self.size*self.size):
            i,j=self.rng.integers(self.size,size=2); spin=self.spins[i,j]
            neighbours=self.spins[(i-1)%self.size,j]+self.spins[(i+1)%self.size,j]+self.spins[i,(j-1)%self.size]+self.spins[i,(j+1)%self.size]
            delta=2*spin*(self.coupling*neighbours+self.field)
            if delta<=0 or self.rng.random()<np.exp(-delta/self.temperature): self.spins[i,j]=-spin
    def run(self, sweeps: int):
        values=np.empty((sweeps,2))
        for k in range(sweeps): self.sweep(); values[k]=(self.magnetization(),self.energy()/(self.size*self.size))
        return values

@dataclass
class Ising3DModel(LatticeSystem):
    """3-D cubic Ising model with periodic boundaries and Metropolis dynamics."""
    size: int=16; coupling: float=1.; field: float=0.; temperature: float=4.5115; seed: int|None=None; name: str="ising-3d"
    spins: np.ndarray=dataclass_field(init=False); rng: np.random.Generator=dataclass_field(init=False)
    def __post_init__(self):
        self.rng=np.random.default_rng(self.seed); self.spins=self.rng.choice((-1,1),size=(self.size,self.size,self.size))
    def energy(self):
        s=self.spins
        bonds=np.roll(s,-1,0)+np.roll(s,-1,1)+np.roll(s,-1,2)
        return float(-self.coupling*np.sum(s*bonds)-self.field*np.sum(s))
    def magnetization(self): return float(self.spins.mean())
    def sweep(self):
        for _ in range(self.size**3):
            i,j,k=self.rng.integers(self.size,size=3); spin=self.spins[i,j,k]
            neighbours=(self.spins[(i-1)%self.size,j,k]+self.spins[(i+1)%self.size,j,k]+
                        self.spins[i,(j-1)%self.size,k]+self.spins[i,(j+1)%self.size,k]+
                        self.spins[i,j,(k-1)%self.size]+self.spins[i,j,(k+1)%self.size])
            delta=2*spin*(self.coupling*neighbours+self.field)
            if delta<=0 or self.rng.random()<np.exp(-delta/self.temperature): self.spins[i,j,k]=-spin
    def run(self,sweeps:int):
        values=np.empty((sweeps,2))
        for index in range(sweeps): self.sweep(); values[index]=(self.magnetization(),self.energy()/(self.size**3))
        return values

@dataclass
class PottsModel(LatticeSystem):
    """2-D q-state Potts model with periodic boundaries and Metropolis dynamics."""
    size: int=32; states: int=3; coupling: float=1.; field: float=0.; temperature: float=1.; seed: int|None=None; name: str="potts"
    spins: np.ndarray=dataclass_field(init=False); rng: np.random.Generator=dataclass_field(init=False)
    def __post_init__(self):
        if self.states<2: raise ValueError("Potts model requires at least two states.")
        self.rng=np.random.default_rng(self.seed); self.spins=self.rng.integers(self.states,size=(self.size,self.size))
    def energy(self):
        s=self.spins
        # Cast before addition: NumPy boolean addition is logical rather than
        # integer-valued and would under-count sites aligned in both directions.
        aligned=(s==np.roll(s,-1,0)).astype(int)+(s==np.roll(s,-1,1)).astype(int)
        return float(-self.coupling*np.sum(aligned)-self.field*np.sum(s==0))
    def order_parameter(self):
        fractions=np.bincount(self.spins.ravel(),minlength=self.states)/(self.size*self.size)
        return float((self.states*fractions.max()-1)/(self.states-1))
    def sweep(self):
        for _ in range(self.size*self.size):
            i,j=self.rng.integers(self.size,size=2); old=self.spins[i,j]; new=self.rng.integers(self.states-1); new+=new>=old
            neighbours=np.array((self.spins[(i-1)%self.size,j],self.spins[(i+1)%self.size,j],self.spins[i,(j-1)%self.size],self.spins[i,(j+1)%self.size]))
            field_change=int(new==0)-int(old==0)
            delta=-self.coupling*(np.sum(neighbours==new)-np.sum(neighbours==old))-self.field*field_change
            if delta<=0 or self.rng.random()<np.exp(-delta/self.temperature): self.spins[i,j]=new
    def run(self,sweeps:int):
        values=np.empty((sweeps,2))
        for k in range(sweeps): self.sweep(); values[k]=(self.order_parameter(),self.energy()/(self.size*self.size))
        return values
