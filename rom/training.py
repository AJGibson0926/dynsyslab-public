"""Declarative state-domain coverage plans for ROM training campaigns."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Iterable

import numpy as np


@dataclass(frozen=True)
class AxisAlignedStateDomain:
    minimum: tuple[float,...]
    maximum: tuple[float,...]
    coordinate_system: str="model"

    def __post_init__(self):
        minimum=np.asarray(self.minimum,dtype=float); maximum=np.asarray(self.maximum,dtype=float)
        if minimum.ndim!=1 or minimum.size==0 or minimum.shape!=maximum.shape:
            raise ValueError("State-domain bounds must be equally sized nonempty vectors.")
        if np.any(~np.isfinite(minimum)) or np.any(~np.isfinite(maximum)) or np.any(maximum<minimum):
            raise ValueError("State-domain bounds must be finite and ordered.")
        object.__setattr__(self,"minimum",tuple(float(value) for value in minimum))
        object.__setattr__(self,"maximum",tuple(float(value) for value in maximum))

    @property
    def dimension(self) -> int: return len(self.minimum)

    def expanded_to_include(self, other: "AxisAlignedStateDomain", *, margin_fraction: float=0.) -> "AxisAlignedStateDomain":
        if other.dimension!=self.dimension or other.coordinate_system!=self.coordinate_system:
            raise ValueError("State domains must use the same dimension and coordinate system.")
        if margin_fraction<0. or not np.isfinite(margin_fraction):
            raise ValueError("margin_fraction must be finite and nonnegative.")
        current_low=np.asarray(self.minimum); current_high=np.asarray(self.maximum)
        required_low=np.asarray(other.minimum); required_high=np.asarray(other.maximum)
        low=np.minimum(current_low,required_low); high=np.maximum(current_high,required_high)
        margin=margin_fraction*np.maximum(high-low,np.finfo(float).eps)
        low=low-np.where(required_low<current_low,margin,0.)
        high=high+np.where(required_high>current_high,margin,0.)
        return AxisAlignedStateDomain(tuple(low),tuple(high),self.coordinate_system)

    def contains(self, values: Iterable[Iterable[float]]|np.ndarray) -> np.ndarray:
        states=np.asarray(values,dtype=float)
        if states.ndim==1: states=states.reshape(1,-1)
        if states.ndim!=2 or states.shape[1]!=self.dimension:
            raise ValueError("States do not match the state-domain dimension.")
        return np.all((states>=np.asarray(self.minimum))&(states<=np.asarray(self.maximum)),axis=1)

    def summary(self) -> dict[str,Any]:
        return {"encoding":"axis-aligned-state-domain/v1","minimum":list(self.minimum),
            "maximum":list(self.maximum),"coordinate_system":self.coordinate_system}

    @classmethod
    def from_krom(cls, krom) -> "AxisAlignedStateDomain":
        domain=krom.training_domain
        if domain is None: raise ValueError("KROM has no archived training domain.")
        return cls(tuple(domain["minimum"]),tuple(domain["maximum"]),
            str(domain.get("coordinate_system") or krom.coordinate_system))


@dataclass(frozen=True)
class StateCoveragePlan:
    current: AxisAlignedStateDomain
    required: AxisAlignedStateDomain
    planned: AxisAlignedStateDomain
    initial_conditions: tuple[tuple[float,...],...]
    points_per_axis: tuple[int,...]
    outside_current_only: bool

    def summary(self) -> dict[str,Any]:
        points=np.asarray(self.initial_conditions,dtype=float)
        return {"method":"cartesian-state-domain-expansion/v1","current":self.current.summary(),
            "required":self.required.summary(),"planned":self.planned.summary(),
            "points_per_axis":list(self.points_per_axis),"outside_current_only":self.outside_current_only,
            "initial_condition_count":len(self.initial_conditions),
            "sampled_minimum":None if not len(points) else np.min(points,axis=0).tolist(),
            "sampled_maximum":None if not len(points) else np.max(points,axis=0).tolist(),
            "required_covered":bool(np.all(self.planned.contains((self.required.minimum,self.required.maximum))))}


def plan_state_domain_expansion(current: AxisAlignedStateDomain, required: AxisAlignedStateDomain, *,
                                points_per_axis: Iterable[int], margin_fraction: float=.1,
                                outside_current_only: bool=True,
                                anchors: Iterable[Iterable[float]]=()) -> StateCoveragePlan:
    """Create a bounded Cartesian IC plan; never schedules or solves it."""
    planned=current.expanded_to_include(required,margin_fraction=margin_fraction)
    counts=tuple(int(value) for value in points_per_axis)
    if len(counts)!=planned.dimension or any(value<2 for value in counts):
        raise ValueError("points_per_axis must provide at least two points for every state axis.")
    total=int(np.prod(counts,dtype=np.int64))
    if total>100000: raise ValueError("Cartesian coverage plans are capped at 100,000 points.")
    axes=[np.linspace(low,high,count) for low,high,count in zip(planned.minimum,planned.maximum,counts)]
    samples=np.asarray(list(product(*axes)),dtype=float)
    if outside_current_only: samples=samples[~current.contains(samples)]
    anchor_values=np.asarray(list(anchors),dtype=float)
    if anchor_values.size:
        if anchor_values.ndim!=2 or anchor_values.shape[1]!=planned.dimension:
            raise ValueError("Coverage anchors do not match the state dimension.")
        if not np.all(planned.contains(anchor_values)):
            raise ValueError("Coverage anchors must lie inside the planned domain.")
        samples=np.vstack((samples,anchor_values))
    if not samples.size: raise ValueError("Coverage plan produced no new initial conditions.")
    samples=np.unique(samples,axis=0)
    return StateCoveragePlan(current,required,planned,tuple(tuple(float(value) for value in row) for row in samples),
        counts,bool(outside_current_only))


__all__=["AxisAlignedStateDomain","StateCoveragePlan","plan_state_domain_expansion"]
