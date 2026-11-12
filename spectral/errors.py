"""Neutral free-trajectory validation for Koopman eigenfunctions.

The fitted generator equation and the finite-time Koopman semigroup are
different validation questions.  This module records both without coupling
the result to Viewer, a particular solver, or a particular fitting method.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from typing import Any, Callable, Literal, Sequence
from uuid import uuid4

import numpy as np

from time_series import TimeSeries


DatasetRole = Literal["training", "test"]


def _hash_array(digest, value: np.ndarray) -> None:
    array=np.ascontiguousarray(value)
    digest.update(str(array.dtype).encode()); digest.update(str(array.shape).encode())
    digest.update(array.view(np.uint8))


@dataclass(frozen=True)
class FreeTrajectoryDataset:
    """Named free trajectories used to validate, but never refit, eigenfunctions."""

    trajectories: tuple[TimeSeries,...]
    label: str
    role: DatasetRole
    coordinate_system: str = "model"
    provenance: dict[str,Any] = field(default_factory=dict)

    def __init__(self, trajectories: Sequence[TimeSeries], *, label: str,
                 role: DatasetRole, coordinate_system: str="model",
                 provenance: dict[str,Any]|None=None):
        expanded=[]
        for trajectory in trajectories:
            if not isinstance(trajectory,TimeSeries):
                raise TypeError("FreeTrajectoryDataset trajectories must be TimeSeries objects.")
            expanded.extend(trajectory.members())
        items=tuple(expanded)
        if not items or not isinstance(label,str) or not label.strip():
            raise ValueError("A labeled free-trajectory dataset requires at least one trajectory.")
        if role not in {"training","test"}:
            raise ValueError("Free-trajectory dataset role must be 'training' or 'test'.")
        if coordinate_system not in {"model","physical"}:
            raise ValueError("coordinate_system must be 'model' or 'physical'.")
        dimensions=set()
        for item in items:
            values=np.asarray(item.values)
            if values.ndim!=2 or values.shape[0]<2 or np.any(~np.isfinite(values)):
                raise ValueError("Free validation trajectories must be finite sample-by-state matrices.")
            dimensions.add(values.shape[1])
        if len(dimensions)!=1:
            raise ValueError("All free validation trajectories must share a state dimension.")
        object.__setattr__(self,"trajectories",items)
        object.__setattr__(self,"label",label.strip())
        object.__setattr__(self,"role",role)
        object.__setattr__(self,"coordinate_system",coordinate_system)
        object.__setattr__(self,"provenance",dict(provenance or {}))

    @property
    def state_dimension(self) -> int:
        return int(self.trajectories[0].values.shape[1])

    @property
    def sample_count(self) -> int:
        return int(sum(item.sample_count for item in self.trajectories))

    @property
    def fingerprint(self) -> str:
        digest=sha256()
        for trajectory in self.trajectories:
            _hash_array(digest,trajectory.time); _hash_array(digest,trajectory.values)
        digest.update(self.coordinate_system.encode())
        return digest.hexdigest()

    def summary(self) -> dict[str,Any]:
        return {"encoding":"koopman-free-validation-dataset/v1","label":self.label,
            "role":self.role,"coordinate_system":self.coordinate_system,
            "trajectory_count":len(self.trajectories),"sample_count":self.sample_count,
            "state_dimension":self.state_dimension,"fingerprint_sha256":self.fingerprint,
            "provenance":dict(self.provenance)}


@dataclass(frozen=True)
class ErrorStatistics:
    sample_pairs: int
    absolute_rms: float
    maximum_absolute: float
    relative_l2: float
    median_pointwise_relative: float
    p95_pointwise_relative: float

    def summary(self) -> dict[str,Any]:
        return {"sample_pairs":self.sample_pairs,"absolute_rms":self.absolute_rms,
            "maximum_absolute":self.maximum_absolute,"relative_l2":self.relative_l2,
            "median_pointwise_relative":self.median_pointwise_relative,
            "p95_pointwise_relative":self.p95_pointwise_relative}

    @classmethod
    def from_summary(cls, value: dict[str,Any]) -> "ErrorStatistics":
        return cls(int(value["sample_pairs"]),float(value["absolute_rms"]),
            float(value["maximum_absolute"]),float(value["relative_l2"]),
            float(value["median_pointwise_relative"]),float(value["p95_pointwise_relative"]))


@dataclass(frozen=True)
class SemigroupHorizonError:
    horizon: float
    statistics: ErrorStatistics

    def summary(self) -> dict[str,Any]:
        return {"horizon":self.horizon,"statistics":self.statistics.summary()}


@dataclass(frozen=True)
class EigenfunctionErrorEvaluation:
    """One immutable validation record attached to one eigenfunction."""

    evaluation_id: str
    dataset: dict[str,Any]
    eigenvalue: complex
    one_step: ErrorStatistics
    horizons: tuple[SemigroupHorizonError,...]
    anchored_long_term: ErrorStatistics
    generator: ErrorStatistics
    coverage: dict[str,Any]
    metadata: dict[str,Any] = field(default_factory=dict)

    def summary(self) -> dict[str,Any]:
        return {"encoding":"koopman-eigenfunction-error-evaluation/v1",
            "evaluation_id":self.evaluation_id,"dataset":dict(self.dataset),
            "eigenvalue":{"real":float(self.eigenvalue.real),"imag":float(self.eigenvalue.imag)},
            "one_step":self.one_step.summary(),
            "horizons":[item.summary() for item in self.horizons],
            "anchored_long_term":self.anchored_long_term.summary(),
            "generator":self.generator.summary(),"coverage":dict(self.coverage),
            "metadata":dict(self.metadata)}

    @classmethod
    def from_summary(cls, value: dict[str,Any]) -> "EigenfunctionErrorEvaluation":
        if value.get("encoding")!="koopman-eigenfunction-error-evaluation/v1":
            raise ValueError("Unsupported eigenfunction-error evaluation encoding.")
        eigen=value["eigenvalue"]
        return cls(str(value["evaluation_id"]),dict(value["dataset"]),
            complex(float(eigen["real"]),float(eigen["imag"])),
            ErrorStatistics.from_summary(value["one_step"]),
            tuple(SemigroupHorizonError(float(item["horizon"]),
                ErrorStatistics.from_summary(item["statistics"])) for item in value["horizons"]),
            ErrorStatistics.from_summary(value["anchored_long_term"]),
            ErrorStatistics.from_summary(value["generator"]),dict(value["coverage"]),
            dict(value.get("metadata") or {}))


def _statistics(actual: np.ndarray, predicted: np.ndarray) -> ErrorStatistics:
    actual=np.asarray(actual,dtype=complex).reshape(-1)
    predicted=np.asarray(predicted,dtype=complex).reshape(-1)
    if actual.size!=predicted.size or not actual.size:
        raise ValueError("Error statistics require equal, nonempty actual and predicted arrays.")
    error=actual-predicted; magnitude=np.abs(error)
    scale=max(float(np.max(np.abs(actual))),np.finfo(float).tiny)
    floor=np.sqrt(np.finfo(float).eps)*scale
    relative=magnitude/np.maximum(np.abs(actual),floor)
    return ErrorStatistics(int(actual.size),float(np.sqrt(np.mean(magnitude**2))),
        float(np.max(magnitude)),float(np.linalg.norm(error)/max(np.linalg.norm(actual),np.finfo(float).tiny)),
        float(np.median(relative)),float(np.quantile(relative,.95)))


def _default_horizons(trajectories: Sequence[TimeSeries]) -> tuple[float,...]:
    minimum_step=min(float(np.min(np.diff(item.time))) for item in trajectories)
    shortest=min(float(item.time[-1]-item.time[0]) for item in trajectories)
    candidates=(10*minimum_step,100*minimum_step,.01*shortest,.05*shortest,
        .1*shortest,.25*shortest,.5*shortest)
    return tuple(sorted({float(value) for value in candidates if minimum_step<value<=shortest}))


def _coverage(dataset: FreeTrajectoryDataset, training_domain: dict[str,Any]|None) -> dict[str,Any]:
    if not training_domain:
        return {"training_domain_available":False,"outside_sample_count":None,
            "outside_fraction":None,"maximum_relative_excursion":None}
    minimum=np.asarray(training_domain.get("minimum"),dtype=float)
    maximum=np.asarray(training_domain.get("maximum"),dtype=float)
    if minimum.shape!=(dataset.state_dimension,) or maximum.shape!=minimum.shape or np.any(maximum<minimum):
        raise ValueError("Training-domain bounds do not match the free validation data.")
    states=np.vstack([item.values for item in dataset.trajectories])
    scale=np.maximum(maximum-minimum,np.finfo(float).eps)
    excursion=np.maximum(np.maximum((minimum-states)/scale,(states-maximum)/scale),0.)
    outside=np.any(excursion>0.,axis=1)
    return {"training_domain_available":True,"outside_sample_count":int(np.sum(outside)),
        "outside_fraction":float(np.mean(outside)),
        "maximum_relative_excursion":float(np.max(excursion)),
        "training_minimum":minimum.tolist(),"training_maximum":maximum.tolist()}


def evaluate_eigenfunction_errors(
    eigenfunction,
    dataset: FreeTrajectoryDataset,
    *,
    horizons: Sequence[float]|None=None,
    coordinate_map: Callable[[np.ndarray],np.ndarray]|None=None,
    training_domain: dict[str,Any]|None=None,
    metadata: dict[str,Any]|None=None,
) -> EigenfunctionErrorEvaluation:
    """Evaluate generator, one-step, horizon, and anchored semigroup defects.

    Long-horizon values are interpolated at exactly ``t+horizon``.  This makes
    the result valid for nonuniform test trajectories as well as regular
    training grids.  No coefficient is refit during evaluation.
    """
    if coordinate_map is None and getattr(eigenfunction,"coordinate_observable",None) is not None:
        transform=eigenfunction.coordinate_observable.evaluate
    else:
        transform=(lambda state:state) if coordinate_map is None else coordinate_map
    lifted=[]
    for trajectory in dataset.trajectories:
        lifted.append(np.asarray([eigenfunction.psi(transform(row)) for row in trajectory.values],dtype=complex))
    eigenvalue=complex(eigenfunction.eigenvalue.continuous_time)

    step_actual=[]; step_predicted=[]; anchored_actual=[]; anchored_predicted=[]
    generator_actual=[]; generator_predicted=[]
    for trajectory,psi in zip(dataset.trajectories,lifted):
        dt=np.diff(trajectory.time)
        step_actual.append(psi[1:]); step_predicted.append(np.exp(eigenvalue*dt)*psi[:-1])
        elapsed=trajectory.time-trajectory.time[0]
        anchored_actual.append(psi[1:]); anchored_predicted.append(np.exp(eigenvalue*elapsed[1:])*psi[0])
        edge_order=2 if psi.size>=3 else 1
        dpsi=np.gradient(psi,trajectory.time,edge_order=edge_order)
        generator_actual.append(dpsi); generator_predicted.append(eigenvalue*psi)
    one_step=_statistics(np.concatenate(step_actual),np.concatenate(step_predicted))
    anchored=_statistics(np.concatenate(anchored_actual),np.concatenate(anchored_predicted))
    generator_actual_values=np.concatenate(generator_actual)
    generator_predicted_values=np.concatenate(generator_predicted)
    generator=_statistics(generator_actual_values,generator_predicted_values)
    generator_residual=generator_actual_values-generator_predicted_values
    generator_relative_to_lambda_psi=float(np.linalg.norm(generator_residual)/max(
        np.linalg.norm(generator_predicted_values),np.finfo(float).tiny))

    requested=_default_horizons(dataset.trajectories) if horizons is None else tuple(float(item) for item in horizons)
    if any(not np.isfinite(item) or item<=0. for item in requested):
        raise ValueError("Semigroup horizons must be finite and positive.")
    horizon_results=[]
    for horizon in sorted(set(requested)):
        actual=[]; predicted=[]
        for trajectory,psi in zip(dataset.trajectories,lifted):
            targets=trajectory.time+horizon
            valid=targets<=trajectory.time[-1]+32*np.finfo(float).eps*max(1.,abs(trajectory.time[-1]))
            if not np.any(valid): continue
            target=targets[valid]
            interpolated=np.interp(target,trajectory.time,psi.real)+1j*np.interp(target,trajectory.time,psi.imag)
            actual.append(interpolated); predicted.append(np.exp(eigenvalue*horizon)*psi[valid])
        if actual:
            horizon_results.append(SemigroupHorizonError(horizon,
                _statistics(np.concatenate(actual),np.concatenate(predicted))))
    return EigenfunctionErrorEvaluation(str(uuid4()),dataset.summary(),eigenvalue,one_step,
        tuple(horizon_results),anchored,generator,_coverage(dataset,training_domain),
        {"generator_derivative":"trajectory-numerical-gradient",
         "generator_relative_l2_denominator":"observed trajectory dpsi/dt",
         "generator_relative_to_lambda_psi":generator_relative_to_lambda_psi,
         "coordinate_observable":None if getattr(eigenfunction,"coordinate_observable",None) is None else eigenfunction.coordinate_observable.summary(),
         "semigroup_convention":"psi(x(t+h)) versus exp(lambda*h) psi(x(t))",
         **dict(metadata or {})})


__all__=["DatasetRole","EigenfunctionErrorEvaluation","ErrorStatistics",
    "FreeTrajectoryDataset","SemigroupHorizonError","evaluate_eigenfunction_errors"]
