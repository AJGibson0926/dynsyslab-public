"""Small, deterministic mode-selection policies for executable KROMs."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ObservableModeSelectionPolicy:
    """Grouped greedy state-observability policy for polynomial eigenfunctions.

    A complex eigenfunction is an indivisible two-column real/imaginary group.
    Selection is performed in dictionary-feature space, so cost depends on the
    small dictionary dimension rather than the number of trajectory samples
    times the number of candidate modes.
    """

    r2_threshold: float=0.99
    max_modes: int=30
    minimum_relative_improvement: float=1e-10
    control_authority_weight: float=0.05

    def __post_init__(self):
        if not 0.<self.r2_threshold<=1.:
            raise ValueError("r2_threshold must lie in (0, 1].")
        if self.max_modes<1:
            raise ValueError("max_modes must be positive.")
        if self.minimum_relative_improvement<0.:
            raise ValueError("minimum_relative_improvement cannot be negative.")
        if self.control_authority_weight<0.:
            raise ValueError("control_authority_weight cannot be negative.")


@dataclass(frozen=True)
class ObservableSelectionStep:
    step: int
    retained_index: int
    full_spectrum_index: int
    continuous_eigenvalue: complex
    raw_residual_reduction: float
    weighted_score: float
    relative_control_authority: float
    r2_by_state: tuple[float,...]

    def summary(self) -> dict[str,Any]:
        payload=asdict(self)
        value=self.continuous_eigenvalue
        payload["continuous_eigenvalue"]={"real":float(value.real),"imag":float(value.imag)}
        return payload


@dataclass(frozen=True)
class ObservableModeSelectionReport:
    policy: ObservableModeSelectionPolicy
    selected_retained_indices: tuple[int,...]
    selected_full_spectrum_indices: tuple[int,...]
    final_r2_by_state: tuple[float,...]
    steps: tuple[ObservableSelectionStep,...]
    sample_count: int
    state_dimension: int
    candidate_count: int
    dictionary_feature_count: int
    stop_reason: str

    def summary(self) -> dict[str,Any]:
        return {
            "method":"grouped-polynomial-state-observability/v1",
            "policy":asdict(self.policy),
            "selected_retained_indices":list(self.selected_retained_indices),
            "selected_full_spectrum_indices":list(self.selected_full_spectrum_indices),
            "final_r2_by_state":list(self.final_r2_by_state),
            "steps":[item.summary() for item in self.steps],
            "sample_count":self.sample_count,
            "state_dimension":self.state_dimension,
            "candidate_count":self.candidate_count,
            "dictionary_feature_count":self.dictionary_feature_count,
            "stop_reason":self.stop_reason,
        }


def _polynomial_features(krom, states: np.ndarray) -> tuple[np.ndarray,list[np.ndarray]]:
    dictionary_metadata=dict(krom.spectrum.metadata.get("dictionary") or {})
    preprocessing=dict(krom.spectrum.metadata.get("state_preprocessing") or {})
    if dictionary_metadata.get("encoding")!="polynomial-dictionary/v1":
        raise ValueError("Observable polynomial selection requires polynomial-dictionary/v1 metadata.")
    required={"center","scale","multiplier"}
    if not required.issubset(preprocessing):
        raise ValueError("Observable polynomial selection requires complete state preprocessing.")
    x=np.asarray(states,dtype=float)
    if x.ndim!=2 or x.shape[1]!=np.asarray(krom.system.initial_state).size:
        raise ValueError("Selection states must have shape (samples, system state dimension).")
    center=np.asarray(preprocessing["center"],dtype=float)
    scale=np.asarray(preprocessing["scale"],dtype=float)
    multiplier=np.asarray(preprocessing["multiplier"],dtype=float)
    transformed=(x-center)/scale*multiplier
    from spectral.dictionary import PolynomialDictionary
    dictionary=PolynomialDictionary(degree=int(dictionary_metadata["degree"]),
        include_constant=bool(dictionary_metadata.get("include_constant",True)))
    features=np.asarray(dictionary.transform(transformed),dtype=float)
    groups=[]
    for function in krom.eigenfunctions:
        if function.coefficients is None:
            raise ValueError("Every candidate eigenfunction must provide polynomial coefficients.")
        coefficient=np.asarray(function.coefficients,dtype=complex).reshape(-1)
        if coefficient.size!=features.shape[1]:
            raise ValueError("Eigenfunction coefficient width does not match the polynomial dictionary.")
        columns=[coefficient.real]
        if not function.is_real: columns.append(coefficient.imag)
        groups.append(np.column_stack(columns))
    return features,groups


def select_observable_modes(krom, states: np.ndarray, policy: ObservableModeSelectionPolicy|None=None) -> ObservableModeSelectionReport:
    """Select nonredundant KROM modes that greedily reconstruct model state."""
    selected_policy=policy or ObservableModeSelectionPolicy()
    x=np.asarray(states,dtype=float)
    features,groups=_polynomial_features(krom,x)
    centered_features=features-np.mean(features,axis=0,keepdims=True)
    centered_state=x-np.mean(x,axis=0,keepdims=True)
    total_by_state=np.sum(centered_state*centered_state,axis=0)
    if np.any(total_by_state<=np.finfo(float).eps):
        raise ValueError("Every selected state component must vary over the training data.")
    feature_gram=centered_features.T@centered_features
    residual=centered_state.copy()
    total_residual=float(np.sum(residual*residual))
    if not np.isfinite(total_residual) or total_residual<=0.:
        raise ValueError("Selection training state variance must be finite and positive.")

    try:
        authority=np.linalg.norm(np.asarray(krom.B_psi(np.asarray(krom.system.initial_state,dtype=float),
            output="complex",state_coordinates="physical",time_coordinates="model"),dtype=complex),axis=1)
    except (AttributeError,NotImplementedError,TypeError,ValueError,np.linalg.LinAlgError):
        authority=np.zeros(len(groups),dtype=float)
    authority=np.asarray(authority,dtype=float).reshape(-1)
    authority_scale=float(np.max(authority)) if authority.size and np.any(np.isfinite(authority)) else 0.
    relative_authority=np.zeros_like(authority) if authority_scale<=0. else np.nan_to_num(authority/authority_scale)

    selected=[]; selected_columns=[]; steps=[]; r2=np.zeros(x.shape[1],dtype=float)
    stop_reason="max_modes"
    for step in range(1,selected_policy.max_modes+1):
        feature_residual=centered_features.T@residual
        best=None
        for index,columns in enumerate(groups):
            if index in selected: continue
            gram=columns.T@feature_gram@columns
            cross=columns.T@feature_residual
            reduction=float(np.trace(cross.T@np.linalg.pinv(gram,rcond=1e-12)@cross))
            score=reduction*(1.+selected_policy.control_authority_weight*float(relative_authority[index]))
            candidate=(score,reduction,-index,index)
            if best is None or candidate>best: best=candidate
        if best is None or best[1]/total_residual<selected_policy.minimum_relative_improvement:
            stop_reason="minimum_relative_improvement"; break
        _,reduction,_,index=best
        selected.append(index); selected_columns.append(groups[index])
        design=centered_features@np.column_stack(selected_columns)
        coefficients=np.linalg.lstsq(design,centered_state,rcond=1e-12)[0]
        residual=centered_state-design@coefficients
        r2=1.-np.sum(residual*residual,axis=0)/total_by_state
        full_index=int(krom._retained_indices[index])
        steps.append(ObservableSelectionStep(step,index,full_index,complex(krom.eigenvalues[index].continuous_time),
            reduction,float(best[0]),float(relative_authority[index]),tuple(float(value) for value in r2)))
        if np.all(r2>=selected_policy.r2_threshold):
            stop_reason="r2_threshold"; break
    return ObservableModeSelectionReport(selected_policy,tuple(selected),
        tuple(int(krom._retained_indices[index]) for index in selected),tuple(float(value) for value in r2),
        tuple(steps),int(x.shape[0]),int(x.shape[1]),len(groups),int(features.shape[1]),stop_reason)


__all__=["ObservableModeSelectionPolicy","ObservableModeSelectionReport","ObservableSelectionStep","select_observable_modes"]
