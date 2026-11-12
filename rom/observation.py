"""Validated state-observation maps for Koopman reduced-order models."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Literal

import numpy as np

from dynsyslib.observables import Observable


ObservationSplit=Literal["interleaved","tail"]


def split_lifted_coordinates(values: np.ndarray, eigenfunctions) -> np.ndarray:
    """Encode nonredundant complex Koopman coordinates as real columns."""
    coordinates=np.asarray(values,dtype=complex)
    if coordinates.ndim==1: coordinates=coordinates.reshape(1,-1)
    if coordinates.ndim!=2 or coordinates.shape[1]!=len(eigenfunctions):
        raise ValueError("Lifted coordinates must have one column per nonredundant eigenfunction.")
    columns=[]
    for index,function in enumerate(eigenfunctions):
        columns.append(coordinates[:,index].real[:,None])
        if not function.is_real: columns.append(coordinates[:,index].imag[:,None])
    return np.hstack(columns) if columns else np.empty((coordinates.shape[0],0),dtype=float)


@dataclass(frozen=True)
class ObservationMapQualityPolicy:
    """Promotion thresholds for a learned state decoder."""

    minimum_r2: float=0.99
    maximum_normalized_rmse: float=0.10
    maximum_condition_number: float=1e12

    def __post_init__(self):
        if not -np.inf<self.minimum_r2<=1.:
            raise ValueError("minimum_r2 must not exceed 1.")
        if self.maximum_normalized_rmse<0. or not np.isfinite(self.maximum_normalized_rmse):
            raise ValueError("maximum_normalized_rmse must be finite and nonnegative.")
        if self.maximum_condition_number<1. or not np.isfinite(self.maximum_condition_number):
            raise ValueError("maximum_condition_number must be finite and at least 1.")


@dataclass(frozen=True)
class StateObservationMap:
    """Affine map from real/imaginary split Koopman coordinates to state."""

    coefficients: np.ndarray
    intercept: np.ndarray
    quality: dict[str,Any]
    coordinate_system: str="model"
    encoding: str="real-imaginary-split-affine/v1"

    def __post_init__(self):
        coefficients=np.array(self.coefficients,dtype=float,copy=True)
        intercept=np.array(self.intercept,dtype=float,copy=True).reshape(-1)
        if coefficients.ndim!=2 or coefficients.shape[1]!=intercept.size:
            raise ValueError("Observation coefficients must have shape (lifted_features, state_dimension).")
        if not coefficients.size or np.any(~np.isfinite(coefficients)) or np.any(~np.isfinite(intercept)):
            raise ValueError("Observation-map arrays must be finite and nonempty.")
        if self.encoding!="real-imaginary-split-affine/v1":
            raise ValueError("Unsupported observation-map encoding.")
        coefficients.setflags(write=False); intercept.setflags(write=False)
        object.__setattr__(self,"coefficients",coefficients)
        object.__setattr__(self,"intercept",intercept)
        object.__setattr__(self,"quality",dict(self.quality))

    @property
    def lifted_feature_dimension(self) -> int: return int(self.coefficients.shape[0])

    @property
    def state_dimension(self) -> int: return int(self.coefficients.shape[1])

    @property
    def validated(self) -> bool: return bool(self.quality.get("passed",False))

    def reconstruct(self, lifted_state: np.ndarray, eigenfunctions) -> np.ndarray:
        values=np.asarray(lifted_state,dtype=complex)
        one=values.ndim==1
        design=split_lifted_coordinates(values,eigenfunctions)
        if design.shape[1]!=self.lifted_feature_dimension:
            raise ValueError("Observation map is not aligned with this KROM's lifted coordinates.")
        result=design@self.coefficients+self.intercept
        return result[0] if one else result

    def summary(self) -> dict[str,Any]:
        return {"encoding":self.encoding,"coordinate_system":self.coordinate_system,
            "lifted_feature_dimension":self.lifted_feature_dimension,
            "state_dimension":self.state_dimension,"validated":self.validated,
            "quality":dict(self.quality)}


@dataclass(frozen=True)
class ObservableReconstruction:
    """Koopman mode matrix ``C`` for a declared physical/model observable.

    ``coefficients`` has shape ``(observable_dimension, retained_modes)`` and
    uses the nonredundant complex convention.  A retained complex coordinate
    implicitly includes its conjugate contribution, so reconstruction is
    ``c psi`` for a real mode and ``2 Re(c psi)`` for a complex pair.
    """

    observable: Observable
    coefficients: np.ndarray
    intercept: np.ndarray
    quality: dict[str, Any]
    provenance: dict[str, Any]
    encoding: str = "nonredundant-complex-koopman-modes/v1"

    def __post_init__(self):
        coefficients = np.array(self.coefficients, dtype=complex, copy=True)
        intercept = np.array(self.intercept, dtype=float, copy=True).reshape(-1)
        if coefficients.ndim != 2 or coefficients.shape[0] != self.observable.output_dimension:
            raise ValueError("Observable reconstruction C must have one row per observable output.")
        if coefficients.shape[1] < 1 or intercept.size != coefficients.shape[0]:
            raise ValueError("Observable reconstruction arrays have incompatible dimensions.")
        if np.any(~np.isfinite(coefficients)) or np.any(~np.isfinite(intercept)):
            raise ValueError("Observable reconstruction arrays must be finite.")
        if self.encoding != "nonredundant-complex-koopman-modes/v1":
            raise ValueError("Unsupported observable reconstruction encoding.")
        coefficients.setflags(write=False); intercept.setflags(write=False)
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "intercept", intercept)
        object.__setattr__(self, "quality", dict(self.quality))
        object.__setattr__(self, "provenance", dict(self.provenance))

    @property
    def observable_dimension(self) -> int:
        return int(self.coefficients.shape[0])

    @property
    def lifted_dimension(self) -> int:
        return int(self.coefficients.shape[1])

    @property
    def validated(self) -> bool:
        return bool(self.quality.get("passed", False))

    def split_matrix(self, eigenfunctions) -> np.ndarray:
        """Return real ``C`` aligned with interleaved real/imag coordinates."""
        if len(eigenfunctions) != self.lifted_dimension:
            raise ValueError("Observable reconstruction is not aligned with these eigenfunctions.")
        columns = []
        for index, function in enumerate(eigenfunctions):
            coefficient = self.coefficients[:, index]
            if function.is_real:
                columns.append(coefficient.real[:, None])
            else:
                columns.extend((2.0 * coefficient.real[:, None], -2.0 * coefficient.imag[:, None]))
        return np.hstack(columns)

    def reconstruct(self, lifted_state: np.ndarray, eigenfunctions) -> np.ndarray:
        values = np.asarray(lifted_state, dtype=complex)
        one = values.ndim == 1
        design = split_lifted_coordinates(values, eigenfunctions)
        result = design @ self.split_matrix(eigenfunctions).T + self.intercept
        return result[0] if one else result

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": self.encoding,
            "observable": self.observable.summary(),
            "observable_dimension": self.observable_dimension,
            "lifted_dimension": self.lifted_dimension,
            "validated": self.validated,
            "quality": dict(self.quality),
            "provenance": dict(self.provenance),
        }

    def subset(self, positions: list[int] | tuple[int, ...]) -> "ObservableReconstruction":
        selected = np.asarray(positions, dtype=int)
        return ObservableReconstruction(
            self.observable,
            self.coefficients[:, selected],
            self.intercept,
            {**self.quality, "subset_of_validated_map": self.validated, "passed": False},
            {**self.provenance, "subset_positions": selected.tolist()},
        )


def _metrics(actual: np.ndarray, predicted: np.ndarray) -> dict[str,Any]:
    residual=np.asarray(actual,dtype=float)-np.asarray(predicted,dtype=float)
    centered=actual-np.mean(actual,axis=0,keepdims=True)
    variance=np.sum(centered*centered,axis=0)
    squared=np.sum(residual*residual,axis=0)
    r2=np.where(variance>np.finfo(float).eps,1.-squared/variance,np.where(squared<=np.finfo(float).eps,1.,-np.inf))
    rmse=np.sqrt(np.mean(residual*residual,axis=0))
    scale=np.ptp(actual,axis=0)
    normalized=np.where(scale>np.finfo(float).eps,rmse/scale,np.where(rmse<=np.finfo(float).eps,0.,np.inf))
    return {"sample_count":int(actual.shape[0]),"r2_by_state":r2.tolist(),
        "rmse_by_state":rmse.tolist(),"normalized_rmse_by_state":normalized.tolist(),
        "maximum_absolute_error_by_state":np.max(np.abs(residual),axis=0).tolist()}


def _lift(krom, states: np.ndarray) -> np.ndarray:
    values=np.asarray(states,dtype=float)
    if values.ndim!=2 or values.shape[1]!=np.asarray(krom.system.initial_state).size:
        raise ValueError("Observation-map states must have shape (samples, system state dimension).")
    if values.shape[0]<3 or np.any(~np.isfinite(values)):
        raise ValueError("Observation-map states require at least three finite samples.")
    dictionary_metadata=dict(krom.spectrum.metadata.get("dictionary") or {})
    preprocessing=dict(krom.spectrum.metadata.get("state_preprocessing") or {})
    coefficient_rows=[function.coefficients for function in krom.eigenfunctions]
    if (dictionary_metadata.get("encoding")=="polynomial-dictionary/v1" and
            {"center","scale","multiplier"}.issubset(preprocessing) and
            coefficient_rows and all(item is not None for item in coefficient_rows)):
        from spectral.dictionary import PolynomialDictionary
        coordinates=np.asarray(krom.coordinate_observable.evaluate(values),dtype=float)
        transformed=(coordinates-np.asarray(preprocessing["center"],dtype=float))/np.asarray(preprocessing["scale"],dtype=float)
        transformed*=np.asarray(preprocessing["multiplier"],dtype=float)
        dictionary=PolynomialDictionary(degree=int(dictionary_metadata["degree"]),
            include_constant=bool(dictionary_metadata.get("include_constant",True)))
        features=np.asarray(dictionary.transform(transformed),dtype=float)
        coefficients=np.column_stack([np.asarray(item,dtype=complex).reshape(-1) for item in coefficient_rows])
        if features.shape[1]!=coefficients.shape[0]:
            raise ValueError("Archived polynomial coefficients do not match the dictionary width.")
        complex_values=features@coefficients
    else:
        complex_values=np.vstack([krom.psi(row,output="complex",state_coordinates="model") for row in values])
    return split_lifted_coordinates(complex_values,krom.eigenfunctions)


def fit_state_observation_map(krom, states: np.ndarray, *, validation_states: np.ndarray|None=None,
                              validation_fraction: float=0.2, split: ObservationSplit="interleaved",
                              ridge: float=0., quality_policy: ObservationMapQualityPolicy|None=None,
                              require_quality: bool=True) -> StateObservationMap:
    """Fit and held-out validate ``x = intercept + C psi``.

    Supplying ``validation_states`` is preferred for independent validation.
    Otherwise a deterministic held-out split is used and recorded explicitly.
    """
    policy=quality_policy or ObservationMapQualityPolicy()
    values=np.asarray(states,dtype=float)
    if ridge<0. or not np.isfinite(ridge): raise ValueError("ridge must be finite and nonnegative.")
    if validation_states is None:
        if not 0.<validation_fraction<.5: raise ValueError("validation_fraction must lie in (0, 0.5).")
        count=values.shape[0]
        validation_count=max(1,int(np.ceil(count*validation_fraction)))
        if count-validation_count<2: raise ValueError("Observation-map fitting needs at least two fit samples.")
        if split=="tail": validation_indices=np.arange(count-validation_count,count)
        elif split=="interleaved":
            stride=max(2,int(round(1./validation_fraction)))
            validation_indices=np.arange(0,count,stride)[:validation_count]
        else: raise ValueError("split must be 'interleaved' or 'tail'.")
        mask=np.ones(count,dtype=bool); mask[validation_indices]=False
        fit_states=values[mask]; held_out=values[validation_indices]
        validation_source={"kind":"internal-held-out","split":split,"fraction":float(validation_fraction)}
    else:
        fit_states=values; held_out=np.asarray(validation_states,dtype=float)
        validation_source={"kind":"independent-states","split":None,"fraction":None}
    fit_design=_lift(krom,fit_states); validation_design=_lift(krom,held_out)
    center=np.mean(fit_design,axis=0); state_center=np.mean(fit_states,axis=0)
    centered=fit_design-center; target=fit_states-state_center
    feature_scale=np.std(centered,axis=0)
    feature_scale=np.where(feature_scale>64.*np.finfo(float).eps,feature_scale,1.)
    normalized=centered/feature_scale
    gram=normalized.T@normalized
    regularized=gram+ridge*np.eye(gram.shape[0])
    normalized_coefficients=np.linalg.pinv(regularized,rcond=1e-12)@normalized.T@target
    coefficients=normalized_coefficients/feature_scale[:,None]
    intercept=state_center-center@coefficients
    train_metrics=_metrics(fit_states,fit_design@coefficients+intercept)
    validation_metrics=_metrics(held_out,validation_design@coefficients+intercept)
    condition=float(np.linalg.cond(regularized))
    r2=np.asarray(validation_metrics["r2_by_state"],dtype=float)
    nrmse=np.asarray(validation_metrics["normalized_rmse_by_state"],dtype=float)
    passed=bool(np.all(r2>=policy.minimum_r2) and np.all(nrmse<=policy.maximum_normalized_rmse)
                and condition<=policy.maximum_condition_number)
    quality={"method":"held-out-affine-state-observation/v1","passed":passed,
        "policy":asdict(policy),"ridge":float(ridge),"condition_number":condition,
        "validation_source":validation_source,"fit":train_metrics,"validation":validation_metrics}
    if require_quality and not passed:
        raise ValueError("Fitted KROM observation map failed held-out quality thresholds; rerun with require_quality=False only for a labeled diagnostic artifact.")
    return StateObservationMap(coefficients,intercept,quality,coordinate_system="model")


def fit_observable_reconstruction(
    krom,
    states: np.ndarray,
    observable: Observable,
    *,
    validation_states: np.ndarray | None = None,
    ridge: float = 0.0,
) -> ObservableReconstruction:
    """Fit ``h(x) = c0 + C Psi(x)`` with an independent validation report."""
    values = np.asarray(states, dtype=float)
    if values.ndim != 2 or values.shape[1] != observable.state_dimension:
        raise ValueError("Observable reconstruction states have the wrong shape.")
    if ridge < 0.0 or not np.isfinite(ridge):
        raise ValueError("ridge must be finite and nonnegative.")
    if validation_states is None:
        indices = np.arange(values.shape[0])
        validation_mask = indices % 5 == 0
        if np.sum(~validation_mask) < 2 or np.sum(validation_mask) < 1:
            raise ValueError("Observable reconstruction needs at least three samples.")
        fit_states = values[~validation_mask]; held_out = values[validation_mask]
        validation_source = "deterministic-interleaved-20-percent"
    else:
        fit_states = values; held_out = np.asarray(validation_states, dtype=float)
        validation_source = "independent-states"
    fit_complex = np.vstack([krom.psi(row, output="complex") for row in fit_states])
    fit_design = split_lifted_coordinates(fit_complex, krom.eigenfunctions)
    target = observable.evaluate(fit_states)
    augmented = np.column_stack((np.ones(fit_design.shape[0]), fit_design))
    gram = augmented.T @ augmented + ridge * np.diag(
        np.r_[0.0, np.ones(fit_design.shape[1])]
    )
    fit = np.linalg.pinv(gram, rcond=1e-12) @ augmented.T @ target
    intercept = fit[0]; split_coefficients = fit[1:].T
    complex_coefficients = np.empty(
        (observable.output_dimension, len(krom.eigenfunctions)), dtype=complex
    )
    offset = 0
    for index, function in enumerate(krom.eigenfunctions):
        if function.is_real:
            complex_coefficients[:, index] = split_coefficients[:, offset]
            offset += 1
        else:
            complex_coefficients[:, index] = 0.5 * (
                split_coefficients[:, offset] - 1j * split_coefficients[:, offset + 1]
            )
            offset += 2
    provisional = ObservableReconstruction(
        observable,
        complex_coefficients,
        intercept,
        {"passed": False},
        {"method": "least-squares-observable-reconstruction/v1"},
    )
    fit_prediction = provisional.reconstruct(fit_complex, krom.eigenfunctions)
    validation_complex = np.vstack([krom.psi(row, output="complex") for row in held_out])
    validation_prediction = provisional.reconstruct(validation_complex, krom.eigenfunctions)
    quality = {
        "passed": True,
        "method": "held-out-observable-reconstruction/v1",
        "ridge": float(ridge),
        "condition_number": float(np.linalg.cond(gram)),
        "validation_source": validation_source,
        "fit": _metrics(target, fit_prediction),
        "validation": _metrics(observable.evaluate(held_out), validation_prediction),
    }
    return ObservableReconstruction(
        observable,
        complex_coefficients,
        intercept,
        quality,
        {"method": "least-squares-observable-reconstruction/v1"},
    )


def transfer_observable_reconstruction(source_krom, target_krom, states: np.ndarray) -> ObservableReconstruction:
    """Transfer H-DMD ``C`` after correcting refitted eigenfunction scales.

    Koopman eigenfunctions are only defined up to a nonzero scalar.  SINDy
    therefore cannot safely inherit H-DMD modes verbatim: each target
    eigenfunction is aligned to its source counterpart on common state samples
    before the modal coefficient is transferred.
    """
    source=source_krom.observable_reconstruction
    if source is None: raise ValueError("Source KROM has no observable reconstruction C.")
    if source_krom.lifted_dimension!=target_krom.lifted_dimension:
        raise ValueError("Source and target KROMs must retain the same number of nonredundant modes.")
    samples=np.asarray(states,dtype=float)
    if samples.ndim!=2 or samples.shape[0]<3: raise ValueError("C transfer requires common state samples.")
    source_values=[]; target_values=[]; scales=[]; residuals=[]
    for source_value,target_value,source_function,target_function in zip(
            source_krom.eigenvalues,target_krom.eigenvalues,source_krom.eigenfunctions,target_krom.eigenfunctions):
        if not np.isclose(source_value.continuous_time,target_value.continuous_time,rtol=1e-7,atol=1e-9):
            raise ValueError("Source and target C-transfer spectra are not eigenvalue aligned.")
        left=np.asarray([source_function.psi(source_krom.coordinate_observable.evaluate(row)) for row in samples],dtype=complex)
        right=np.asarray([target_function.psi(target_krom.coordinate_observable.evaluate(row)) for row in samples],dtype=complex)
        denominator=np.vdot(right,right)
        scale=0j if abs(denominator)<=np.finfo(float).eps else np.vdot(right,left)/denominator
        approximation=scale*right
        residual=np.linalg.norm(left-approximation)/max(np.linalg.norm(left),np.finfo(float).eps)
        source_values.append(left); target_values.append(right); scales.append(scale); residuals.append(float(residual))
    coefficients=source.coefficients*np.asarray(scales,dtype=complex)[None,:]
    quality={"passed":bool(np.all(np.isfinite(residuals)) and max(residuals)<=.2),"method":"eigenfunction-scale-aligned-C-transfer/v1",
        "maximum_accepted_relative_residual":.2,
        "relative_alignment_residual_by_mode":residuals,"maximum_relative_alignment_residual":float(max(residuals))}
    return ObservableReconstruction(source.observable,coefficients,source.intercept,quality,
        {"method":"H-DMD-to-SINDy C transfer with complex least-squares eigenfunction scaling",
         "source":source.summary(),"sample_count":int(samples.shape[0]),
         "complex_scales":[{"real":float(value.real),"imag":float(value.imag)} for value in scales]})


__all__=["ObservationMapQualityPolicy","ObservationSplit","StateObservationMap",
         "ObservableReconstruction","fit_observable_reconstruction",
         "transfer_observable_reconstruction",
         "fit_state_observation_map","split_lifted_coordinates"]
