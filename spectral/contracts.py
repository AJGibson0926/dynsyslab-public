"""Typed, provenance-preserving representations of Koopman spectra."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Callable, Literal
from uuid import uuid4
import numpy as np

TimeDomain = Literal["discrete", "continuous"]
EigenfunctionRepresentationKind = Literal[
    "analytic",
    "scattered-interpolation",
    "sampled-values",
    "unknown",
]
GradientBackend = Literal["auto", "analytic", "finite_difference"]

def _complex_record(value: complex) -> dict[str, float]:
    return {"real":float(np.real(value)),"imag":float(np.imag(value))}

class ScatteredEigenfunctionInterpolator:
    """Fast inverse-distance interpolation over scattered state samples.

    H-DMD and fitted-boundary eigenfunctions are authoritative at their stored
    trajectory samples.  A KD tree supplies the state-space continuation needed
    by KROM evaluation without turning every control call into a full scan of
    the training ensemble.
    """

    def __init__(self, states: np.ndarray, values: np.ndarray, *, neighbors: int=8):
        from scipy.spatial import cKDTree
        samples=np.asarray(states,dtype=float)
        coordinates=np.asarray(values,dtype=complex).reshape(-1)
        if samples.ndim==1: samples=samples[:,None]
        if samples.ndim!=2 or samples.shape[0]!=coordinates.size:
            raise ValueError("Interpolation states and eigenfunction values must have equal sample counts.")
        if samples.shape[0]<1 or np.any(~np.isfinite(samples)) or np.any(~np.isfinite(coordinates)):
            raise ValueError("Scattered eigenfunction samples must be finite and nonempty.")
        if int(neighbors)!=neighbors or neighbors<1:
            raise ValueError("neighbors must be a positive integer.")
        self.samples=np.array(samples,copy=True)
        self.values=np.array(coordinates,copy=True)
        self.neighbors=min(int(neighbors),coordinates.size)
        self._tree=cKDTree(self.samples)

    @property
    def state_dimension(self) -> int:
        return int(self.samples.shape[1])

    def __call__(self, state: np.ndarray) -> complex:
        point=np.asarray(state,dtype=float).reshape(-1)
        if point.size!=self.state_dimension or np.any(~np.isfinite(point)):
            raise ValueError("State dimension does not match the eigenfunction interpolation data.")
        distance,indices=self._tree.query(point,k=self.neighbors)
        distance=np.atleast_1d(np.asarray(distance,dtype=float))
        indices=np.atleast_1d(np.asarray(indices,dtype=int))
        if distance[0]<=1e-13: return complex(self.values[indices[0]])
        weights=1./np.maximum(distance,1e-13)**2
        return complex(np.sum(weights*self.values[indices])/np.sum(weights))

    def summary(self) -> dict[str, Any]:
        return {"encoding":"kd-tree-inverse-distance-scattered-interpolant/v1",
                "sample_count":int(self.samples.shape[0]),"state_dimension":self.state_dimension,
                "neighbors":self.neighbors,"extrapolation":"nearest-neighbor weighted continuation"}


def _trajectory_interpolant(states: np.ndarray, values: np.ndarray) -> ScatteredEigenfunctionInterpolator:
    return ScatteredEigenfunctionInterpolator(states,values)

@dataclass(frozen=True)
class DynamicalSystemReference:
    """Mandatory archive-backed pointer and parameter snapshot for a source system.

    The reference is deliberately serializable rather than a live object so
    archive-only postprocessing remains independent of dynsyslib.  KROM uses
    it to reconstruct the live library object when it needs ``B(x)``.
    """
    source_archive: str | None
    system_name: str | None
    physical_parameters: dict[str, Any]

    @classmethod
    def from_archive_metadata(cls, metadata: dict[str, Any] | None, source_archive: str | None=None):
        metadata=metadata or {}
        return cls(source_archive,metadata.get("system"),dict(metadata.get("system_parameters") or {}))

@dataclass
class KoopmanEigenvalue:
    """A paired discrete/continuous eigenvalue with explicit derivation origin."""
    discrete_time: complex
    continuous_time: complex
    delta_t: float
    original_domain: TimeDomain
    system: DynamicalSystemReference
    metadata: dict[str, Any] = field(default_factory=dict)
    dual: KoopmanEigenvalue | None = field(default=None, repr=False, compare=False)

    @classmethod
    def from_discrete(cls, value: complex, delta_t: float, system: DynamicalSystemReference, metadata: dict[str, Any] | None=None):
        if delta_t<=0: raise ValueError("Koopman eigenvalue delta_t must be positive.")
        discrete=complex(value)
        return cls(discrete,np.log(discrete)/delta_t,float(delta_t),"discrete",system,dict(metadata or {}))

    @classmethod
    def from_continuous(cls, value: complex, delta_t: float, system: DynamicalSystemReference, metadata: dict[str, Any] | None=None):
        if delta_t<=0: raise ValueError("Koopman eigenvalue delta_t must be positive.")
        continuous=complex(value)
        return cls(np.exp(continuous*delta_t),continuous,float(delta_t),"continuous",system,dict(metadata or {}))

    @property
    def computed_domain(self) -> TimeDomain:
        return "continuous" if self.original_domain=="discrete" else "discrete"

    @property
    def is_real(self) -> bool:
        return bool(np.isclose(self.continuous_time.imag, 0., atol=1e-10, rtol=1e-8))

    @property
    def is_complex(self) -> bool:
        return not self.is_real

    @property
    def growth_decay_rate(self) -> float:
        """Continuous-time exponential rate, positive for growth and negative for decay."""
        return float(self.continuous_time.real)

    @property
    def angular_frequency(self) -> float:
        """Signed continuous-time angular frequency in radians per time unit."""
        return float(self.continuous_time.imag)

    @property
    def frequency(self) -> float:
        """Signed continuous-time cyclic frequency in cycles per time unit."""
        return float(self.angular_frequency/(2*np.pi))

    @property
    def q_factor(self) -> float:
        """Modal quality factor ``|omega|/(2|growth_decay_rate|)``."""
        omega=abs(self.angular_frequency); rate=abs(self.growth_decay_rate)
        if omega==0.: return 0.
        return float(np.inf if rate==0. else omega/(2*rate))

    def summary(self) -> dict[str, Any]:
        return {"discrete_time":_complex_record(self.discrete_time),"continuous_time":_complex_record(self.continuous_time),
                "delta_t":self.delta_t,"original_domain":self.original_domain,"computed_domain":self.computed_domain,
                "is_real":self.is_real,"is_complex":self.is_complex,"growth_decay_rate":self.growth_decay_rate,
                "angular_frequency":self.angular_frequency,"frequency":self.frequency,"q_factor":self.q_factor,
                "dual_index":self.metadata.get("dual_index"),
                "metadata":self.metadata}

@dataclass
class KoopmanEigenfunction:
    """One eigenfunction's stored representation and optional callable state map."""
    eigenvalue: KoopmanEigenvalue
    representation: str
    values: np.ndarray | None = None
    coefficients: np.ndarray | None = None
    state_callable: Callable[[np.ndarray], complex] | None = None
    gradient_callable: Callable[[np.ndarray], np.ndarray] | None = None
    residual_operator: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    representation_kind: EigenfunctionRepresentationKind | None = None
    dual: KoopmanEigenfunction | None = field(default=None, repr=False, compare=False)
    error_evaluations: dict[str, Any] = field(default_factory=dict, repr=False, compare=False)
    coordinate_observable: Any | None = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.representation_kind is None:
            if self.gradient_callable is not None or self.coefficients is not None:
                self.representation_kind="analytic"
            elif isinstance(self.state_callable,ScatteredEigenfunctionInterpolator):
                self.representation_kind="scattered-interpolation"
            elif self.values is not None:
                self.representation_kind="sampled-values"
            else:
                self.representation_kind="unknown"
        if self.representation_kind not in {"analytic","scattered-interpolation","sampled-values","unknown"}:
            raise ValueError("Unsupported Koopman eigenfunction representation_kind.")
        self.error_evaluations=dict(self.error_evaluations)

    @property
    def has_state_callable(self) -> bool:
        return self.state_callable is not None

    @property
    def gradient_backend(self) -> str | None:
        """The default derivative backend implied by this representation."""
        if self.gradient_callable is not None: return "analytic"
        if self.state_callable is not None: return "finite_difference"
        return None

    @property
    def has_gradient(self) -> bool:
        return self.gradient_backend is not None

    @property
    def is_real(self) -> bool:
        return self.eigenvalue.is_real

    @property
    def is_complex(self) -> bool:
        return self.eigenvalue.is_complex

    def psi(self, x: np.ndarray, *, output: Literal["complex","split"]="complex") -> complex | np.ndarray:
        """Evaluate the eigenfunction in the spectrum's declared state coordinates.

        ``output='complex'`` returns the conventional complex Koopman value.
        ``output='split'`` returns ``[real(psi), imag(psi)]`` as a float vector.
        """
        if self.state_callable is None:
            raise NotImplementedError("This eigenfunction has no state-space evaluation callable.")
        value=complex(self.state_callable(np.asarray(x,dtype=float)))
        if output=="complex": return value
        if output=="split": return np.array((value.real,value.imag),dtype=float)
        raise ValueError("output must be 'complex' or 'split'.")

    def psi_from_state(self, state: np.ndarray, *, output: Literal["complex","split"]="complex"):
        """Evaluate after applying this eigenfunction's observable coordinates."""
        x=np.asarray(state,dtype=float)
        coordinates=x if self.coordinate_observable is None else self.coordinate_observable.evaluate(x)
        return self.psi(coordinates,output=output)

    def nabla_psi(self, x: np.ndarray, *, output: Literal["complex","split"]="complex",
                  backend: GradientBackend="auto", relative_step: float=1e-5) -> np.ndarray:
        """Gradient in eigenfunction coordinates with representation-aware dispatch.

        ``auto`` uses an analytic derivative when the representation provides
        one and otherwise differentiates its state callable with a fourth-order
        centered stencil.  Explicit backend selection remains available for
        validation and parity experiments.
        """
        state=np.asarray(x,dtype=float).reshape(-1)
        if backend not in {"auto","analytic","finite_difference"}:
            raise ValueError("backend must be 'auto', 'analytic', or 'finite_difference'.")
        if relative_step<=0. or not np.isfinite(relative_step):
            raise ValueError("relative_step must be finite and positive.")
        selected=self.gradient_backend if backend=="auto" else backend
        if selected=="analytic":
            if self.gradient_callable is None:
                raise NotImplementedError("This eigenfunction representation has no analytic derivative.")
            gradient=np.asarray(self.gradient_callable(state),dtype=complex).reshape(-1)
        elif selected=="finite_difference" and self.state_callable is not None:
            h=relative_step*np.maximum(1.,np.abs(state)); gradient=np.empty(state.size,dtype=complex)
            for axis in range(state.size):
                delta=np.zeros_like(state); delta[axis]=h[axis]
                gradient[axis]=(-self.psi(state+2*delta)+8*self.psi(state+delta)-8*self.psi(state-delta)+self.psi(state-2*delta))/(12*h[axis])
        else: raise NotImplementedError("This eigenfunction has no state-space derivative backend.")
        if output=="complex": return gradient
        if output=="split": return np.vstack((gradient.real,gradient.imag))
        raise ValueError("output must be 'complex' or 'split'.")

    def nabla_psi_from_state(self, state: np.ndarray, *, output: Literal["complex","split"]="complex",
                             backend: GradientBackend="auto", observable_backend: str="auto",
                             relative_step: float=1e-5) -> np.ndarray:
        """Compose ``nabla_y psi`` with ``J_h=dh/dx`` for ``y=h(x)``."""
        x=np.asarray(state,dtype=float).reshape(-1)
        if self.coordinate_observable is None:
            return self.nabla_psi(x,output=output,backend=backend,relative_step=relative_step)
        observable=self.coordinate_observable
        y=np.asarray(observable.evaluate(x),dtype=float)
        selected=("analytic" if observable_backend=="auto" and observable.has_analytic_jacobian
                  else "finite_difference" if observable_backend=="auto" else observable_backend)
        jacobian=np.asarray(observable.jacobian(x,backend=selected,relative_step=relative_step),dtype=complex)
        gradient=np.asarray(self.nabla_psi(y,backend=backend,relative_step=relative_step),dtype=complex)@jacobian
        if output=="complex": return gradient
        if output=="split": return np.vstack((gradient.real,gradient.imag))
        raise ValueError("output must be 'complex' or 'split'.")

    def residuals(self, states: np.ndarray, *, derivatives: np.ndarray | None=None, time: np.ndarray | None=None,
                  system: Any | None=None, system_time: float=0.) -> dict[str, Any]:
        """Evaluate generator-PDE residuals on a test trajectory or state set."""
        if self.residual_operator is None:
            raise NotImplementedError("This eigenfunction does not provide a generator residual operator.")
        x=np.asarray(states,dtype=float)
        if x.ndim==1: x=x[:,None]
        if derivatives is None:
            if system is not None:
                derivatives=np.asarray([system.rhs(system_time,row) for row in x],dtype=float)
            elif time is not None:
                t=np.asarray(time,dtype=float).reshape(-1)
                if t.size!=x.shape[0] or np.any(np.diff(t)<=0): raise ValueError("Residual time values must be strictly increasing and match states.")
                derivatives=np.gradient(x,t,axis=0,edge_order=2)
            else: raise ValueError("Provide derivatives, time, or a system to compute eigenfunction residuals.")
        residual=np.asarray(self.residual_operator(x,np.asarray(derivatives,dtype=float)))
        magnitude=np.abs(residual)
        return {"residual":residual,"rms":float(np.sqrt(np.mean(magnitude**2))),"max_abs":float(np.max(magnitude)),"samples":int(x.shape[0])}

    def attach_error_evaluation(self, evaluation, *, replace: bool=False):
        """Attach one immutable training/test validation record to this function."""
        from .errors import EigenfunctionErrorEvaluation
        if not isinstance(evaluation,EigenfunctionErrorEvaluation):
            raise TypeError("evaluation must be an EigenfunctionErrorEvaluation.")
        if not np.isclose(evaluation.eigenvalue,self.eigenvalue.continuous_time,rtol=1e-10,atol=1e-12):
            raise ValueError("An error evaluation must match the eigenfunction eigenvalue.")
        if evaluation.evaluation_id in self.error_evaluations and not replace:
            raise KeyError(f"Eigenfunction error evaluation {evaluation.evaluation_id!r} already exists.")
        self.error_evaluations[evaluation.evaluation_id]=evaluation
        return evaluation

    def evaluate_errors(self, dataset, *, horizons=None, coordinate_map=None,
                        training_domain=None, metadata=None, attach: bool=True):
        """Evaluate this function against a named free training or test dataset."""
        from .errors import evaluate_eigenfunction_errors
        evaluation=evaluate_eigenfunction_errors(self,dataset,horizons=horizons,
            coordinate_map=coordinate_map,training_domain=training_domain,metadata=metadata)
        if attach: self.attach_error_evaluation(evaluation)
        return evaluation

    def summary(self) -> dict[str, Any]:
        interpolation=(self.state_callable.summary() if isinstance(self.state_callable,ScatteredEigenfunctionInterpolator) else None)
        return {"representation":self.representation,"representation_kind":self.representation_kind,
                "gradient_backend":self.gradient_backend,"interpolation":interpolation,"has_values":self.values is not None,
                "has_coefficients":self.coefficients is not None,"has_psi":self.has_state_callable,"has_nabla_psi":self.has_gradient,"has_residual_operator":self.residual_operator is not None,
                "is_real":self.is_real,"is_complex":self.is_complex,"dual_index":self.metadata.get("dual_index"),
                "coordinate_observable":None if self.coordinate_observable is None else self.coordinate_observable.summary(),
                "error_evaluations":[item.summary() for item in self.error_evaluations.values()],
                "metadata":self.metadata}

@dataclass
class KoopmanSpectrum:
    """Everything produced by one Koopman/DMD extraction run."""
    system: DynamicalSystemReference
    eigenvalues: list[KoopmanEigenvalue]
    eigenfunctions: list[KoopmanEigenfunction]
    modes: list[np.ndarray] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    dynamical_system: Any | None = field(default=None, repr=False, compare=False)
    eigenvalue_source: KoopmanSpectrum | None = field(default=None, repr=False, compare=False)

    @property
    def dynamical_system_reference(self) -> DynamicalSystemReference:
        """Archive-safe pointer to the originating dynsyslib system and parameters."""
        return self.system

    def __post_init__(self):
        # Eigenvalue-only spectra are valid inputs to later construction methods
        # such as Koopman-with-SINDy. Once functions are supplied, pairing must
        # remain exact and ordered.
        if self.eigenfunctions and len(self.eigenvalues)!=len(self.eigenfunctions): raise ValueError("A populated spectrum requires one eigenfunction per eigenvalue.")
        if self.eigenfunctions and any(function.eigenvalue is not value for function,value in zip(self.eigenfunctions,self.eigenvalues)):
            raise ValueError("Each eigenfunction must point to its corresponding spectrum eigenvalue instance.")
        if self.dynamical_system is None:
            try:
                from dynsyslib import resolve_archived_system
                self.dynamical_system=resolve_archived_system(self.system.system_name,self.system.physical_parameters)
            except Exception as error:
                raise ValueError("KoopmanSpectrum requires a resolvable dynamical-system reference or an explicit dynamical_system= object.") from error
        self.metadata.setdefault("spectrum_id",str(uuid4()))
        if self.eigenvalue_source is not None:
            if len(self.eigenvalues)!=len(self.eigenvalue_source.eigenvalues) or any(left is not right for left,right in zip(self.eigenvalues,self.eigenvalue_source.eigenvalues)):
                raise ValueError("A derived spectrum must reuse the exact eigenvalue instances from its eigenvalue_source.")
            self.metadata["eigenvalue_source_spectrum_id"]=self.eigenvalue_source.metadata["spectrum_id"]

    def pair_conjugates(self, *, atol: float=1e-9, rtol: float=1e-7) -> None:
        """Link numerically conjugate complex modes without discarding spectrum data."""
        for index,value in enumerate(self.eigenvalues):
            value.dual=None; value.metadata.pop("dual_index",None)
            if index < len(self.eigenfunctions):
                self.eigenfunctions[index].dual=None; self.eigenfunctions[index].metadata.pop("dual_index",None)
        for index,value in enumerate(self.eigenvalues):
            if value.is_real or value.dual is not None: continue
            target=np.conjugate(value.continuous_time)
            candidate=next((other_index for other_index,other in enumerate(self.eigenvalues)
                            if other_index != index and other.dual is None and other.is_complex
                            and np.isclose(other.continuous_time,target,atol=atol,rtol=rtol)),None)
            if candidate is None: continue
            other=self.eigenvalues[candidate]; value.dual=other; other.dual=value
            value.metadata["dual_index"]=candidate; other.metadata["dual_index"]=index
            if index < len(self.eigenfunctions) and candidate < len(self.eigenfunctions):
                first,second=self.eigenfunctions[index],self.eigenfunctions[candidate]
                first.dual=second; second.dual=first
                first.metadata["dual_index"]=candidate; second.metadata["dual_index"]=index

    def summary(self, *, include_eigenvalues: bool=True) -> dict[str, Any]:
        return {"system":{"source_archive":self.system.source_archive,"system_name":self.system.system_name,
                           "physical_parameters":self.system.physical_parameters},"eigenvalues":[value.summary() for value in self.eigenvalues],
                "eigenfunctions":[function.summary() for function in self.eigenfunctions],"mode_count":len(self.modes),
                "eigenvalue_source_spectrum_id":None if self.eigenvalue_source is None else self.eigenvalue_source.metadata["spectrum_id"],"metadata":self.metadata} if include_eigenvalues else {"system":{"source_archive":self.system.source_archive,"system_name":self.system.system_name,
                "physical_parameters":self.system.physical_parameters},"eigenvalue_count":len(self.eigenvalues),"eigenfunctions":[function.summary() for function in self.eigenfunctions],"mode_count":len(self.modes),
                "eigenvalue_source_spectrum_id":self.eigenvalue_source.metadata["spectrum_id"],"metadata":self.metadata}

def discrete_spectrum(eigenvalues: np.ndarray, eigenfunction_values: np.ndarray | None, modes: np.ndarray | None, *,
                      time: np.ndarray, source_metadata: dict[str, Any] | None, source_archive: str | None,
                      method: str, representation: str, coefficients: np.ndarray | None=None, method_metadata: dict[str, Any] | None=None,
                      interpolation_states: np.ndarray | None=None) -> KoopmanSpectrum:
    """Build a typed spectrum from a discrete-time fitted Koopman operator."""
    samples=np.asarray(time,dtype=float)
    delta_t=float(samples[1]-samples[0])
    system=DynamicalSystemReference.from_archive_metadata(source_metadata,source_archive)
    shared={"method":method,**dict(method_metadata or {})}
    values=[KoopmanEigenvalue.from_discrete(value,delta_t,system,{**shared,"index":index}) for index,value in enumerate(np.asarray(eigenvalues))]
    functions=[]
    coordinates=None if eigenfunction_values is None else np.asarray(eigenfunction_values)
    coefficient_matrix=None if coefficients is None else np.asarray(coefficients)
    for index,value in enumerate(values):
        trajectory_values=None if coordinates is None else coordinates[:,index]
        state_callable=None if interpolation_states is None or trajectory_values is None else _trajectory_interpolant(interpolation_states,trajectory_values)
        functions.append(KoopmanEigenfunction(value,representation,
            values=trajectory_values,
            coefficients=None if coefficient_matrix is None else coefficient_matrix[:,index],
            state_callable=state_callable,
            metadata={**shared,"index":index,"value_samples":0 if coordinates is None else int(coordinates.shape[0]),
                      "state_representation":"scattered-interpolation" if state_callable is not None else "stored_values"},
            representation_kind="scattered-interpolation" if state_callable is not None else "sampled-values"))
    mode_list=[] if modes is None else [np.asarray(modes)[:,index] for index in range(np.asarray(modes).shape[1])]
    training_domain=None
    if interpolation_states is not None:
        samples_x=np.asarray(interpolation_states,dtype=float)
        training_domain={"encoding":"axis-aligned-training-domain/v1","minimum":np.min(samples_x,axis=0).tolist(),
            "maximum":np.max(samples_x,axis=0).tolist(),"sample_count":int(samples_x.shape[0]),
            "coordinate_system":str(((source_metadata or {}).get("postprocessing_preprocessing") or {}).get("coordinate_system","dimensional")),
            "source":"eigenfunction interpolation states"}
    metadata={**shared,"delta_t":delta_t,"eigenvalue_origin":"discrete operator eigendecomposition"}
    if training_domain is not None: metadata["training_domain"]=training_domain
    return KoopmanSpectrum(system,values,functions,mode_list,metadata)


__all__=["DynamicalSystemReference","EigenfunctionRepresentationKind",
    "GradientBackend","KoopmanEigenfunction","KoopmanEigenvalue",
    "KoopmanSpectrum","ScatteredEigenfunctionInterpolator","TimeDomain",
    "discrete_spectrum"]
