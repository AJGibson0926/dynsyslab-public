"""Koopman-with-SINDy: sparse generator-PDE eigenfunction construction.

Archive-only trajectory mode depends only on archives and NumPy/SciPy.  The
optional ``gamma_method='system_rhs'`` path and domain mode each lazily import
``dynsyslib`` to evaluate a dynsyslib ODE's RHS.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path
from typing import Any, Literal, Sequence
import numpy as np
from time import perf_counter
from ..core import Trajectory, load_solution
from time_series import TimeSeries
from .dictionary import ObservableDictionary, PolynomialDictionary
from .spectrum import DynamicalSystemReference, KoopmanEigenfunction, KoopmanEigenvalue, KoopmanSpectrum
from postprocessing.resources import require_memory_budget

EigenvalueDomain=Literal["discrete","continuous"]
GammaMethod=Literal["system_rhs","finite_difference"]
NullSpaceBackend=Literal["svd","gram-eigh"]
DirectionalDerivativeBackend=Literal["directional-finite-difference","dictionary-gradient"]
SparsificationPolicy=Literal["iterative-relative","matlab-single-absolute"]

@dataclass
class KoopmanSINDyResult:
    spectrum: KoopmanSpectrum
    training_states: np.ndarray
    training_coordinates: np.ndarray
    training_eigenfunction_values: np.ndarray | None
    coefficients: np.ndarray
    metadata: dict[str, Any]
    krom: Any | None = None

    def evaluate_training_eigenfunctions(self, columns: Sequence[int] | np.ndarray | None=None) -> np.ndarray:
        """Return materialized or exactly reconstructed training values."""
        if self.training_eigenfunction_values is not None:
            values=np.asarray(self.training_eigenfunction_values,dtype=complex)
            if columns is None: return values
            indices=np.asarray(columns,dtype=int).reshape(-1)
            if np.any(indices<0) or np.any(indices>=values.shape[1]):
                raise IndexError("Training-eigenfunction column index is out of range.")
            return values[:,indices]
        from .virtual_arrays import polynomial_eigenfunction_values
        return polynomial_eigenfunction_values(self.training_coordinates,self.coefficients,
            dictionary=self.metadata.get("dictionary") or {},
            state_preprocessing=self.metadata.get("state_preprocessing") or {},columns=columns)

@dataclass(frozen=True)
class _StateTransform:
    """Fitted physical-state to dictionary-coordinate transformation."""
    center: np.ndarray
    scale: np.ndarray
    multiplier: np.ndarray

    def __call__(self, states: np.ndarray) -> np.ndarray:
        return (np.asarray(states,dtype=float)-self.center)/self.scale*self.multiplier

    def metadata(self) -> dict[str, Any]:
        return {"center":self.center.tolist(),"scale":self.scale.tolist(),"multiplier":self.multiplier.tolist()}

def _state_transform(states: np.ndarray, *, normalize: bool, rescale: float | Sequence[float] | None) -> _StateTransform:
    x=np.asarray(states,dtype=float)
    if x.ndim!=2 or x.shape[0]<1: raise ValueError("State preprocessing requires a nonempty two-dimensional state array.")
    center=np.mean(x,axis=0) if normalize else np.zeros(x.shape[1])
    scale=np.std(x,axis=0) if normalize else np.ones(x.shape[1])
    scale=np.where(scale>np.finfo(float).eps,scale,1.)
    multiplier=np.ones(x.shape[1]) if rescale is None else np.broadcast_to(np.asarray(rescale,dtype=float),x.shape[1]).copy()
    if not np.all(np.isfinite(multiplier)) or np.any(multiplier<=0):
        raise ValueError("rescale must be a positive finite scalar or one positive factor per state component.")
    return _StateTransform(center,scale,multiplier)

def _features(dictionary: ObservableDictionary, transform: _StateTransform, states: np.ndarray) -> np.ndarray:
    return np.asarray(dictionary.transform(transform(states)),dtype=complex)

def _dictionary_metadata(dictionary: ObservableDictionary) -> dict[str,Any]:
    if isinstance(dictionary,PolynomialDictionary):
        return {"encoding":"polynomial-dictionary/v1","class":"PolynomialDictionary",
                "degree":int(dictionary.degree),"include_constant":bool(dictionary.include_constant)}
    return {"encoding":"unserializable-dictionary/v1","class":type(dictionary).__name__}

def _directional_derivative(dictionary: ObservableDictionary, states: np.ndarray, velocities: np.ndarray, *, transform: _StateTransform,
                            relative_step: float=1e-6, backend: DirectionalDerivativeBackend="directional-finite-difference") -> np.ndarray:
    """Evaluate ``∇Theta(x)·f(x)`` with an explicitly selected backend.

    The established directional finite-difference route remains the default.
    ``dictionary-gradient`` uses the dictionary Jacobian contract; it is exact
    for :class:`PolynomialDictionary` and is useful for parity experiments
    against implementations that construct polynomial derivatives directly.
    """
    x=np.asarray(states,dtype=float); velocity=np.asarray(velocities,dtype=float)
    if backend=="dictionary-gradient":
        transformed=transform(x)
        jacobian=np.asarray(dictionary.gradient(transformed),dtype=complex)
        transformed_velocity=velocity*(transform.multiplier/transform.scale)
        return np.einsum("sfd,sd->sf",jacobian,transformed_velocity)
    if backend!="directional-finite-difference":
        raise ValueError("directional_derivative_backend must be 'directional-finite-difference' or 'dictionary-gradient'.")
    scale=relative_step*np.maximum(1.,np.max(np.abs(x),axis=1,keepdims=True))
    return (_features(dictionary,transform,x+scale*velocity)-_features(dictionary,transform,x-scale*velocity))/(2.*scale)

def _rhs_gamma(dictionary: ObservableDictionary, transform: _StateTransform, system: Any, time: np.ndarray, states: np.ndarray, *,
               coordinate_system: str="dimensional", coordinate_observable=None,
               directional_derivative_backend: DirectionalDerivativeBackend="directional-finite-difference") -> np.ndarray:
    """Use a dynsyslib ODE RHS and a directional dictionary derivative for Gamma."""
    try: from dynsyslib.interfaces import ODESystem
    except ImportError as error: raise ImportError("gamma_method='system_rhs' requires the optional dynsyslib package.") from error
    if not isinstance(system,ODESystem):
        raise TypeError("gamma_method='system_rhs' requires a dynsyslib ODESystem instance supplied as system=.")
    rhs=system.dimensionless_rhs if coordinate_system=="dimensionless" and hasattr(system,"dimensionless_rhs") else system.rhs
    velocity=np.asarray([rhs(float(t),row) for t,row in zip(np.asarray(time,dtype=float),np.asarray(states,dtype=float))],dtype=float)
    if coordinate_observable is None:
        modeled_states=np.asarray(states,dtype=float); modeled_velocity=velocity
    else:
        modeled_states=coordinate_observable.evaluate(states)
        modeled_velocity=np.vstack([coordinate_observable.jacobian(row)@direction for row,direction in zip(states,velocity)])
    return _directional_derivative(dictionary,modeled_states,modeled_velocity,transform=transform,
        backend=directional_derivative_backend)


from .sindy_nullspace import (
    _conjugate_fit_plan, _generator_gram, _generator_gram_terms,
    _gram_smallest_vector, _matlab_single_absolute_null_vector,
    _sparse_null_vector, _sparse_null_vector_from_gram,
)







def _coerce_trajectories(trajectories: Sequence[str | Path | Trajectory] | TimeSeries) -> list[Trajectory]:
    if isinstance(trajectories,TimeSeries): trajectories=(trajectories,)
    result=[]
    for item in trajectories:
        if isinstance(item,TimeSeries):
            members=item.members()
            for index,member in enumerate(members):
                result.append(Trajectory(member,{"time_series_format":"dynsyslab.TimeSeries/v1","pooled_input":item.is_pooled,
                    "pooled_member_index":index if item.is_pooled else None,"pooled_member_count":item.ensemble_size},None))
            continue
        source=item if isinstance(item,Trajectory) else load_solution(item)
        time=np.asarray(source.time,dtype=float); states=np.asarray(source.values,dtype=float)
        if states.ndim!=2 or states.shape[0]<3 or time.size!=states.shape[0] or np.any(np.diff(time)<=0):
            raise ValueError("Each Koopman-with-SINDy trajectory must have at least three increasing-time state samples.")
        result.append(source)
    if not result: raise ValueError("Provide at least one archived trajectory or a domain.")
    dimension=result[0].values.shape[1]
    if any(item.values.ndim!=2 or item.values.shape[1]!=dimension for item in result):
        raise ValueError("All Koopman-with-SINDy trajectories must have the same state dimension.")
    return result

def _domain_data(domain: tuple[np.ndarray, np.ndarray], system: Any, dictionary: ObservableDictionary, density: int, domain_time: float,
                 *, normalize: bool, rescale: float | Sequence[float] | None, gamma_method: GammaMethod, coordinate_observable=None,
                 coordinate_system: str="dimensional",
                 directional_derivative_backend: DirectionalDerivativeBackend="directional-finite-difference"):
    # The sole optional cross-package dependency, loaded only in this domain branch.
    try: from dynsyslib.interfaces import ODESystem
    except ImportError as error: raise ImportError("Domain-based Koopman-with-SINDy requires the optional dynsyslib package.") from error
    if not isinstance(system,ODESystem): raise TypeError("Domain mode requires a dynsyslib.ODESystem instance so its equations of motion are available.")
    lower,upper=(np.asarray(part,dtype=float).reshape(-1) for part in domain)
    if lower.shape!=upper.shape or lower.size<1 or np.any(upper<=lower): raise ValueError("domain must be (lower_bounds, upper_bounds) with upper > lower componentwise.")
    if density<3: raise ValueError("domain_density must be at least 3 samples per state dimension.")
    count=density**lower.size
    if count>250_000: raise ValueError("Requested domain grid exceeds 250,000 points; lower domain_density explicitly.")
    mesh=np.meshgrid(*[np.linspace(lo,hi,density) for lo,hi in zip(lower,upper)],indexing="ij")
    states=np.column_stack([axis.ravel() for axis in mesh])
    rhs=system.dimensionless_rhs if coordinate_system=="dimensionless" and hasattr(system,"dimensionless_rhs") else system.rhs
    derivatives=np.asarray([rhs(domain_time,row) for row in states],dtype=float)
    modeled_states=states if coordinate_observable is None else coordinate_observable.evaluate(states)
    modeled_derivatives=derivatives if coordinate_observable is None else np.vstack([coordinate_observable.jacobian(row)@direction for row,direction in zip(states,derivatives)])
    transform=_state_transform(modeled_states,normalize=normalize,rescale=rescale)
    theta=_features(dictionary,transform,modeled_states)
    # A domain has no temporal samples to differentiate.  Both modes therefore
    # use the supplied vector field; the directional derivative remains a
    # finite-difference approximation to the dictionary gradient.
    gamma=_directional_derivative(dictionary,modeled_states,modeled_derivatives,transform=transform,
        backend=directional_derivative_backend)
    parameters=asdict(system) if is_dataclass(system) else {key:repr(value) for key,value in vars(system).items()}
    reference=DynamicalSystemReference(None,getattr(system,"name",type(system).__name__),parameters)
    return states,theta,gamma,reference,{"input":"domain","domain_lower":lower.tolist(),"domain_upper":upper.tolist(),"domain_density":density,"domain_time":domain_time,"gamma_method":gamma_method,"optional_library_dependency_used":True,"preprocessing":{"sample_count":int(states.shape[0]),"state_dimension":int(states.shape[1]),"modeled_coordinate_dimension":int(modeled_states.shape[1]),"input":"domain","state_normalization":bool(normalize),"state_rescale":transform.metadata(),"coordinate_system":coordinate_system},"delay_embedding":{"type":"none","delays":1,"window_samples":1,"window_duration":0.}},transform

def _domain_point_data(domain_points: np.ndarray, system: Any, dictionary: ObservableDictionary, domain_time: float,
                       *, normalize: bool, rescale: float | Sequence[float] | None, gamma_method: GammaMethod,
                       coordinate_observable=None,
                       coordinate_system: str="dimensional",
                       directional_derivative_backend: DirectionalDerivativeBackend="directional-finite-difference"):
    """Prepare an explicit, possibly nonuniform set of state-space samples.

    This keeps externally generated domains as data rather than attempting to
    infer and regenerate their mesh.  It is especially useful for parity
    experiments where endpoint and colon-grid behavior are part of the
    numerical provenance.
    """
    try: from dynsyslib.interfaces import ODESystem
    except ImportError as error: raise ImportError("Domain-point Koopman-with-SINDy requires the optional dynsyslib package.") from error
    if not isinstance(system,ODESystem):
        raise TypeError("domain_points mode requires a dynsyslib.ODESystem instance so its equations of motion are available.")
    states=np.asarray(domain_points,dtype=float)
    expected_dimension=int(getattr(system,"state_dimension",np.asarray(system.initial_state).size))
    if states.ndim!=2 or states.shape[0]<3 or states.shape[1]!=expected_dimension:
        raise ValueError("domain_points must contain at least three sample-major states aligned with the system dimension.")
    if np.any(~np.isfinite(states)):
        raise ValueError("domain_points must be finite.")
    rhs=system.dimensionless_rhs if coordinate_system=="dimensionless" and hasattr(system,"dimensionless_rhs") else system.rhs
    derivatives=np.asarray([rhs(float(domain_time),row) for row in states],dtype=float)
    modeled_states=states if coordinate_observable is None else coordinate_observable.evaluate(states)
    modeled_derivatives=derivatives if coordinate_observable is None else np.vstack([
        coordinate_observable.jacobian(row)@direction for row,direction in zip(states,derivatives)
    ])
    transform=_state_transform(modeled_states,normalize=normalize,rescale=rescale)
    theta=_features(dictionary,transform,modeled_states)
    gamma=_directional_derivative(dictionary,modeled_states,modeled_derivatives,transform=transform,
        backend=directional_derivative_backend)
    parameters=asdict(system) if is_dataclass(system) else {key:repr(value) for key,value in vars(system).items()}
    reference=DynamicalSystemReference(None,getattr(system,"name",type(system).__name__),parameters)
    metadata={"input":"domain_points","domain_time":float(domain_time),"gamma_method":gamma_method,
        "optional_library_dependency_used":True,
        "preprocessing":{"sample_count":int(states.shape[0]),"state_dimension":int(states.shape[1]),
            "modeled_coordinate_dimension":int(modeled_states.shape[1]),"input":"domain_points",
            "state_normalization":bool(normalize),"state_rescale":transform.metadata(),
            "coordinate_system":coordinate_system},
        "delay_embedding":{"type":"none","delays":1,"window_samples":1,"window_duration":0.}}
    return states,theta,gamma,reference,metadata,transform

def _system_reference(system: Any) -> DynamicalSystemReference:
    parameters=asdict(system) if is_dataclass(system) else {key:repr(value) for key,value in vars(system).items()}
    return DynamicalSystemReference(None,getattr(system,"name",type(system).__name__),parameters)

def _trajectory_data(trajectories: Sequence[str | Path | Trajectory] | TimeSeries, dictionary: ObservableDictionary, system: Any | None=None,
                     *, normalize: bool, rescale: float | Sequence[float] | None, gamma_method: GammaMethod, coordinate_observable=None,
                     directional_derivative_backend: DirectionalDerivativeBackend="directional-finite-difference"):
    sources=_coerce_trajectories(trajectories); states=[]; theta=[]; gamma=[]
    expected_dimension=None if system is None else int(getattr(system,"state_dimension",np.asarray(system.initial_state).size))
    raw_states=[np.asarray(source.values,dtype=float) for source in sources]
    modeled_states=[x if coordinate_observable is None else coordinate_observable.evaluate(x) for x in raw_states]
    transform=_state_transform(np.vstack(modeled_states),normalize=normalize,rescale=rescale)
    for source,x,y in zip(sources,raw_states,modeled_states):
        t=np.asarray(source.time,dtype=float); features=_features(dictionary,transform,y)
        if x.ndim!=2 or t.size<3 or t.size!=x.shape[0] or np.any(np.diff(t)<=0):
            raise ValueError("Each Koopman-with-SINDy trajectory must have at least three increasing-time samples and a two-dimensional (sample, state) value array.")
        if expected_dimension is not None and x.shape[1]!=expected_dimension:
            raise ValueError(f"Pooled/individual trajectory state dimension {x.shape[1]} does not match {getattr(system,'name','system')} state dimension {expected_dimension}.")
        states.append(x); theta.append(features)
        coordinate_system=str((source.metadata.get("postprocessing_preprocessing") or {}).get("coordinate_system","dimensional"))
        gamma.append(np.gradient(features,t,axis=0,edge_order=2) if gamma_method=="finite_difference" else _rhs_gamma(
            dictionary,transform,system,t,x,coordinate_system=coordinate_system,
            coordinate_observable=coordinate_observable,
            directional_derivative_backend=directional_derivative_backend))
    source_name=None if sources[0].source is None else str(sources[0].source)
    reference=_system_reference(system) if system is not None else DynamicalSystemReference.from_archive_metadata(sources[0].metadata,source_name)
    first_steps=np.diff(sources[0].time)
    analysis_delta_t=float(first_steps[0]) if first_steps.size and np.allclose(first_steps,first_steps[0],rtol=1e-9,atol=1e-12) else None
    preprocessing={"sample_count":int(sum(item.time.size for item in sources)),"state_dimension":int(sources[0].values.shape[1]),"modeled_coordinate_dimension":int(modeled_states[0].shape[1]),"analysis_delta_t":analysis_delta_t,"input":"trajectories","pooled_input":bool(any(item.metadata.get('pooled_input',False) for item in sources)),"trajectory_boundaries_preserved":True,"state_normalization":bool(normalize),"state_rescale":transform.metadata()}
    # Preserve explicit transformations supplied by an upstream analysis stage,
    # such as decimation or state normalization before Hankel-DMD.
    preprocessing.update(sources[0].metadata.get("postprocessing_preprocessing",{}))
    return np.vstack(states),np.vstack(theta),np.vstack(gamma),reference,{"input":"trajectories","trajectory_archives":[None if item.source is None else str(item.source) for item in sources],"trajectory_count":len(sources),"gamma_method":gamma_method,"gamma_rhs_coordinate_system":str((sources[0].metadata.get("postprocessing_preprocessing") or {}).get("coordinate_system","dimensional")),"optional_library_dependency_used":gamma_method=="system_rhs","preprocessing":preprocessing,"delay_embedding":{"type":"none","delays":1,"window_samples":1,"window_duration":0.}},sources[0],transform

def _eigenvalues(spectrum: KoopmanSpectrum | None, eigenvalues: np.ndarray | None, eigenvalue_domain: EigenvalueDomain | None, delta_t: float | None, reference: DynamicalSystemReference):
    if spectrum is not None:
        if eigenvalues is not None or eigenvalue_domain is not None or delta_t is not None: raise ValueError("Provide either spectrum or numerical eigenvalues, not both.")
        return list(spectrum.eigenvalues),spectrum.system
    if eigenvalues is None or eigenvalue_domain not in {"discrete","continuous"}:
        raise ValueError("Provide a KoopmanSpectrum or numerical eigenvalues with eigenvalue_domain='discrete' or 'continuous'.")
    if delta_t is None: raise ValueError("Numerical eigenvalues require delta_t to construct the paired continuous/discrete Koopman eigenvalues.")
    build=KoopmanEigenvalue.from_discrete if eigenvalue_domain=="discrete" else KoopmanEigenvalue.from_continuous
    return [build(value,delta_t,reference,{"method":"koopman-with-sindy input","index":index}) for index,value in enumerate(np.asarray(eigenvalues,dtype=complex).reshape(-1))],reference

def _eigenvalue_transfer(values: list[KoopmanEigenvalue]) -> dict:
    """Validate paired spectral domains and declare the generator input domain."""
    if not values: raise ValueError("Koopman-with-SINDy needs at least one eigenvalue.")
    discrete_count=continuous_count=0; forward_errors=[]; inverse_errors=[]
    for value in values:
        if value.original_domain=="discrete":
            discrete_count+=1; expected=np.log(value.discrete_time)/value.delta_t
            if not np.isfinite(expected) or not np.isclose(value.continuous_time,expected,rtol=1e-12,atol=1e-12):
                raise ValueError("Inconsistent Koopman eigenvalue pair: continuous_time must equal log(discrete_time)/delta_t for a discrete H-DMD source.")
        else:
            continuous_count+=1; expected=np.exp(value.continuous_time*value.delta_t)
            if not np.isfinite(expected) or not np.isclose(value.discrete_time,expected,rtol=1e-12,atol=1e-12):
                raise ValueError("Inconsistent Koopman eigenvalue pair: discrete_time must equal exp(continuous_time*delta_t) for a continuous source.")
        forward_errors.append(abs(np.exp(value.continuous_time*value.delta_t)-value.discrete_time))
        if value.original_domain=="discrete": inverse_errors.append(abs(np.log(value.discrete_time)/value.delta_t-value.continuous_time))
    return {"source_original_domains":{"discrete":discrete_count,"continuous":continuous_count},
            "generator_eigenvalue_domain":"continuous","generator_conversion":"lambda=log(mu)/delta_t for discrete-source eigenvalues",
            "display_eigenvalue_domain":"discrete","delta_t_values":sorted(set(float(value.delta_t) for value in values)),
            "max_exp_roundtrip_error":float(max(forward_errors)),"max_log_conversion_error":float(max(inverse_errors,default=0.))}

def koopman_with_sindy(*, spectrum: KoopmanSpectrum | None=None, eigenvalues: np.ndarray | None=None, eigenvalue_domain: EigenvalueDomain | None=None,
                       delta_t: float | None=None, trajectories: Sequence[str | Path | Trajectory] | TimeSeries | None=None,
                       domain: tuple[np.ndarray, np.ndarray] | None=None, domain_points: np.ndarray | None=None,
                       system: Any | None=None, domain_density: int=20, domain_time: float=0.,
                       dictionary: ObservableDictionary | None=None, polynomial_degree: int=4, threshold: float=.05, max_iterations: int=20,
                       normalize: bool=True, rescale: float | Sequence[float] | None=None, gamma_method: GammaMethod="system_rhs",
                       compute_training_errors: bool=False, reuse_conjugate_fits: bool=True, null_space_backend: NullSpaceBackend="svd",
                       directional_derivative_backend: DirectionalDerivativeBackend="directional-finite-difference",
                       sparsification_policy: SparsificationPolicy="iterative-relative",
                       compute_training_values: bool=True,
                       max_memory_mb: float|None=1024.0, coordinate_observable=None,
                       coordinate_system: str|None=None) -> KoopmanSINDyResult:
    """Construct sparse generator-PDE eigenfunctions from supplied eigenvalues.

    By default, states are z-score normalized before a fourth-order polynomial
    dictionary is evaluated. ``rescale`` can add a positive scalar or one
    positive multiplier per state component after normalization.
    """
    started=perf_counter()
    if sum(item is not None for item in (trajectories,domain,domain_points))!=1:
        raise ValueError("Provide exactly one of trajectories, domain, or domain_points.")
    if system is None and spectrum is not None: system=spectrum.dynamical_system
    resolved_coordinate_system=coordinate_system
    if resolved_coordinate_system is None and system is not None and (domain is not None or domain_points is not None):
        resolved_coordinate_system=("dimensionless" if bool(getattr(system,"prefer_dimensionless_form",False)) else "dimensional")
    if resolved_coordinate_system not in {None,"dimensionless","dimensional"}:
        raise ValueError("coordinate_system must be 'dimensionless', 'dimensional', or None.")
    if system is not None:
        from dynsyslib.observables import coerce_coordinate_observable
        state_dimension=int(getattr(system,"state_dimension",np.asarray(system.initial_state).size))
        coordinate_observable=coerce_coordinate_observable(coordinate_observable,state_dimension,
            system=system,coordinate_system=resolved_coordinate_system)
    elif coordinate_observable is not None:
        from dynsyslib.observables import Observable
        if not isinstance(coordinate_observable,Observable):
            raise ValueError("Named or stacked observable coordinates require system= for resolution.")
    if gamma_method not in {"system_rhs","finite_difference"}:
        raise ValueError("gamma_method must be 'system_rhs' or 'finite_difference'.")
    if null_space_backend not in {"svd","gram-eigh"}: raise ValueError("null_space_backend must be 'svd' or 'gram-eigh'.")
    if directional_derivative_backend not in {"directional-finite-difference","dictionary-gradient"}:
        raise ValueError("directional_derivative_backend must be 'directional-finite-difference' or 'dictionary-gradient'.")
    if sparsification_policy not in {"iterative-relative","matlab-single-absolute"}:
        raise ValueError("sparsification_policy must be 'iterative-relative' or 'matlab-single-absolute'.")
    if sparsification_policy=="matlab-single-absolute" and null_space_backend!="svd":
        raise ValueError("matlab-single-absolute sparsification requires null_space_backend='svd'.")
    if trajectories is not None and gamma_method=="system_rhs" and system is None:
        raise ValueError("gamma_method='system_rhs' requires the originating dynsyslib ODE passed as system=. Use gamma_method='finite_difference' for archive-only analysis.")
    if polynomial_degree<1: raise ValueError("polynomial_degree must be at least one.")
    if max_iterations<1: raise ValueError("max_iterations must be at least one.")
    features=dictionary or PolynomialDictionary(degree=polynomial_degree)
    dictionary_metadata=_dictionary_metadata(features)
    primary=None
    if trajectories is not None:
        if isinstance(trajectories,TimeSeries) and trajectories.is_pooled and system is None:
            raise ValueError("Pooled TimeSeries input requires the originating dynsyslib system so Koopman-with-SINDy can validate state-space dimension.")
        states,theta,gamma,reference,metadata,primary,transform=_trajectory_data(
            trajectories,features,system,normalize=normalize,rescale=rescale,
            gamma_method=gamma_method,coordinate_observable=coordinate_observable,
            directional_derivative_backend=directional_derivative_backend)
    elif domain is not None:
        if system is None: raise ValueError("Domain-based Koopman-with-SINDy requires a dynsyslib system instance.")
        states,theta,gamma,reference,metadata,transform=_domain_data(
            domain,system,features,domain_density,domain_time,normalize=normalize,
            rescale=rescale,gamma_method=gamma_method,coordinate_observable=coordinate_observable,
            coordinate_system=resolved_coordinate_system or "dimensional",
            directional_derivative_backend=directional_derivative_backend)
    else:
        if system is None: raise ValueError("Domain-point Koopman-with-SINDy requires a dynsyslib system instance.")
        states,theta,gamma,reference,metadata,transform=_domain_point_data(
            domain_points,system,features,domain_time,normalize=normalize,rescale=rescale,
            gamma_method=gamma_method,coordinate_observable=coordinate_observable,
            coordinate_system=resolved_coordinate_system or "dimensional",
            directional_derivative_backend=directional_derivative_backend)
    data_prepared=perf_counter()
    training_coordinates=(np.asarray(states,dtype=float) if coordinate_observable is None else
                          np.asarray(coordinate_observable.evaluate(states),dtype=float))
    values,reference=_eigenvalues(spectrum,eigenvalues,eigenvalue_domain,delta_t,reference)
    eigenvalue_transfer=_eigenvalue_transfer(values)
    # Each eigenfunction fit factors a dense (samples x features) complex
    # generator matrix.  Check before starting the expensive SVD loop.
    require_memory_budget("Koopman-with-SINDy",theta.shape,gamma.shape,(theta.shape[0],theta.shape[1]),
                          max_memory_mb=max_memory_mb,workspace_factor=6.0,itemsize=np.dtype(complex).itemsize)
    fit_plan=_conjugate_fit_plan(values,theta,gamma,enabled=reuse_conjugate_fits)
    coefficient_rows=[None]*len(values); direct_factorizations=0
    backend_diagnostics={"requested":null_space_backend,"svd_subproblems":0,"gram_eigh_subproblems":0,"gram_eigh_fallbacks":0}
    gram_terms=_generator_gram_terms(theta,gamma) if null_space_backend=="gram-eigh" else None
    fit_prepared=perf_counter()
    for index,value in enumerate(values):
        source_index,conjugated=fit_plan[index]
        if conjugated:
            coefficient_rows[index]=np.conjugate(coefficient_rows[source_index])
        else:
            if null_space_backend=="svd":
                matrix=gamma-value.continuous_time*theta
                coefficient_rows[index]=(_matlab_single_absolute_null_vector(
                    matrix,threshold,diagnostics=backend_diagnostics)
                    if sparsification_policy=="matlab-single-absolute" else
                    _sparse_null_vector(matrix,threshold,max_iterations,diagnostics=backend_diagnostics))
            else:
                eigenvalue=value.continuous_time
                coefficient_rows[index]=_sparse_null_vector_from_gram(_generator_gram(gram_terms,eigenvalue),threshold,max_iterations,
                    matrix_factory=lambda eigenvalue=eigenvalue: gamma-eigenvalue*theta,diagnostics=backend_diagnostics)
            direct_factorizations+=1
    fitted=perf_counter()
    coefficients=np.vstack(coefficient_rows)
    if compute_training_values:
        require_memory_budget("Koopman-with-SINDy training values",(theta.shape[0],coefficients.shape[0]),
                              max_memory_mb=max_memory_mb,workspace_factor=1.,itemsize=np.dtype(complex).itemsize)
        function_values=theta@coefficients.T
    else: function_values=None
    evaluated=perf_counter()
    if null_space_backend=="gram-eigh":
        gram_count=backend_diagnostics["gram_eigh_subproblems"]; fallback_count=backend_diagnostics["gram_eigh_fallbacks"]
        backend_diagnostics["effective"]="svd-fallback" if gram_count==0 else ("gram-eigh-with-fallback" if fallback_count else "gram-eigh")
        if not np.isfinite(backend_diagnostics.get("minimum_relative_gap",np.inf)): backend_diagnostics["minimum_relative_gap"]=None
    else: backend_diagnostics["effective"]="svd"
    acceleration={"conjugate_fit_reuse_requested":bool(reuse_conjugate_fits),
                  "conjugate_fit_reuse_active":bool(direct_factorizations<len(values)),
                  "direct_eigenvalue_fits":int(direct_factorizations),
                  "reused_conjugate_fits":int(len(values)-direct_factorizations),
                  "eigenfunction_evaluation":"batched-matrix-multiply" if compute_training_values else "deferred-polynomial-evaluation",
                  "training_values_materialized":bool(compute_training_values),"null_space_backend":backend_diagnostics}
    functions=[]
    for index,(value,coefficient) in enumerate(zip(values,coefficients)):
        psi_values=None if function_values is None else function_values[:,index]
        def psi(x, coeff=coefficient, basis=features, state_transform=transform): return complex(_features(basis,state_transform,np.asarray(x,dtype=float).reshape(1,-1))[0]@coeff)
        def residual_operator(x, velocity, coeff=coefficient, basis=features, eigenvalue=value, state_transform=transform,
                              derivative_backend=directional_derivative_backend,coordinates=coordinate_observable):
            raw=np.asarray(x,dtype=float); raw_velocity=np.asarray(velocity,dtype=float)
            if coordinates is None:
                modeled=raw; modeled_velocity=raw_velocity
            else:
                modeled=np.asarray(coordinates.evaluate(raw),dtype=float)
                modeled_velocity=np.vstack([coordinates.jacobian(row)@direction
                    for row,direction in zip(raw,raw_velocity)])
            return _directional_derivative(basis,modeled,modeled_velocity,transform=state_transform,
                backend=derivative_backend)@coeff-eigenvalue.continuous_time*(_features(basis,state_transform,modeled)@coeff)
        scale_factor=transform.multiplier/transform.scale
        def gradient(x, coeff=coefficient, basis=features, state_transform=transform, factor=scale_factor):
            transformed=state_transform(np.asarray(x,dtype=float).reshape(1,-1))
            return (np.asarray(coeff)@np.asarray(basis.gradient(transformed))[0])*factor
        function=KoopmanEigenfunction(value,"sparse dictionary coefficients",values=psi_values,coefficients=coefficient,state_callable=psi,gradient_callable=gradient,residual_operator=residual_operator,
            metadata={"method":"koopman-with-sindy","index":index,"dictionary_class":type(features).__name__,"dictionary":dictionary_metadata,"polynomial_degree":polynomial_degree if dictionary is None else None,"threshold":threshold,"max_iterations":max_iterations,"feature_count":int(theta.shape[1]),"state_preprocessing":transform.metadata(),"generator_eigenvalue_domain":"continuous",
                      "coordinate_observable":None if coordinate_observable is None else coordinate_observable.summary(),
                      "eigenfunction_coordinate_dimension":int(training_coordinates.shape[1])},
            representation_kind="analytic",coordinate_observable=coordinate_observable)
        if compute_training_errors:
            residual=(gamma-value.continuous_time*theta)@coefficient; magnitude=np.abs(residual)
            function.metadata["training_residual"]={"rms":float(np.sqrt(np.mean(magnitude**2))),"max_abs":float(np.max(magnitude))}
        functions.append(function)
    functions_built=perf_counter()
    stage_timing={"data_and_dictionary_seconds":data_prepared-started,"fit_preparation_seconds":fit_prepared-data_prepared,
                  "sparse_nullspace_fit_seconds":fitted-fit_prepared,"batched_eigenfunction_evaluation_seconds":evaluated-fitted,
                  "eigenfunction_object_build_seconds":functions_built-evaluated,"pre_krom_total_seconds":functions_built-started,
                  "clock":"time.perf_counter"}
    training_domain={"encoding":"axis-aligned-training-domain/v1","minimum":np.min(states,axis=0).tolist(),
        "maximum":np.max(states,axis=0).tolist(),"sample_count":int(states.shape[0]),
        "coordinate_system":str((metadata.get("preprocessing") or {}).get("coordinate_system","dimensional")),
        "source":"Koopman-with-SINDy training states"}
    coordinate_training_domain={"encoding":"axis-aligned-observable-coordinate-domain/v1",
        "minimum":np.min(training_coordinates,axis=0).tolist(),"maximum":np.max(training_coordinates,axis=0).tolist(),
        "sample_count":int(training_coordinates.shape[0]),"coordinate_observable":None if coordinate_observable is None else coordinate_observable.summary(),
        "source":"Koopman-with-SINDy modeled observable coordinates"}
    output=KoopmanSpectrum(reference,values,functions,[],{"method":"koopman-with-sindy","dictionary_class":type(features).__name__,"dictionary":dictionary_metadata,"polynomial_degree":polynomial_degree if dictionary is None else None,"threshold":threshold,"max_iterations":max_iterations,"feature_count":int(theta.shape[1]),"normalize":bool(normalize),"rescale":None if rescale is None else np.asarray(rescale,dtype=float).tolist(),"gamma_method":gamma_method,"null_space_backend":null_space_backend,"directional_derivative_backend":directional_derivative_backend,"sparsification_policy":sparsification_policy,"compute_training_values":bool(compute_training_values),"max_memory_mb":max_memory_mb,"state_preprocessing":transform.metadata(),"coordinate_observable":None if coordinate_observable is None else coordinate_observable.summary(),"eigenvalue_transfer":eigenvalue_transfer,"fit_acceleration":acceleration,"training_domain":training_domain,"coordinate_training_domain":coordinate_training_domain,**metadata},dynamical_system=system,eigenvalue_source=spectrum)
    from rom import KROM
    krom=KROM(output,system,coordinate_observable=coordinate_observable)
    stage_timing["krom_construction_seconds"]=perf_counter()-functions_built
    stage_timing["total_seconds"]=perf_counter()-started
    output.metadata["pipeline_timing"]=stage_timing
    return KoopmanSINDyResult(output,states,training_coordinates,function_values,coefficients,{"primary_source":None if primary is None or primary.source is None else str(primary.source),**output.metadata,"krom":krom.metadata},krom)
