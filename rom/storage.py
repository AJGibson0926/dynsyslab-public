"""Load and reserve executable coefficient-based KROM archive subsets."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import numpy as np


def _metadata(path: Path) -> tuple[dict[str,Any],dict[str,Any]]:
    with np.load(path,allow_pickle=False) as data:
        outer=json.loads(str(data["metadata_json"]))
    return outer,dict(outer.get("analysis_metadata") or {})


def _eigenvalues(path: Path, analysis: dict[str,Any]) -> np.ndarray:
    with np.load(path,allow_pickle=False) as data:
        if "eigenvalues" in data.files: return np.asarray(data["eigenvalues"],dtype=complex)
    source=analysis.get("eigenvalue_source_archive")
    if source is None: raise KeyError("KROM archive has neither eigenvalues nor an exact eigenvalue-source archive.")
    source_path=Path(source).resolve(); _,source_analysis=_metadata(source_path)
    return _eigenvalues(source_path,source_analysis)


def _control_representation_arrays(krom) -> dict[str,np.ndarray]:
    """Serialize self-contained fitted input maps."""
    arrays={}
    constant=(krom.control_representations or {}).get("constant")
    if constant is not None:
        arrays["constant_control_matrix"]=constant.matrix
        if constant.discrete_matrix is not None:
            arrays["constant_control_matrix_discrete"]=constant.discrete_matrix
    bilinear=(krom.control_representations or {}).get("bilinear")
    if bilinear is not None:
        arrays["bilinear_control_operators"]=bilinear.operators
        arrays["bilinear_control_offsets"]=bilinear.offsets
        if bilinear.conjugate_operators is not None:
            arrays["bilinear_control_conjugate_operators"]=bilinear.conjugate_operators
    return arrays


def _observable_reconstruction_arrays(krom) -> dict[str,np.ndarray]:
    """Serialize every named observable map while retaining legacy default keys."""
    arrays={}
    if krom.observable_reconstruction is not None:
        arrays["observable_reconstruction_coefficients"]=krom.observable_reconstruction.coefficients
        arrays["observable_reconstruction_intercept"]=krom.observable_reconstruction.intercept
    for index,name in enumerate(sorted(krom.observable_reconstructions)):
        reconstruction=krom.observable_reconstructions[name]
        arrays[f"observable_reconstruction_{index:03d}_coefficients"]=reconstruction.coefficients
        arrays[f"observable_reconstruction_{index:03d}_intercept"]=reconstruction.intercept
    return arrays


def load_krom_archive(path: str|Path):
    """Reconstruct an executable polynomial-eigenfunction KROM from an archive."""
    source_path=Path(path).resolve(); outer,analysis=_metadata(source_path)
    krom_metadata=dict(analysis.get("krom") or {})
    registry_metadata=(analysis.get("observable_reconstructions") or
                       krom_metadata.get("observable_reconstructions") or {})
    with np.load(source_path,allow_pickle=False) as data:
        coefficients=np.asarray(data["coefficients"],dtype=complex)
        training_states=None if "training_states" not in data.files else np.asarray(data["training_states"],dtype=float)
        training_coordinates=None if "training_coordinates" not in data.files else np.asarray(data["training_coordinates"],dtype=float)
        observation_coefficients=None if "observation_coefficients" not in data.files else np.asarray(data["observation_coefficients"],dtype=float)
        observation_intercept=None if "observation_intercept" not in data.files else np.asarray(data["observation_intercept"],dtype=float)
        reconstruction_coefficients=None if "observable_reconstruction_coefficients" not in data.files else np.asarray(data["observable_reconstruction_coefficients"],dtype=complex)
        reconstruction_intercept=None if "observable_reconstruction_intercept" not in data.files else np.asarray(data["observable_reconstruction_intercept"],dtype=float)
        constant_matrix=None if "constant_control_matrix" not in data.files else np.asarray(data["constant_control_matrix"],dtype=complex)
        constant_discrete=None if "constant_control_matrix_discrete" not in data.files else np.asarray(data["constant_control_matrix_discrete"],dtype=complex)
        bilinear_operators=None if "bilinear_control_operators" not in data.files else np.asarray(data["bilinear_control_operators"],dtype=complex)
        bilinear_offsets=None if "bilinear_control_offsets" not in data.files else np.asarray(data["bilinear_control_offsets"],dtype=complex)
        bilinear_conjugate=None if "bilinear_control_conjugate_operators" not in data.files else np.asarray(data["bilinear_control_conjugate_operators"],dtype=complex)
        reconstruction_registry_arrays={}
        for index,name in enumerate(sorted(registry_metadata)):
            coefficient_key=f"observable_reconstruction_{index:03d}_coefficients"
            intercept_key=f"observable_reconstruction_{index:03d}_intercept"
            if coefficient_key in data.files and intercept_key in data.files:
                reconstruction_registry_arrays[name]=(
                    np.asarray(data[coefficient_key],dtype=complex),
                    np.asarray(data[intercept_key],dtype=float))
    if bilinear_conjugate is not None and bilinear_operators is None:
        raise ValueError("KROM archive has conjugate-bilinear operators without a primary bilinear map.")
    discrete=_eigenvalues(source_path,analysis)
    if coefficients.ndim!=2 or coefficients.shape[0]!=discrete.size:
        raise ValueError("Archived KROM coefficients and eigenvalues are not aligned.")
    primary=analysis.get("primary_source") or outer.get("source_archive")
    if primary is None: raise KeyError("Archived KROM has no primary solution lineage.")
    with np.load(Path(primary),allow_pickle=False) as data:
        solution_metadata=json.loads(str(data["metadata_json"]))
    from dynsyslib import Observable, observable_from_summary, resolve_archived_system
    from spectral.dictionary import PolynomialDictionary
    from spectral.contracts import DynamicalSystemReference, KoopmanEigenfunction, KoopmanEigenvalue, KoopmanSpectrum
    from .koopman import KROM
    from .control_inputs import BilinearControlMatrix, ConstantControlMatrix
    from .observation import ObservableReconstruction, StateObservationMap
    system=resolve_archived_system(solution_metadata.get("system"),solution_metadata.get("system_parameters"))
    observable_metadata=analysis.get("coordinate_observable") or krom_metadata.get("coordinate_observable")
    coordinate_observable=(Observable.identity(int(np.asarray(system.initial_state).size))
        if observable_metadata is None else observable_from_summary(observable_metadata))
    if training_coordinates is None and training_states is not None:
        training_coordinates=np.asarray(coordinate_observable.evaluate(training_states),dtype=float)
    reference=DynamicalSystemReference.from_archive_metadata(solution_metadata,str(Path(primary).resolve()))
    delta_t=float((analysis.get("krom") or {}).get("delta_t") or solution_metadata.get("delta_t"))
    values=[KoopmanEigenvalue.from_discrete(value,delta_t,reference,{"archive_index":index}) for index,value in enumerate(discrete)]
    dictionary_metadata=dict(analysis.get("dictionary") or {})
    degree=int(dictionary_metadata.get("degree") or analysis.get("polynomial_degree") or 4)
    include_constant=bool(dictionary_metadata.get("include_constant",True))
    dictionary_metadata={"encoding":"polynomial-dictionary/v1","degree":degree,"include_constant":include_constant}
    dictionary=PolynomialDictionary(degree=degree,include_constant=include_constant)
    preprocessing=dict(analysis.get("state_preprocessing") or {})
    center=np.asarray(preprocessing["center"],dtype=float); scale=np.asarray(preprocessing["scale"],dtype=float); multiplier=np.asarray(preprocessing["multiplier"],dtype=float)
    if center.shape!=scale.shape or center.shape!=multiplier.shape or np.any(scale<=0) or np.any(multiplier<=0):
        raise ValueError("Archived KROM state preprocessing is invalid.")
    def transform(state): return (np.asarray(state,dtype=float)-center)/scale*multiplier
    def eigenfunction(eigenvalue,coefficient,index,*,role="retained"):
        def psi(state,coeff=coefficient): return complex(dictionary.transform(transform(state).reshape(1,-1))[0]@coeff)
        def gradient(state,coeff=coefficient):
            transformed=transform(state).reshape(1,-1)
            return (np.asarray(coeff)@np.asarray(dictionary.gradient(transformed))[0])*(multiplier/scale)
        return KoopmanEigenfunction(eigenvalue,"archived sparse polynomial coefficients",
            coefficients=np.array(coefficient,copy=True),state_callable=psi,gradient_callable=gradient,
            metadata={"method":"koopman-with-sindy","archive_index":index,"archive_role":role,
                      "dictionary":dictionary_metadata,"state_preprocessing":preprocessing,
                      "coordinate_observable":coordinate_observable.summary()},
            coordinate_observable=coordinate_observable)
    functions=[eigenfunction(value,coefficient,index)
        for index,(value,coefficient) in enumerate(zip(values,coefficients))]
    coordinate_metadata=dict(analysis.get("preprocessing") or {})
    if not coordinate_metadata and solution_metadata.get("integration_form") == "dimensionless":
        coordinate_metadata={"coordinate_system":"dimensionless","source_integration_form":"dimensionless"}
    spectrum_metadata={**analysis,"method":str(analysis.get("method") or "koopman-with-sindy"),"dictionary":dictionary_metadata,"state_preprocessing":preprocessing,
                       "preprocessing":coordinate_metadata,
                       "loaded_from_archive":str(source_path)}
    if "training_domain" not in spectrum_metadata and training_states is not None:
        spectrum_metadata["training_domain"]={"encoding":"axis-aligned-training-domain/v1",
            "minimum":np.min(training_states,axis=0).tolist(),"maximum":np.max(training_states,axis=0).tolist(),
            "sample_count":int(training_states.shape[0]),"coordinate_system":coordinate_metadata.get("coordinate_system","dimensional"),
            "source":"archived-training-states"}
    spectrum=KoopmanSpectrum(reference,values,functions,[],spectrum_metadata,dynamical_system=system)
    observation=None
    observation_metadata=analysis.get("observation_map")
    if observation_metadata is not None:
        if observation_coefficients is None or observation_intercept is None:
            raise ValueError("KROM archive declares an observation map but omits its arrays.")
        observation=StateObservationMap(observation_coefficients,observation_intercept,
            dict(observation_metadata.get("quality") or {}),
            coordinate_system=str(observation_metadata.get("coordinate_system") or "model"),
            encoding=str(observation_metadata.get("encoding") or "real-imaginary-split-affine/v1"))
    elif observation_coefficients is not None or observation_intercept is not None:
        raise ValueError("KROM archive has incomplete observation-map metadata.")
    reconstruction=None; reconstructions={}
    reconstruction_metadata=analysis.get("observable_reconstruction") or krom_metadata.get("observable_reconstruction")
    if reconstruction_metadata is not None:
        if reconstruction_coefficients is None or reconstruction_intercept is None:
            raise ValueError("KROM archive declares observable reconstruction C but omits its arrays.")
        reconstruction_observable=observable_from_summary(reconstruction_metadata["observable"])
        reconstruction=ObservableReconstruction(reconstruction_observable,reconstruction_coefficients,
            reconstruction_intercept,dict(reconstruction_metadata.get("quality") or {}),
            dict(reconstruction_metadata.get("provenance") or {}),
            encoding=str(reconstruction_metadata.get("encoding") or "nonredundant-complex-koopman-modes/v1"))
    for name,summary in registry_metadata.items():
        arrays=reconstruction_registry_arrays.get(name)
        if arrays is None:
            if reconstruction is not None and name==str(reconstruction.observable.name):
                reconstructions[name]=reconstruction; continue
            raise ValueError(f"KROM archive declares observable reconstruction {name!r} but omits its arrays.")
        observable=observable_from_summary(summary["observable"])
        reconstructions[name]=ObservableReconstruction(observable,arrays[0],arrays[1],
            dict(summary.get("quality") or {}),dict(summary.get("provenance") or {}),
            encoding=str(summary.get("encoding") or "nonredundant-complex-koopman-modes/v1"))
    if reconstruction is not None and not reconstructions:
        reconstructions.setdefault(str(reconstruction.observable.name),reconstruction)
    representations={}; representation_metadata=dict(krom_metadata.get("control_representations") or {})
    if constant_matrix is not None:
        details=dict(representation_metadata.get("constant") or {})
        representations["constant"]=ConstantControlMatrix(constant_matrix,"continuous",details.get("delta_t"),
            dict(details.get("provenance") or {}),constant_discrete)
    if bilinear_operators is not None or bilinear_offsets is not None:
        if bilinear_operators is None or bilinear_offsets is None:
            raise ValueError("KROM archive has incomplete bilinear control arrays.")
        details=dict(representation_metadata.get("bilinear") or {})
        representations["bilinear"]=BilinearControlMatrix(bilinear_operators,bilinear_offsets,"continuous",
            details.get("delta_t"),dict(details.get("provenance") or {}),bilinear_conjugate)
    default_representation=str(krom_metadata.get("default_control_representation") or "exact")
    if default_representation!="exact" and default_representation not in representations:
        default_representation="exact"
    default_observable=(krom_metadata.get("default_observable_reconstruction") or
                        analysis.get("default_observable_reconstruction"))
    if default_observable in reconstructions:
        reconstruction=reconstructions[default_observable]
    return KROM(spectrum,system,observation_map=observation,
        coordinate_observable=coordinate_observable,observable_reconstruction=reconstruction,
        observable_reconstructions=reconstructions,
        default_observable_reconstruction=default_observable,
        control_representations=representations,default_control_representation=default_representation)


def reserve_krom_archive(path: str|Path, selection, *, include_conjugate_duals: bool=True,
                         label: str="reserved-krom", selection_metadata: dict[str,Any]|None=None) -> Path:
    """Extract and immutably publish a reusable KROM subset archive."""
    source_path=Path(path).resolve(); outer,analysis=_metadata(source_path)
    source=load_krom_archive(source_path)
    derived=source.extract(selection,include_conjugate_duals=include_conjugate_duals)
    source_indices=[int(value.metadata["extraction_source_index"]) for value in derived.spectrum.eigenvalues]
    coefficients=np.vstack([np.asarray(function.coefficients,dtype=complex) for function in derived.spectrum.eigenfunctions])
    with np.load(source_path,allow_pickle=False) as data:
        training_states=np.asarray(data["training_states"],dtype=float)
        training_coordinates=(np.asarray(data["training_coordinates"],dtype=float)
            if "training_coordinates" in data.files else np.asarray(derived.coordinate_observable.evaluate(training_states),dtype=float))
    from analysis_data import load_derived_array, load_solution, save_derived
    dictionary=derived.spectrum.metadata.get("dictionary") or {}
    preprocessing=derived.spectrum.metadata.get("state_preprocessing") or {}
    try:
        from spectral.evaluation import polynomial_eigenfunction_values
        sampled=polynomial_eigenfunction_values(training_coordinates,coefficients,
            dictionary=dictionary,state_preprocessing=preprocessing)
        sampled_source="coefficient-evaluation"
    except (KeyError,TypeError,ValueError):
        try:
            sampled=load_derived_array(source_path,"eigenfunction_values",columns=source_indices)
            sampled_source="source-array"
        except KeyError:
            sampled=np.empty((training_states.shape[0],0),dtype=complex); sampled_source="unavailable"
    primary=analysis.get("primary_source") or outer.get("source_archive")
    if primary is None: raise KeyError("Source KROM has no primary solution lineage.")
    extraction=dict(derived.spectrum.metadata.get("extraction") or {})
    metadata={"method":"koopman-with-sindy-krom-extract","primary_source":str(Path(primary).resolve()),
        "source_krom_archive":str(source_path),"reserved":True,"reservation_label":label,
        "dictionary":derived.spectrum.metadata.get("dictionary"),"polynomial_degree":derived.spectrum.metadata.get("polynomial_degree"),
        "preprocessing":derived.spectrum.metadata.get("preprocessing"),"delay_embedding":derived.spectrum.metadata.get("delay_embedding"),
        "coordinate_observable":derived.coordinate_observable.summary(),
        "observable_reconstruction":None if derived.observable_reconstruction is None else derived.observable_reconstruction.summary(),
        "state_preprocessing":derived.spectrum.metadata.get("state_preprocessing"),"extraction":extraction,
        "selection":dict(selection_metadata or {}),
        "training_domain":derived.spectrum.metadata.get("training_domain") or {"encoding":"axis-aligned-training-domain/v1",
            "minimum":np.min(training_states,axis=0).tolist(),"maximum":np.max(training_states,axis=0).tolist(),
            "sample_count":int(training_states.shape[0]),"coordinate_system":derived.coordinate_system,"source":"archived-training-states"},
        "krom":derived.metadata,"storage":{"format":"koopman-krom-subset/v1","coefficient_vectors_saved":True,
        "training_states_saved":True,"eigenfunction_values_saved":bool(sampled.shape[1]),"eigenfunction_values_source":sampled_source}}
    arrays={"eigenvalues":np.asarray([value.discrete_time for value in derived.spectrum.eigenvalues]),
            "coefficients":coefficients,"training_states":training_states,
            "training_coordinates":training_coordinates}
    if sampled.shape[1]: arrays["eigenfunction_values"]=sampled
    arrays.update(_observable_reconstruction_arrays(derived))
    arrays.update(_control_representation_arrays(derived))
    return save_derived(load_solution(primary),"koopman-krom-extract",analysis_metadata=metadata,**arrays)


def reserve_control_representations_archive(
    path: str|Path,
    representations: dict[str,Any],
    *,
    default: str="exact",
    label: str="control-equipped-krom",
) -> Path:
    """Immutably publish fitted input maps on an existing KROM archive.

    The autonomous spectrum and coefficient representations are copied without
    refitting.  A later :func:`reserve_krom_archive` call can therefore exercise
    the same independent closed control-map restriction as an in-memory
    extraction.
    """
    if not representations:
        raise ValueError("At least one fitted control representation is required.")
    if set(representations)-{"constant","bilinear"}:
        raise ValueError("Stored control representations must be named 'constant' or 'bilinear'.")
    if default not in {"exact",*representations}:
        raise ValueError("The requested default control representation is unavailable.")
    source_path=Path(path).resolve(); outer,analysis=_metadata(source_path)
    model=load_krom_archive(source_path)
    for name,representation in representations.items():
        model=model.with_control_representation(name,representation)
    model.default_control_representation=default
    model._validate_control_representation()
    primary=analysis.get("primary_source") or outer.get("source_archive")
    if primary is None: raise KeyError("Source KROM has no primary solution lineage.")
    with np.load(source_path,allow_pickle=False) as data:
        if "coefficients" not in data.files or "training_states" not in data.files:
            raise KeyError("Source KROM must retain coefficients and training_states.")
        coefficients=np.asarray(data["coefficients"],dtype=complex)
        training_states=np.asarray(data["training_states"],dtype=float)
        training_coordinates=(np.asarray(data["training_coordinates"],dtype=float)
            if "training_coordinates" in data.files else np.asarray(model.coordinate_observable.evaluate(training_states),dtype=float))
    from analysis_data import load_solution, save_derived
    metadata={"method":"koopman-krom-control-representations",
        "primary_source":str(Path(primary).resolve()),"source_krom_archive":str(source_path),
        "reserved":True,"reservation_label":label,
        "dictionary":model.spectrum.metadata.get("dictionary"),
        "polynomial_degree":model.spectrum.metadata.get("polynomial_degree"),
        "preprocessing":model.spectrum.metadata.get("preprocessing"),
        "delay_embedding":model.spectrum.metadata.get("delay_embedding"),
        "coordinate_observable":model.coordinate_observable.summary(),
        "observable_reconstruction":None if model.observable_reconstruction is None else model.observable_reconstruction.summary(),
        "state_preprocessing":model.spectrum.metadata.get("state_preprocessing"),
        "training_domain":model.training_domain,"krom":model.metadata,
        "storage":{"format":"koopman-krom-control-representations/v1",
            "source_spectrum_unchanged":True,"coefficient_vectors_saved":True,
            "training_states_saved":True,"control_representations_saved":sorted(representations)}}
    arrays={"eigenvalues":_eigenvalues(source_path,analysis),"coefficients":coefficients,
        "training_states":training_states,"training_coordinates":training_coordinates}
    arrays.update(_observable_reconstruction_arrays(model))
    if model.observation_map is not None:
        arrays["observation_coefficients"]=model.observation_map.coefficients
        arrays["observation_intercept"]=model.observation_map.intercept
        metadata["observation_map"]=model.observation_map.summary()
    arrays.update(_control_representation_arrays(model))
    return save_derived(load_solution(primary),"koopman-krom-control-representations",
        analysis_metadata=metadata,**arrays)


def reserve_observation_map_archive(path: str|Path, observation_map, *,
                                    label: str="observed-krom") -> Path:
    """Immutably publish a KROM derivative with a fitted state decoder."""
    from .observation import StateObservationMap
    if not isinstance(observation_map,StateObservationMap):
        raise TypeError("observation_map must be a StateObservationMap.")
    source_path=Path(path).resolve(); outer,analysis=_metadata(source_path)
    source=load_krom_archive(source_path).with_observation_map(observation_map)
    primary=analysis.get("primary_source") or outer.get("source_archive")
    if primary is None: raise KeyError("Source KROM has no primary solution lineage.")
    with np.load(source_path,allow_pickle=False) as data:
        if "training_states" not in data.files: raise KeyError("Source KROM has no training_states array.")
        training_states=np.asarray(data["training_states"],dtype=float)
        coefficients=np.asarray(data["coefficients"],dtype=complex)
        training_coordinates=(np.asarray(data["training_coordinates"],dtype=float)
            if "training_coordinates" in data.files else np.asarray(source.coordinate_observable.evaluate(training_states),dtype=float))
    from analysis_data import load_solution, save_derived
    metadata={"method":"koopman-krom-observation-map","primary_source":str(Path(primary).resolve()),
        "source_krom_archive":str(source_path),"reserved":True,"reservation_label":label,
        "dictionary":source.spectrum.metadata.get("dictionary"),
        "polynomial_degree":source.spectrum.metadata.get("polynomial_degree"),
        "preprocessing":source.spectrum.metadata.get("preprocessing"),
        "delay_embedding":source.spectrum.metadata.get("delay_embedding"),
        "coordinate_observable":source.coordinate_observable.summary(),
        "observable_reconstruction":None if source.observable_reconstruction is None else source.observable_reconstruction.summary(),
        "state_preprocessing":source.spectrum.metadata.get("state_preprocessing"),
        "training_domain":source.training_domain,"observation_map":observation_map.summary(),
        "krom":source.metadata,"storage":{"format":"koopman-krom-observation-map/v1",
            "coefficient_vectors_saved":True,"training_states_saved":True,
            "observation_map_arrays_saved":True,"source_arrays_unchanged":True}}
    arrays={"eigenvalues":_eigenvalues(source_path,analysis),"coefficients":coefficients,
        "training_states":training_states,"training_coordinates":training_coordinates,
        "observation_coefficients":observation_map.coefficients,
        "observation_intercept":observation_map.intercept}
    arrays.update(_observable_reconstruction_arrays(source))
    arrays.update(_control_representation_arrays(source))
    return save_derived(load_solution(primary),"koopman-krom-observation-map",analysis_metadata=metadata,**arrays)


__all__=["load_krom_archive","reserve_control_representations_archive",
    "reserve_krom_archive","reserve_observation_map_archive"]
