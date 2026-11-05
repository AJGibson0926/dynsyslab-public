"""Versioned archive representations for Koopman-with-SINDy results."""
from __future__ import annotations

from typing import Any, Literal

import numpy as np


SINDyArchiveStorage=Literal["coefficient-first","full"]


def _storage_policy(policy: str) -> SINDyArchiveStorage:
    if policy not in {"coefficient-first","full"}:
        raise ValueError("archive_storage must be 'coefficient-first' or 'full'.")
    return policy


def sindy_archive_payload(result: Any, *, storage_policy: str) -> tuple[dict[str,np.ndarray],dict[str,Any]]:
    """Return nonredundant arrays and explicit virtual-array metadata."""
    policy=_storage_policy(storage_policy)
    coefficients=np.asarray(result.coefficients)
    states=np.asarray(result.training_states,dtype=float)
    coordinates=np.asarray(result.training_coordinates,dtype=float)
    raw_values=result.training_eigenfunction_values
    values=None if raw_values is None else np.asarray(raw_values)
    expected_shape=(states.shape[0],coefficients.shape[0])
    if coefficients.ndim!=2 or states.ndim!=2 or coordinates.ndim!=2 or coordinates.shape[0]!=states.shape[0] or values is not None and values.shape!=expected_shape:
        raise ValueError("Koopman-with-SINDy archive arrays have inconsistent shapes.")
    reused=result.spectrum.eigenvalue_source is not None
    arrays={"coefficients":coefficients,"training_states":states,"training_coordinates":coordinates}
    if not reused:
        arrays["eigenvalues"]=np.asarray([value.discrete_time for value in result.spectrum.eigenvalues])
    if result.krom is not None:
        reconstruction=result.krom.observable_reconstruction
        if reconstruction is not None:
            arrays["observable_reconstruction_coefficients"]=reconstruction.coefficients
            arrays["observable_reconstruction_intercept"]=reconstruction.intercept
        constant=(result.krom.control_representations or {}).get("constant")
        if constant is not None:
            arrays["constant_control_matrix"]=constant.matrix
            if constant.discrete_matrix is not None:
                arrays["constant_control_matrix_discrete"]=constant.discrete_matrix
        bilinear=(result.krom.control_representations or {}).get("bilinear")
        if bilinear is not None:
            arrays["bilinear_control_operators"]=bilinear.operators
            arrays["bilinear_control_offsets"]=bilinear.offsets

    if policy=="full":
        if values is None:
            evaluator=getattr(result,"evaluate_training_eigenfunctions",None)
            if callable(evaluator): values=np.asarray(evaluator(),dtype=complex)
            else:
                from .virtual_arrays import polynomial_eigenfunction_values
                values=polynomial_eigenfunction_values(coordinates,coefficients,
                    dictionary=result.metadata.get("dictionary") or {},
                    state_preprocessing=result.metadata.get("state_preprocessing") or {})
        arrays["eigenfunction_values"]=values
        storage={"format":"koopman-sindy/full-values/v1","policy":policy,
                 "full_coordinate_series_saved":True,"coefficient_vectors_saved":True,
                 "virtual_arrays":{}}
    else:
        dictionary=result.metadata.get("dictionary") or {}
        if dictionary.get("encoding")!="polynomial-dictionary/v1":
            raise ValueError("coefficient-first storage currently requires a serializable PolynomialDictionary.")
        transform=result.metadata.get("state_preprocessing") or {}
        required={"center","scale","multiplier"}
        if not required.issubset(transform):
            raise ValueError("coefficient-first storage requires complete state-preprocessing metadata.")
        storage={"format":"koopman-sindy/coefficient-first/v1","policy":policy,
                 "full_coordinate_series_saved":False,"coefficient_vectors_saved":True,
                 "virtual_arrays":{"eigenfunction_values":{
                     "encoding":"polynomial-dictionary-evaluation/v1",
                     "states":"training_coordinates","coefficients":"coefficients",
                     "dictionary":dict(dictionary),"state_preprocessing":dict(transform),
                     "shape":[int(expected_shape[0]),int(expected_shape[1])],"dtype":"complex128"}}}

    metadata={**result.metadata,
        "spectrum":result.spectrum.summary(include_eigenvalues=not reused),
        "krom":None if result.krom is None else result.krom.metadata,
        "eigenvalue_storage":{"reused_from_spectrum":reused,
            "source_spectrum_id":None if not reused else result.spectrum.eigenvalue_source.metadata["spectrum_id"],
            "eigenvalues_saved":not reused},
        "storage":storage}
    return arrays,metadata


__all__=["SINDyArchiveStorage","sindy_archive_payload"]
