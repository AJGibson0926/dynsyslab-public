"""Exact and fitted control representations plus attached error evidence."""
from __future__ import annotations

from copy import deepcopy
from typing import Any
import numpy as np

from .control_inputs import (BilinearControlMatrix, ConstantControlMatrix,
                             complex_to_split)

class _KROMControlMixin:
    """Exact and fitted control representations plus attached error evidence."""

    def _validate_control_representation(self) -> None:
        if self.control_mode not in {"analytic","bilinear","linear"}:
            raise ValueError("control_mode must be 'analytic', 'bilinear', or 'linear'.")
        # Preserve the original one-active-representation constructor while
        # normalizing it into a simultaneous named collection.
        modes=self.lifted_dimension
        representations=dict(self.control_representations or {})
        inputs=None
        try: inputs=self.input_dimension
        except (AttributeError,NotImplementedError,TypeError,ValueError): pass
        if self.control_mode=="linear":
            if self.B1 is None:
                if "constant" not in representations:
                    raise ValueError("control_mode='linear' requires B1 or a named constant representation.")
            else:
                matrix=np.asarray(self.B1,dtype=complex)
                if matrix.ndim!=2 or matrix.shape[0]!=modes:
                    raise ValueError("Linear B1 must have shape (retained_modes, input_dimension).")
                if inputs is not None and matrix.shape[1]!=inputs:
                    raise ValueError("Linear B1 input dimension does not match the system.")
                if self.B2 is not None: raise ValueError("Linear KROM control uses B1 only; do not provide B2.")
                representations.setdefault("constant",ConstantControlMatrix(matrix,provenance={"source":"legacy-B1"}))
        elif self.control_mode=="bilinear":
            if self.B1 is None or self.B2 is None:
                if "bilinear" not in representations:
                    raise ValueError("control_mode='bilinear' requires B1/B2 or a named bilinear representation.")
            else:
                operator=np.asarray(self.B1,dtype=complex); offsets=np.asarray(self.B2,dtype=complex)
                if operator.shape!=(modes,modes):
                    raise ValueError("Legacy bilinear B1 must have shape (retained_modes, retained_modes).")
                if offsets.shape==(modes,): offsets=offsets[:,None]
                if offsets.ndim!=2 or offsets.shape[0]!=modes:
                    raise ValueError("Legacy bilinear B2 must have shape (retained_modes, input_dimension).")
                if inputs is not None and offsets.shape[1]!=inputs:
                    raise ValueError("Bilinear B2 input dimension does not match the system.")
                if offsets.shape[1]!=1:
                    raise ValueError("Legacy B1/B2 syntax supports one input; use BilinearControlMatrix for multiple inputs.")
                representations.setdefault("bilinear",BilinearControlMatrix(operator[None,:,:],offsets,provenance={"source":"legacy-B1-B2"}))
        for name,value in representations.items():
            if name=="constant" and not isinstance(value,ConstantControlMatrix):
                raise TypeError("The 'constant' representation must be ConstantControlMatrix.")
            if name=="bilinear" and not isinstance(value,BilinearControlMatrix):
                raise TypeError("The 'bilinear' representation must be BilinearControlMatrix.")
            if name not in {"constant","bilinear"}:
                raise ValueError(f"Unknown stored control representation {name!r}.")
            if value.lifted_dimension!=modes:
                raise ValueError(f"Stored {name} representation does not match KROM lifted dimension.")
            if inputs is not None and value.input_dimension!=inputs:
                raise ValueError(f"Stored {name} representation does not match system input dimension.")
            if value.time_domain!="continuous":
                raise ValueError("Executable KROM control representations must be continuous-time; retain raw discrete fits in provenance/discrete_matrix.")
        self.control_representations=representations
        legacy_default={"analytic":"exact","linear":"constant","bilinear":"bilinear"}[self.control_mode]
        # A caller may still mutate the legacy public fields and then request
        # validation.  In that compatibility path, a non-analytic control_mode
        # is an explicit request to activate the corresponding fitted model.
        if self.default_control_representation is None or (
                self.control_mode!="analytic" and self.default_control_representation=="exact"):
            self.default_control_representation=legacy_default
        if (self.default_control_representation!="exact" and
                self.default_control_representation not in self.available_control_representations()):
            raise ValueError("Default control representation is not available on this KROM.")
        # Keep legacy fields synchronized for old archive/tests without making
        # them the authoritative storage boundary.
        if self.default_control_representation=="constant":
            self.control_mode="linear"; self.B1=np.array(representations["constant"].matrix,copy=True); self.B2=None
        elif self.default_control_representation=="bilinear":
            self.control_mode="bilinear"
            active=representations["bilinear"]
            if active.input_dimension==1:
                self.B1=np.array(active.operators[0],copy=True); self.B2=np.array(active.offsets,copy=True)
            else:
                # Legacy B1 represents one input only; named representations
                # remain authoritative for multi-input systems.
                self.B1=None; self.B2=None
        else:
            self.control_mode="analytic"; self.B1=None; self.B2=None

    def available_control_representations(self) -> list[str]:
        result=[]
        try:
            state=self._model_state(np.asarray(self.system.initial_state,dtype=float),"physical")
            # Exact coupling is conceptually available whenever the physical
            # input map exists.  Analytic gradients are a stricter capability;
            # B_psi emits the established NotImplementedError unless the caller
            # explicitly opts into finite differences.
            exact=self._input_matrix(state,"model") is not None
        except (AttributeError,NotImplementedError,TypeError,ValueError): exact=False
        if exact: result.append("exact")
        result.extend(name for name in ("bilinear","constant") if name in (self.control_representations or {}))
        return result

    def with_control_representation(self, name: str, representation, *, make_default: bool=False) -> "KROM":
        """Return a derivative retaining all existing input models plus one named fit."""
        if name not in {"constant","bilinear"}:
            raise ValueError("Stored representation name must be 'constant' or 'bilinear'.")
        derived=self.extract(range(len(self.eigenvalues)))
        derived.control_error_evaluations=deepcopy(self.control_error_evaluations)
        derived.control_representations={**derived.control_representations,name:representation}
        if make_default: derived.default_control_representation=name
        derived._validate_control_representation()
        derived.spectrum.metadata["control_representations"]={key:value.summary() for key,value in derived.control_representations.items()}
        return derived

    def evaluate_eigenfunction_errors(self, dataset, *, horizons=None, metadata=None,
                                      attach: bool=True) -> dict[int,Any]:
        """Evaluate every retained eigenfunction on one free training/test set."""
        from spectral.errors import FreeTrajectoryDataset
        from time_series import TimeSeries
        if not isinstance(dataset,FreeTrajectoryDataset):
            raise TypeError("dataset must be a FreeTrajectoryDataset.")
        normalized=dataset
        if dataset.coordinate_system=="physical" and self.coordinate_system=="dimensionless":
            trajectories=[]
            for item in dataset.trajectories:
                time=self._model_time(item.time,"physical")
                states=np.asarray([self._model_state(row,"physical") for row in item.values])
                trajectories.append(TimeSeries(time,states))
            normalized=FreeTrajectoryDataset(trajectories,label=dataset.label,role=dataset.role,
                coordinate_system="model",provenance={**dataset.provenance,
                "physical_dataset_fingerprint_sha256":dataset.fingerprint,
                "coordinate_transform":"physical-to-krom-model"})
        elif dataset.coordinate_system!="model" and self.coordinate_system=="dimensional":
            normalized=FreeTrajectoryDataset(dataset.trajectories,label=dataset.label,role=dataset.role,
                coordinate_system="model",provenance={**dataset.provenance,
                "physical_dataset_fingerprint_sha256":dataset.fingerprint,
                "coordinate_transform":"identity-physical-to-model"})
        coordinate_map=lambda state:np.asarray(self.coordinate_observable.evaluate(state),dtype=float)
        result={}
        for index,function in enumerate(self.eigenfunctions):
            result[index]=function.evaluate_errors(normalized,horizons=horizons,
                coordinate_map=coordinate_map,training_domain=self.training_domain,
                metadata={"krom_method":self.method,"krom_retained_index":index,
                **dict(metadata or {})},attach=attach)
        return result

    def attach_control_error_evaluation(self, evaluation, *, replace: bool=False):
        """Attach a control-map fit and its phase-space error field by unique id."""
        from .control_error import ControlApproximationEvaluation
        if not isinstance(evaluation,ControlApproximationEvaluation):
            raise TypeError("evaluation must be a ControlApproximationEvaluation.")
        if evaluation.representation.lifted_dimension!=self.lifted_dimension:
            raise ValueError("Control-error evaluation does not match KROM lifted dimension.")
        if evaluation.representation.input_dimension!=self.input_dimension:
            raise ValueError("Control-error evaluation does not match KROM input dimension.")
        if evaluation.evaluation_id in self.control_error_evaluations and not replace:
            raise KeyError(f"Control error evaluation {evaluation.evaluation_id!r} already exists.")
        self.control_error_evaluations[evaluation.evaluation_id]=evaluation
        return evaluation

    def evaluate_control_approximation(self, representation, controlled_dataset, bounds, *,
                                       label: str, attach: bool=True, **options):
        """Evaluate and optionally attach a candidate B_psi approximation."""
        from .control_error import evaluate_control_approximation
        evaluation=evaluate_control_approximation(self,representation,controlled_dataset,bounds,
            label=label,**options)
        if attach: self.attach_control_error_evaluation(evaluation)
        return evaluation

    def compare_control_error_evaluations(self, evaluation_ids=None, *, metric: str="relative.mean"):
        """Compare retained regressions without changing the active control map."""
        from .control_error import compare_control_approximation_evaluations
        selected=list(self.control_error_evaluations.values()) if evaluation_ids is None else [
            self.control_error_evaluations[str(item)] for item in evaluation_ids]
        return compare_control_approximation_evaluations(selected,metric=metric)

    def B_psi(self,state: np.ndarray, *, output: str="split", state_coordinates: str="model", time_coordinates: str="model",
              gradient_backend: str="auto", representation: str|None=None) -> np.ndarray:
        x=self._model_state(state,state_coordinates)
        selected=self.default_control_representation if representation in {None,"default"} else representation
        if selected not in self.available_control_representations():
            raise ValueError(f"Control representation {selected!r} is unavailable; available={self.available_control_representations()}.")
        if selected=="exact":
            if gradient_backend not in {"auto","analytic","finite_difference"}: raise ValueError("gradient_backend must be 'auto', 'analytic', or 'finite_difference'.")
            original=self._input_matrix(x,time_coordinates)
            if original is None: raise ValueError("This KROM's system has no compatible input matrix.")
            original=np.asarray(original,dtype=float)
            if original.ndim!=2 or original.shape[0]!=x.size: raise ValueError("System B(x) must have shape (state_dimension, input_dimension).")
            complex_result=self.nabla_psi(x,output="complex",state_coordinates="model",
                gradient_backend=gradient_backend)@original
        elif selected=="constant": complex_result=self.control_representations["constant"].evaluate()
        else:
            bilinear=self.control_representations["bilinear"]
            complex_result=bilinear.evaluate(self.psi(x,output="complex"))
        if output=="complex": return complex_result
        if output=="split": return complex_to_split(complex_result,self.eigenfunctions)
        raise ValueError("output must be 'complex' or 'split'.")
