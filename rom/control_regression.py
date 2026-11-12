"""Controlled-data identification of approximate Koopman input operators."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from scipy.linalg import expm

from infrastructure.controlled_data import ControlledDataset, ControlledTrajectory

from .control_inputs import BilinearControlMatrix, ConstantControlMatrix


@dataclass(frozen=True)
class ControlRegressionResult:
    representation: ConstantControlMatrix | BilinearControlMatrix
    diagnostics: dict[str,Any]

    def attach(self,krom,*,make_default: bool=False):
        name="constant" if isinstance(self.representation,ConstantControlMatrix) else "bilinear"
        return krom.with_control_representation(name,self.representation,make_default=make_default)

    def evaluate(self,krom,controlled_dataset,bounds,*,label: str,attach: bool=True,**options):
        """Create the persistent error record associated with this regression."""
        return krom.evaluate_control_approximation(self.representation,controlled_dataset,bounds,
            label=label,attach=attach,fit_diagnostics=self.diagnostics,**options)


@dataclass(frozen=True)
class SymmetricPulseControlMapSamples:
    """Pointwise control-map estimates recovered from opposite pulse pairs."""

    states: np.ndarray
    matrices: np.ndarray
    diagnostics: dict[str,Any]

    def __post_init__(self):
        states=np.array(self.states,dtype=float,copy=True)
        matrices=np.array(self.matrices,dtype=complex,copy=True)
        if states.ndim!=2 or matrices.ndim!=3 or matrices.shape[0]!=states.shape[0]:
            raise ValueError("Symmetric-pulse state and control-map samples do not align.")
        if states.shape[0]<2 or np.any(~np.isfinite(states)) or np.any(~np.isfinite(matrices)):
            raise ValueError("Symmetric-pulse samples must be finite and nonempty.")
        states.setflags(write=False); matrices.setflags(write=False)
        object.__setattr__(self,"states",states); object.__setattr__(self,"matrices",matrices)
        object.__setattr__(self,"diagnostics",dict(self.diagnostics))


def estimate_symmetric_pulse_control_map(krom,dataset: ControlledDataset) -> SymmetricPulseControlMapSamples:
    """Estimate ``B_psi(x)`` from paired opposite-command first responses.

    Pairing is read from each trajectory's ``pair_id`` and ``polarity``
    metadata.  No analytic gradient or exact KROM input map is consulted.
    """
    groups: dict[int,dict[int,ControlledTrajectory]]={}
    for trajectory in dataset.trajectories:
        metadata=trajectory.metadata or {}
        if "pair_id" not in metadata or "polarity" not in metadata:
            raise ValueError("Symmetric-pulse trajectories require pair_id and polarity metadata.")
        pair_id=int(metadata["pair_id"]); polarity=int(metadata["polarity"])
        if polarity not in {-1,1} or polarity in groups.setdefault(pair_id,{}):
            raise ValueError("Each symmetric-pulse pair needs unique -1 and +1 trajectories.")
        groups[pair_id][polarity]=trajectory
    states=[]; matrices=[]; response_norm=[]
    for pair_id in sorted(groups):
        pair=groups[pair_id]
        if set(pair)!={-1,1}: raise ValueError(f"Symmetric-pulse pair {pair_id} is incomplete.")
        negative,positive=pair[-1],pair[1]
        if not np.array_equal(negative.states[0],positive.states[0]):
            raise ValueError(f"Symmetric-pulse pair {pair_id} does not share an exact anchor.")
        dt=positive.delta_t
        if not np.isclose(negative.delta_t,dt,rtol=0.,atol=0.):
            raise ValueError(f"Symmetric-pulse pair {pair_id} has inconsistent time steps.")
        input_difference=positive.inputs[0]-negative.inputs[0]
        if np.any(input_difference==0.):
            raise ValueError(f"Symmetric-pulse pair {pair_id} has no opposing input difference.")
        if input_difference.size!=krom.input_dimension:
            raise ValueError("Symmetric-pulse inputs do not align with the KROM input dimension.")
        psi_positive=krom.psi(positive.states[1],output="complex",state_coordinates="model")
        psi_negative=krom.psi(negative.states[1],output="complex",state_coordinates="model")
        difference=np.asarray(psi_positive-psi_negative,dtype=complex)
        # The standard paired design excites each input separately.  Division
        # is therefore column-wise and remains explicit for future MIMO data.
        active=np.flatnonzero(input_difference)
        if active.size!=1:
            raise ValueError("A symmetric-pulse pair must excite exactly one input channel.")
        matrix=np.zeros((difference.size,input_difference.size),dtype=complex)
        matrix[:,active[0]]=difference/(float(input_difference[active[0]])*dt)
        states.append(positive.states[0]); matrices.append(matrix); response_norm.append(float(np.linalg.norm(difference)))
    values=np.asarray(response_norm,dtype=float)
    diagnostics={"method":"symmetric-opposite-pulse-control-map-estimator/v1",
        "pair_count":len(states),"analytic_Bpsi_used":False,"sample_time":dataset.delta_t,
        "lifted_response_norm_minimum":float(np.min(values)),"lifted_response_norm_median":float(np.median(values)),
        "lifted_response_norm_maximum":float(np.max(values)),"dataset":dataset.summary()}
    return SymmetricPulseControlMapSamples(np.asarray(states),np.asarray(matrices),diagnostics)


def fit_constant_control_map_samples(krom,samples: SymmetricPulseControlMapSamples) -> ControlRegressionResult:
    """Fit a constant input matrix to statistically estimated map samples."""
    if samples.matrices.shape[1:]!=(krom.lifted_dimension,krom.input_dimension):
        raise ValueError("Control-map samples do not align with the KROM.")
    matrix=np.mean(samples.matrices,axis=0)
    residual=samples.matrices-matrix[None,:,:]
    diagnostics={"method":"sampled-pointwise-constant-control-map-regression/v1",
        "sample_count":int(samples.states.shape[0]),"complex_rmse":float(np.sqrt(np.mean(np.abs(residual)**2))),
        "relative_residual":float(np.linalg.norm(residual)/max(np.linalg.norm(samples.matrices),np.finfo(float).tiny)),
        "source_estimator":samples.diagnostics}
    return ControlRegressionResult(ConstantControlMatrix(matrix,"continuous",None,diagnostics),diagnostics)


def fit_bilinear_control_map_samples(
    krom,
    samples: SymmetricPulseControlMapSamples,
    *,
    ridge: float=0.,
    rcond: float=1e-12,
    conjugate_complete: bool=False,
) -> ControlRegressionResult:
    """Fit ``M Psi+b`` (optionally widely linear) to measured map samples."""
    if samples.matrices.shape[1:]!=(krom.lifted_dimension,krom.input_dimension):
        raise ValueError("Control-map samples do not align with the KROM.")
    if ridge<0. or not np.isfinite(ridge) or rcond<=0. or not np.isfinite(rcond):
        raise ValueError("ridge and rcond settings are invalid.")
    lifted=_complex_lift(krom,samples.states)
    parts=[lifted]
    if conjugate_complete: parts.append(np.conjugate(lifted))
    design=np.column_stack((*parts,np.ones(lifted.shape[0],dtype=complex)))
    penalty=np.eye(design.shape[1],dtype=float); penalty[-1,-1]=0.
    if ridge==0.: inverse=np.linalg.pinv(design,rcond=rcond); solve_method="svd-pseudoinverse"
    else:
        inverse=np.linalg.pinv(design.conj().T@design+ridge*penalty,rcond=rcond)@design.conj().T
        solve_method="ridge-normal-equations"
    modes=krom.lifted_dimension; inputs=krom.input_dimension
    operators=np.empty((inputs,modes,modes),dtype=complex)
    conjugate=np.empty_like(operators) if conjugate_complete else None
    offsets=np.empty((modes,inputs),dtype=complex); predicted=np.empty_like(samples.matrices)
    for input_index in range(inputs):
        coefficients=inverse@samples.matrices[:,:,input_index]
        operators[input_index]=coefficients[:modes].T
        cursor=modes
        if conjugate_complete:
            conjugate[input_index]=coefficients[cursor:cursor+modes].T; cursor+=modes
        offsets[:,input_index]=coefficients[cursor]
        predicted[:,:,input_index]=design@coefficients
    residual=predicted-samples.matrices
    diagnostics={"method":"sampled-pointwise-bilinear-control-map-regression/v1",
        "conjugate_complete":bool(conjugate_complete),"sample_count":int(samples.states.shape[0]),
        "design_shape":list(design.shape),"design_rank":int(np.linalg.matrix_rank(design)),
        "design_condition_number":float(np.linalg.cond(design)),"solve_method":solve_method,
        "ridge":float(ridge),"rcond":float(rcond),"complex_rmse":float(np.sqrt(np.mean(np.abs(residual)**2))),
        "relative_residual":float(np.linalg.norm(residual)/max(np.linalg.norm(samples.matrices),np.finfo(float).tiny)),
        "source_estimator":samples.diagnostics}
    representation=BilinearControlMatrix(operators,offsets,"continuous",None,diagnostics,conjugate)
    return ControlRegressionResult(representation,diagnostics)


def _complex_lift(krom, states: np.ndarray) -> np.ndarray:
    return np.vstack([krom.psi(row,output="complex",state_coordinates="model") for row in states])


def _split_to_complex_rows(matrix: np.ndarray, eigenfunctions) -> np.ndarray:
    values=np.asarray(matrix,dtype=float); rows=[]; offset=0
    for function in eigenfunctions:
        if function.is_real:
            rows.append(values[offset].astype(complex)); offset+=1
        else:
            rows.append(values[offset]+1j*values[offset+1]); offset+=2
    if offset!=values.shape[0]: raise ValueError("Split matrix does not align with KROM eigenfunctions.")
    return np.vstack(rows)


def _discrete_pair(krom,delta_t: float) -> tuple[np.ndarray,np.ndarray]:
    generator=krom.A(output="split",time_coordinates="model")
    dimension=generator.shape[0]
    augmented=np.zeros((2*dimension,2*dimension),dtype=float)
    augmented[:dimension,:dimension]=generator
    augmented[:dimension,dimension:]=np.eye(dimension)
    transition=expm(augmented*delta_t)
    return transition[:dimension,:dimension],transition[:dimension,dimension:]


def fit_constant_control_matrix(
    krom,
    dataset: ControlledDataset,
    *,
    reconstruction=None,
    method: Literal["multistep-output","continuous-residual"]="multistep-output",
    ridge: float=0.0,
) -> ControlRegressionResult:
    """Fit a constant lifted ``B`` without replacing the exact KROM input map.

    ``multistep-output`` implements the Korda--Mezic Duffing code path: fixed
    discrete ``A`` and observable reconstruction ``C`` are held constant while
    all multi-step output errors are regressed against interval-held controls.
    The raw discrete fit is retained and converted consistently for continuous
    KLQR use. ``continuous-residual`` is a local derivative-domain alternative.
    """
    if ridge<0. or not np.isfinite(ridge): raise ValueError("ridge must be finite and nonnegative.")
    if dataset.input_dimension!=krom.input_dimension:
        raise ValueError("Controlled dataset input dimension does not match the KROM system.")
    if method=="multistep-output":
        if reconstruction is None: reconstruction=krom.observable_reconstruction
        if reconstruction is None:
            raise ValueError("multistep-output constant-B regression requires an observable reconstruction C.")
        if reconstruction.lifted_dimension!=krom.lifted_dimension:
            raise ValueError("Reconstruction C is not aligned with the KROM.")
        dt=dataset.delta_t; transition,integral=_discrete_pair(krom,dt)
        C=reconstruction.split_matrix(krom.eigenfunctions); intercept=reconstruction.intercept
        lifted_dimension=transition.shape[0]; inputs=dataset.input_dimension
        blocks=[]; targets=[]
        for trajectory in dataset.trajectories:
            observed=reconstruction.observable.evaluate(trajectory.states)
            # Only the initial lifted state enters the autonomous multi-step
            # prediction.  Avoid evaluating Psi on every controlled sample.
            free_state=np.asarray(krom.psi(trajectory.states[0],output="split"),dtype=float)
            control_response=np.zeros((C.shape[0],inputs,lifted_dimension),dtype=float)
            for step in range(1,trajectory.states.shape[0]):
                # The original Korda--Mezic design is a control convolution,
                #   S_k = sum_j kron(u_j, C A^(k-1-j)).
                # Advancing that convolution recursively is algebraically
                # identical and changes assembly from O(horizon^2) to
                # O(horizon), which is material for the retained 20k-sample
                # controlled datasets used by campaign searches.
                control_response=np.einsum(
                    "oim,mn->oin",control_response,transition,optimize=True
                )+C[:,None,:]*trajectory.inputs[step-1][None,:,None]
                free_state=transition@free_state
                blocks.append(control_response.reshape(C.shape[0],inputs*lifted_dimension))
                targets.append(observed[step]-intercept-C@free_state)
        design=np.vstack(blocks); target=np.concatenate(targets)
        gram=design.T@design+ridge*np.eye(design.shape[1])
        vector=np.linalg.pinv(gram,rcond=1e-12)@design.T@target
        discrete_split=vector.reshape((lifted_dimension,inputs),order="F")
        continuous_split=np.linalg.lstsq(integral,discrete_split,rcond=1e-12)[0]
        predicted=design@vector; residual=target-predicted
        complex_matrix=_split_to_complex_rows(continuous_split,krom.eigenfunctions)
        discrete_complex=_split_to_complex_rows(discrete_split,krom.eigenfunctions)
        diagnostics={"method":"korda-mezic-multistep-output/v1","sample_time":dt,
            "design_assembly":"recursive-control-convolution/O(horizon)",
            "design_shape":list(design.shape),"design_rank":int(np.linalg.matrix_rank(design)),
            "condition_number":float(np.linalg.cond(gram)),"ridge":float(ridge),
            "output_rmse":float(np.sqrt(np.mean(residual*residual))),
            "maximum_absolute_output_residual":float(np.max(np.abs(residual))),
            "discrete_to_continuous":"Bd = integral_0^dt exp(A s) ds Bc; least-squares solve",
            "dataset":dataset.summary(),"reconstruction":reconstruction.summary()}
        representation=ConstantControlMatrix(complex_matrix,"continuous",dt,
            {**diagnostics,"raw_fit_time_domain":"discrete"},discrete_complex)
        return ControlRegressionResult(representation,diagnostics)
    if method!="continuous-residual":
        raise ValueError("method must be 'multistep-output' or 'continuous-residual'.")
    lift=[]; inputs=[]; derivatives=[]
    for trajectory in dataset.trajectories:
        psi=_complex_lift(krom,trajectory.states)
        dpsi=np.gradient(psi,trajectory.time,axis=0,edge_order=2)
        lift.append(psi[:-1]); inputs.append(trajectory.inputs); derivatives.append(dpsi[:-1])
    psi=np.vstack(lift); U=np.vstack(inputs); dpsi=np.vstack(derivatives)
    residual=dpsi-psi@krom.A(output="complex").T
    coefficient=np.linalg.pinv(U.T@U+ridge*np.eye(U.shape[1]),rcond=1e-12)@U.T@residual
    matrix=coefficient.T
    error=residual-U@coefficient
    diagnostics={"method":"continuous-lifted-residual/v1","design_shape":list(U.shape),
        "design_rank":int(np.linalg.matrix_rank(U)),"ridge":float(ridge),
        "lifted_rmse":float(np.sqrt(np.mean(np.abs(error)**2))),"dataset":dataset.summary()}
    return ControlRegressionResult(ConstantControlMatrix(matrix,"continuous",dataset.delta_t,diagnostics),diagnostics)


def fit_bilinear_control_matrix(
    krom,
    dataset: ControlledDataset,
    *,
    ridge: float=0.0,
) -> ControlRegressionResult:
    """Fit ``dPsi-A Psi = sum_j (M_j Psi+b_j) u_j`` in continuous time."""
    if ridge<0. or not np.isfinite(ridge): raise ValueError("ridge must be finite and nonnegative.")
    rows=[]; residuals=[]
    for trajectory in dataset.trajectories:
        psi=_complex_lift(krom,trajectory.states)
        derivative=np.gradient(psi,trajectory.time,axis=0,edge_order=2)
        active=psi[:-1]; controls=trajectory.inputs
        features=np.hstack([
            np.hstack((controls[:,index,None]*active,controls[:,index,None]))
            for index in range(dataset.input_dimension)
        ])
        rows.append(features); residuals.append(derivative[:-1]-active@krom.A(output="complex").T)
    design=np.vstack(rows); target=np.vstack(residuals)
    gram=design.conj().T@design+ridge*np.eye(design.shape[1])
    coefficients=np.linalg.pinv(gram,rcond=1e-12)@design.conj().T@target
    modes=krom.lifted_dimension; operators=[]; offsets=np.empty((modes,dataset.input_dimension),dtype=complex)
    for index in range(dataset.input_dimension):
        block=coefficients[index*(modes+1):(index+1)*(modes+1)]
        operators.append(block[:modes].T); offsets[:,index]=block[modes]
    error=target-design@coefficients
    diagnostics={"method":"continuous-bilinear-lifted-residual/v1","design_shape":list(design.shape),
        "design_rank":int(np.linalg.matrix_rank(design)),"condition_number":float(np.linalg.cond(gram)),
        "ridge":float(ridge),"lifted_rmse":float(np.sqrt(np.mean(np.abs(error)**2))),
        "dataset":dataset.summary()}
    representation=BilinearControlMatrix(np.asarray(operators),offsets,"continuous",dataset.delta_t,diagnostics)
    return ControlRegressionResult(representation,diagnostics)


def fit_bilinear_control_map(
    krom,
    states: np.ndarray,
    *,
    ridge: float=0.0,
    rcond: float=1e-12,
    gradient_backend: str="auto",
) -> ControlRegressionResult:
    """Regress ``B_psi(x) = B1 Psi(x) + B2`` from exact map samples.

    Unlike the controlled-dynamics residual fit, this projection needs state
    samples and an executable exact ``B_psi`` backend but no forced trajectory.
    It therefore applies uniformly to analytic, H-DMD-interpolated, and
    fitted-boundary KROMs.
    """
    values=np.asarray(states,dtype=float)
    state_dimension=int(np.asarray(krom.system.initial_state).size)
    if values.ndim!=2 or values.shape[1]!=state_dimension or values.shape[0]<3 or np.any(~np.isfinite(values)):
        raise ValueError("states must be a finite sample-major state matrix with at least three rows.")
    if ridge<0. or not np.isfinite(ridge) or rcond<=0. or not np.isfinite(rcond):
        raise ValueError("ridge and rcond settings are invalid.")
    lifted=_complex_lift(krom,values)
    exact=np.asarray([krom.B_psi(row,output="complex",representation="exact",
        gradient_backend=gradient_backend) for row in values],dtype=complex)
    design=np.column_stack((lifted,np.ones(lifted.shape[0],dtype=complex)))
    modes=krom.lifted_dimension; inputs=krom.input_dimension
    penalty=np.diag(np.r_[np.ones(modes),0.])
    if ridge==0.:
        inverse=np.linalg.pinv(design,rcond=rcond); solve_method="svd-pseudoinverse"
    else:
        inverse=np.linalg.pinv(design.conj().T@design+ridge*penalty,rcond=rcond)@design.conj().T
        solve_method="ridge-normal-equations"
    operators=np.empty((inputs,modes,modes),dtype=complex); offsets=np.empty((modes,inputs),dtype=complex)
    predicted=np.empty_like(exact)
    for input_index in range(inputs):
        coefficients=inverse@exact[:,:,input_index]
        operators[input_index]=coefficients[:modes].T; offsets[:,input_index]=coefficients[modes]
        predicted[:,:,input_index]=design@coefficients
    residual=predicted-exact
    diagnostics={"method":"pointwise-bilinear-exact-control-map-regression/v1",
        "purpose":"isolate B1*Psi+B2 representation error from dynamical derivative error",
        "gradient_backend":gradient_backend,"eigenfunction_representation_kinds":
        [item.representation_kind for item in krom.eigenfunctions],"solve_method":solve_method,
        "design_shape":list(design.shape),"design_rank":int(np.linalg.matrix_rank(design)),
        "condition_number":float(np.linalg.cond(design)),"ridge":float(ridge),"rcond":float(rcond),
        "lifted_rmse":float(np.sqrt(np.mean(np.abs(residual)**2))),
        "relative_residual":float(np.linalg.norm(residual)/max(np.linalg.norm(exact),np.finfo(float).tiny)),
        "sample_count":int(values.shape[0])}
    representation=BilinearControlMatrix(operators,offsets,"continuous",krom.delta_t,diagnostics)
    return ControlRegressionResult(representation,diagnostics)


def fit_conjugate_bilinear_control_map(
    krom,
    states: np.ndarray,
    *,
    ridge: float=0.0,
    rcond: float=1e-12,
    gradient_backend: str="auto",
) -> ControlRegressionResult:
    """Fit ``B_psi=M psi+N conjugate(psi)+b`` from exact map samples.

    This candidate is the conjugate-complete equivalent of fitting a general
    real-linear bilinear map on the KROM's real/imaginary split lift.  It is
    additive and opt-in; the existing complex-linear bilinear fit is unchanged.
    """
    values=np.asarray(states,dtype=float)
    state_dimension=int(np.asarray(krom.system.initial_state).size)
    if values.ndim!=2 or values.shape[1]!=state_dimension or values.shape[0]<3 or np.any(~np.isfinite(values)):
        raise ValueError("states must be a finite sample-major state matrix with at least three rows.")
    if ridge<0. or not np.isfinite(ridge) or rcond<=0. or not np.isfinite(rcond):
        raise ValueError("ridge and rcond settings are invalid.")
    lifted=_complex_lift(krom,values)
    exact=np.asarray([krom.B_psi(row,output="complex",representation="exact",
        gradient_backend=gradient_backend) for row in values],dtype=complex)
    design=np.column_stack((lifted,np.conjugate(lifted),np.ones(lifted.shape[0],dtype=complex)))
    modes=krom.lifted_dimension; inputs=krom.input_dimension
    penalty=np.diag(np.r_[np.ones(2*modes),0.])
    if ridge==0.:
        inverse=np.linalg.pinv(design,rcond=rcond); solve_method="svd-pseudoinverse"
    else:
        inverse=np.linalg.pinv(design.conj().T@design+ridge*penalty,rcond=rcond)@design.conj().T
        solve_method="ridge-normal-equations"
    operators=np.empty((inputs,modes,modes),dtype=complex)
    conjugate=np.empty_like(operators)
    offsets=np.empty((modes,inputs),dtype=complex); predicted=np.empty_like(exact)
    for input_index in range(inputs):
        coefficients=inverse@exact[:,:,input_index]
        operators[input_index]=coefficients[:modes].T
        conjugate[input_index]=coefficients[modes:2*modes].T
        offsets[:,input_index]=coefficients[2*modes]
        predicted[:,:,input_index]=design@coefficients
    residual=predicted-exact
    diagnostics={"method":"pointwise-conjugate-complete-bilinear-control-map-regression/v1",
        "purpose":"general real-linear map on split coordinates; equivalent to explicit conjugate dual features",
        "candidate_only":True,"gradient_backend":gradient_backend,"solve_method":solve_method,
        "design_shape":list(design.shape),"design_rank":int(np.linalg.matrix_rank(design)),
        "condition_number":float(np.linalg.cond(design)),"ridge":float(ridge),"rcond":float(rcond),
        "lifted_rmse":float(np.sqrt(np.mean(np.abs(residual)**2))),
        "relative_residual":float(np.linalg.norm(residual)/max(np.linalg.norm(exact),np.finfo(float).tiny)),
        "sample_count":int(values.shape[0])}
    representation=BilinearControlMatrix(operators,offsets,"continuous",krom.delta_t,
        diagnostics,conjugate)
    return ControlRegressionResult(representation,diagnostics)


__all__=["ControlledDataset","ControlledTrajectory","ControlRegressionResult",
         "fit_bilinear_control_map","fit_bilinear_control_matrix",
         "fit_conjugate_bilinear_control_map","fit_constant_control_matrix"]
