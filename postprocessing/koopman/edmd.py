"""Extended DMD with user-supplied observable dictionaries."""
from __future__ import annotations
import numpy as np
from time_series import TimeSeries
from .dictionary import ObservableDictionary
from .result import KoopmanResult
from .validation import require_uniform_time_steps
from .spectrum import discrete_spectrum

def edmd(states: np.ndarray | TimeSeries, dictionary: ObservableDictionary, regularization: float=0., *, time: np.ndarray|None=None, metadata: dict|None=None, source_archive: str|None=None) -> KoopmanResult:
    if isinstance(states,TimeSeries):
        if time is not None: raise ValueError("Pass either a TimeSeries or separate states and time, not both.")
        if states.is_pooled: raise ValueError("This EDMD call received a pooled TimeSeries. Use series.members() to preserve trajectory boundaries.")
        time=states.time; states=states.values
    if time is None: raise ValueError("Array EDMD input requires an explicit time array.")
    x=np.asarray(states,dtype=float)
    if x.ndim==1: x=x[:,None]
    validation=require_uniform_time_steps(time,metadata)
    if x.shape[0]!=np.asarray(time).size: raise ValueError("EDMD states and time arrays must have equal sample counts.")
    if x.ndim!=2 or x.shape[0]<3: raise ValueError("EDMD requires at least three time-ordered state samples.")
    psi_x,psi_y=dictionary.transform(x[:-1]),dictionary.transform(x[1:])
    gram=psi_x.T@psi_x; cross=psi_x.T@psi_y
    K=np.linalg.solve(gram+regularization*np.eye(gram.shape[0]),cross)
    eigvals,eigvecs=np.linalg.eig(K); eigenfunctions=dictionary.transform(x)@eigvecs
    modes=np.linalg.lstsq(eigenfunctions,x,rcond=None)[0].T
    result=KoopmanResult("edmd",eigvals,modes,eigenfunctions,K,{"observables":psi_x.shape[1],"regularization":regularization})
    result.metadata.update(validation)
    result.metadata["preprocessing"]={"state_dimension":int(x.shape[1]),"sample_count":int(x.shape[0]),"analysis_delta_t":validation.get("delta_t"),"dictionary_class":type(dictionary).__name__}
    result.metadata["delay_embedding"]={"type":"none","delays":1,"window_samples":1,"window_duration":0.}
    result.spectrum=discrete_spectrum(result.eigenvalues,result.eigenfunction_values,result.modes,time=time,source_metadata=metadata,source_archive=source_archive,
                                      method="edmd",representation="EDMD dictionary eigenvector coefficients and trajectory values",coefficients=eigvecs,method_metadata={**result.metadata,"dictionary_class":type(dictionary).__name__})
    for function,coefficient in zip(result.spectrum.eigenfunctions,eigvecs.T):
        function.state_callable=lambda state, coeff=coefficient, basis=dictionary: complex(np.asarray(basis.transform(np.asarray(state,dtype=float).reshape(1,-1)))[0]@coeff)
        function.gradient_callable=lambda state, coeff=coefficient, basis=dictionary: np.asarray(coeff)@np.asarray(basis.gradient(np.asarray(state,dtype=float).reshape(1,-1)))[0]
    from rom import KROM
    result.krom=KROM(result.spectrum)
    result.metadata["krom"]=result.krom.metadata
    result.metadata["spectrum"]=result.spectrum.summary(); return result
