"""Self-verifying storage for KLQE designs and result time series."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np

from time_series import TimeSeries


def _sha256_file(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()

def save_klqe_result(path: str | Path, result: KLQEResult) -> Path:
    target = Path(path).expanduser().resolve(); target.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, Any] = {
        "time": result.lifted_estimate.time,
        "lifted_estimate": result.lifted_estimate.values,
        "state_estimate": result.state_estimate.values,
        "measurement": result.measurement.values,
        "measurement_estimate": result.measurement_estimate.values,
        "innovation": result.innovation.values,
        "metadata_json": np.asarray(json.dumps(result.summary(), sort_keys=True)),
    }
    if result.true_state is not None: arrays["true_state"] = result.true_state.values
    if result.state_error is not None: arrays["state_error"] = result.state_error.values
    np.savez_compressed(target, **arrays); return target

def save_klqe_design(
    path: str | Path,
    design: KLQEDesign,
    *,
    krom_archive: str | Path,
    measurement_name: str,
    state_output_name: str,
) -> Path:
    """Publish a self-verifying numerical design tied to one KROM archive."""
    source=Path(krom_archive).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    fingerprint=_sha256_file(source)
    metadata={
        "encoding":"koopman-klqe-design-artifact/v1",
        "krom_archive":str(source),"krom_sha256":fingerprint,
        "measurement_name":str(measurement_name),"state_output_name":str(state_output_name),
        "design":design.summary(),
    }
    target=Path(path).expanduser().resolve(); target.parent.mkdir(parents=True,exist_ok=True)
    model=design.observation_model
    np.savez_compressed(target,A=model.A,C_measurement=model.C_measurement,
        measurement_intercept=model.measurement_intercept,C_state=model.C_state,
        state_intercept=model.state_intercept,Q=design.Q,R=design.R,covariance=design.covariance,
        gain=design.gain,observer_poles=design.observer_poles,
        metadata_json=np.asarray(json.dumps(metadata,sort_keys=True)))
    return target

def load_klqe_design(path: str | Path) -> tuple[KLQEDesign,dict[str,Any]]:
    """Restore a design only after its KROM fingerprint and matrices agree."""
    from .koopman import design_klqe_from_krom
    source=Path(path).expanduser().resolve()
    with np.load(source,allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        stored={name:np.asarray(data[name]) for name in (
            "A","C_measurement","measurement_intercept","C_state","state_intercept",
            "Q","R","covariance","gain","observer_poles")}
    krom_path=Path(metadata["krom_archive"]).resolve()
    if not krom_path.is_file(): raise FileNotFoundError(krom_path)
    actual=_sha256_file(krom_path)
    if actual!=metadata["krom_sha256"]:
        raise ValueError("The KLQE design's KROM archive fingerprint no longer matches.")
    from rom import load_krom_archive
    krom=load_krom_archive(krom_path)
    design=design_klqe_from_krom(krom,metadata["measurement_name"],metadata["state_output_name"],
                                 stored["Q"],stored["R"])
    current={"A":design.observation_model.A,"C_measurement":design.observation_model.C_measurement,
        "measurement_intercept":design.observation_model.measurement_intercept,
        "C_state":design.observation_model.C_state,"state_intercept":design.observation_model.state_intercept,
        "covariance":design.covariance,"gain":design.gain,"observer_poles":design.observer_poles}
    for name,value in current.items():
        if not np.allclose(value,stored[name],rtol=1e-10,atol=1e-12):
            raise ValueError(f"Restored KLQE design disagrees with stored {name}.")
    return design,metadata

def load_klqe_result(path: str | Path) -> KLQEResult:
    from .koopman import KLQEResult
    source = Path(path).expanduser().resolve()
    with np.load(source, allow_pickle=False) as data:
        summary = json.loads(str(data["metadata_json"])); time = np.asarray(data["time"], dtype=float)
        truth = None if "true_state" not in data.files else TimeSeries(time, data["true_state"])
        error = None if "state_error" not in data.files else TimeSeries(time, data["state_error"])
        return KLQEResult(
            TimeSeries(time, data["lifted_estimate"]), TimeSeries(time, data["state_estimate"]),
            TimeSeries(time, data["measurement"]), TimeSeries(time, data["measurement_estimate"]),
            TimeSeries(time, data["innovation"]), summary.get("metadata") or {}, truth, error,
        )
