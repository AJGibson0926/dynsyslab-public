"""Canonical persistence for aligned controlled-trajectory ensembles."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .control_regression import ControlledDataset, ControlledTrajectory


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def save_controlled_dataset(dataset: ControlledDataset,destination: str|Path) -> tuple[Path,Path]:
    """Save a shared-grid controlled ensemble and its audit manifest."""
    path=Path(destination).expanduser().resolve()
    if path.suffix.lower()!=".npz": path=path.with_suffix(".npz")
    path.parent.mkdir(parents=True,exist_ok=True)
    reference=dataset.trajectories[0].time
    if any(not np.array_equal(item.time,reference) for item in dataset.trajectories[1:]):
        raise ValueError("Canonical controlled-dataset persistence requires one shared time grid.")
    arrays={
        "time":reference,
        "states":np.stack([item.states for item in dataset.trajectories]),
        "inputs":np.stack([item.inputs for item in dataset.trajectories]),
        "metadata_json":np.asarray(json.dumps(dict(dataset.metadata or {}),sort_keys=True)),
        "trajectory_metadata_json":np.asarray(json.dumps([dict(item.metadata or {}) for item in dataset.trajectories],sort_keys=True)),
    }
    np.savez_compressed(path,**arrays)
    manifest=path.with_suffix(".json")
    payload={**dataset.summary(),"archive":str(path),"archive_sha256":_sha256(path),
        "array_keys":sorted(arrays),"array_shapes":{name:list(value.shape) for name,value in arrays.items()}}
    manifest.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return path,manifest


def load_controlled_dataset(source: str|Path) -> ControlledDataset:
    """Load and validate a canonical controlled ensemble."""
    path=Path(source).expanduser().resolve()
    with np.load(path,allow_pickle=False) as archive:
        required={"time","states","inputs","metadata_json","trajectory_metadata_json"}
        if not required.issubset(archive.files):
            raise ValueError(f"Controlled dataset is missing arrays: {sorted(required-set(archive.files))}")
        time=np.asarray(archive["time"],dtype=float)
        states=np.asarray(archive["states"],dtype=float)
        inputs=np.asarray(archive["inputs"],dtype=float)
        metadata=json.loads(str(archive["metadata_json"].item()))
        trajectory_metadata=json.loads(str(archive["trajectory_metadata_json"].item()))
    if states.ndim!=3 or inputs.ndim!=3 or states.shape[0]!=inputs.shape[0] or len(trajectory_metadata)!=states.shape[0]:
        raise ValueError("Controlled-dataset stacked arrays have incompatible shapes.")
    trajectories=[ControlledTrajectory(time,states[index],inputs[index],trajectory_metadata[index])
        for index in range(states.shape[0])]
    return ControlledDataset(trajectories,metadata)


__all__=["load_controlled_dataset","save_controlled_dataset"]
