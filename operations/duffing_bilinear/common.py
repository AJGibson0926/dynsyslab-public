"""Shared contract, paths, serialization, and metrics for the bounded search."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = ROOT / "examples" / "duffing-minimal-exact-klqr-transfer-release" / "manifest.json"
DEFAULT_CONTROL_SOURCE = ROOT / "archives" / "experiments" / "korda_duffing_candidate" / "korda_mezic_duffing_random_forced_full_v3.mat"
DEFAULT_OUTPUT = ROOT / "archives" / "experiments" / "duffing_bilinear_krom_search_20260829"

INITIAL_STATE = np.asarray((0.5, 0.0))
TARGET_STATE = np.asarray((-0.25, 0.0))
LEFT_EQUILIBRIUM = np.asarray((-0.5, 0.0))
ACTIVE_WINDOW = (1.0, 10.0)
FINAL_TIME = 20.0

def _sha(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imaginary": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, complex):
        return {"real": float(value.real), "imaginary": float(value.imag)}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)

def _write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_jsonable(payload), indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path

def _control_metrics(result) -> dict[str, Any]:
    time = np.asarray(result.time, dtype=float)
    state = np.asarray(result.values, dtype=float)
    command = np.asarray(result.control_time_series.values[:, 0], dtype=float)
    release = min(int(np.searchsorted(time, ACTIVE_WINDOW[1])), time.size - 1)
    settled = (time >= 2.0) & (time <= ACTIVE_WINDOW[1])
    active = (time >= ACTIVE_WINDOW[0]) & (time <= ACTIVE_WINDOW[1])
    call_values = getattr(result, "control_call_values", None)
    call_peak = None if call_values is None else float(np.max(np.abs(call_values)))
    return {
        "state_at_release": state[release],
        "distance_to_target_at_release": float(np.linalg.norm(state[release] - TARGET_STATE)),
        "final_state": state[-1],
        "final_distance_to_left_equilibrium": float(
            np.linalg.norm(state[-1] - LEFT_EQUILIBRIUM)
        ),
        "controlled_state_rmse_after_settling": float(
            np.sqrt(np.mean(np.sum((state[settled] - TARGET_STATE) ** 2, axis=1)))
        ),
        "control_peak_abs": float(np.max(np.abs(command))),
        "control_call_peak_abs": call_peak,
        "control_99th_percentile_abs_during_active": float(
            np.quantile(np.abs(command[active]), 0.99)
        ),
        "control_rms_during_active": float(np.sqrt(np.mean(command[active] ** 2))),
        "controller_call_count": int(result.metadata.get("controller_call_count") or 0),
    }
