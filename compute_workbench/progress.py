"""Low-rate, atomic progress snapshots for external compute clients."""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
from time import monotonic
from typing import Any


PROGRESS_SCHEMA = "dynsyslab.compute-progress/v1"


class ProgressReporter:
    """Reduce solver events to one cheap, pollable JSON status document.

    Numerical state is never retained. Writes occur only for already-throttled
    heartbeats, trajectory boundaries, and lifecycle changes.
    """

    def __init__(self, path: str | Path, *, kind: str,
                 planned_trajectories: int = 0) -> None:
        self.path = Path(path).expanduser().resolve()
        self.kind = str(kind)
        self.planned = max(0, int(planned_trajectories))
        self.completed = 0
        self.current_fraction = 0.0
        self.sequence = 0
        self.started = monotonic()
        self.determinate = False
        self._publish("running", "Starting compute job", force=True)

    @property
    def fraction(self) -> float:
        if self.planned <= 0:
            return 0.0
        value = (self.completed + self.current_fraction) / self.planned
        return min(1.0, max(0.0, value))

    def event(self, payload: dict[str, Any]) -> None:
        """Consume one ordinary solver event without changing solver behavior."""
        event = str(payload.get("event", "progress"))
        if event == "started":
            self.current_fraction = 0.0
            message = "Integrating trajectory"
        elif event == "heartbeat":
            self.current_fraction = min(1.0, max(
                0.0, float(payload.get("progress_fraction", 0.0))))
            self.determinate = True
            message = "Integrating trajectory"
        elif event == "hierarchy_progress":
            self.completed = max(self.completed, int(
                payload.get("completed_trajectories", self.completed)))
            self.planned = max(self.planned, int(
                payload.get("total_trajectories", self.planned)))
            self.current_fraction = 0.0
            self.determinate = True
            message = f"Archived trajectory {self.completed} of {self.planned}"
        elif event == "completed":
            self.current_fraction = 1.0
            message = "Finalizing trajectory archive"
        elif event == "failed":
            self._publish("failed", str(payload.get("message", "Compute failed")),
                          event=payload, force=True)
            return
        else:
            message = event.replace("_", " ").capitalize()
        self._publish("running", message, event=payload)

    def complete(self, message: str = "Compute job complete") -> None:
        """Publish a terminal complete snapshot."""
        if self.planned:
            self.completed = self.planned
        self.current_fraction = 0.0
        self.determinate = True
        self._publish("complete", message, force=True)

    def fail(self, message: str) -> None:
        """Publish a terminal failure snapshot."""
        self._publish("failed", message, force=True)

    def _publish(self, state: str, message: str, *, event: dict[str, Any] | None = None,
                 force: bool = False) -> None:
        del force  # documents that lifecycle writes intentionally bypass no gate
        self.sequence += 1
        elapsed = max(0.0, monotonic() - self.started)
        fraction = 1.0 if state == "complete" else self.fraction
        eta = None
        if state == "running" and fraction > 0.0:
            eta = max(0.0, elapsed * (1.0 - fraction) / fraction)
        payload: dict[str, Any] = {
            "schema": PROGRESS_SCHEMA,
            "state": state,
            "kind": self.kind,
            "sequence": self.sequence,
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "determinate": bool(self.determinate or state == "complete"),
            "fraction": fraction,
            "percent": 100.0 * fraction,
            "completed_trajectories": self.completed,
            "planned_trajectories": self.planned,
            "elapsed_seconds": elapsed,
            "estimated_remaining_seconds": eta,
            "message": message,
        }
        if event:
            for key in ("simulation_time", "accepted_steps", "system_index",
                        "cluster_index", "trajectory_index"):
                if key in event:
                    payload[key] = event[key]
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_name(self.path.name + ".tmp")
        temporary.write_text(json.dumps(payload, separators=(",", ":")),
                             encoding="utf-8")
        os.replace(temporary, self.path)


__all__ = ["PROGRESS_SCHEMA", "ProgressReporter"]
