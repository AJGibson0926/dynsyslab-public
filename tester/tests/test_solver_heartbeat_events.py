"""Characterize the low-rate heartbeat event exposed to progress clients."""
from __future__ import annotations

from types import SimpleNamespace
import unittest
from unittest.mock import patch

from solver.models import HistoryOptions
from solver.ode import _TelemetryHeartbeat


class _Session:
    def __init__(self) -> None:
        self.events: list[tuple[str, dict]] = []

    def telemetry(self, event: str, **details) -> None:
        self.events.append((event, details))


class SolverHeartbeatEventTests(unittest.TestCase):
    def test_throttled_telemetry_heartbeat_reuses_on_event_channel(self):
        events: list[dict] = []
        request = SimpleNamespace(
            history=HistoryOptions(telemetry_interval_seconds=1.0),
            on_event=events.append,
        )
        session = _Session()
        with patch("solver.ode.monotonic", side_effect=(10.0, 10.2, 11.1)):
            heartbeat = _TelemetryHeartbeat(request, session, 0.0, 20.0)
            heartbeat.report(2.0)
            heartbeat.report(10.0, accepted_steps=17)
        self.assertEqual(len(session.events), 1)
        self.assertEqual(session.events[0][0], "heartbeat")
        self.assertEqual(events, [{
            "event": "heartbeat", "elapsed_wall_seconds": 1.1,
            "simulation_time": 10.0, "progress_fraction": 0.5,
            "accepted_steps": 17,
        }])


if __name__ == "__main__":
    unittest.main()
