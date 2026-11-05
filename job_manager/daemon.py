"""Single-instance local scheduler service for polling and watchdog refresh."""
from __future__ import annotations
from pathlib import Path
from time import sleep
from typing import Any
import json, os

class JobManagerDaemon:
    """Small local polling loop that publishes watchdog health."""
    def __init__(self, manager: Any, *, interval_seconds: float=10.):
        self.manager=manager; self.interval_seconds=max(.1,float(interval_seconds)); self.path=Path(manager.root)/"daemon.json"
    def tick(self) -> dict[str,Any]:
        report={"pid":os.getpid(),"watchdog":self.manager.watch_paused_workers()}
        self.path.write_text(json.dumps(report,default=str),encoding="utf-8")
        return report
    def serve(self, *, max_ticks: int|None=None) -> None:
        count=0
        while max_ticks is None or count<max_ticks:
            self.tick(); count+=1; sleep(self.interval_seconds)
