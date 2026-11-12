"""Private viewing-session lifecycle and passive refresh coordination.

This service deliberately knows only source descriptors and the secondary
Viewer subscription registry.  It does not read, derive, or mutate primary
archive data; rendering remains in :mod:`viewer.core`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from ._compute import ViewerCompute
from ._sessions import poll_events, register_sources


class ViewingSession:
    """Own a disposable compute session and passive source subscriptions."""

    def __init__(self) -> None:
        self.identifier=ViewerCompute.start_viewing_session()
        self.compute=ViewerCompute(session_id=self.identifier)
        self._event_cursor=0

    def synchronize(self, sources: Iterable[Any]) -> None:
        records=[]
        for entry in sources:
            source=entry.source; path=None
            # Source classes are intentionally duck-typed here to avoid a
            # dependency back on the public Viewer data model.
            for name in ("path", "campaign_path", "manifest_path"):
                value=getattr(source,name,None)
                if value is not None:
                    path=value; break
            nested=getattr(source,"source",None)
            if path is None and isinstance(nested,(str,Path)):
                path=nested
            if path is not None:
                records.append({"name":entry.name,"data_type":entry.data_type,
                                "path":str(Path(path).expanduser().resolve())})
        register_sources(self.identifier,records)

    def poll(self) -> list[dict[str,Any]]:
        events,self._event_cursor=poll_events(self.identifier,after=self._event_cursor)
        return events
