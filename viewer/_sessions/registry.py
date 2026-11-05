"""Compatibility imports for the former private Viewer event registry."""

from infrastructure.source_events import (
    notify_source_update,
    poll_events,
    register_sources,
)

__all__ = ["notify_source_update", "poll_events", "register_sources"]
