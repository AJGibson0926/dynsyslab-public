"""Closed registry of built-in post-processing method handlers."""
from .base import MethodHandler
from .decomposition import DECOMPOSITION_HANDLERS
from .koopman import KOOPMAN_HANDLERS

METHOD_HANDLERS={**KOOPMAN_HANDLERS,**DECOMPOSITION_HANDLERS}
SUPPORTED_METHODS=frozenset(METHOD_HANDLERS)


def handler_for(name: str) -> MethodHandler:
    try:
        return METHOD_HANDLERS[name]
    except KeyError as error:
        raise ValueError(f"Unsupported post-processing method {name!r}.") from error


__all__=["METHOD_HANDLERS","SUPPORTED_METHODS","MethodHandler","handler_for"]
