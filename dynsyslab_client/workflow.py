"""Dependency-free ctypes binding to the canonical DynSysLab workflow ABI."""
from __future__ import annotations

from collections.abc import Callable, Mapping
import ctypes
from dataclasses import dataclass
import json
import os
from pathlib import Path
from typing import Any

WORKFLOW_ABI_VERSION = 4
STATUS_OK, STATUS_BUFFER_TOO_SMALL = 0, 2
STATUS_CANCELLED, STATUS_RESOURCE_LIMIT = 4, 6
_CANCEL = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p)
_PROGRESS = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double,
                             ctypes.c_size_t, ctypes.c_size_t, ctypes.c_void_p)


class _Callbacks(ctypes.Structure):
    _fields_ = [("cancellation", _CANCEL), ("progress", _PROGRESS),
                ("user_data", ctypes.c_void_p), ("maximum_threads", ctypes.c_size_t),
                ("maximum_memory_bytes", ctypes.c_size_t)]


@dataclass(frozen=True)
class NativeProgress:
    """One progress event emitted by compiled execution."""

    stage: str
    message: str
    fraction: float
    completed: int
    planned: int


class NativeWorkflowError(RuntimeError):
    """A status and diagnostic returned by the compiled workflow."""

    def __init__(self, status: int, message: str):
        self.status = status
        super().__init__(message or f"native workflow failed with status {status}")


class NativeWorkflowCancelled(NativeWorkflowError):
    """Compiled execution observed a cancellation request."""


class NativeWorkflowResourceLimit(NativeWorkflowError):
    """Compiled execution rejected or exceeded a resource limit."""


def _names() -> tuple[str, ...]:
    if os.name == "nt":
        return ("dynsyslab_workflow_c_api.dll",)
    if os.uname().sysname == "Darwin":
        return ("libdynsyslab_workflow_c_api.dylib",)
    return ("libdynsyslab_workflow_c_api.so",)


def installation_root() -> Path:
    """Return the explicit or conventional compiled installation root."""
    override = os.environ.get("DYNSYSLAB_INSTALL_DIR")
    if override:
        return Path(override).expanduser().resolve()
    if os.name == "nt":
        return (Path(os.environ.get("LOCALAPPDATA", Path.home())) /
                "Programs" / "DynSysLab").resolve()
    return (Path.home() / ".local" / "lib" / "dynsyslab").resolve()


def library_candidates(project_root: str | Path | None = None) -> tuple[Path, ...]:
    """Return deterministic library locations without searching ``PATH``."""
    override = os.environ.get("DYNSYSLAB_WORKFLOW_LIBRARY")
    if override:
        return (Path(override).expanduser().resolve(),)
    roots = [installation_root() / "bin"]
    if project_root:
        project = Path(project_root).expanduser().resolve()
        roots.extend((project / "build" / "native-suite-fast" / "bin",
                      project / "build" / "native-suite" / "bin" / "Release"))
    return tuple(dict.fromkeys((root / name).resolve()
                               for root in roots for name in _names()))


def find_library(project_root: str | Path | None = None) -> Path:
    """Locate the compiled workflow library or raise an actionable error."""
    candidates = library_candidates(project_root)
    for path in candidates:
        if path.is_file():
            return path
    locations = "\n  ".join(map(str, candidates))
    raise FileNotFoundError(
        "The compiled DynSysLab workflow is not installed. Install DynSysLab "
        "or set DYNSYSLAB_WORKFLOW_LIBRARY.\n  " + locations)


def _text(value: bytes | None) -> str:
    return value.decode("utf-8", errors="replace") if value else ""


class NativeWorkflow:
    """ABI-checked handle to the canonical compiled workflow service."""

    def __init__(self, library: str | Path | None = None,
                 project_root: str | Path | None = None):
        self.path = Path(library).resolve() if library else find_library(project_root)
        self._search = (os.add_dll_directory(str(self.path.parent))
                        if os.name == "nt" and hasattr(os, "add_dll_directory") else None)
        self._library = ctypes.CDLL(str(self.path))
        self._declare()
        found = int(self._library.dsl_workflow_abi_version())
        if found != WORKFLOW_ABI_VERSION:
            raise RuntimeError(f"workflow ABI mismatch: expected {WORKFLOW_ABI_VERSION}, got {found}")

    def _declare(self) -> None:
        api = self._library
        api.dsl_workflow_abi_version.argtypes = []
        api.dsl_workflow_abi_version.restype = ctypes.c_uint32
        api.dsl_workflow_plan_json.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p)]
        api.dsl_workflow_plan_json.restype = ctypes.c_int
        api.dsl_workflow_run_json.argtypes = [ctypes.c_char_p, ctypes.c_char_p,
            ctypes.c_char_p, ctypes.POINTER(_Callbacks), ctypes.POINTER(ctypes.c_void_p)]
        api.dsl_workflow_run_json.restype = ctypes.c_int
        api.dsl_workflow_result_destroy.argtypes = [ctypes.c_void_p]
        api.dsl_workflow_result_json.argtypes = [ctypes.c_void_p, ctypes.c_void_p,
                                                  ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
        api.dsl_workflow_result_json.restype = ctypes.c_int
        api.dsl_workflow_last_error.argtypes = [ctypes.c_void_p, ctypes.c_size_t,
                                                 ctypes.POINTER(ctypes.c_size_t)]
        api.dsl_workflow_last_error.restype = ctypes.c_int

    def _read(self, function: Any, *prefix: Any) -> str:
        required = ctypes.c_size_t()
        status = int(function(*prefix, None, 0, ctypes.byref(required)))
        if status not in (STATUS_OK, STATUS_BUFFER_TOO_SMALL):
            return ""
        buffer = ctypes.create_string_buffer(max(1, required.value))
        if int(function(*prefix, buffer, len(buffer), ctypes.byref(required))) != STATUS_OK:
            return ""
        return _text(buffer.value)

    def last_error(self) -> str:
        """Return the current thread's diagnostic from the compiled library."""
        return self._read(self._library.dsl_workflow_last_error)

    def _error(self, status: int, output: ctypes.c_void_p) -> None:
        if output:
            self._library.dsl_workflow_result_destroy(output)
        kind = (NativeWorkflowCancelled if status == STATUS_CANCELLED else
                NativeWorkflowResourceLimit if status == STATUS_RESOURCE_LIMIT else
                NativeWorkflowError)
        raise kind(status, self.last_error())

    def _result(self, output: ctypes.c_void_p) -> dict[str, Any]:
        try:
            return json.loads(self._read(self._library.dsl_workflow_result_json, output))
        finally:
            self._library.dsl_workflow_result_destroy(output)

    def observable(self, query: Mapping[str, Any]) -> dict[str, Any]:
        """Evaluate values and analytic derivatives through the compiled engine."""
        function = getattr(self._library, "dsl_workflow_observable_json", None)
        if function is None:
            raise RuntimeError("This native installation predates observable queries; update DynSysLab.")
        function.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p)]
        function.restype = ctypes.c_int
        request = {"schema": "dynsyslab.observable-query/v1", **dict(query)}
        output = ctypes.c_void_p()
        status = int(function(json.dumps(request, separators=(",", ":")).encode(), ctypes.byref(output)))
        if status != STATUS_OK:
            self._error(status, output)
        return self._result(output)

    def koopman_spectral(self, query: Mapping[str, Any]) -> dict[str, Any]:
        """Fit, inspect, evaluate or roll out a spectral KROM in the compiled engine.

        Arrays are row-major; complex values use [real, imag]. See SPECTRAL_KROM.md.
        """
        function = getattr(self._library, "dsl_workflow_koopman_spectral_json", None)
        if function is None:
            raise RuntimeError("This native installation predates spectral KROM queries.")
        function.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p)]
        function.restype = ctypes.c_int
        request = {"schema": "dynsyslab.koopman-spectral-query/v1", **dict(query)}
        output = ctypes.c_void_p()
        status = int(function(json.dumps(request, separators=(",", ":")).encode(), ctypes.byref(output)))
        if status != STATUS_OK:
            self._error(status, output)
        return self._result(output)

    def plan(self, job: Mapping[str, Any]) -> dict[str, Any]:
        """Validate and plan one job without numerical execution."""
        encoded = json.dumps(dict(job), separators=(",", ":")).encode()
        output = ctypes.c_void_p()
        status = int(self._library.dsl_workflow_plan_json(encoded, ctypes.byref(output)))
        if status != STATUS_OK:
            self._error(status, output)
        return self._result(output)

    def run(self, job: Mapping[str, Any], archive_root: str | Path, *,
            progress_path: str | Path | None = None,
            progress: Callable[[NativeProgress], None] | None = None,
            cancelled: Callable[[], bool] | None = None, maximum_threads: int = 0,
            maximum_memory_bytes: int = 0) -> dict[str, Any]:
        """Execute one job entirely inside the compiled service."""
        callback_error: list[BaseException] = []
        def cancel(_: object) -> int:
            try:
                return int(bool(callback_error or (cancelled and cancelled())))
            except BaseException as error:
                callback_error.append(error); return 1
        def report(stage: bytes | None, message: bytes | None, fraction: float,
                   completed: int, planned: int, _: object) -> None:
            if progress:
                try:
                    progress(NativeProgress(_text(stage), _text(message), fraction,
                                            int(completed), int(planned)))
                except BaseException as error:
                    callback_error.append(error)
        cancel_cb, progress_cb = _CANCEL(cancel), _PROGRESS(report)
        callbacks = _Callbacks(cancel_cb, progress_cb, None, max(0, maximum_threads),
                               max(0, maximum_memory_bytes))
        output = ctypes.c_void_p()
        status = int(self._library.dsl_workflow_run_json(
            json.dumps(dict(job), separators=(",", ":")).encode(),
            str(Path(archive_root).expanduser().resolve()).encode(),
            str(Path(progress_path).expanduser().resolve()).encode() if progress_path else None,
            ctypes.byref(callbacks), ctypes.byref(output)))
        if callback_error:
            if output:
                self._library.dsl_workflow_result_destroy(output)
            raise callback_error[0]
        if status != STATUS_OK:
            self._error(status, output)
        result = self._result(output)
        return _compatibility_fields(result)


def _compatibility_fields(result: dict[str, Any]) -> dict[str, Any]:
    """Add stable convenience aliases without reproducing scientific work."""
    result.update({"status": "complete", "execution_backend": "native-cpp",
                   "python_required": False})
    archives = list(result.get("archives", ()))
    if result.get("kind") == "postprocess":
        result["derived_archives"] = archives
    else:
        result["solution_archives"] = archives
        result["trajectory_count"] = len(archives)
    if result.get("kind") == "solve":
        result["system_count"] = 1
        result["hierarchy_archive"] = result.get("hierarchy_manifest")
    elif result.get("kind") == "campaign":
        path = Path(result.get("manifest", ""))
        manifest = (json.loads(path.read_text(encoding="utf-8"))
                    if path.is_file() else {})
        points = list(manifest.get("points", ()))
        result["campaign_manifest"] = str(path)
        result["campaign_id"] = manifest.get("id")
        result["system_count"] = len(points)
        result["hierarchy_archives"] = [
            point["hierarchy_manifest"] for point in points
            if point.get("hierarchy_manifest")
        ]
    return result


def plan_job(job: Mapping[str, Any], **options: Any) -> dict[str, Any]:
    """Plan through the compiled engine using a short functional API."""
    library, project = options.pop("library", None), options.pop("project_root", None)
    if options:
        raise TypeError(f"unknown native plan options: {sorted(options)}")
    return NativeWorkflow(library, project).plan(job)


def run_job(job: Mapping[str, Any], archive_root: str | Path,
            **options: Any) -> dict[str, Any]:
    """Run through the compiled engine using a short functional API."""
    return NativeWorkflow(options.pop("library", None),
                          options.pop("project_root", None)).run(job, archive_root, **options)
