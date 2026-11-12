"""Small benchmark harness with explicit warmup and distribution reporting."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from statistics import median
from time import perf_counter
from typing import Any, Callable, Iterable

from .runtime import runtime_fingerprint


@dataclass(frozen=True)
class BenchmarkCase:
    """One repeatable benchmark callable with warmup and sample counts."""
    name: str
    function: Callable[[], Any]
    category: str
    iterations: int = 5
    warmup: int = 1

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("BenchmarkCase.name cannot be empty.")
        if self.iterations < 1 or self.warmup < 0:
            raise ValueError("Benchmark iterations must be positive and warmup cannot be negative.")


@dataclass(frozen=True)
class BenchmarkSample:
    """Timing distribution recorded for one benchmark case."""
    name: str
    category: str
    seconds: tuple[float, ...]
    minimum_seconds: float
    median_seconds: float
    maximum_seconds: float


@dataclass(frozen=True)
class BenchmarkReport:
    """Versioned benchmark result with runtime identity and timing samples."""
    schema: str
    created_utc: str
    runtime: dict[str, Any]
    samples: tuple[BenchmarkSample, ...]

    def payload(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "created_utc": self.created_utc,
            "runtime": self.runtime,
            "samples": [asdict(sample) for sample in self.samples],
        }


def run_benchmarks(cases: Iterable[BenchmarkCase]) -> BenchmarkReport:
    """Run benchmark cases serially after their requested warmup iterations."""
    samples = []
    for case in cases:
        for _ in range(case.warmup):
            case.function()
        timings = []
        for _ in range(case.iterations):
            started = perf_counter()
            case.function()
            timings.append(perf_counter() - started)
        values = tuple(timings)
        samples.append(
            BenchmarkSample(
                case.name,
                case.category,
                values,
                min(values),
                median(values),
                max(values),
            )
        )
    return BenchmarkReport(
        "dynsyslab.benchmark-report/v1",
        datetime.now(timezone.utc).isoformat(),
        runtime_fingerprint(),
        tuple(samples),
    )
