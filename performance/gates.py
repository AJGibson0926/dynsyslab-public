"""Performance admission/regression gates for versioned benchmark reports."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class PerformanceGateResult:
    """Admission result for one required candidate speedup."""
    name: str
    reference_seconds: float
    candidate_seconds: float
    speedup: float
    required_speedup: float
    passed: bool


@dataclass(frozen=True)
class PerformanceGateReport:
    """Versioned collection of performance admissions for comparable runs."""
    schema: str
    comparable_runtime: bool
    runtime_differences: dict[str, tuple[Any, Any]]
    results: tuple[PerformanceGateResult, ...]

    @property
    def passed(self) -> bool:
        return self.comparable_runtime and all(result.passed for result in self.results)

    def payload(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "comparable_runtime": self.comparable_runtime,
            "runtime_differences": self.runtime_differences,
            "passed": self.passed,
            "results": [asdict(result) for result in self.results],
        }


def evaluate_speedups(
    reference: Mapping[str, Any],
    candidate: Mapping[str, Any],
    required_speedups: Mapping[str, float],
    *,
    require_same_runtime: bool = True,
) -> PerformanceGateReport:
    """Compare matching benchmark medians across two versioned reports."""
    if reference.get("schema") != "dynsyslab.benchmark-report/v1" or candidate.get("schema") != "dynsyslab.benchmark-report/v1":
        raise ValueError("Performance gates require dynsyslab.benchmark-report/v1 payloads.")
    keys = ("machine", "platform", "python", "numpy", "scipy")
    reference_runtime = reference.get("runtime") or {}
    candidate_runtime = candidate.get("runtime") or {}
    differences = {
        key: (reference_runtime.get(key), candidate_runtime.get(key))
        for key in keys
        if reference_runtime.get(key) != candidate_runtime.get(key)
    }
    reference_samples = {sample["name"]: sample for sample in reference.get("samples", ())}
    candidate_samples = {sample["name"]: sample for sample in candidate.get("samples", ())}
    results = []
    for name, required in required_speedups.items():
        if required <= 0:
            raise ValueError("Required speedups must be positive.")
        if name not in reference_samples or name not in candidate_samples:
            raise KeyError(f"Benchmark {name!r} must be present in both reports.")
        before = float(reference_samples[name]["median_seconds"])
        after = float(candidate_samples[name]["median_seconds"])
        speedup = float("inf") if after == 0.0 else before / after
        results.append(PerformanceGateResult(name, before, after, speedup, float(required), speedup >= required))
    return PerformanceGateReport(
        "dynsyslab.performance-gate/v1",
        not differences or not require_same_runtime,
        differences,
        tuple(results),
    )


def evaluate_relative_speedups(
    report: Mapping[str, Any],
    requirements: Mapping[tuple[str, str], float],
) -> PerformanceGateReport:
    """Compare candidate/reference cases measured in the same report."""
    if report.get("schema") != "dynsyslab.benchmark-report/v1":
        raise ValueError("Performance gates require a dynsyslab.benchmark-report/v1 payload.")
    samples = {sample["name"]: sample for sample in report.get("samples", ())}
    results = []
    for (candidate_name, reference_name), required in requirements.items():
        if required <= 0:
            raise ValueError("Required speedups must be positive.")
        missing = [name for name in (candidate_name, reference_name) if name not in samples]
        if missing:
            raise KeyError(f"Benchmarks missing from report: {missing}.")
        before = float(samples[reference_name]["median_seconds"])
        after = float(samples[candidate_name]["median_seconds"])
        speedup = float("inf") if after == 0.0 else before / after
        results.append(PerformanceGateResult(
            f"{candidate_name}/{reference_name}", before, after, speedup,
            float(required), speedup >= required,
        ))
    return PerformanceGateReport(
        "dynsyslab.relative-performance-gate/v1", True, {}, tuple(results),
    )


__all__ = ["PerformanceGateReport", "PerformanceGateResult", "evaluate_speedups", "evaluate_relative_speedups"]
