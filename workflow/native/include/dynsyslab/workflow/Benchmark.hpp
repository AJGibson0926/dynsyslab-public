#pragma once
// Native benchmark harness and performance gates (`dynsyslab-cli benchmark`).
// Each case runs `warmup` untimed iterations, then `iterations` timed ones; the report keeps
// every sample with min/median/max and the runtime fingerprint (build, compiler, machine).
// Medians are compared because they resist the outliers of a busy machine (Chen & Revels
// 2016, "Robust benchmarking in noisy environments", arXiv:1608.04295).

#include <dynsyslab/workflow/Json.hpp>

#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct BenchmarkOptions {
    std::size_t iterations{5};
    std::size_t warmup{1};
    std::vector<std::string> cases;  // empty: the whole kernel suite
};

// Names of the built-in kernel cases.
[[nodiscard]] std::vector<std::string> benchmarkCaseNames();
// dynsyslab.benchmark-report/v1
[[nodiscard]] json::Value runBenchmarks(const BenchmarkOptions& options);
// dynsyslab.performance-gate/v1: candidate median speedup over the reference per case
// (reference / candidate) must reach `required` (by case name, or "*" for every shared case).
// Runtimes must match (machine, build type, compiler) unless requireSameRuntime is false.
[[nodiscard]] json::Value benchmarkGate(const json::Value& reference, const json::Value& candidate,
                                        const json::Value& required, bool requireSameRuntime = true);

}  // namespace dynsyslab::workflow
