#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace dynsyslab::core {

struct ExecutionLimits {
    std::size_t maximumThreads{};
    std::size_t maximumMemoryBytes{};
};

struct ExecutionProgress {
    std::string stage;
    std::string message;
    double fraction{};
    std::size_t completed{};
    std::size_t planned{};
};

// Process-wide worker-thread cap set by the queue and CLI (`--max-threads`, job key
// `execution.max_threads`); zero means no cap. It also sets OMP_NUM_THREADS and the BLAS
// thread variables for libraries started afterwards. Kernels that spawn threads ask
// workerThreads (defined in Execution.cpp so every module shares one cap).
void setProcessThreadLimit(std::size_t threads);
[[nodiscard]] std::size_t processThreadLimit();
[[nodiscard]] std::size_t workerThreads(std::size_t wanted);

using CancellationProbe = std::function<bool()>;
using ExecutionProgressCallback = std::function<void(const ExecutionProgress&)>;

class ExecutionCancelled final : public std::runtime_error {
public:
    ExecutionCancelled() : std::runtime_error("DynSysLab execution cancelled") {}
};

class ExecutionLimitExceeded final : public std::runtime_error {
public:
    explicit ExecutionLimitExceeded(const std::string& message)
        : std::runtime_error(message) {}
};

class ExecutionContext {
public:
    ExecutionLimits limits;
    CancellationProbe cancellation;
    ExecutionProgressCallback progress;

    [[nodiscard]] bool cancelled() const {
        return cancellation && cancellation();
    }
    void checkpoint() const {
        if (cancelled()) throw ExecutionCancelled{};
    }
    void report(ExecutionProgress event) const {
        checkpoint();
        if (progress) progress(event);
    }
    void requireMemory(std::size_t bytes, std::string_view purpose) const {
        if (limits.maximumMemoryBytes && bytes > limits.maximumMemoryBytes)
            throw ExecutionLimitExceeded(
                "DynSysLab memory limit exceeded by " + std::string(purpose));
    }
};

}  // namespace dynsyslab::core
