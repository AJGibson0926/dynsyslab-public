#include <dynsyslab/core/Execution.hpp>

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <string>
#include <thread>

namespace dynsyslab::core {
namespace {
std::atomic<std::size_t>& limit() {
    static std::atomic<std::size_t> value{0};
    return value;
}

void setEnvironment(const char* name, const std::string& value) {
#ifdef _WIN32
    _putenv_s(name, value.c_str());
#else
    setenv(name, value.c_str(), 1);
#endif
}
}  // namespace

void setProcessThreadLimit(std::size_t threads) {
    limit() = threads;
    if (!threads) return;
    const auto text = std::to_string(threads);
    for (const char* name : {"OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                             "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"})
        setEnvironment(name, text);
}

std::size_t processThreadLimit() { return limit().load(); }

std::size_t workerThreads(std::size_t wanted) {
    const std::size_t hardware = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    const std::size_t cap = limit().load();
    return std::max<std::size_t>(1, std::min({wanted, hardware, cap ? cap : hardware}));
}

}  // namespace dynsyslab::core
