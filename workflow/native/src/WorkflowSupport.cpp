#include "WorkflowSupport.hpp"
#include <dynsyslab/core/OutputPath.hpp>

#include <chrono>
#include <cctype>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow::detail {

std::string jsonEscape(const std::string& value) {
    std::string result;
    for (const char item : value) {
        if (item == '"' || item == '\\') result.push_back('\\');
        if (item == '\n') result += "\\n";
        else if (item != '\r') result.push_back(item);
    }
    return result;
}

std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t raw = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y%m%dT%H%M%S") << 'Z';
    return out.str();
}

std::string safeName(std::string value, const std::string& fallback) {
    for (char& item : value)
        if (!std::isalnum(static_cast<unsigned char>(item))
            && item != '-' && item != '_') item = '-';
    return value.empty() ? fallback : value;
}

namespace {
// Start time of each live progress file, for the heartbeat and remaining-time estimate.
std::map<std::wstring, std::chrono::steady_clock::time_point>& progressStarts() {
    static std::map<std::wstring, std::chrono::steady_clock::time_point> starts;
    return starts;
}
std::mutex& progressLock() { static std::mutex lock; return lock; }

std::string utcNow() {
    const std::time_t raw = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out; out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}
}  // namespace

void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned,
                    const std::string& message) {
    if (path.empty()) return;
    (void)dynsyslab::core::prepareOutputFile(path, "progress file");
    // Elapsed time since this file's first report (a report at zero restarts the clock);
    // the remaining time extrapolates the mean rate so far, -1 until there is progress.
    double elapsed = 0.0;
    {
        std::lock_guard<std::mutex> guard(progressLock());
        auto& starts = progressStarts();
        const auto now = std::chrono::steady_clock::now();
        auto found = starts.find(path.wstring());
        if (found == starts.end() || (fraction <= 0.0 && complete == 0))
            found = starts.insert_or_assign(path.wstring(), now).first;
        elapsed = std::chrono::duration<double>(now - found->second).count();
        if (fraction >= 1.0) starts.erase(found);
    }
    const double remaining = fraction > 0.0 && fraction < 1.0
        ? elapsed * (1.0 - fraction) / fraction : fraction >= 1.0 ? 0.0 : -1.0;
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    out << "{\"schema\":\"dynsyslab.compute-progress/v1\",\"determinate\":true,"
        << "\"percent\":" << fraction * 100.0
        << ",\"heartbeat_utc\":\"" << utcNow() << "\",\"elapsed_seconds\":" << elapsed
        << ",\"estimated_remaining_seconds\":" << remaining
        << ",\"seconds_per_trajectory\":" << (complete ? elapsed / complete : -1.0)
        << ",\"completed_trajectories\":"
        << complete << ",\"planned_trajectories\":" << planned
        << ",\"state\":\"" << (fraction >= 1.0 ? "complete" : "running")
        << "\",\"message\":\"" << jsonEscape(message) << "\"}";
    out.close();
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish workflow progress");
#else
    std::filesystem::rename(temporary, path);
#endif
}

}  // namespace dynsyslab::workflow::detail
