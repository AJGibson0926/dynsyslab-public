#include "WorkflowSupport.hpp"

#include <chrono>
#include <cctype>
#include <ctime>
#include <fstream>
#include <iomanip>
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

void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned,
                    const std::string& message) {
    if (path.empty()) return;
    std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    out << "{\"schema\":\"dynsyslab.compute-progress/v1\",\"determinate\":true,"
        << "\"percent\":" << fraction * 100.0
        << ",\"estimated_remaining_seconds\":-1,\"completed_trajectories\":"
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
