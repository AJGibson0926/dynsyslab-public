// Codec selection for written archives (Archive.hpp ArchiveCompression).
#include <dynsyslab/core/Archive.hpp>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

// -1: not set in this process; otherwise the ArchiveCompression value.
std::atomic<int>& processSetting() {
    static std::atomic<int> value{-1};
    return value;
}

}  // namespace

ArchiveCompression parseArchiveCompression(const std::string& name) {
    std::string value = name;
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (value == "deflate") return ArchiveCompression::deflate;
    if (value == "lzma") return ArchiveCompression::lzma;
    if (value == "stored" || value == "store") return ArchiveCompression::stored;
    throw std::invalid_argument("archive compression must be deflate, lzma or stored, not '" + name + "'");
}

std::string archiveCompressionName(ArchiveCompression compression) {
    switch (compression) {
    case ArchiveCompression::lzma: return "lzma";
    case ArchiveCompression::stored: return "stored";
    case ArchiveCompression::deflate: break;
    }
    return "deflate";
}

ArchiveCompression archiveCompression() {
    if (const int set = processSetting().load(); set >= 0) return static_cast<ArchiveCompression>(set);
    if (const char* environment = std::getenv("DYNSYSLAB_COMPRESSION"); environment && *environment)
        return parseArchiveCompression(environment);
    return ArchiveCompression::deflate;
}

void setArchiveCompression(ArchiveCompression compression) {
    processSetting().store(static_cast<int>(compression));
}

}  // namespace dynsyslab::core
