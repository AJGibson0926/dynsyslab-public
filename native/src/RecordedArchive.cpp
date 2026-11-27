// readNpz / readNpzIndex with recorded archives (Archive.hpp): a derived archive that stores how to
// recompute its arrays is replayed by the reader its writer registered for its schema. Time arrays
// kept in a shared grid (TimeGrid.hpp) are restored first.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TimeGrid.hpp>

#include <map>
#include <mutex>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::mutex& registryMutex() {
    static std::mutex mutex;
    return mutex;
}
std::map<std::string, RecordedArchiveReader>& registry() {
    static std::map<std::string, RecordedArchiveReader> readers;
    return readers;
}
constexpr const char* kKey = "\"recorded_archive\"";

}  // namespace

void registerRecordedArchiveReader(const std::string& schema, RecordedArchiveReader reader) {
    std::lock_guard<std::mutex> lock(registryMutex());
    registry()[schema] = std::move(reader);
}

bool isRecordedArchive(const std::string& metadataJson) {
    if (metadataJson.find(kKey) == std::string::npos) return false;
    return !json::parse(metadataJson).at("recorded_archive").isNull();
}

NpzRecord readNpz(const std::filesystem::path& path) {
    auto stored = readNpzStored(path);
    attachTimeGrids(path, stored);  // archives that reference a shared time grid (TimeGrid.hpp)
    const auto text = stored.metadataJson();
    if (!isRecordedArchive(text)) return stored;
    const auto schema = json::parse(text).at("recorded_archive").at("schema").string();
    RecordedArchiveReader reader;
    {
        std::lock_guard<std::mutex> lock(registryMutex());
        if (const auto found = registry().find(schema); found != registry().end()) reader = found->second;
    }
    if (!reader)
        throw std::runtime_error(path.u8string() + " is a recorded archive (" + schema +
                                 "): its arrays are recomputed by the DynSysLab workflow library, which is not loaded here. "
                                 "Use `dynsyslab-cli materialize` to write a self-contained copy.");
    return reader(path, std::move(stored));
}

NpzIndex readNpzIndex(const std::filesystem::path& path) {
    auto index = readNpzIndexStored(path);
    addTimeGridNames(index);
    if (!isRecordedArchive(index.metadataJson)) return index;
    const auto metadata = json::parse(index.metadataJson);  // kept alive: the loop reads into it
    for (const auto& name : metadata.at("recorded_archive").at("arrays").array())
        if (!index.contains(name.string())) index.arrays.push_back(name.string());
    return index;
}

}  // namespace dynsyslab::core
