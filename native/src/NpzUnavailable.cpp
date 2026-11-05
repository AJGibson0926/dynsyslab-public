#include <dynsyslab/core/Archive.hpp>

#include <algorithm>
#include <stdexcept>

namespace dynsyslab::core {

NpzRecord readNpz(const std::filesystem::path&) {
    throw std::runtime_error("This DynSysLab build does not include native NPZ support");
}

bool NpzIndex::contains(const std::string& name) const {
    return std::find(arrays.begin(),arrays.end(),name)!=arrays.end();
}

NpzIndex readNpzIndex(const std::filesystem::path&) {
    throw std::runtime_error("This DynSysLab build does not include native NPZ support");
}

void writeNpzAtomic(const std::filesystem::path&, const NpzRecord&) {
    throw std::runtime_error("This DynSysLab build does not include native NPZ support");
}

void writeNpzAtomicNew(const std::filesystem::path& path, const NpzRecord& record) {
    writeNpzAtomic(path,record);
}

}  // namespace dynsyslab::core
