#pragma once
#include <dynsyslab/core/KoopmanArchive.hpp>

namespace dynsyslab::core::detail {
SystemConfiguration systemFromMetadata(const NpzRecord& archive);
struct ArchiveBinding {
    SystemConfiguration system;
    std::string evidence;
};
ArchiveBinding resolveArchiveBinding(const NpzRecord& archive,
    const SystemConfiguration& requested, const ModelBindingOptions& options);
}
