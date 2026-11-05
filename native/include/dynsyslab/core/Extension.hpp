#pragma once

#include <dynsyslab/core/System.hpp>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ExtensionDescriptor {
    std::string identifier;
    std::string version;
    std::string capabilitiesJson;
    std::filesystem::path library;
};

class ExtensionRegistry {
public:
    static ExtensionRegistry& instance();
    std::string load(const std::filesystem::path& library);
    [[nodiscard]] std::vector<ExtensionDescriptor> extensions() const;
    [[nodiscard]] std::vector<SystemDescriptor> systems() const;
    [[nodiscard]] std::unique_ptr<ODESystem> createSystem(
        const SystemConfiguration& configuration) const;
    [[nodiscard]] std::string execute(const std::string& extension,
                                      const std::string& capability,
                                      const std::string& requestJson) const;

private:
    ExtensionRegistry();
    ~ExtensionRegistry();
    ExtensionRegistry(const ExtensionRegistry&) = delete;
    ExtensionRegistry& operator=(const ExtensionRegistry&) = delete;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace dynsyslab::core
