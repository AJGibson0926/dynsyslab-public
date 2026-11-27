#pragma once

#include <dynsyslab/core/System.hpp>

#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
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

// A non-ODE model supplied by an extension (extension_api.h delay, field and lattice tables).
// Callbacks throw std::runtime_error with the extension's message when it reports a failure.
struct ExtensionModel {
    std::string family;          // "dde", "pde" or "lattice"
    std::string identifier;
    std::string label;
    std::string descriptorJson;
    std::function<void(const double* window, std::size_t windowCount, double* derivative,
                       std::size_t stateCount)> delayRhs;
    std::function<std::vector<double>(const std::string& configurationJson, std::size_t count)> initialField;
    std::function<void(double time, const double* field, std::size_t count, double* derivative)> fieldRhs;
    std::function<void(std::uint64_t seed, std::int32_t* sites, std::size_t count)> latticeInitial;
    std::function<void(std::uint64_t seed, std::uint64_t sweep, std::int32_t* sites, std::size_t count,
                       double* observables, std::size_t observableCount)> latticeSweep;
};

class ExtensionRegistry {
public:
    static ExtensionRegistry& instance();
    std::string load(const std::filesystem::path& library);
    [[nodiscard]] std::vector<ExtensionDescriptor> extensions() const;
    [[nodiscard]] std::vector<SystemDescriptor> systems() const;
    [[nodiscard]] std::unique_ptr<ODESystem> createSystem(
        const SystemConfiguration& configuration) const;
    [[nodiscard]] std::vector<ExtensionModel> models() const;
    [[nodiscard]] std::optional<ExtensionModel> model(const std::string& identifier) const;
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
