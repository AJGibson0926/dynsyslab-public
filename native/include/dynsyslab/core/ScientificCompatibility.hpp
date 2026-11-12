#pragma once

#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::core {

// Compare resolved physics, not preset labels, initial conditions or forcing.
// Returns every differing contract field, for useful preflight diagnostics.
[[nodiscard]] std::vector<std::string> systemCompatibilityDifferences(
    const SystemConfiguration& expected, const SystemConfiguration& actual);
void requireCompatibleSystems(const SystemConfiguration& expected,
                              const SystemConfiguration& actual);

// Content identity includes polynomial coordinates/normalization, ordered modes,
// spectrum and time convention, optional physical system, and input dimension.
// It excludes storage paths, data pointers, fit errors and attached input maps.
[[nodiscard]] std::string koopmanModelIdentity(const KoopmanModel& model);

struct ModelBindingOptions {
    // Legacy explicit-system API remains readable, but labels missing provenance.
    // Workflow consumers disable this unless the user supplies an explicit reason.
    bool allowUnverifiedLegacy{true};
    std::string unverifiedLegacyReason;
    std::string systemTransferReason;
};

// No conversion of a sampled input response into continuous B-psi is implied.
void requireContinuousControlModel(const KoopmanModel& model,
                                  const std::string& representation);
} // namespace dynsyslab::core
