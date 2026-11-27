#pragma once
#include <filesystem>
#include <vector>

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
    std::string systemTransferReason;
    // Scattered eigenfunction objects are interpolated on the one domain their method built them on.
    // Non-empty: that domain's name, checked (a mismatch is refused; the loader never re-domains).
    // For spectral K-ROMs, "state" selects the interpolation instead of the frozen expansion.
    std::string interpolationDomain;
    // Scattered K-ROMs: interpolant method object {"name", "neighbors", "power", "width"}; null keeps the archive's.
    json::Value interpolationMethod;
    // The representation the model evaluates, by name ("polynomial", "interpolation", ...); empty
    // keeps the krom's first.
    std::string representation;
};

// No conversion of a sampled input response into continuous B-psi is implied.
void requireContinuousControlModel(const KoopmanModel& model,
                                  const std::string& representation);
} // namespace dynsyslab::core
