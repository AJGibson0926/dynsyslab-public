#pragma once

#include <dynsyslab/core/System.hpp>

namespace dynsyslab::workflow {

struct ResolvedSystemDefinition {
    dynsyslab::core::SystemDescriptor descriptor;
    dynsyslab::core::SystemConfiguration configuration;
    dynsyslab::core::CoordinateScaling scaling;
};

// Resolves presets, defaults, units, and scientific identity across every
// native system family. UI, planning, and adapters use this one boundary.
[[nodiscard]] ResolvedSystemDefinition resolveSystemDefinition(
    const dynsyslab::core::SystemConfiguration& requested);

}  // namespace dynsyslab::workflow
