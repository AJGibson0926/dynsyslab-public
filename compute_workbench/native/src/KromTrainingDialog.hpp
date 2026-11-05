#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace dsl::compute {

struct KromTrainingSetup {
    std::string label;
    std::string system;
    std::string preset;
    std::string parametersJson;
    std::string coordinateSystem;
    std::string initialConditionsJson;
    std::size_t trajectoryCount{};
    std::string finalTime;
    std::string outputStep;
    std::string preflightJson;
    std::size_t selectedKrom{};
};

// Plans a system-agnostic KROM training-domain expansion directly from the
// selected KROM's retained training states. The returned setup is complete:
// system identity, exact IC grid, integration settings, and safety preflight.
std::optional<KromTrainingSetup> showKromTrainingDialog(
    HWND owner, const std::vector<KromDefinition>& compatibleKroms,
    std::size_t initialSelection, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
