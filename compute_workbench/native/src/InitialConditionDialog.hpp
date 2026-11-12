#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>
#include <vector>

namespace dsl::compute {

// Edits state-vector explicit values or random-sampling bounds. The solver's
// deterministic generator remains the sole owner of random draws.
std::optional<std::string> showInitialConditionDialog(
    HWND owner, const std::string& mode, std::size_t trajectoryCount,
    const std::vector<std::string>& stateLabels,
    const std::vector<double>& defaultState, const std::string& current,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
