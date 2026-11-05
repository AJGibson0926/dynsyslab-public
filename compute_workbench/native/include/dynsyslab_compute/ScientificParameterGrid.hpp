#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <functional>
#include <map>
#include <string>

namespace dsl::compute {

// Compact, read-only scientific notation view. Parameter identifiers remain
// canonical in the model; this control only translates them for presentation.
using ScientificParameterEditHandler =
    std::function<void(const std::string& key, double value)>;

enum class ScientificParameterSection {
    physical,
    derived,
    dimensionless
};

HWND makeScientificParameterGrid(
    HWND parent, const dsl::theme::Fonts& fonts,
    ScientificParameterEditHandler onEdit = {});
void setScientificParameterGridSection(
    HWND grid, ScientificParameterSection section,
    const std::map<std::string, double>& values,
    const std::map<std::string, std::string>& units);
void setScientificParameterGridSectionMessage(
    HWND grid, ScientificParameterSection section, const std::wstring& message);
void setScientificParameterGridLongPrecision(HWND grid, bool enabled);

}  // namespace dsl::compute
