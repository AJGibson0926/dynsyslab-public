#include "dynsyslab_compute/SolverPanel.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_compute/ScientificParameterGrid.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/workflow/SystemResolution.hpp>

#include <map>
#include <stdexcept>

namespace dsl::compute {
void SolverPanel::updateParameterSummaries() {
    try {
        const std::size_t index = systemIndex();
        if (!registry_ || index >= registry_->systems.size())
            throw std::runtime_error("No dynamical system is selected.");
        const auto& definition = registry_->systems[index];
        dynsyslab::core::SystemConfiguration requested;
        requested.identifier = definition.name;
        const int presetIndex = static_cast<int>(
            SendMessageW(preset_, CB_GETCURSEL, 0, 0));
        if (presetIndex > 0
            && static_cast<std::size_t>(presetIndex - 1) < definition.presets.size())
            requested.preset = definition.presets[static_cast<std::size_t>(presetIndex - 1)].key;
        const auto overrides = dsl::json::parse(
            jsonObject(utf8Text(parameters_)));
        for (const auto& entry : overrides.object())
            requested.parameters[entry.first] = entry.second.number();
        const auto resolved = dynsyslab::workflow::resolveSystemDefinition(requested);
        setScientificParameterGridSection(parameterSummary_,
            ScientificParameterSection::physical,
            resolved.configuration.parameters,
            resolved.configuration.parameterUnits);
        std::map<std::string, double> groups, derived;
        std::map<std::string, std::string> groupUnits, derivedUnits;
        for (const auto& entry : resolved.configuration.dimensionlessGroups) {
            const auto unit = resolved.configuration.dimensionlessGroupUnits.find(entry.first);
            const bool dimensional = unit != resolved.configuration.dimensionlessGroupUnits.end()
                && !unit->second.empty() && unit->second != "dimensionless";
            (dimensional ? derived : groups)[entry.first] = entry.second;
            if (unit != resolved.configuration.dimensionlessGroupUnits.end())
                (dimensional ? derivedUnits : groupUnits)[entry.first] = unit->second;
        }
        if (!groups.empty()) setScientificParameterGridSection(parameterSummary_,
            ScientificParameterSection::dimensionless, groups, groupUnits);
        else setScientificParameterGridSectionMessage(parameterSummary_,
            ScientificParameterSection::dimensionless,
            resolved.descriptor.preferredCoordinates == "dimensionless"
                ? L"Dimensionless coordinates; no named groups declared."
                : L"No nondimensionalization declared.");
        if (!derived.empty()) setScientificParameterGridSection(parameterSummary_,
            ScientificParameterSection::derived, derived, derivedUnits);
        else setScientificParameterGridSectionMessage(parameterSummary_,
            ScientificParameterSection::derived,
            L"No dimensional scales derived.");
    } catch (const std::exception& error) {
        const std::wstring message = L"Unable to resolve parameters: "
            + dsl::wide(error.what());
        setScientificParameterGridSectionMessage(parameterSummary_,
            ScientificParameterSection::physical, message);
        setScientificParameterGridSectionMessage(parameterSummary_,
            ScientificParameterSection::dimensionless,
            L"Unavailable until parameters are valid.");
        setScientificParameterGridSectionMessage(parameterSummary_,
            ScientificParameterSection::derived,
            L"Unavailable until parameters are valid.");
    }
}

void SolverPanel::setParameterOverride(const std::string& key, double value) {
    auto overrides = dsl::json::parse(jsonObject(utf8Text(parameters_)));
    overrides.object()[key] = value;
    SetWindowTextW(parameters_, dsl::wide(dsl::json::dump(overrides)).c_str());
    clearTrainingPreflight();
    updateParameterSummaries();
}

}  // namespace dsl::compute
