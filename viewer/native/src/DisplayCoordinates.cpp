#include "DisplayCoordinates.hpp"

#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/SystemLabels.hpp>
#include <dynsyslab_viewer/Model.hpp>
#include <dynsyslab_viewer/MathText.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>

namespace dsl {
namespace {

std::string token(std::string value) {
    std::string result;
    for (const unsigned char item : value)
        if (std::isalnum(item)) result.push_back(
            static_cast<char>(std::tolower(item)));
    return result;
}

std::string identifier(const json::Value& metadata) {
    auto value = metadata.at("system").string(
        metadata.at("system_configuration").at("identifier").string());
    const auto system = token(value);
    if (system == "rpe" || system == "rayleighplesset") return "rpe";
    if (system == "reddyszeri") return "reddy_szeri";
    return value;
}

double positive(const json::Value& value, const char* key) {
    const double result = value.at(key).number();
    return std::isfinite(result) && result > 0.0 ? result : 0.0;
}

std::string storedCoordinates(const json::Value& metadata,
                              const std::string& desired) {
    auto stored = metadata.at("archive_coordinate_system").string(
        metadata.at("stored_coordinate_system").string());
    if (!stored.empty()) return stored;
    const auto backend = metadata.at("execution_backend").at("name").string();
    // Historical Python solver archives explicitly preserve dimensional
    // TimeSeries arrays even when their equations were advanced in model
    // coordinates. The old Python Viewer applied this transform at display.
    if (backend == "python-reference" ||
        metadata.at("time_series_format").string() == "dynsyslab.TimeSeries/v1")
        return "dimensional";
    return desired;
}

std::vector<std::string> labels(const std::string& system,
                                std::size_t dimension) {
    std::vector<std::string> result;
    for (const auto& descriptor : dynsyslab::core::availableSystems())
        if (descriptor.identifier == system) { result = dynsyslab::core::systemStateLabels(descriptor, dimension); break; }
    while (result.size() < dimension)
        result.push_back("state " + std::to_string(result.size()));
    result.resize(dimension); return result;
}

}  // namespace

bool hasTransformedValues(const json::Value& metadata) {
    const auto space = metadata.at("value_space").string();
    return space == "observation-coordinates" || space == "history-coordinates" ||
        metadata.at("method").string() == "coordinate-transform";
}

DisplayCoordinateTransform displayCoordinateTransform(
    const json::Value& metadata, std::size_t stateDimension, const std::string& requested) {
    DisplayCoordinateTransform result;
    const auto system = identifier(metadata);
    const auto archivedForm = metadata.at("integration_form").string(
        metadata.at("coordinate_system").string("dimensional"));
    const auto desired = requested.empty() ? archivedForm : requested;
    const auto stored = storedCoordinates(metadata, archivedForm);
    result.coordinates = desired;
    result.timeLabel = desired == "dimensionless" ? "\\tau" : "t";
    result.stateMultipliers.assign(stateDimension, 1.0);
    if (hasTransformedValues(metadata)) {
        result.coordinates = stored;
        result.timeLabel = stored == "dimensionless" ? "\\tau" : "time";
        // Display labels (TeX) when the coordinate change recorded them, else identifiers.
        const auto& labels = metadata.at("coordinate_labels").array();
        const auto& names = labels.size() == stateDimension ? labels : metadata.at("coordinate_names").array();
        for (std::size_t i = 0; i < stateDimension; ++i) {
            const auto name = i < names.size() ? names[i].string() : "";
            result.stateLabels.push_back(name.empty() ? "q_{" + std::to_string(i+1) + "}" :
                scientificMathLabel(name));
        }
        return result;
    }
    result.stateLabels = labels(system, stateDimension);
    for (auto& label : result.stateLabels) label = scientificMathLabel(label);
    if (system == "rpe" || system == "marmottant") {
        const auto& legacy = metadata.at("system_parameters");
        const auto& resolved = metadata.at("resolved_system_parameters");
        double radius = positive(legacy, "equilibrium_radius_m");
        if (!radius) radius = positive(resolved, "equilibrium_radius_m");
        if (!radius) radius = positive(resolved, "R0");
        const auto& scaling = metadata.at("nondimensionalization");
        double timeScale = positive(scaling, "time_scale_s");
        if (!timeScale) timeScale = positive(metadata, "time_scale_to_seconds");
        if (!timeScale) timeScale = positive(resolved, "characteristic_time_s");
        if (stateDimension >= 2) {
            result.stateLabels[0] = desired == "dimensionless" ? "r" : "R\\; (m)";
            result.stateLabels[1] = desired == "dimensionless" ? "\\frac{dr}{d\\tau}"
                : "\\frac{dR}{dt}\\; (m/s)";
            result.timeLabel = desired == "dimensionless" ? "\\tau" : "t\\; (s)";
        }
        if (desired == stored) return result;
        if (radius && timeScale) {
            const bool dimensionless = desired == "dimensionless";
            result.timeMultiplier = dimensionless ? 1.0 / timeScale : timeScale;
            if (stateDimension >= 1) result.stateMultipliers[0] = dimensionless ? 1.0 / radius : radius;
            if (stateDimension >= 2) result.stateMultipliers[1] = dimensionless ? timeScale / radius : radius / timeScale;
            return result;
        }
    }
    if (desired == stored) return result;
    result.coordinates = stored; result.suffix.clear(); result.timeLabel = stored == "dimensionless" ? "\\tau" : "t";
    result.stateLabels = labels(system, stateDimension);
    for (auto& label : result.stateLabels) label = scientificMathLabel(label);
    if ((system == "rpe" || system == "marmottant") && stateDimension >= 2) {
        result.stateLabels[0] = stored == "dimensionless" ? "r" : "R\\; (m)";
        result.stateLabels[1] = stored == "dimensionless" ? "\\frac{dr}{d\\tau}" : "\\frac{dR}{dt}\\; (m/s)";
        result.timeLabel = stored == "dimensionless" ? "\\tau" : "t\\; (s)";
    }
    result.warning = "The requested coordinate conversion is unavailable because archived scaling metadata is incomplete.";
    return result;
}

std::vector<double> transformDisplayTime(
    std::vector<double> values, const DisplayCoordinateTransform& transform) {
    if (transform.timeMultiplier != 1.0)
        for (auto& value : values) value *= transform.timeMultiplier;
    return values;
}

dynsyslab::core::NpyArray transformDisplayStates(
    const dynsyslab::core::NpyArray& values,
    const DisplayCoordinateTransform& transform) {
    if (values.shape.size() != 2 || transform.stateMultipliers.empty()) return values;
    auto scaled = values.doubles();
    const std::size_t columns = values.shape[1];
    bool changed = false;
    for (std::size_t column = 0;
         column < std::min(columns, transform.stateMultipliers.size()); ++column) {
        const double multiplier = transform.stateMultipliers[column];
        if (multiplier == 1.0) continue;
        changed = true;
        for (std::size_t row = 0; row < values.shape[0]; ++row)
            scaled[row * columns + column] *= multiplier;
    }
    return changed ? dynsyslab::core::NpyArray::fromDoubles(values.shape, scaled) : values;
}

}  // namespace dsl
