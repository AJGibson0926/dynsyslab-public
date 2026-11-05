#include <dynsyslab/workflow/ControlTargets.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Targets.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::workflow {
namespace {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

void finiteVector(const std::vector<double>& values, const char* name) {
    if (std::any_of(values.begin(), values.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument(std::string(name) + " must be finite");
}

dynsyslab::core::StateTargetFunction fixedTarget(
    const json::Value& source, std::size_t dimension) {
    auto state = numbers(source.at("state"));
    if (state.size() != dimension)
        throw std::invalid_argument("control target dimension mismatch");
    finiteVector(state, "control target");
    return [state = std::move(state)](double) { return state; };
}

dynsyslab::core::StateTargetFunction waveformTarget(
    const json::Value& source, std::size_t dimension) {
    auto offset = numbers(source.at("offset"));
    if (offset.empty()) offset = numbers(source.at("state"));
    if (offset.size() != dimension)
        throw std::invalid_argument("waveform target offset dimension mismatch");
    auto amplitudes = numbers(source.at("amplitudes"));
    auto frequencies = numbers(source.at("frequencies"));
    auto angular = numbers(source.at("angular_frequencies"));
    if (!frequencies.empty() && !angular.empty())
        throw std::invalid_argument(
            "waveform target cannot mix cyclic and angular frequencies");
    constexpr double twoPi = 6.283185307179586476925286766559;
    if (angular.empty()) {
        angular = std::move(frequencies);
        for (double& value : angular) value *= twoPi;
    }
    auto phases = numbers(source.at("phases_rad"));
    if (phases.empty()) phases.assign(amplitudes.size(), 0.0);
    if (amplitudes.empty() || amplitudes.size() != angular.size()
        || amplitudes.size() != phases.size())
        throw std::invalid_argument(
            "waveform target amplitudes, frequencies, and phases must match");
    finiteVector(offset, "waveform target offset");
    finiteVector(amplitudes, "waveform target amplitudes");
    finiteVector(angular, "waveform target frequencies");
    finiteVector(phases, "waveform target phases");
    return [offset = std::move(offset), amplitudes = std::move(amplitudes),
            angular = std::move(angular), phases = std::move(phases), dimension]
        (double time) {
        auto result = offset;
        for (std::size_t index = 0; index < amplitudes.size(); ++index) {
            const double phase = angular[index] * time + phases[index];
            result[0] += amplitudes[index] * std::sin(phase);
            if (dimension > 1)
                result[1] += amplitudes[index] * angular[index] * std::cos(phase);
        }
        return result;
    };
}

dynsyslab::core::StateTargetFunction trajectoryTarget(
    const json::Value& source, std::size_t dimension) {
    const std::string path = source.at("source").string();
    if (path.empty()) throw std::invalid_argument("trajectory target source is empty");
    const auto archive = dynsyslab::core::readNpz(std::filesystem::u8path(path));
    const std::string timeName = source.at("time_array").string("time");
    std::string valueName = source.at("value_array").string();
    if (valueName.empty()) valueName = archive.contains("control_target_values")
        ? "control_target_values" : "values";
    if (!archive.contains(timeName) || !archive.contains(valueName))
        throw std::invalid_argument("trajectory target archive lacks time or state values");
    auto time = archive.at(timeName).doubles();
    auto values = archive.at(valueName).doubles();
    const auto& shape = archive.at(valueName).shape;
    if (time.empty() || shape.size() != 2 || shape[0] != time.size()
        || shape[1] != dimension || values.size() != time.size() * dimension)
        throw std::invalid_argument("trajectory target archive shape is incompatible");
    finiteVector(time, "trajectory target time");
    finiteVector(values, "trajectory target values");
    if (!std::is_sorted(time.begin(), time.end())
        || std::adjacent_find(time.begin(), time.end()) != time.end())
        throw std::invalid_argument("trajectory target time must be strictly increasing");
    return [time = std::move(time), values = std::move(values), dimension]
        (double query) {
        const auto copyRow = [&](std::size_t row) {
            return std::vector<double>(values.begin() + row * dimension,
                values.begin() + (row + 1) * dimension);
        };
        if (query <= time.front()) return copyRow(0);
        if (query >= time.back()) return copyRow(time.size() - 1);
        const auto upper = std::upper_bound(time.begin(), time.end(), query);
        const std::size_t right = static_cast<std::size_t>(upper - time.begin());
        const std::size_t left = right - 1;
        const double weight = (query - time[left]) / (time[right] - time[left]);
        std::vector<double> result(dimension);
        for (std::size_t state = 0; state < dimension; ++state)
            result[state] = values[left * dimension + state] * (1.0 - weight)
                + values[right * dimension + state] * weight;
        return result;
    };
}

}  // namespace

dynsyslab::core::StateTargetFunction resolveControlTarget(
    const json::Value& source, std::size_t dimension) {
    const std::string kind = source.at("kind").string("constant");
    if (kind == "constant" || kind == "fixed_point")
        return fixedTarget(source, dimension);
    if (kind == "waveform") return waveformTarget(source, dimension);
    if (kind == "trajectory") return trajectoryTarget(source, dimension);
    if (kind == "rpe_canonical") {
        if (dimension != 2)
            throw std::invalid_argument("RPE target requires a two-state system");
        const auto target = dynsyslab::core::rpeControlTarget(
            source.at("name").string());
        return [target](double time) { return target.evaluate(time); };
    }
    throw std::invalid_argument("unknown native control target kind: " + kind);
}

}  // namespace dynsyslab::workflow
