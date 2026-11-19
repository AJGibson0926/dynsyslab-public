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

// Optional smooth onset for any state target (the ramp toggle). A ramp block
//   {"duration": T, "from": [...], "derivative_pairs": [[position, velocity], ...]}
// gives x_ref = from + s(t)(w(t) - from), with the quintic smoothstep s rising
// 0 -> 1 over [0, T]. For each derivative pair, the velocity component also gets
// s'(t)(w_position - from_position), so it stays the exact derivative of the
// ramped position. The default pairs are [[0, 1]] for a state of dimension 2 or
// more (r and r' for the RPE), and none otherwise. "ramp" absent, null or false,
// or {"enabled": false}, disables the ramp: the target is used as given.
bool rampRequested(const json::Value& ramp) {
    if (ramp.isNull()) return false;
    if (std::holds_alternative<bool>(ramp.storage())) return ramp.boolean();
    return ramp.at("enabled").boolean(true);
}

dynsyslab::core::StateTargetFunction withRamp(
    dynsyslab::core::StateTargetFunction base, const json::Value& ramp, std::size_t dimension) {
    if (!rampRequested(ramp)) return base;
    if (!std::holds_alternative<json::Value::Object>(ramp.storage()))
        throw std::invalid_argument("target ramp must be an object, null or false");
    const double duration = ramp.at("duration").number(0.0);
    auto from = numbers(ramp.at("from"));
    if (!(duration > 0.0) || !std::isfinite(duration))
        throw std::invalid_argument("target ramp.duration must be positive");
    if (from.size() != dimension)
        throw std::invalid_argument("target ramp.from dimension mismatch");
    finiteVector(from, "target ramp.from");
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    if (ramp.at("derivative_pairs").isNull()) {
        if (dimension > 1) pairs.emplace_back(0, 1);
    } else {
        for (const auto& pair : ramp.at("derivative_pairs").array()) {
            const auto ij = numbers(pair);
            if (ij.size() != 2 || ij[0] < 0 || ij[1] < 0 || ij[0] != std::floor(ij[0])
                || ij[1] != std::floor(ij[1]) || ij[0] >= double(dimension) || ij[1] >= double(dimension))
                throw std::invalid_argument("target ramp.derivative_pairs must hold [position, velocity] state indices");
            pairs.emplace_back(std::size_t(ij[0]), std::size_t(ij[1]));
        }
    }
    return [base = std::move(base), duration, from = std::move(from), pairs = std::move(pairs)]
        (double time) {
        auto result = base(time);
        const double x = std::clamp(time / duration, 0.0, 1.0);
        const double s = x * x * x * (10.0 - 15.0 * x + 6.0 * x * x);
        const double ds = (time > 0.0 && time < duration)
            ? 30.0 * x * x * (1.0 - x) * (1.0 - x) / duration : 0.0;
        std::vector<double> lead(pairs.size());
        for (std::size_t p = 0; p < pairs.size(); ++p) lead[p] = result[pairs[p].first] - from[pairs[p].first];
        for (std::size_t k = 0; k < result.size(); ++k) result[k] = from[k] + s * (result[k] - from[k]);
        for (std::size_t p = 0; p < pairs.size(); ++p) result[pairs[p].second] += ds * lead[p];
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

namespace {
dynsyslab::core::StateTargetFunction baseTarget(const json::Value& source, std::size_t dimension) {
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
}  // namespace

dynsyslab::core::StateTargetFunction resolveControlTarget(
    const json::Value& source, std::size_t dimension) {
    return withRamp(baseTarget(source, dimension), source.at("ramp"), dimension);
}

}  // namespace dynsyslab::workflow
