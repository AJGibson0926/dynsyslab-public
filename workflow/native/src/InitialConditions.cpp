#include <dynsyslab/workflow/Solve.hpp>

#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

dynsyslab::core::SystemConfiguration configuration(const json::Value& job) {
    const auto& source = job.at("system");
    dynsyslab::core::SystemConfiguration result;
    result.identifier = source.at("name").string();
    result.preset = source.at("preset").string();
    for (const auto& parameter : source.at("parameters").object())
        result.parameters.emplace(parameter.first, parameter.second.number());
    return result;
}

std::size_t assignedCluster(std::size_t trajectory, std::size_t trajectories,
                            std::size_t clusters) {
    return std::min(clusters - 1, trajectory * clusters / trajectories);
}

std::vector<double> generate(const json::Value& spec,
                             const dynsyslab::core::ODESystem& system,
                             bool controlled, std::size_t count) {
    const std::size_t dimension = system.descriptor().stateDimension;
    const std::string mode = spec.at("mode").string("default");
    if (mode == "default") {
        const auto seed = controlled ? system.controlledInitialState() : system.initialState();
        std::vector<double> result; result.reserve(count * dimension);
        for (std::size_t index = 0; index < count; ++index)
            result.insert(result.end(), seed.begin(), seed.end());
        return result;
    }
    if (mode == "explicit") {
        const auto& rows = spec.at("values").array();
        std::vector<double> result;
        if (rows.size() == dimension && !rows.empty() && rows.front().array().empty())
            result = numbers(spec.at("values"));
        else for (const auto& row : rows) {
            const auto values = numbers(row);
            if (values.size() != dimension)
                throw std::invalid_argument("explicit IC dimension mismatch");
            result.insert(result.end(), values.begin(), values.end());
        }
        if (result.size() != count * dimension)
            throw std::invalid_argument("explicit IC count mismatch");
        return result;
    }
    if (mode == "sphere" || mode == "ellipsoid") {
        const auto& distribution = spec.at("distribution");
        auto center = numbers(distribution.at("center"));
        if (center.empty()) center.assign(dimension, 0.0);
        if (center.size() != dimension)
            throw std::invalid_argument("IC distribution center must cover every state");
        std::vector<double> axes;
        if (mode == "sphere") {
            const double radius = distribution.at("radius").number(1.0);
            if (!(radius > 0.0))
                throw std::invalid_argument("sphere IC radius must be positive");
            axes.assign(dimension, radius);
        } else {
            axes = numbers(distribution.at("axes"));
            if (axes.size() != dimension || std::any_of(axes.begin(), axes.end(),
                    [](double value) { return !(value > 0.0); }))
                throw std::invalid_argument(
                    "ellipsoid IC axes must be positive and cover every state");
        }
        const double inner = distribution.at("inner_radius").number(1.0);
        const double outer = distribution.at("outer_radius").number(inner);
        if (inner < 0.0 || outer < inner || !(outer > 0.0))
            throw std::invalid_argument(
                "IC shell radii require 0 <= inner <= outer and outer > 0");
        std::mt19937_64 engine(static_cast<std::uint64_t>(spec.at("seed").number()));
        std::normal_distribution<double> direction(0.0, 1.0);
        std::uniform_real_distribution<double> shell(0.0, 1.0);
        std::vector<double> result(count * dimension), vector(dimension);
        for (std::size_t row = 0; row < count; ++row) {
            double norm = 0.0;
            do {
                norm = 0.0;
                for (double& value : vector) {
                    value = direction(engine); norm += value * value;
                }
            } while (!(norm > 0.0));
            norm = std::sqrt(norm);
            const double radial = inner == outer ? inner : std::pow(
                std::pow(inner, static_cast<double>(dimension))
                + shell(engine) * (std::pow(outer, static_cast<double>(dimension))
                - std::pow(inner, static_cast<double>(dimension))),
                1.0 / static_cast<double>(dimension));
            for (std::size_t column = 0; column < dimension; ++column)
                result[row * dimension + column] = center[column]
                    + axes[column] * radial * vector[column] / norm;
        }
        return result;
    }
    if (mode != "uniform" && mode != "normal")
        throw std::invalid_argument("unknown IC mode");
    const auto& ranges = spec.at("ranges").array();
    if (ranges.size() != dimension)
        throw std::invalid_argument("IC ranges must cover every state");
    std::mt19937_64 engine(static_cast<std::uint64_t>(spec.at("seed").number()));
    std::vector<double> result(count * dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
        const auto bounds = numbers(ranges[column]);
        if (bounds.size() != 2 || bounds[0] > bounds[1])
            throw std::invalid_argument("invalid IC range");
        std::uniform_real_distribution<double> uniform(bounds[0], bounds[1]);
        std::normal_distribution<double> normal(0.5 * (bounds[0] + bounds[1]),
                                                (bounds[1] - bounds[0]) / 6.0);
        for (std::size_t row = 0; row < count; ++row) {
            const double value = mode == "uniform" ? uniform(engine) : normal(engine);
            result[row * dimension + column] = std::clamp(value, bounds[0], bounds[1]);
        }
    }
    return result;
}

}  // namespace

std::vector<std::string> supportedOdeInitialConditionModes() {
    return {"default", "explicit", "uniform", "normal", "sphere", "ellipsoid"};
}

OdeInitialConditionSet resolveOdeInitialConditions(const json::Value& job) {
    auto system = dynsyslab::core::createSystem(configuration(job));
    const bool controlled = job.at("controller").at("kind").string("none") != "none";
    const auto& root = job.at("initial_conditions");
    OdeInitialConditionSet result;
    result.trajectoryCount = static_cast<std::size_t>(
        std::max(1.0, root.at("count").number(1.0)));
    result.stateDimension = system->descriptor().stateDimension;
    result.stateLabels = system->descriptor().stateLabels;
    if (result.stateLabels.size() != result.stateDimension)
        result.stateLabels = system->descriptor().stateNames;
    const std::size_t clusters = std::clamp<std::size_t>(
        static_cast<std::size_t>(std::max(1.0, root.at("clusters").number(1.0))),
        1, result.trajectoryCount);
    const auto& clusterSpecs = root.at("cluster_specs").array();
    if (clusterSpecs.empty()) {
        result.values = generate(root, *system, controlled, result.trajectoryCount);
        return result;
    }
    if (clusterSpecs.size() != clusters)
        throw std::invalid_argument("IC cluster_specs must contain one specification per cluster");
    // A null entry inherits the system-wide generator. This lets one Cluster
    // become specialized without perturbing every other proposed IC.
    result.values = generate(root, *system, controlled, result.trajectoryCount);
    for (std::size_t cluster = 0; cluster < clusters; ++cluster) {
        if (clusterSpecs[cluster].isNull()) continue;
        std::size_t members = 0;
        for (std::size_t trajectory = 0; trajectory < result.trajectoryCount; ++trajectory)
            if (assignedCluster(trajectory, result.trajectoryCount, clusters) == cluster)
                ++members;
        const auto generated = generate(clusterSpecs[cluster], *system, controlled, members);
        std::size_t source = 0;
        for (std::size_t trajectory = 0; trajectory < result.trajectoryCount; ++trajectory) {
            if (assignedCluster(trajectory, result.trajectoryCount, clusters) != cluster) continue;
            std::copy_n(generated.begin() + static_cast<std::ptrdiff_t>(
                source * result.stateDimension), result.stateDimension,
                result.values.begin() + static_cast<std::ptrdiff_t>(
                    trajectory * result.stateDimension));
            ++source;
        }
    }
    return result;
}

json::Value odeInitialConditionSetJson(const OdeInitialConditionSet& set) {
    json::Value::Object root;
    root["schema"] = "dynsyslab.initial-conditions/v1";
    root["trajectory_count"] = static_cast<double>(set.trajectoryCount);
    root["state_dimension"] = static_cast<double>(set.stateDimension);
    json::Value::Array labels;
    for (const auto& label : set.stateLabels) labels.emplace_back(label);
    root["state_labels"] = json::Value(std::move(labels));
    json::Value::Array rows;
    for (std::size_t row = 0; row < set.trajectoryCount; ++row) {
        json::Value::Array values;
        for (std::size_t column = 0; column < set.stateDimension; ++column)
            values.emplace_back(set.values[row * set.stateDimension + column]);
        rows.emplace_back(std::move(values));
    }
    root["values"] = json::Value(std::move(rows));
    return json::Value(std::move(root));
}

}  // namespace dynsyslab::workflow
