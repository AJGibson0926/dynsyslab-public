#include <dynsyslab/workflow/Solve.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/DimensionalSystem.hpp>
#include <dynsyslab/core/System.hpp>
#include "SolveContinuation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
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

// initial_conditions.clip: one [low, high] or null per state, applied after generation.
void clip(const json::Value& source, OdeInitialConditionSet& set) {
    if (source.isNull()) return;
    const auto& rows = source.array();
    if (rows.size() != set.stateDimension) throw std::invalid_argument("IC clip must cover every state");
    for (std::size_t column = 0; column < set.stateDimension; ++column) {
        if (rows[column].isNull()) continue;
        const auto bounds = numbers(rows[column]);
        if (bounds.size() != 2 || !(bounds[0] <= bounds[1])) throw std::invalid_argument("invalid IC clip range");
        for (std::size_t row = 0; row < set.trajectoryCount; ++row) {
            auto& value = set.values[row * set.stateDimension + column];
            value = std::clamp(value, bounds[0], bounds[1]);
        }
    }
}

std::size_t assignedCluster(std::size_t trajectory, std::size_t trajectories,
                            std::size_t clusters) {
    return std::min(clusters - 1, trajectory * clusters / trajectories);
}

// {"mode": "file", "path": "table.csv", "sha256": "...", "values_sha256": "...", "layout": "rows" |
//  "columns", "first": 0, "norm": 1.0, "norm_tolerance": 1e-12}: a table of initial states (one state
// per row, or per column), refused unless its bytes (sha256) or its parsed states (values_sha256)
// have the given SHA-256 (when given) and, with "norm",
// every state has that Euclidean norm. The first `count` states from `first` are used.
std::vector<double> fileTable(const json::Value& spec, std::size_t dimension, std::size_t count) {
    const auto path = std::filesystem::u8path(spec.at("path").string());
    if (path.empty() || !std::filesystem::is_regular_file(path))
        throw std::invalid_argument("IC file not found: " + path.u8string());
    const auto expected = spec.at("sha256").string();
    if (!expected.empty() && dynsyslab::core::sha256File(path) != expected)
        throw std::invalid_argument("IC file " + path.u8string() + " failed its SHA-256 check");
    std::ifstream in(path);
    std::vector<std::vector<double>> table;
    for (std::string line; std::getline(in, line);) {
        std::vector<double> row;
        std::stringstream cells(line);
        for (std::string cell; std::getline(cells, cell, ',');)
            if (cell.find_first_not_of(" \t\r") != std::string::npos) row.push_back(std::stod(cell));
        if (!row.empty()) table.push_back(std::move(row));
    }
    const bool columns = spec.at("layout").string("rows") == "columns";
    if (!columns && spec.at("layout").string("rows") != "rows") throw std::invalid_argument("IC file layout must be rows or columns");
    const std::size_t states = columns ? (table.empty() ? 0 : table.front().size()) : table.size();
    for (const auto& row : table)
        if (row.size() != (columns ? states : dimension) || (columns && table.size() != dimension))
            throw std::invalid_argument("IC file table does not hold states of the system dimension");
    // values_sha256: SHA-256 of every state in the table as little-endian float64, one state after
    // another (independent of the file's text formatting and line endings).
    if (const auto digest = spec.at("values_sha256").string(); !digest.empty()) {
        std::string bytes;
        for (std::size_t k = 0; k < states; ++k)
            for (std::size_t i = 0; i < dimension; ++i) {
                const double value = columns ? table[i][k] : table[k][i];
                unsigned char raw[sizeof(double)];
                std::memcpy(raw, &value, sizeof(double));
                bytes.append(reinterpret_cast<const char*>(raw), sizeof(double));  // the native order is little-endian here
            }
        if (dynsyslab::core::sha256Text(bytes) != digest)
            throw std::invalid_argument("IC file " + path.u8string() + " failed its values_sha256 check");
    }
    const auto first = std::size_t(spec.at("first").number(0));
    if (first + count > states) throw std::invalid_argument("IC file holds fewer states than requested");
    std::vector<double> result;
    const double norm = spec.at("norm").number(0.0), tolerance = spec.at("norm_tolerance").number(1e-12);
    for (std::size_t k = first; k < first + count; ++k) {
        double squared = 0;
        for (std::size_t i = 0; i < dimension; ++i) {
            const double value = columns ? table[i][k] : table[k][i];
            result.push_back(value);
            squared += value * value;
        }
        if (norm > 0 && std::abs(std::sqrt(squared) - norm) > tolerance)
            throw std::invalid_argument("IC file state " + std::to_string(k) + " is not on the stated norm");
    }
    return result;
}

std::vector<double> generate(const json::Value& spec,
                             const dynsyslab::core::ODESystem& system,
                             bool controlled, std::size_t count);

std::vector<std::pair<double, double>> rangePairs(const json::Value& ranges, std::size_t dimension) {
    const auto& rows = ranges.array();
    if (rows.size() != dimension) throw std::invalid_argument("IC ranges must cover every state");
    std::vector<std::pair<double, double>> result;
    for (const auto& row : rows) {
        const auto bounds = numbers(row);
        if (bounds.size() != 2 || !(bounds[0] <= bounds[1]))
            throw std::invalid_argument("invalid IC range");
        result.emplace_back(bounds[0], bounds[1]);
    }
    return result;
}

// Stratified (Latin-hypercube) sampling: every state's range is cut into `count` equal strata,
// each stratum is used once, strata are permuted independently per state and each point is
// jittered inside its stratum by U(jitter[0], jitter[1]) (McKay, Beckman & Conover,
// Technometrics 21 (1979) 239-245).
std::vector<double> stratified(const json::Value& spec, std::size_t dimension, std::size_t count) {
    const auto ranges = rangePairs(spec.at("ranges"), dimension);
    auto jitter = numbers(spec.at("jitter"));
    if (jitter.empty()) jitter = {0.1, 0.9};
    if (jitter.size() != 2 || !(0.0 <= jitter[0] && jitter[0] <= jitter[1] && jitter[1] <= 1.0))
        throw std::invalid_argument("stratified IC jitter must be [low, high] inside [0, 1]");
    std::mt19937_64 engine(static_cast<std::uint64_t>(spec.at("seed").number()));
    std::uniform_real_distribution<double> offset(jitter[0], jitter[1]);
    std::vector<double> result(count * dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
        std::vector<std::size_t> order(count);
        for (std::size_t index = 0; index < count; ++index) order[index] = index;
        std::shuffle(order.begin(), order.end(), engine);
        for (std::size_t row = 0; row < count; ++row) {
            const double fraction = (static_cast<double>(order[row]) + offset(engine)) / static_cast<double>(count);
            result[row * dimension + column] = ranges[column].first
                + (ranges[column].second - ranges[column].first) * fraction;
        }
    }
    return result;
}

// Points on a multisine target tube, x(tau) = offset + sum a_k sin(w_k tau + phi_k) in the
// first state and its derivative in the second (when present), at `count` times spread over
// tau_span (end excluded), plus Gaussian jitter of standard deviation jitter[i] per state.
std::vector<double> targetTube(const json::Value& spec, std::size_t dimension, std::size_t count) {
    const auto amplitudes = numbers(spec.at("amplitudes"));
    const auto frequencies = numbers(spec.at("angular_frequencies"));
    auto phases = numbers(spec.at("phases"));
    if (phases.empty()) phases.assign(amplitudes.size(), 0.0);
    if (amplitudes.empty() || amplitudes.size() != frequencies.size() || phases.size() != amplitudes.size())
        throw std::invalid_argument("target-tube IC needs matching amplitudes, angular_frequencies and phases");
    auto span = numbers(spec.at("tau_span"));
    if (span.empty()) span = {0.0, 1.0};
    if (span.size() != 2 || !(span[0] < span[1])) throw std::invalid_argument("target-tube tau_span must be [start, stop]");
    auto spread = numbers(spec.at("jitter"));
    spread.resize(dimension, 0.0);
    const double offset = spec.at("offset").number();
    std::mt19937_64 engine(static_cast<std::uint64_t>(spec.at("seed").number()));
    std::normal_distribution<double> normal(0.0, 1.0);
    std::vector<double> result(count * dimension, 0.0);
    for (std::size_t row = 0; row < count; ++row) {
        const double tau = span[0] + (span[1] - span[0]) * static_cast<double>(row) / static_cast<double>(count);
        double value = offset, rate = 0.0;
        for (std::size_t k = 0; k < amplitudes.size(); ++k) {
            const double angle = frequencies[k] * tau + phases[k];
            value += amplitudes[k] * std::sin(angle);
            rate += amplitudes[k] * frequencies[k] * std::cos(angle);
        }
        result[row * dimension] = value;
        if (dimension > 1) result[row * dimension + 1] = rate;
    }
    for (std::size_t column = 0; column < dimension; ++column)
        for (std::size_t row = 0; row < count; ++row)
            if (spread[column] > 0.0) result[row * dimension + column] += spread[column] * normal(engine);
    return result;
}

// The last computed state of each parent archive (one archive per trajectory).
std::vector<double> continuation(const json::Value& spec, std::size_t dimension, std::size_t count) {
    const auto paths = detail::continuationArchives(spec);
    if (paths.size() != count) throw std::invalid_argument("continue_from needs one archive per trajectory");
    std::vector<double> result;
    for (const auto& path : paths) {
        const auto parent = detail::continuationPoint(path);
        if (parent.state.size() != dimension)
            throw std::invalid_argument("continue_from archive state dimension does not match the system");
        result.insert(result.end(), parent.state.begin(), parent.state.end());
    }
    return result;
}

// {"mode": "archive_states", "archive": "run.npz", "times": [t_1, ...]}: the archived states at
// the samples nearest the given times (each within half a sample interval), one per time, as
// anchors for restarts such as pulse branches (with repeat: one anchor per branch set).
std::vector<double> archiveStates(const json::Value& spec, std::size_t dimension, std::size_t count) {
    const auto record = dynsyslab::core::readNpz(std::filesystem::u8path(spec.at("archive").string()));
    const auto time = record.at("time").doubles(), values = record.at("values").doubles();
    const auto times = numbers(spec.at("times"));
    if (time.size() < 2 || values.size() != time.size() * dimension)
        throw std::invalid_argument("archive_states archive must hold time and values of the system's states");
    if (times.size() != count) throw std::invalid_argument("archive_states needs one time per distinct initial state");
    std::vector<double> result;
    for (const double t : times) {
        const auto at = std::lower_bound(time.begin(), time.end(), t);
        std::size_t k = std::size_t(at - time.begin());
        if (k == time.size() || (k > 0 && t - time[k - 1] < time[k] - t)) --k;
        const double spacing = time[std::min(k + 1, time.size() - 1)] - time[k > 0 ? k - 1 : 0];
        if (std::abs(time[k] - t) > 0.5 * spacing) throw std::invalid_argument("archive_states time outside the archive");
        result.insert(result.end(), values.begin() + std::ptrdiff_t(k * dimension), values.begin() + std::ptrdiff_t((k + 1) * dimension));
    }
    return result;
}

std::vector<double> generate(const json::Value& spec,
                             const dynsyslab::core::ODESystem& system,
                             bool controlled, std::size_t count) {
    const std::size_t dimension = system.descriptor().stateDimension;
    const std::string mode = spec.at("mode").string("default");
    if (mode == "stratified") return stratified(spec, dimension, count);
    if (mode == "target_tube") return targetTube(spec, dimension, count);
    if (mode == "continue_from") return continuation(spec, dimension, count);
    if (mode == "archive_states") return archiveStates(spec, dimension, count);
    if (mode == "file") return fileTable(spec, dimension, count);
    if (mode == "concatenate") {
        std::vector<double> result;
        for (const auto& part : spec.at("parts").array()) {
            const double members = part.at("count").number(0.0);
            if (!(members >= 1.0) || std::floor(members) != members)
                throw std::invalid_argument("concatenate IC parts need a positive integer count");
            const auto values = generate(part, system, controlled, static_cast<std::size_t>(members));
            result.insert(result.end(), values.begin(), values.end());
        }
        if (result.size() != count * dimension)
            throw std::invalid_argument("concatenate IC part counts must sum to the trajectory count");
        return result;
    }
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
    return {"default", "explicit", "uniform", "normal", "sphere", "ellipsoid",
            "stratified", "target_tube", "concatenate", "continue_from", "archive_states", "file"};
}

OdeInitialConditionSet resolveOdeInitialConditions(const json::Value& job) {
    // Default ICs follow the integration coordinates (a dimensional solve scales them).
    auto system = dynsyslab::core::systemInCoordinates(dynsyslab::core::createSystem(configuration(job)),
        job.at("time").at("coordinate_system").string("automatic"));
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
        const double repeat = root.at("repeat").number(1.0);
        if (!(repeat >= 1.0) || std::floor(repeat) != repeat
            || result.trajectoryCount % static_cast<std::size_t>(repeat))
            throw std::invalid_argument("initial_conditions.repeat must divide the trajectory count");
        const auto copies = static_cast<std::size_t>(repeat);
        const auto distinct = generate(root, *system, controlled, result.trajectoryCount / copies);
        // repeat: every generated state is used by `repeat` consecutive trajectories
        // (for example the -/+ members of a symmetric pulse pair).
        for (std::size_t row = 0; row * result.stateDimension < distinct.size(); ++row)
            for (std::size_t copy = 0; copy < copies; ++copy)
                result.values.insert(result.values.end(),
                    distinct.begin() + static_cast<std::ptrdiff_t>(row * result.stateDimension),
                    distinct.begin() + static_cast<std::ptrdiff_t>((row + 1) * result.stateDimension));
        clip(root.at("clip"), result);
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
