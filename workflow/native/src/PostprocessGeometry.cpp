#include "PostprocessGeometry.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <set>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {

struct SourceShape {
    std::string path;
    bool resolved{};
    std::size_t samples{};
    std::size_t dimension{};
    bool uniform{true};
};

// Methods that factor a Hankel (delay-embedded) matrix of one source (Arbabi & Mezic 2017,
// "Ergodic theory, dynamic mode decomposition, and computation of spectral properties of the
// Koopman operator", SIAM J. Appl. Dyn. Syst. 16:2096-2126).
bool hankelMethod(const std::string& name) {
    static const std::set<std::string> names{"hdmd", "hankel_dmd", "hankel_svd", "hankel_pod",
        "hankel_matrix", "havok", "koopman_spectral", "resdmd", "mpedmd", "hermitian_dmd",
        "tlsdmd", "multdmd", "stochastic_resdmd", "specrkhs"};
    return names.count(name) != 0;
}

std::size_t integer(const json::Value& value, std::size_t fallback) {
    const double number = value.number(static_cast<double>(fallback));
    return number >= 0.0 ? static_cast<std::size_t>(number) : fallback;
}

SourceShape inspect(const json::Value& source) {
    SourceShape shape;
    shape.path = source.string(source.at("archive").string());
    const std::filesystem::path path = dynsyslab::core::wide(shape.path);
    if (shape.path.empty() || shape.path.find("${") != std::string::npos
        || !std::filesystem::is_regular_file(path) || path.extension() != L".npz")
        return shape;
    // Planning stays cheap: archives above 256 MiB are not opened (listed as unresolved).
    if (std::filesystem::file_size(path) > (std::uintmax_t(256) << 20)) return shape;
    const auto archive = dynsyslab::core::readNpz(path);
    const char* timeNames[] = {"time", "controlled_time", "control_time", "free_time"};
    const char* valueNames[] = {"values", "states", "controlled_states", "free_training_values"};
    for (const char* name : timeNames)
        if (archive.contains(name)) {
            const auto time = archive.at(name).doubles();
            shape.samples = time.size();
            if (time.size() > 2) {
                const double step = time[1] - time[0];
                for (std::size_t index = 2; index < time.size(); ++index)
                    if (std::abs((time[index] - time[index - 1]) - step)
                        > 1e-9 * std::max(1.0, std::abs(step))) { shape.uniform = false; break; }
            }
            break;
        }
    for (const char* name : valueNames)
        if (archive.contains(name)) {
            const auto& dims = archive.at(name).shape;
            shape.dimension = dims.size() > 1 ? dims.back() : 1;
            break;
        }
    shape.resolved = shape.samples > 0;
    return shape;
}

double binomial(std::size_t n, std::size_t k) {
    double result = 1.0;
    for (std::size_t index = 1; index <= k; ++index)
        result = result * static_cast<double>(n - k + index) / static_cast<double>(index);
    return result;
}

}  // namespace

json::Value postprocessGeometry(const json::Value& job) {
    json::Value::Array sources, stages, issues;
    std::vector<SourceShape> shapes;
    for (const auto& source : job.at("sources").array()) {
        auto shape = inspect(source);
        sources.emplace_back(json::Value::Object{{"archive", shape.path}, {"resolved", shape.resolved},
            {"samples", static_cast<double>(shape.samples)},
            {"state_dimension", static_cast<double>(shape.dimension)},
            {"uniform_time_steps", shape.uniform}});
        if (shape.resolved && !shape.uniform)
            issues.emplace_back("source " + shape.path + " is not uniformly sampled; resample it first");
        shapes.push_back(std::move(shape));
    }
    double peak = 0.0;
    for (const auto& stage : job.at("methods").array()) {
        if (!stage.at("enabled").boolean(true)) continue;
        const auto name = stage.at("name").string();
        const auto& settings = stage.at("settings");
        json::Value::Object item{{"method", name}};
        const SourceShape* shape = nullptr;
        for (const auto& candidate : shapes) if (candidate.resolved) { shape = &candidate; break; }
        if (hankelMethod(name) && shape) {
            const auto& window = settings.at("window");
            const bool modern = !window.isNull();
            const std::size_t stride = integer(modern ? window.at("sample_stride") : settings.at("sample_stride"), 1);
            const std::size_t samples = (shape->samples - 1) / std::max<std::size_t>(1, stride) + 1;
            const std::size_t delays = modern ? integer(window.at("delays"), 0)
                : integer(settings.at("sampling_window"), integer(settings.at("matrix_size"), 10));
            const std::size_t delayStride = integer(modern ? window.at("delay_stride") : settings.at("delay_stride"), 1);
            const std::size_t span = delays ? (delays - 1) * std::max<std::size_t>(1, delayStride) + 1 : 0;
            const std::size_t available = samples > span ? samples - span + 1 : 0;
            std::size_t columns = modern ? integer(window.at("columns"), 0)
                : (settings.at("sampling_window").isNull() ? 0 : integer(settings.at("matrix_size"), 0));
            columns = columns ? std::min(columns, available) : available;
            const std::size_t observables = std::max<std::size_t>(1, settings.at("observables").array().size());
            const std::size_t rows = delays * observables;
            const std::size_t rankLimit = std::min(rows, columns);
            const std::size_t requestedRank = integer(settings.at("rank"), 0);
            // Dense SVD working set: the matrix, its factors and LAPACK workspace (factor 4).
            const double bytes = 4.0 * static_cast<double>(rows) * columns * sizeof(double);
            peak = std::max(peak, bytes);
            item["hankel_shape"] = json::Value(json::Value::Array{
                json::Value(static_cast<double>(rows)), json::Value(static_cast<double>(columns))});
            item["retained_rank_upper_bound"] = static_cast<double>(
                requestedRank ? std::min(requestedRank, rankLimit) : rankLimit);
            item["estimated_dense_working_set_bytes"] = bytes;
            if (columns < 2) issues.emplace_back(name + ": the source is too short for two Hankel columns");
        }
        if ((name == "sindy" || name == "sindy_dynamics") && shape) {
            const std::size_t degree = integer(settings.at("polynomial_degree"), 4);
            const double features = binomial(shape->dimension + degree, degree);
            const double bytes = 2.0 * static_cast<double>(shape->samples) * features * sizeof(double);
            peak = std::max(peak, bytes);
            item["polynomial_degree"] = static_cast<double>(degree);
            item["feature_count"] = features;
            item["estimated_dense_working_set_bytes"] = bytes;
        }
        stages.emplace_back(std::move(item));
    }
    return json::Value::Object{{"schema", "dynsyslab.postprocess-geometry/v1"},
        {"valid", issues.empty()}, {"issues", json::Value(std::move(issues))},
        {"sources", json::Value(std::move(sources))}, {"stages", json::Value(std::move(stages))},
        {"estimated_peak_bytes", peak}};
}

void requirePostprocessMemory(const json::Value& job, const dynsyslab::core::ExecutionLimits& limits) {
    double budget = limits.maximumMemoryBytes ? static_cast<double>(limits.maximumMemoryBytes) : 0.0;
    const auto take = [&budget](double bytes) {
        if (bytes > 0.0) budget = budget > 0.0 ? std::min(budget, bytes) : bytes;
    };
    take(job.at("execution").at("max_memory_mib").number(0.0) * 1024.0 * 1024.0);
    for (const auto& stage : job.at("methods").array())
        take(stage.at("settings").at("max_memory_mb").number(0.0) * 1024.0 * 1024.0);
    if (budget <= 0.0) return;
    const auto geometry = postprocessGeometry(job);
    const double peak = geometry.at("estimated_peak_bytes").number();
    if (peak > budget)
        throw dynsyslab::core::ExecutionLimitExceeded(
            "post-processing preflight: estimated dense working set "
            + std::to_string(static_cast<long long>(peak / (1024.0 * 1024.0)))
            + " MiB exceeds the memory budget of "
            + std::to_string(static_cast<long long>(budget / (1024.0 * 1024.0))) + " MiB");
}

}  // namespace dynsyslab::workflow::detail
