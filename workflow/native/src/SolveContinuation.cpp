#include "SolveContinuation.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/LongPath.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <variant>

namespace dynsyslab::workflow::detail {

ContinuationPoint continuationPoint(const std::filesystem::path& archive) {
    const auto record = core::readNpz(archive);
    if (!record.contains("time") || !record.contains("values"))
        throw std::invalid_argument("continuation archive has no time and values arrays");
    const auto time = record.at("time").doubles();
    const auto& valuesArray = record.at("values");
    const auto values = valuesArray.doubles();
    if (time.empty() || values.size() % time.size())
        throw std::invalid_argument("continuation archive arrays are inconsistent");
    const std::size_t width = values.size() / time.size();
    // The last row whose state is finite (a truncated or paused parent ends early).
    std::size_t rows = time.size();
    while (rows) {
        bool finite = std::isfinite(time[rows - 1]);
        for (std::size_t column = 0; column < width && finite; ++column)
            finite = std::isfinite(values[(rows - 1) * width + column]);
        if (finite) break;
        --rows;
    }
    if (!rows) throw std::invalid_argument("continuation archive has no computed sample");
    ContinuationPoint result;
    result.time = time[rows - 1];
    result.state.assign(values.begin() + static_cast<std::ptrdiff_t>((rows - 1) * width),
                        values.begin() + static_cast<std::ptrdiff_t>(rows * width));
    result.completedSamples = rows;
    return result;
}

std::vector<std::string> continuationArchives(const json::Value& spec) {
    std::vector<std::string> paths;
    if (!spec.at("archive").string().empty()) paths.push_back(spec.at("archive").string());
    for (const auto& item : spec.at("archives").array()) paths.push_back(item.string());
    return paths;
}

json::Value continuationMetadata(const json::Value& job, std::size_t index) {
    const auto& initial = job.at("initial_conditions");
    if (initial.at("mode").string() != "continue_from") return json::Value::Object{};
    const auto paths = continuationArchives(initial);
    if (index >= paths.size()) return json::Value::Object{};
    const auto path = std::filesystem::u8path(paths[index]);
    const auto parent = continuationPoint(path);
    return json::Value::Object{
        {"parent_archive", paths[index]},
        {"parent_archive_sha256", core::sha256File(core::extendedLengthPath(path))},
        {"parent_time", parent.time},
        {"parent_completed_samples", static_cast<double>(parent.completedSamples)}};
}

void requireContinuationStart(const json::Value& job, double start) {
    const auto& initial = job.at("initial_conditions");
    if (initial.at("mode").string() != "continue_from") return;
    for (const auto& path : continuationArchives(initial)) {
        const double parent = continuationPoint(std::filesystem::u8path(path)).time;
        if (std::abs(parent - start) > 1.0e-12 * std::max(1.0, std::abs(start)))
            throw std::invalid_argument("continue_from: time.start must equal the parent archive's last time");
    }
}

json::Value extensionJob(const std::filesystem::path& parent, double stop) {
    const auto record = core::readNpz(parent);
    const auto metadata = json::parse(record.metadataJson());
    const auto point = continuationPoint(parent);
    if (!(stop > point.time)) throw std::invalid_argument("extension stop must follow the parent's last time");
    json::Value::Object system{{"name", metadata.at("system").string()},
                               {"preset", metadata.at("system_preset").string()},
                               {"parameters", metadata.at("system_parameters")}};
    auto time = metadata.at("time_configuration");
    auto& timeObject = time.object();
    timeObject.erase("output_times");
    timeObject.erase("sample_times");
    timeObject["start"] = point.time;
    timeObject["stop"] = stop;
    json::Value::Object job{
        {"schema", std::string("dynsyslab.compute-job/v1")},
        {"kind", std::string("solve")},
        {"label", metadata.at("label").string("native-solve") + "-continued"},
        {"system", json::Value(std::move(system))},
        {"time", std::move(time)},
        {"solver", metadata.at("solver_configuration")},
        {"controller", metadata.at("controller_configuration")},
        {"initial_conditions", json::Value::Object{
            {"mode", std::string("continue_from")}, {"count", 1.0},
            {"archive", parent.generic_u8string()}}}};
    return json::Value(std::move(job));
}

void writeSolveOutcome(const std::filesystem::path& manifest, const json::Value& job,
                       const core::NativeSolution& solution,
                       const std::vector<std::filesystem::path>& archives,
                       const std::vector<std::size_t>& kept, const std::vector<double>& initials) {
    const std::size_t count = solution.trajectoryCount, n = solution.stateDimension;
    if (kept.size() == count && !solution.paused) return;
    json::Value::Array mask, rejected, resume, notStarted;
    std::vector<double> waiting;
    for (std::size_t trajectory = 0; trajectory < count; ++trajectory) {
        const auto found = std::find(kept.begin(), kept.end(), trajectory);
        mask.emplace_back(found != kept.end());
        const auto& reason = solution.terminationReasons[trajectory];
        if (reason == "paused" && solution.completedSamples[trajectory] == 0) {
            notStarted.emplace_back(static_cast<double>(trajectory));
            waiting.insert(waiting.end(), initials.begin() + static_cast<std::ptrdiff_t>(trajectory * n),
                           initials.begin() + static_cast<std::ptrdiff_t>((trajectory + 1) * n));
        } else if (reason == "paused" && found != kept.end()) {
            resume.push_back(extensionJob(archives[static_cast<std::size_t>(found - kept.begin())],
                                          job.at("time").at("stop").number()));
        } else if (found == kept.end()) {
            rejected.push_back(json::Value::Object{{"trajectory_index", static_cast<double>(trajectory)},
                {"reason", reason},
                {"completed_samples", static_cast<double>(solution.completedSamples[trajectory])}});
        }
    }
    if (!waiting.empty()) {
        auto remaining = job;
        json::Value::Array rows;
        for (std::size_t row = 0; row * n < waiting.size(); ++row) {
            json::Value::Array values;
            for (std::size_t column = 0; column < n; ++column) values.emplace_back(waiting[row * n + column]);
            rows.emplace_back(std::move(values));
        }
        remaining.object()["initial_conditions"] = json::Value::Object{
            {"mode", std::string("explicit")}, {"count", static_cast<double>(rows.size())},
            {"values", json::Value(std::move(rows))}};
        remaining.object()["label"] = job.at("label").string("native-solve") + "-remaining";
        if (auto found = remaining.object().find("history");
            found != remaining.object().end() && std::holds_alternative<json::Value::Object>(found->second.storage()))
            found->second.object().erase("pause_file");
        resume.push_back(std::move(remaining));
    }
    json::Value::Object record{{"schema", std::string("dynsyslab.solve-outcome/v1")},
        {"candidate_trajectory_count", static_cast<double>(count)},
        {"accepted_trajectory_count", static_cast<double>(kept.size())},
        {"accepted_member_mask", json::Value(std::move(mask))},
        {"rejected", json::Value(std::move(rejected))}, {"paused", solution.paused},
        {"not_started", json::Value(std::move(notStarted))},
        {"resume_jobs", json::Value(std::move(resume))}};
    auto path = manifest;
    path += L".outcome.json";
    std::ofstream stream(core::extendedLengthPath(path), std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to write the solve outcome record");
    stream << json::dump(json::Value(std::move(record)));
}

}  // namespace dynsyslab::workflow::detail
