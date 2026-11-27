// Controller call log and partial history on failure for ODE solves (SolveHistory.hpp).
#include "SolveHistory.hpp"

#include "SolveArchive.hpp"

#include <stdexcept>

namespace dynsyslab::workflow::detail {

void applyHistoryOptions(const json::Value& job, core::IntegrationOptions& options,
                         core::ControllerCallLog& log, std::size_t stateDimension,
                         std::size_t inputDimension) {
    const auto& history = job.at("history");
    options.keepPartialHistoryOnFailure = history.at("keep_partial_history_on_failure").boolean(false);
    if (!history.at("controller_call_log").boolean(false)) return;
    const double limit = history.at("controller_call_log_limit").number(0.0);
    if (!(limit >= 0.0)) throw std::invalid_argument("history.controller_call_log_limit must be nonnegative");
    log.stateDimension = stateDimension;
    log.inputDimension = inputDimension;
    log.maximumEntries = static_cast<std::size_t>(limit);
    options.controllerCallLog = &log;
}

void appendControllerCalls(core::NpzRecord& archive, json::Value::Object& metadata,
                           const core::IntegrationOptions& options, std::size_t trajectory) {
    const auto* log = options.controllerCallLog;
    if (!log) return;
    const auto entries = log->entriesOf(trajectory);
    const std::size_t count = entries.size(), n = log->stateDimension, m = log->inputDimension;
    std::vector<double> time, states, inputs, kind, status;
    for (const auto index : entries) {
        time.push_back(log->time[index]);
        states.insert(states.end(), log->states.begin() + index * n, log->states.begin() + (index + 1) * n);
        inputs.insert(inputs.end(), log->inputs.begin() + index * m, log->inputs.begin() + (index + 1) * m);
        kind.push_back(log->kind[index]);
        status.push_back(log->status[index]);
    }
    archive.arrays["controller_call_time"] = core::NpyArray::fromDoubles({count}, time);
    archive.arrays["controller_call_state"] = core::NpyArray::fromDoubles({count, n}, states);
    archive.arrays["controller_call_input"] = core::NpyArray::fromDoubles({count, m}, inputs);
    archive.arrays["controller_call_kind"] = core::NpyArray::fromDoubles({count}, kind);
    archive.arrays["controller_call_status"] = core::NpyArray::fromDoubles({count}, status);
    metadata["controller_call_log"] = json::Value::Object{
        {"entries", static_cast<double>(count)},
        {"dropped_entries_all_trajectories", static_cast<double>(log->droppedEntries)},
        {"kind_codes", json::Value::Object{{"right_hand_side", 0.0}, {"sample_tick", 1.0}, {"output_sample", 2.0}}},
        {"status_codes", json::Value::Object{{"accepted", 1.0}, {"rejected", 0.0}, {"unknown", -1.0}}}};
}

core::NativeSolution integrateKeepingHistory(
    const std::string& label, const std::function<std::filesystem::path(std::size_t)>& archivePath,
    const std::function<core::NativeSolution()>& run) {
    try {
        return run();
    } catch (const core::IntegrationFailure& failure) {
        const auto& partial = failure.partial();
        std::string kept;
        for (std::size_t trajectory = 0; trajectory < partial.trajectoryCount; ++trajectory) {
            const std::size_t rows = partial.completedSamples[trajectory];
            if (!rows) continue;
            core::IntegrationCheckpoint point;
            point.trajectory = trajectory;
            point.sample = rows - 1;
            point.time = partial.time[rows - 1];
            point.stateDimension = partial.stateDimension;
            point.times = partial.time.data();
            point.values = partial.values.data() + trajectory * partial.sampleCount * partial.stateDimension;
            point.state = point.values + (rows - 1) * partial.stateDimension;
            point.inputDimension = partial.inputDimension;
            point.controls = partial.controlValues.empty() ? nullptr
                : partial.controlValues.data() + trajectory * partial.sampleCount * partial.inputDimension;
            json::Value::Object meta{{"schema", std::string("dynsyslab.solution-checkpoint/v1")},
                {"source_job_label", label}, {"trajectory_index", static_cast<double>(trajectory)},
                {"time", point.time}, {"completed_samples", static_cast<double>(rows)},
                {"kept_on_failure", true}, {"termination_reason", partial.terminationReasons[trajectory]}};
            const auto path = partialArchivePath(archivePath(trajectory));
            writePartialArchive(path, point, json::dump(json::Value(std::move(meta))));
            kept += (kept.empty() ? "" : ", ") + path.generic_u8string();
        }
        throw std::runtime_error(std::string(failure.what())
            + (kept.empty() ? "" : "; partial history kept in " + kept));
    }
}

}  // namespace dynsyslab::workflow::detail
