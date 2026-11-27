#include "SolveSettings.hpp"

#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {

std::vector<double> numbers(const json::Value& source, const char* name) {
    std::vector<double> result;
    for (const auto& item : source.array()) {
        const double value = item.number(std::numeric_limits<double>::quiet_NaN());
        if (!std::isfinite(value)) throw std::invalid_argument(std::string(name) + " must be finite");
        result.push_back(value);
    }
    return result;
}

void requireIncreasing(const std::vector<double>& values, const char* name) {
    for (std::size_t index = 1; index < values.size(); ++index)
        if (!(values[index] > values[index - 1]))
            throw std::invalid_argument(std::string(name) + " must be strictly increasing");
}

// solver.state_bounds: one [low, high] (or null for unbounded) per state component.
void stateBounds(const json::Value& source, core::IntegrationOptions& options) {
    if (source.isNull()) return;
    const auto& rows = source.array();
    const double inf = (std::numeric_limits<double>::infinity)();
    for (const auto& row : rows) {
        double low = -inf, high = inf;
        if (!row.isNull()) {
            const auto& pair = row.array();
            if (pair.size() != 2) throw std::invalid_argument("solver.state_bounds entries are [low, high] or null");
            if (!pair[0].isNull()) low = pair[0].number();
            if (!pair[1].isNull()) high = pair[1].number();
            if (!(low <= high)) throw std::invalid_argument("solver.state_bounds requires low <= high");
        }
        options.stateLowerBounds.push_back(low);
        options.stateUpperBounds.push_back(high);
    }
}

}  // namespace

const char* failurePolicyName(FailurePolicy policy) {
    switch (policy) {
    case FailurePolicy::error: return "error";
    case FailurePolicy::truncate: return "truncate";
    case FailurePolicy::reject: return "reject";
    }
    return "error";
}

SolveSettings solveSettings(const json::Value& job) {
    const auto& time = job.at("time");
    const auto& solver = job.at("solver");
    SolveSettings result;
    auto& options = result.integration;
    options.outputTimes = numbers(time.at("output_times"), "time.output_times");
    if (!options.outputTimes.empty()) {
        requireIncreasing(options.outputTimes, "time.output_times");
        if (options.outputTimes.size() < 2) throw std::invalid_argument("time.output_times needs two times");
        options.initialTime = time.at("start").number(options.outputTimes.front());
        options.finalTime = time.at("stop").number(options.outputTimes.back());
        if (options.initialTime != options.outputTimes.front()
            || options.finalTime != options.outputTimes.back())
            throw std::invalid_argument("time.output_times must start at time.start and end at time.stop");
        options.outputStep = time.at("output_step").number(
            (options.finalTime - options.initialTime) / static_cast<double>(options.outputTimes.size() - 1));
    } else {
        options.initialTime = time.at("start").number();
        options.finalTime = time.at("stop").number();
        options.outputStep = time.at("output_step").number();
    }
    if (!(options.initialTime < options.finalTime) || !(options.outputStep > 0.0))
        throw std::invalid_argument("solve time requires start < stop and positive output step");
    options.relativeTolerance = solver.at("rtol").number(1.0e-8);
    options.absoluteTolerance = solver.at("atol").number(1.0e-10);
    options.maximumStep = solver.at("max_step").number(0.0);
    options.initialStep = solver.at("step_size").number(0.0);
    const double budget = solver.at("max_internal_steps").number(10000000.0);
    if (budget < 1.0 || std::floor(budget) != budget || budget > 9.0e15)
        throw std::invalid_argument("max_internal_steps must be a positive integer");
    options.maximumInternalSteps = static_cast<std::size_t>(budget);
    options.method = core::integrationMethodFromName(solver.at("method").string("RK45"));
    stateBounds(solver.at("state_bounds"), options);
    options.maximumStateMagnitude = solver.at("max_state_magnitude").number(0.0);
    if (!(options.maximumStateMagnitude >= 0.0))
        throw std::invalid_argument("solver.max_state_magnitude must be nonnegative");

    // "error" aborts the solve; "truncate" archives each trajectory up to its last computed
    // sample with integration_terminated; "reject" drops failed trajectories (quality gate)
    // and records them in the hierarchy manifest's rejected list.
    const auto failure = solver.at("on_failure").string("error");
    if (failure == "error") result.failure = FailurePolicy::error;
    else if (failure == "truncate") result.failure = FailurePolicy::truncate;
    else if (failure == "reject") result.failure = FailurePolicy::reject;
    else throw std::invalid_argument("solver.on_failure must be error, truncate or reject");
    options.truncateOnFailure = result.failure != FailurePolicy::error;

    // Post-solve resampling: time.sample_times (shorthand) or interpolation {times|time_step}.
    const auto& interpolation = job.at("interpolation");
    auto sampleTimes = numbers(time.at("sample_times"), "time.sample_times");
    auto explicitTimes = numbers(interpolation.at("times"), "interpolation.times");
    if (!sampleTimes.empty() && !explicitTimes.empty())
        throw std::invalid_argument("give resampling times through time.sample_times or interpolation.times, not both");
    if (explicitTimes.empty()) explicitTimes = std::move(sampleTimes);
    const double step = interpolation.at("time_step").number(0.0);
    if (!explicitTimes.empty()) {
        result.resample.times = std::move(explicitTimes);
    } else if (interpolation.at("enabled").boolean(false) || step > 0.0) {
        if (!(step > 0.0)) throw std::invalid_argument("interpolation needs times or a positive time_step");
        const auto count = static_cast<std::size_t>(std::floor(
            (options.finalTime - options.initialTime) / step + 0.5 + 1.0e-9));
        for (std::size_t index = 0; index <= count; ++index)
            result.resample.times.push_back(std::min(options.initialTime + index * step, options.finalTime));
        if (result.resample.times.back() < options.finalTime - 1.0e-12 * std::max(1.0, std::abs(options.finalTime)))
            result.resample.times.push_back(options.finalTime);
    }
    if (!result.resample.times.empty()) {
        requireIncreasing(result.resample.times, "resampling times");
        const double slack = 32.0 * std::numeric_limits<double>::epsilon()
            * std::max({1.0, std::abs(options.initialTime), std::abs(options.finalTime)});
        if (result.resample.times.front() < options.initialTime - slack
            || result.resample.times.back() > options.finalTime + slack)
            throw std::invalid_argument("resampling times must lie within the solved time span");
        for (auto& value : result.resample.times)
            value = std::min(std::max(value, options.initialTime), options.finalTime);
        result.resample.enabled = true;
        result.resample.method = core::interpolationFromName(interpolation.at("method").string("pchip"));
    }

    const auto& history = job.at("history");
    result.checkpoints.interval = history.at("checkpoint_interval").number(0.0);
    const auto pause = history.at("pause_file").string();
    if (!(result.checkpoints.interval >= 0.0))
        throw std::invalid_argument("history.checkpoint_interval must be nonnegative");
    result.checkpoints.pauseFile = std::filesystem::u8path(pause);
    result.checkpoints.wallInterval = history.at("checkpoint_wall_seconds").number(0.0);
    result.checkpoints.deadlineUnix = history.at("deadline_unix_seconds").number(0.0);
    if (!(result.checkpoints.wallInterval >= 0.0) || !(result.checkpoints.deadlineUnix >= 0.0))
        throw std::invalid_argument("history.checkpoint_wall_seconds and deadline_unix_seconds must be nonnegative");
    const bool wallClock = result.checkpoints.wallInterval > 0.0 || result.checkpoints.deadlineUnix > 0.0;
    result.checkpoints.enabled = result.checkpoints.interval > 0.0 || !pause.empty() || wallClock;
    // Wall-clock keys need a call at every sample; CheckpointClock then gates by model time too.
    options.checkpointInterval = wallClock ? 0.0 : result.checkpoints.interval;
    return result;
}

}  // namespace dynsyslab::workflow::detail
