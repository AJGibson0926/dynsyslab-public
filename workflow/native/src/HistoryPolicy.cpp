#include "HistoryPolicy.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace dynsyslab::workflow::detail {
namespace {

json::Value resolution(bool enabled, const std::string& reason, const json::Value& estimate,
                       const json::Value& threshold, const std::string& mode = {}) {
    json::Value::Object result{{"enabled", enabled}, {"reason", reason},
                               {"estimated_seconds", estimate}, {"threshold_seconds", threshold}};
    if (!mode.empty()) result["mode"] = mode;
    return json::Value(std::move(result));
}

}  // namespace

json::Value applyHistoryPolicy(const json::Value& submitted) {
    if (submitted.at("kind").string() != "solve") return submitted;
    if (!submitted.at("history").at("policy_resolution").isNull()) return submitted;
    const auto& explicitPolicy = submitted.at("history_policy");
    const auto& policy = explicitPolicy.isNull()
        ? (submitted.at("campaign").at("enabled").boolean(false)
               ? submitted.at("campaign").at("checkpoint_policy") : explicitPolicy)
        : explicitPolicy;
    if (policy.isNull()) return submitted;
    json::Value job = submitted;
    auto& history = job.at("history");
    const auto& queue = submitted.at("queue");
    json::Value estimate = policy.at("estimated_seconds");
    if (estimate.isNull()) estimate = queue.at("estimated_seconds");
    if (estimate.isNull()) estimate = queue.at("wall_time_seconds");
    const auto& threshold = policy.at("if_estimated_seconds_over");
    const auto unknown = policy.at("when_estimate_unknown").string("disable");
    if (unknown != "enable" && unknown != "disable")
        throw std::invalid_argument("history_policy.when_estimate_unknown must be enable or disable");
    const auto mode = policy.at("mode").string("checkpoint");
    if (mode != "checkpoint" && mode != "telemetry")
        throw std::invalid_argument("history_policy.mode must be checkpoint or telemetry");
    if (!policy.at("enabled").boolean(false)) {
        history.at("policy_resolution") = resolution(false, "no_policy", estimate, threshold);
        return job;
    }
    if (!threshold.isNull() && estimate.isNull() && unknown == "disable") {
        history.at("policy_resolution") = resolution(false, "unknown_estimate", estimate, threshold);
        return job;
    }
    if (!threshold.isNull() && !estimate.isNull() && estimate.number() < threshold.number()) {
        history.at("policy_resolution") = resolution(false, "below_estimate_threshold", estimate, threshold);
        return job;
    }
    if (mode == "checkpoint") {
        if (!policy.at("checkpoint_interval").isNull())
            history.at("checkpoint_interval") = policy.at("checkpoint_interval");
        for (const char* key : {"checkpoint_interval_seconds", "interval_seconds"})
            if (!policy.at(key).isNull() && std::as_const(history).at("checkpoint_wall_seconds").isNull())
                history.at("checkpoint_wall_seconds") = policy.at(key);
        if (!policy.at("pause_at_checkpoint").isNull())
            history.at("pause_at_checkpoint") = policy.at("pause_at_checkpoint");
    }
    auto record = resolution(true, "policy_matched", estimate, threshold, mode);
    record.at("estimate_was_unknown") = estimate.isNull();
    history.at("policy_resolution") = std::move(record);
    return job;
}

}  // namespace dynsyslab::workflow::detail
