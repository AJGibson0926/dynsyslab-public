#include "CampaignDesign.hpp"
#include "HistoryPolicy.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow::detail {
namespace {

struct QualityProfile { std::size_t points; std::size_t systemsPerJob; };
struct CycleProfile { double transient; double retained; };

// Bifurcation resolution levels carried over from the Python campaign designer.
const std::map<std::string, QualityProfile>& qualityProfiles() {
    static const std::map<std::string, QualityProfile> profiles{
        {"low", {31, 10}}, {"medium", {101, 25}}, {"high", {301, 50}},
        {"ultra_high", {1001, 100}}};
    return profiles;
}

// Lyapunov-exponent screens discard a long transient before measuring (Wolf et al. 1985,
// "Determining Lyapunov exponents from a time series", Physica D 16:285-317).
const std::map<std::string, CycleProfile>& cycleProfiles() {
    static const std::map<std::string, CycleProfile> profiles{
        {"lyapunov_screen", {2500, 1024}}, {"lyapunov_verification", {5000, 4096}}};
    return profiles;
}

const std::map<std::string, json::Value>& executionProfiles() {
    static const std::map<std::string, json::Value> profiles{
        {"fixed", json::Value::Object{{"retry_limit", 0.0}, {"remote_allowed", false}}},
        {"adaptive_local", json::Value::Object{{"retry_limit", 1.0}, {"remote_allowed", false}}},
        {"adaptive_portable", json::Value::Object{{"retry_limit", 1.0}, {"remote_allowed", true}}}};
    return profiles;
}

// "recoverable" checkpoints a run estimated to take over 300 s every 600 s of wall time
// (history.checkpoint_wall_seconds, HistoryPolicy.hpp). A single solve pauses at a checkpoint
// and continues from it; a campaign array pauses at chunk boundaries.
const std::map<std::string, json::Value>& checkpointProfiles() {
    static const std::map<std::string, json::Value> profiles{
        {"fast", json::Value::Object{{"enabled", false}}},
        {"recoverable", json::Value::Object{{"enabled", true}, {"checkpoint_interval_seconds", 600.0},
            {"if_estimated_seconds_over", 300.0}, {"pause_at_checkpoint", true},
            {"granularity", "integration"}}}};
    return profiles;
}

std::size_t positive(const json::Value& value, const std::string& field) {
    const double number = value.number(-1.0);
    if (number < 1.0 || number != std::floor(number))
        throw std::invalid_argument(field + " must be a positive integer");
    return static_cast<std::size_t>(number);
}

const QualityProfile* quality(const json::Value& campaign) {
    const auto name = campaign.at("profile").string();
    if (name.empty()) return nullptr;
    const auto found = qualityProfiles().find(name);
    if (found == qualityProfiles().end())
        throw std::invalid_argument("campaign.profile must be low, medium, high or ultra_high");
    return &found->second;
}

json::Value resolvedAxis(const json::Value& source, const json::Value& job,
                         const QualityProfile* profile) {
    std::string path = source.at("path").string();
    if (path.rfind("driver.", 0) == 0) {
        // The native periodic driver is the waveform controller (System.cpp SinusoidalDriver).
        if (job.at("controller").at("kind").string("none") != "waveform")
            throw std::invalid_argument("a driver.* campaign axis needs controller.kind \"waveform\"");
        path = "controller.settings." + path.substr(7);
    }
    json::Value::Object axis;
    axis["path"] = path;
    if (!source.at("values").isNull()) {
        axis["values"] = source.at("values");
        return json::Value(std::move(axis));
    }
    if (source.at("start").isNull() || source.at("stop").isNull())
        throw std::invalid_argument("campaign axis needs values or start/stop: " + path);
    const std::size_t points = !source.at("points").isNull()
        ? positive(source.at("points"), "campaign axis points")
        : profile ? profile->points : 0;
    if (points < 2)
        throw std::invalid_argument("a linspace campaign axis needs points >= 2 or a campaign.profile");
    const double start = source.at("start").number(), stop = source.at("stop").number();
    json::Value::Array values;
    for (std::size_t index = 0; index < points; ++index)
        values.emplace_back(start + (stop - start) * static_cast<double>(index) / (points - 1));
    axis["values"] = json::Value(std::move(values));
    axis["linspace"] = json::Value::Object{{"start", start}, {"stop", stop},
                                           {"points", static_cast<double>(points)}};
    return json::Value(std::move(axis));
}

double driverFrequency(const json::Value& job) {
    const auto& time = job.at("time");
    if (!time.at("frequency_hz").isNull()) return time.at("frequency_hz").number();
    const auto& settings = job.at("controller").at("settings");
    for (const char* key : {"frequency_hz", "frequency"})
        if (!settings.at(key).isNull()) return settings.at(key).number();
    for (const char* key : {"frequencies_hz", "frequencies"})
        if (!settings.at(key).array().empty()) return settings.at(key).array().front().number();
    throw std::invalid_argument(
        "cycle-based time needs time.frequency_hz or controller.settings.frequency_hz");
}

}  // namespace

json::Value resolveCycleTime(const json::Value& job) {
    const auto& time = job.at("time");
    const bool cycles = !time.at("retained_cycles").isNull() || !time.at("cycle_profile").isNull();
    if (!cycles) return job;
    double transient = time.at("transient_cycles").number(0.0);
    double retained = time.at("retained_cycles").number(0.0);
    if (!time.at("cycle_profile").isNull()) {
        const auto found = cycleProfiles().find(time.at("cycle_profile").string());
        if (found == cycleProfiles().end())
            throw std::invalid_argument("time.cycle_profile must be lyapunov_screen or lyapunov_verification");
        if (time.at("transient_cycles").isNull()) transient = found->second.transient;
        if (time.at("retained_cycles").isNull()) retained = found->second.retained;
    }
    const auto samples = static_cast<double>(positive(time.at("samples_per_cycle"),
                                                      "time.samples_per_cycle"));
    const double frequency = driverFrequency(job);
    if (!(frequency > 0.0) || transient < 0.0 || !(retained > 0.0)
        || transient != std::floor(transient) || retained != std::floor(retained))
        throw std::invalid_argument(
            "cycle-based time needs a positive frequency, whole transient cycles >= 0 and retained cycles > 0");
    json::Value result = job;
    auto& target = result.at("time");
    const double start = time.at("start").number(0.0);
    target.at("start") = start;
    target.at("stop") = start + (transient + retained) / frequency;
    target.at("output_step") = 1.0 / (frequency * samples);
    target.at("retained_start") = start + transient / frequency;
    target.at("transient_cycles") = transient;
    target.at("retained_cycles") = retained;
    target.at("cycle_frequency_hz") = frequency;
    return result;
}

json::Value resolveJobDesign(const json::Value& job) {
    if (job.at("kind").string() != "solve") return job;
    json::Value result = job;
    const auto& campaign = job.at("campaign");
    if (campaign.at("enabled").boolean(false)) {
        const auto* profile = quality(campaign);
        json::Value::Array axes;
        for (const auto& axis : campaign.at("axes").array())
            axes.push_back(resolvedAxis(axis, job, profile));
        auto& target = result.at("campaign");
        target.at("axes") = json::Value(std::move(axes));
        if (profile && campaign.at("systems_per_job").isNull())
            target.at("systems_per_job") = static_cast<double>(profile->systemsPerJob);
        if (!campaign.at("execution_profile").isNull()) {
            const auto found = executionProfiles().find(campaign.at("execution_profile").string());
            if (found == executionProfiles().end())
                throw std::invalid_argument(
                    "campaign.execution_profile must be fixed, adaptive_local or adaptive_portable");
            if (job.at("queue").at("retry_limit").isNull())
                result.at("queue").at("retry_limit") = found->second.at("retry_limit");
        }
        if (!campaign.at("checkpoint_profile").isNull()) {
            const auto found = checkpointProfiles().find(campaign.at("checkpoint_profile").string());
            if (found == checkpointProfiles().end())
                throw std::invalid_argument("campaign.checkpoint_profile must be fast or recoverable");
            target.at("checkpoint_policy") = found->second;
        }
    }
    // A swept frequency is re-resolved per point by the campaign; the base job uses the
    // first value so planning sees a concrete time grid.
    json::Value base = result;
    for (const auto& axis : std::as_const(result).at("campaign").at("axes").array()) {
        const auto path = axis.at("path").string();
        if ((path == "controller.settings.frequency_hz" || path == "controller.settings.frequency"
             || path == "time.frequency_hz")
            && campaign.at("enabled").boolean(false) && !axis.at("values").array().empty())
            json::setPath(base, path, axis.at("values").array().front());
    }
    const auto resolved = resolveCycleTime(base);
    result.at("time") = resolved.at("time");
    return applyHistoryPolicy(result);
}

std::size_t campaignSystemsPerJob(const json::Value& job) {
    const auto& campaign = job.at("campaign");
    if (!campaign.at("systems_per_job").isNull())
        return positive(campaign.at("systems_per_job"), "campaign.systems_per_job");
    const auto* profile = quality(campaign);
    return profile ? profile->systemsPerJob : 0;
}

json::Value campaignProfiles() {
    json::Value::Object qualities, cycles, execution, checkpoint;
    for (const auto& [name, item] : qualityProfiles())
        qualities[name] = json::Value::Object{{"points", static_cast<double>(item.points)},
                                              {"systems_per_job", static_cast<double>(item.systemsPerJob)}};
    for (const auto& [name, item] : cycleProfiles())
        cycles[name] = json::Value::Object{{"transient_cycles", item.transient},
                                           {"retained_cycles", item.retained}};
    for (const auto& [name, item] : executionProfiles()) execution[name] = item;
    for (const auto& [name, item] : checkpointProfiles()) checkpoint[name] = item;
    return json::Value::Object{{"schema", "dynsyslab.campaign-profiles/v1"},
        {"quality", json::Value(std::move(qualities))}, {"cycle", json::Value(std::move(cycles))},
        {"execution", json::Value(std::move(execution))},
        {"checkpoint", json::Value(std::move(checkpoint))}};
}

}  // namespace dynsyslab::workflow::detail
