// Queue submission: the `queue` job block (retry, dependencies, resources, side tray, spec
// claims, identical-spec reuse, extension candidates) and campaign job arrays chunked by
// `campaign.systems_per_job`.
#include <dynsyslab/workflow/Queue.hpp>

#include <dynsyslab/core/CampaignRegistry.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/workflow/Campaign.hpp>
#include <dynsyslab/workflow/Plan.hpp>

#include "CampaignDesign.hpp"
#include "WorkflowSupport.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow {
namespace {

using dynsyslab::core::JobStatus;

std::string fnv(const std::string& text) {
    std::uint64_t hash = 1469598103934665603ull;
    for (const unsigned char value : text) { hash ^= value; hash *= 1099511628211ull; }
    std::ostringstream out; out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

json::Value scientificContent(const json::Value& job) {
    json::Value content = job;
    for (const char* key : {"label", "queue", "execution"}) content.object().erase(key);
    if (!std::as_const(content).at("campaign").isNull())
        content.at("campaign").object().erase("id");
    return content;
}

// Same solve except for its final time: a shorter run could be extended, a longer one cut.
std::string extensionKey(const json::Value& job) {
    if (job.at("kind").string() != "solve" || job.at("campaign").at("enabled").boolean(false))
        return {};
    auto content = scientificContent(job);
    content.at("time").object().erase("stop");
    return fnv(json::dump(content));
}

std::size_t count(const json::Value& value, const char* field) {
    const double number = value.number(0.0);
    if (number < 0.0 || number != std::floor(number))
        throw std::invalid_argument(std::string(field) + " must be a nonnegative integer");
    return static_cast<std::size_t>(number);
}

dynsyslab::core::JobSubmission baseSubmission(const json::Value& job,
                                              dynsyslab::core::JobStore& store) {
    const auto& queue = job.at("queue");
    dynsyslab::core::JobSubmission result;
    result.label = job.at("label").string("job");
    result.kind = job.at("kind").string();
    result.manifestJson = json::dump(job);
    // Deadline requeues (queue.checkpoint_before_deadline_seconds) spend the retry budget;
    // queue.max_requeues raises it (job_manager/scheduling.py DeadlinePolicy.max_requeues).
    result.retryLimit = std::max(count(queue.at("retry_limit"), "queue.retry_limit"),
                                 count(queue.at("max_requeues"), "queue.max_requeues"));
    for (const auto& item : queue.at("depends_on").array()) {
        dynsyslab::core::JobDependency dependency;
        dependency.job = item.string(item.at("job").string());
        if (!item.at("condition").isNull())
            dependency.condition = dynsyslab::core::parseDependencyCondition(
                item.at("condition").string());
        if (dependency.job.empty()) throw std::invalid_argument("queue.depends_on needs job ids");
        if (dependency.condition != dynsyslab::core::DependencyCondition::aftercorr
            && !store.find(dependency.job))
            throw std::invalid_argument("queue.depends_on names an unknown job: " + dependency.job);
        result.dependencies.push_back(std::move(dependency));
    }
    const auto limits = jobExecutionLimits(job);
    const auto& resources = queue.at("resources");
    result.resources.threads = resources.at("threads").isNull() ? limits.maximumThreads
        : count(resources.at("threads"), "queue.resources.threads");
    result.resources.memoryBytes = resources.at("memory_mib").isNull() ? limits.maximumMemoryBytes
        : static_cast<std::size_t>(resources.at("memory_mib").number() * 1024.0 * 1024.0);
    result.reserved = queue.at("hold").boolean(false);
    if (queue.at("spec_claim").boolean(true)) result.specKey = jobSpecKey(job);
    return result;
}

json::Value extensionCandidates(const json::Value& job, const dynsyslab::core::JobStore& store) {
    json::Value::Array result;
    const auto key = extensionKey(job);
    if (key.empty()) return json::Value(std::move(result));
    const double stop = job.at("time").at("stop").number();
    for (const auto& record : store.list(1000, JobStatus::completed)) {
        if (record.kind != "solve") continue;
        const auto other = json::parse(record.manifestJson);
        if (extensionKey(other) != key) continue;
        const double otherStop = other.at("time").at("stop").number();
        if (otherStop == stop) continue;
        result.emplace_back(json::Value::Object{{"job", record.id}, {"stop", otherStop},
            {"relation", otherStop < stop ? "shorter" : "longer"},
            // queue.reuse "extend" continues the longest shorter run (extendedJob).
            {"action", otherStop < stop ? "extend_with_reuse_extend" : "run_in_full"}});
    }
    return json::Value(std::move(result));
}

// queue.reuse "extend": the same solve continued from the end of the longest completed shorter
// run (initial_conditions continue_from its archives, time.start = its stop), recorded under
// `extension_of`. Port of job_manager/reuse_services.py _try_single_solution_extension; the
// parent archives stay immutable and each child archive links its parent with a sha256
// (SolveContinuation.hpp). Null when no shorter run covers every trajectory.
json::Value extendedJob(const json::Value& job, const dynsyslab::core::JobStore& store) {
    const auto key = extensionKey(job);
    if (key.empty() || !job.at("time").at("output_times").isNull()) return {};
    const double stop = job.at("time").at("stop").number();
    const double expected = job.at("initial_conditions").at("count").number(1.0);
    std::optional<dynsyslab::core::JobRecord> best;
    double bestStop = 0.0;
    json::Value::Array archives;
    for (const auto& record : store.list(1000, JobStatus::completed)) {
        if (record.kind != "solve") continue;
        const auto other = json::parse(record.manifestJson);
        if (extensionKey(other) != key) continue;
        const double otherStop = other.at("time").at("stop").number();
        if (!(otherStop < stop) || (best && otherStop <= bestStop)) continue;
        json::Value result;
        try { result = json::parse(record.message); } catch (const std::exception&) { continue; }
        const auto& paths = result.at("archives").array();
        if (static_cast<double>(paths.size()) != expected) continue;  // a rejected member
        best = record; bestStop = otherStop; archives = paths;
    }
    if (!best) return {};
    json::Value extended = job;
    const double start = job.at("time").at("start").number(0.0);
    extended.at("time").at("start") = bestStop;
    extended.at("initial_conditions") = json::Value::Object{{"mode", "continue_from"},
        {"count", expected}, {"archives", json::Value(archives)}};
    extended.at("extension_of") = json::Value::Object{{"job", best->id}, {"stop", bestStop},
        {"archives", json::Value(archives)}, {"requested_start", start},
        {"requested_initial_conditions", job.at("initial_conditions")}};
    return extended;
}

}  // namespace

std::string jobSpecKey(const json::Value& job) {
    return fnv(json::dump(scientificContent(job)));
}

json::Value QueueSubmitResult::toJson() const {
    json::Value::Array ids;
    for (const auto& id : jobIds) ids.emplace_back(id);
    json::Value::Object result;
    result["job_id"] = reusedJob.empty() ? (jobIds.empty() ? std::string{} : jobIds.front())
                                         : reusedJob;
    result["job_ids"] = json::Value(std::move(ids));
    if (!arrayId.empty()) result["array_id"] = arrayId;
    if (!campaignId.empty()) result["campaign_id"] = campaignId;
    if (!reusedJob.empty()) result["reused_job"] = reusedJob;
    if (!extendedFrom.empty()) result["extended_from"] = extendedFrom;
    if (!extensionCandidates.array().empty()) result["extension_candidates"] = extensionCandidates;
    return json::Value(std::move(result));
}

QueueSubmitResult JobQueue::submitDetailed(const json::Value& submitted) {
    validateJob(submitted);
    (void)planJob(submitted);
    const auto job = detail::resolveJobDesign(submitted);
    dynsyslab::core::JobStore store(database_);
    QueueSubmitResult result;
    const auto& queue = job.at("queue");
    const std::string reuse = queue.at("reuse").string("none");
    if (reuse != "none" && reuse != "identical" && reuse != "extend")
        throw std::invalid_argument("queue.reuse must be none, identical or extend");

    const bool campaign = job.at("kind").string() == "solve" && hasEnabledCampaign(job);
    const std::size_t chunk = campaign ? detail::campaignSystemsPerJob(job) : 0;
    const bool array = campaign && (queue.at("array").boolean(false) || chunk > 0);
    if (!array) {
        if (reuse == "extend") {
            if (auto extended = extendedJob(job, store); !extended.isNull()) {
                result.extendedFrom = extended.at("extension_of").at("job").string();
                result.jobIds.push_back(store.submit(baseSubmission(extended, store)));
                return result;
            }
        }
        auto submission = baseSubmission(job, store);
        if (reuse == "identical") {
            if (const auto previous = store.findBySpec(submission.specKey)) {
                result.reusedJob = previous->id;
                return result;
            }
        }
        result.extensionCandidates = extensionCandidates(job, store);
        result.jobIds.push_back(store.submit(submission));
        return result;
    }

    const std::size_t total = campaignPointCount(job);
    const std::size_t size = chunk ? chunk : 1;
    const std::string label = job.at("label").string("campaign");
    result.campaignId = job.at("campaign").at("id").string(
        dynsyslab::core::stableIdentifier(label + "|" + detail::timestamp() + "|campaign"));
    result.arrayId = "array-" + result.campaignId;
    dynsyslab::core::CampaignRegistry registry(database_);
    json::Value plan = job;
    plan.at("campaign").at("id") = result.campaignId;
    plan.object().erase("queue");
    registry.registerPlan({result.campaignId, label,
        job.at("campaign").at("kind").string("control_sweep"), total, json::dump(plan), {}});
    const std::size_t maxConcurrent = count(queue.at("max_concurrent"), "queue.max_concurrent");
    for (std::size_t first = 0, index = 0; first < total; first += size, ++index) {
        const std::size_t points = std::min(size, total - first);
        json::Value child = job;
        child.at("label") = label + "-chunk-" + std::to_string(index);
        child.at("campaign").at("id") = result.campaignId;
        child.at("campaign").at("point_range") = json::Value(json::Value::Array{
            json::Value(static_cast<double>(first)), json::Value(static_cast<double>(points))});
        auto submission = baseSubmission(child, store);
        submission.arrayId = result.arrayId;
        submission.arrayIndex = index;
        submission.arrayMaxConcurrent = maxConcurrent;
        const auto id = store.submit(submission);
        registry.recordChunk({result.campaignId, first, points, id, {}, {}});
        result.jobIds.push_back(id);
    }
    return result;
}

}  // namespace dynsyslab::workflow
