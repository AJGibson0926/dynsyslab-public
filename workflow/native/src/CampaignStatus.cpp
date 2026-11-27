// Campaign registry reports: list, status (queue state per chunk) and validation of the
// published chunk manifests against the original campaign plan.
#include <dynsyslab/workflow/Queue.hpp>

#include <dynsyslab/core/CampaignRegistry.hpp>
#include <dynsyslab/core/JobStore.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Campaign.hpp>

#include "CampaignDesign.hpp"

#include <algorithm>
#include <map>
#include <thread>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

std::filesystem::path databaseOf(const std::filesystem::path& archiveRoot) {
    return archiveRoot / L"job_manager" / L"state.sqlite";
}

dynsyslab::core::CampaignPlanRecord requirePlan(const dynsyslab::core::CampaignRegistry& registry,
                                                const std::string& id) {
    auto plan = registry.plan(id);
    if (!plan) throw std::invalid_argument("campaign is not registered: " + id);
    return *plan;
}

// Expected coordinates of every point, by the campaign's Cartesian ordering.
std::vector<json::Value> plannedCoordinates(const json::Value& plan, std::size_t count) {
    const auto job = detail::resolveJobDesign(plan);
    const auto& axes = job.at("campaign").at("axes").array();
    std::vector<json::Value> result;
    for (std::size_t point = 0; point < count; ++point) {
        json::Value::Array coordinates;
        std::size_t divisor = 1;
        for (const auto& axis : axes) {
            const auto& values = axis.at("values").array();
            coordinates.push_back(values.at((point / divisor) % values.size()));
            divisor *= values.size();
        }
        result.emplace_back(std::move(coordinates));
    }
    return result;
}

json::Value readManifest(const std::string& path) {
    if (path.empty() || !std::filesystem::is_regular_file(dynsyslab::core::wide(path)))
        return {};
    return json::load(dynsyslab::core::wide(path));
}

}  // namespace

json::Value campaignList(const std::filesystem::path& archiveRoot) {
    dynsyslab::core::CampaignRegistry registry(databaseOf(archiveRoot));
    json::Value::Array campaigns;
    for (const auto& plan : registry.plans(1000))
        campaigns.emplace_back(json::Value::Object{{"id", plan.id}, {"label", plan.label},
            {"kind", plan.kind}, {"point_count", static_cast<double>(plan.pointCount)},
            {"chunks", static_cast<double>(registry.chunks(plan.id).size())},
            {"created_utc", plan.createdUtc}});
    return json::Value::Object{{"schema", "dynsyslab.campaign-list/v1"},
                               {"campaigns", json::Value(std::move(campaigns))}};
}

json::Value campaignStatus(const std::filesystem::path& archiveRoot, const std::string& id) {
    dynsyslab::core::CampaignRegistry registry(databaseOf(archiveRoot));
    dynsyslab::core::JobStore store(databaseOf(archiveRoot));
    const auto plan = requirePlan(registry, id);
    std::map<std::string, double> byStatus;
    std::size_t published = 0;
    json::Value::Array chunks;
    for (const auto& chunk : registry.chunks(id)) {
        json::Value::Object item{{"first_point", static_cast<double>(chunk.firstPoint)},
            {"point_count", static_cast<double>(chunk.pointCount)}};
        std::string status = chunk.manifestPath.empty() ? "unscheduled" : "completed";
        if (!chunk.jobId.empty()) {
            item["job_id"] = chunk.jobId;
            if (const auto job = store.find(chunk.jobId)) {
                status = dynsyslab::core::jobStatusName(job->status);
                item["attempts"] = static_cast<double>(job->attempts);
                item["progress"] = job->progress;
                if (job->status != dynsyslab::core::JobStatus::completed) item["message"] = job->message;
            } else status = "missing_job";
        }
        if (!chunk.manifestPath.empty()) {
            item["manifest"] = chunk.manifestPath;
            published += chunk.pointCount;
        }
        item["status"] = status;
        byStatus[status] += 1.0;
        chunks.emplace_back(std::move(item));
    }
    json::Value::Object counts;
    for (const auto& [status, total] : byStatus) counts[status] = total;
    return json::Value::Object{{"schema", "dynsyslab.campaign-status/v1"}, {"id", id},
        {"label", plan.label}, {"kind", plan.kind},
        {"point_count", static_cast<double>(plan.pointCount)},
        {"published_points", static_cast<double>(published)},
        {"percent", plan.pointCount ? 100.0 * published / plan.pointCount : 0.0},
        {"chunk_status", json::Value(std::move(counts))}, {"chunks", json::Value(std::move(chunks))}};
}

json::Value campaignValidation(const std::filesystem::path& archiveRoot, const std::string& id) {
    dynsyslab::core::CampaignRegistry registry(databaseOf(archiveRoot));
    const auto plan = requirePlan(registry, id);
    const auto original = json::parse(plan.planJson);
    const auto expected = plannedCoordinates(original, plan.pointCount);
    std::vector<std::size_t> covered(plan.pointCount, 0);
    json::Value::Array mismatched, missingArchives, unpublished, problems;
    for (const auto& chunk : registry.chunks(id)) {
        const auto manifest = readManifest(chunk.manifestPath);
        if (manifest.isNull()) {
            unpublished.emplace_back(json::Value::Object{
                {"first_point", static_cast<double>(chunk.firstPoint)},
                {"point_count", static_cast<double>(chunk.pointCount)}, {"job_id", chunk.jobId}});
            continue;
        }
        if (manifest.at("id").string() != id)
            problems.emplace_back("manifest " + chunk.manifestPath + " names another campaign");
        for (const auto& point : manifest.at("points").array()) {
            const double raw = point.at("index").number(-1.0);
            if (raw < 0.0 || raw >= static_cast<double>(plan.pointCount)) {
                problems.emplace_back("manifest " + chunk.manifestPath + " has a point outside the plan");
                continue;
            }
            const auto index = static_cast<std::size_t>(raw);
            ++covered[index];
            if (json::dump(point.at("coordinates")) != json::dump(expected[index]))
                mismatched.emplace_back(json::Value::Object{{"index", raw},
                    {"planned", expected[index]}, {"published", point.at("coordinates")}});
            const auto& archives = point.at("solution_archives").array();
            if (archives.empty())
                missingArchives.emplace_back(json::Value::Object{{"index", raw}, {"archive", ""}});
            for (const auto& archive : archives)
                if (!std::filesystem::is_regular_file(dynsyslab::core::wide(archive.string())))
                    missingArchives.emplace_back(json::Value::Object{{"index", raw},
                                                                     {"archive", archive}});
        }
    }
    json::Value::Array missing, duplicated;
    for (std::size_t index = 0; index < covered.size(); ++index) {
        if (!covered[index]) missing.emplace_back(static_cast<double>(index));
        if (covered[index] > 1) duplicated.emplace_back(static_cast<double>(index));
    }
    const bool valid = missing.empty() && duplicated.empty() && mismatched.empty()
        && missingArchives.empty() && problems.empty();
    return json::Value::Object{{"schema", "dynsyslab.campaign-validation/v1"}, {"id", id},
        {"label", plan.label}, {"valid", valid},
        {"point_count", static_cast<double>(plan.pointCount)},
        {"missing_points", json::Value(std::move(missing))},
        {"duplicated_points", json::Value(std::move(duplicated))},
        {"mismatched_points", json::Value(std::move(mismatched))},
        {"missing_archives", json::Value(std::move(missingArchives))},
        {"unpublished_chunks", json::Value(std::move(unpublished))},
        {"problems", json::Value(std::move(problems))}};
}

json::Value campaignRepair(const std::filesystem::path& archiveRoot, const std::string& id) {
    dynsyslab::core::CampaignRegistry registry(databaseOf(archiveRoot));
    const auto plan = requirePlan(registry, id);
    JobQueue queue(archiveRoot);
    json::Value::Array resubmitted;
    for (const auto& chunk : registry.chunks(id)) {
        if (!chunk.manifestPath.empty()) continue;
        if (!chunk.jobId.empty()) {
            const auto job = queue.find(chunk.jobId);
            if (job && job->status != dynsyslab::core::JobStatus::failed
                && job->status != dynsyslab::core::JobStatus::cancelled
                && job->status != dynsyslab::core::JobStatus::quarantined) continue;
        }
        auto child = json::parse(plan.planJson);
        child.at("label") = plan.label + "-repair-" + std::to_string(chunk.firstPoint);
        child.at("campaign").at("id") = id;
        child.at("campaign").at("point_range") = json::Value(json::Value::Array{
            json::Value(static_cast<double>(chunk.firstPoint)),
            json::Value(static_cast<double>(chunk.pointCount))});
        child.at("campaign").object().erase("systems_per_job");  // one task for this chunk
        child.at("campaign").object().erase("profile");
        const auto jobId = queue.submit(child);
        registry.recordChunk({id, chunk.firstPoint, chunk.pointCount, jobId, {}, {}});
        resubmitted.emplace_back(json::Value::Object{
            {"first_point", static_cast<double>(chunk.firstPoint)},
            {"point_count", static_cast<double>(chunk.pointCount)}, {"job_id", jobId}});
    }
    return json::Value::Object{{"schema", "dynsyslab.campaign-repair/v1"}, {"id", id},
                               {"resubmitted", json::Value(std::move(resubmitted))}};
}

json::Value queueResources(const std::filesystem::path& archiveRoot) {
    const std::size_t cores = std::max<unsigned>(1, std::thread::hardware_concurrency());
    std::error_code ignored;
    const auto space = std::filesystem::space(archiveRoot, ignored);
    dynsyslab::core::JobStore store(databaseOf(archiveRoot));
    const auto running = store.list(1000, dynsyslab::core::JobStatus::running).size();
    // Leave one core for the desktop, at most four local workers, halved while jobs run.
    std::size_t workers = std::max<std::size_t>(1, std::min<std::size_t>(4, cores - 1));
    if (running) workers = std::max<std::size_t>(1, std::min(workers, cores / 2));
    return json::Value::Object{{"schema", "dynsyslab.queue-resources/v1"},
        {"logical_cores", static_cast<double>(cores)},
        {"archive_free_bytes", ignored ? -1.0 : static_cast<double>(space.available)},
        {"running_jobs", static_cast<double>(running)},
        {"recommended_local_workers", static_cast<double>(workers)},
        {"backends", json::Value(json::Value::Array{json::Value("in_process"),
                                                    json::Value("command_template")})}};
}

}  // namespace dynsyslab::workflow
