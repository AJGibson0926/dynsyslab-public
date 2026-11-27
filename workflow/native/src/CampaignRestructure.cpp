// Campaign restructuring and coverage review: the native form of the Python JobManager's
// restructure, controlled parallel transition and post-transition review
// (job_manager/execution_services.py, transition_services.py). Native campaign arrays run in
// parallel from submission, so a "transition" is a re-chunking of unstarted chunks.
#include <dynsyslab/workflow/Queue.hpp>

#include <dynsyslab/core/CampaignRegistry.hpp>
#include <dynsyslab/core/JobStore.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Text.hpp>

#include "WorkflowSupport.hpp"

#include <algorithm>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow {
namespace {

using dynsyslab::core::JobStatus;

std::filesystem::path databaseOf(const std::filesystem::path& root) {
    return root / L"job_manager" / L"state.sqlite";
}

std::filesystem::path writeReport(const std::filesystem::path& directory, const std::string& name,
                                  const json::Value& report) {
    std::filesystem::create_directories(core::extendedLengthPath(directory));
    const auto path = directory / core::wide(detail::safeName(name, "report") + ".json");
    std::ofstream out(core::extendedLengthPath(path), std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("unable to write " + path.generic_u8string());
    out << json::dump(report);
    return path;
}

bool unstarted(JobStatus status) {
    return status == JobStatus::queued || status == JobStatus::reserved || status == JobStatus::paused;
}

}  // namespace

json::Value campaignRestructure(const std::filesystem::path& root, const std::string& id,
                                std::size_t pointsPerChunk, const std::string& reason) {
    if (!pointsPerChunk) throw std::invalid_argument("campaign restructure needs points per chunk >= 1");
    core::CampaignRegistry registry(databaseOf(root));
    if (!registry.plan(id)) throw std::invalid_argument("campaign is not registered: " + id);
    core::JobStore store(databaseOf(root));
    JobQueue queue(root);
    json::Value::Array changes;
    for (const auto& chunk : registry.chunks(id)) {
        if (!chunk.manifestPath.empty() || chunk.jobId.empty() || chunk.pointCount <= pointsPerChunk) continue;
        const auto parent = store.find(chunk.jobId);
        if (!parent || !unstarted(parent->status)) continue;
        auto base = json::parse(parent->manifestJson);
        base.at("campaign").object().erase("systems_per_job");  // one task per new chunk
        base.at("campaign").object().erase("profile");
        base.at("queue").at("array") = false;
        json::Value::Array children;
        const auto label = parent->label;
        for (std::size_t first = chunk.firstPoint, part = 0; first < chunk.firstPoint + chunk.pointCount;
             first += pointsPerChunk, ++part) {
            const std::size_t points = std::min(pointsPerChunk, chunk.firstPoint + chunk.pointCount - first);
            json::Value child = base;
            child.at("label") = label + "-part-" + std::to_string(part);
            child.at("campaign").at("point_range") = json::Value(json::Value::Array{
                json::Value(static_cast<double>(first)), json::Value(static_cast<double>(points))});
            child.at("queue").at("logical_parent_job") = parent->id;
            const auto jobId = queue.submit(child);
            registry.recordChunk({id, first, points, jobId, {}, {}});
            children.emplace_back(json::Value::Object{{"job_id", jobId},
                {"first_point", static_cast<double>(first)}, {"point_count", static_cast<double>(points)}});
        }
        store.cancel(parent->id, "restructured into chunks of " + std::to_string(pointsPerChunk) + " points");
        changes.emplace_back(json::Value::Object{{"parent_job_id", parent->id},
            {"first_point", static_cast<double>(chunk.firstPoint)},
            {"point_count", static_cast<double>(chunk.pointCount)},
            {"children", json::Value(std::move(children))}});
    }
    const auto stamp = detail::timestamp();
    json::Value report = json::Value::Object{{"schema", "dynsyslab.job-restructure/v1"},
        {"created_utc", stamp}, {"campaign_id", id},
        {"points_per_chunk", static_cast<double>(pointsPerChunk)},
        {"reason", reason.empty() ? std::string("manual restructuring") : reason},
        {"restructured", json::Value(std::move(changes))}};
    if (!report.at("restructured").array().empty())
        report.at("manifest") = writeReport(root / L"job_manager" / L"restructure",
                                            id + "__" + stamp, report).generic_u8string();
    return report;
}

json::Value campaignReview(const std::filesystem::path& root, const std::string& id) {
    core::CampaignRegistry registry(databaseOf(root));
    const auto plan = registry.plan(id);
    if (!plan) throw std::invalid_argument("campaign is not registered: " + id);
    core::JobStore store(databaseOf(root));
    std::vector<std::size_t> covered(plan->pointCount, 0);
    json::Value::Array issues, repair;
    std::size_t published = 0, pending = 0;
    for (const auto& chunk : registry.chunks(id)) {
        const json::Value where = json::Value::Object{{"first_point", static_cast<double>(chunk.firstPoint)},
            {"point_count", static_cast<double>(chunk.pointCount)}, {"job_id", chunk.jobId}};
        if (chunk.firstPoint + chunk.pointCount > plan->pointCount)
            issues.emplace_back(json::Value::Object{{"issue", "chunk_outside_plan"}, {"chunk", where}});
        for (std::size_t point = chunk.firstPoint;
             point < std::min(chunk.firstPoint + chunk.pointCount, plan->pointCount); ++point) ++covered[point];
        if (!chunk.manifestPath.empty()) {
            if (!std::filesystem::is_regular_file(core::extendedLengthPath(std::filesystem::u8path(chunk.manifestPath))))
                issues.emplace_back(json::Value::Object{{"issue", "manifest_missing"}, {"chunk", where},
                                                        {"manifest", chunk.manifestPath}});
            else published += chunk.pointCount;
            continue;
        }
        std::optional<core::JobRecord> job;
        if (!chunk.jobId.empty()) job = store.find(chunk.jobId);
        if (!job) { repair.push_back(where); continue; }
        if (job->status == JobStatus::completed)
            issues.emplace_back(json::Value::Object{{"issue", "completed_without_manifest"}, {"chunk", where}});
        else if (job->status == JobStatus::failed || job->status == JobStatus::cancelled
                 || job->status == JobStatus::quarantined) repair.push_back(where);
        else ++pending;
    }
    json::Value::Array uncovered, duplicated;
    for (std::size_t point = 0; point < covered.size(); ++point) {
        if (!covered[point]) uncovered.emplace_back(static_cast<double>(point));
        if (covered[point] > 1) duplicated.emplace_back(static_cast<double>(point));
    }
    if (!uncovered.empty()) issues.emplace_back(json::Value::Object{{"issue", "uncovered_points"},
                                                                    {"points", json::Value(std::move(uncovered))}});
    if (!duplicated.empty()) issues.emplace_back(json::Value::Object{{"issue", "duplicated_points"},
                                                                     {"points", json::Value(std::move(duplicated))}});
    const auto stamp = detail::timestamp();
    json::Value report = json::Value::Object{{"schema", "dynsyslab.campaign-review/v1"},
        {"reviewed_utc", stamp}, {"campaign_id", id},
        {"point_count", static_cast<double>(plan->pointCount)},
        {"published_points", static_cast<double>(published)},
        {"pending_chunks", static_cast<double>(pending)},
        {"needs_repair", json::Value(std::move(repair))},
        {"safe_to_execute", issues.empty()}, {"issues", json::Value(std::move(issues))}};
    report.at("report_path") = writeReport(root / L"job_manager" / L"reviews",
                                           "campaign-review__" + id + "__" + stamp, report).generic_u8string();
    return report;
}

}  // namespace dynsyslab::workflow
