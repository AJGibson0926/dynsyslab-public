#include "SolveSet.hpp"

#include "SolveDispatch.hpp"
#include "WorkflowSupport.hpp"

#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Campaign.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow::detail {
namespace {

void atomicJson(const std::filesystem::path& path, const json::Value& value) {
    std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    { std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to create solve-set manifest");
      stream << json::dump(value) << '\n'; }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish solve-set manifest");
#else
    std::filesystem::rename(temporary, path);
#endif
}

std::size_t systemCount(const json::Value& job) {
    std::size_t result = 0;
    for (const auto& batch : job.at("batches").array())
        result += batch.at("systems").array().size();
    return result;
}

}  // namespace

RunResult executeSolveSet(const json::Value& job,
                          const std::filesystem::path& archiveRoot,
                          const std::filesystem::path& progressPath,
                          const dynsyslab::core::ExecutionContext* execution) {
    const std::size_t total = systemCount(job);
    if (!total) throw std::invalid_argument("solve_set contains no configured systems");
    const std::string stamp = timestamp();
    const std::string label = job.at("label").string("solve-set");
    const auto manifest = archiveRoot / L"job_sets"
        / dynsyslab::core::wide(safeName(label, "solve-set") + "__" + stamp
            + ".solve-set.json");
    json::Value::Array hierarchyManifests;
    std::vector<std::filesystem::path> archives;
    std::wostringstream summary;
    summary << L"Native solve job completed\r\n\r\nJob: "
            << dynsyslab::core::wide(label) << L"\r\nBatches: "
            << job.at("batches").array().size() << L"\r\nConfigured systems: "
            << total << L"\r\nPython required: no";
    std::size_t completed = 0;
    atomicProgress(progressPath, 0.0, 0, total, "Preparing solve job");
    for (std::size_t batchIndex = 0; batchIndex < job.at("batches").array().size();
         ++batchIndex) {
        const auto& batch = job.at("batches").array()[batchIndex];
        if (batch.at("systems").array().empty())
            throw std::invalid_argument("Every solve_set batch must contain a system");
        const std::string batchLabel = batch.at("label").string(
            "Batch " + std::to_string(batchIndex + 1));
        const std::string batchId = dynsyslab::core::stableIdentifier(
            label + "|" + stamp + "|batch|" + std::to_string(batchIndex));
        const auto hierarchy = archiveRoot / L"hierarchies" / L"job-sets"
            / dynsyslab::core::wide(safeName(label + "-" + batchLabel, "batch")
                + "__" + stamp + ".hierarchy.json");
        hierarchyManifests.emplace_back(hierarchy.generic_u8string());
        for (std::size_t systemIndex = 0;
             systemIndex < batch.at("systems").array().size(); ++systemIndex) {
            json::Value child = batch.at("systems").array()[systemIndex];
            child.at("schema") = "dynsyslab.compute-job/v1";
            child.at("kind") = "solve";
            json::Value::Object batchMetadata;
            batchMetadata["id"] = batchId;
            batchMetadata["label"] = batchLabel;
            batchMetadata["manifest_path"] = hierarchy.generic_u8string();
            batchMetadata["merge"] = true;
            batchMetadata["system_run_id"] = dynsyslab::core::stableIdentifier(
                batchId + "|system|" + std::to_string(systemIndex));
            child.at("batch") = json::Value(std::move(batchMetadata));
            const std::size_t slot = completed;
            dynsyslab::core::ExecutionContext childExecution;
            if (execution) {
                childExecution.limits = execution->limits;
                childExecution.cancellation = [execution] { return execution->cancelled(); };
            }
            childExecution.progress = [execution, progressPath, slot, total](
                dynsyslab::core::ExecutionProgress event) {
                const double local = std::clamp(event.fraction, 0.0, 1.0);
                const double global = (static_cast<double>(slot) + local) / total;
                const std::string message = "System " + std::to_string(slot + 1)
                    + " of " + std::to_string(total) + ": " + event.message;
                atomicProgress(progressPath, global, slot
                    + static_cast<std::size_t>(local >= 1.0), total, message);
                if (execution) execution->report({"solve_set/" + event.stage,
                    message, global, slot + static_cast<std::size_t>(local >= 1.0), total});
            };
            const auto childProgress = progressPath.empty() ? std::filesystem::path{}
                : progressPath.parent_path() / (progressPath.stem().wstring()
                    + L"-system-" + std::to_wstring(completed) + L".json");
            if (hasEnabledCampaign(child)) {
                auto run = executeCampaign(child, archiveRoot, childProgress, &childExecution);
                for (const auto& point : run.points)
                    archives.insert(archives.end(), point.solve.archives.begin(),
                                    point.solve.archives.end());
            } else {
                auto run = executeSingleSolve(child, archiveRoot, childProgress, &childExecution);
                archives.insert(archives.end(), run.archives.begin(), run.archives.end());
            }
            if (!childProgress.empty()) {
                std::error_code ignored; std::filesystem::remove(childProgress, ignored);
            }
            ++completed;
            atomicProgress(progressPath, static_cast<double>(completed) / total,
                completed, total, "Completed system " + std::to_string(completed)
                    + " of " + std::to_string(total));
        }
    }
    json::Value::Object record;
    record["format"] = "dynsyslab.solve-set/v1";
    record["label"] = label;
    record["created_utc"] = stamp;
    record["batch_count"] = static_cast<double>(job.at("batches").array().size());
    record["system_count"] = static_cast<double>(total);
    record["hierarchy_manifests"] = json::Value(std::move(hierarchyManifests));
    json::Value::Array paths;
    for (const auto& archive : archives) paths.emplace_back(archive.generic_u8string());
    record["solution_archives"] = json::Value(std::move(paths));
    atomicJson(manifest, json::Value(std::move(record)));
    summary << L"\r\n\r\nManifest:\r\n" << manifest.wstring();
    return {"solve_set", summary.str(), std::move(archives), manifest};
}

}  // namespace dynsyslab::workflow::detail
