#include <dynsyslab/workflow/Run.hpp>

#include <dynsyslab/workflow/Campaign.hpp>
#include <dynsyslab/workflow/KromControlSearch.hpp>
#include <dynsyslab/workflow/Pipeline.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>

#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/BuildInfo.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>

#include "SolveDispatch.hpp"
#include "SolveSet.hpp"
#include "ViewerFigureJob.hpp"

#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow {
namespace {

RunResult indexed(RunResult result, const std::filesystem::path& archiveRoot) {
    auto paths = result.archives;
    if (!result.manifest.empty()) paths.push_back(result.manifest);
    // Index maintenance is part of archive publication. A later explicit
    // reconcile can repair the summary if an external writer was interrupted.
    (void)dynsyslab::core::indexArchiveFiles(archiveRoot, paths);
    return result;
}

}  // namespace

void validateJob(const json::Value& job) {
    if (job.at("schema").string() != "dynsyslab.compute-job/v1")
        throw std::invalid_argument("job must use dynsyslab.compute-job/v1");
    const std::string kind = job.at("kind").string();
    if (kind != "solve" && kind != "solve_set" && kind != "postprocess"
        && kind != "krom_control_search" && kind != "pipeline"
        && kind != "job_group" && kind != "viewer_figure")
        throw std::invalid_argument(
            "job kind must be solve, solve_set, postprocess, krom_control_search, pipeline, job_group, or viewer_figure");
}

RunResult executeJob(const json::Value& job,
                     const std::filesystem::path& requestedArchiveRoot,
                     const std::filesystem::path& progressPath,
                     const dynsyslab::core::ExecutionContext* execution) {
    // Resolve once for writers and indexers: relative published paths must not
    // be interpreted a second time beneath the archive root.
    const auto archiveRoot=std::filesystem::absolute(requestedArchiveRoot).lexically_normal();
    if(job.at("kind").string()=="pipeline")  // stages index their own outputs
        return executePipeline(job,archiveRoot,progressPath,execution);
    if(job.at("kind").string()=="job_group")  // its leaf jobs index and record themselves
        return executeJobGroup(job,archiveRoot,progressPath,execution);
    (void)planJob(job);
    // Every archive this job publishes records the resolved job and the build that ran it.
    const dynsyslab::core::ScopedRunProvenance provenance(json::dump(json::Value(json::Value::Object{
        {"schema",json::Value("dynsyslab.run-provenance/v1")},
        {"build",json::parse(dynsyslab::core::buildIdentityJson())},{"resolved_job",job}})));
    if(job.at("kind").string()=="viewer_figure")
        return detail::executeViewerFigure(job,archiveRoot);
    if(job.at("kind").string()=="solve_set")
        return indexed(detail::executeSolveSet(job,archiveRoot,progressPath,execution),
                       archiveRoot);
    if(job.at("kind").string()=="krom_control_search")
        return indexed(executeKromControlSearch(job,archiveRoot,progressPath,execution),
                       archiveRoot);
    if (job.at("kind").string() == "postprocess") {
        auto result = executePostprocess(job, archiveRoot, progressPath, execution);
        return indexed({"postprocess", result.summary(), std::move(result.archives), {},
                        std::move(result.methods)}, archiveRoot);
    }
    if (hasEnabledCampaign(job)) {
        auto result = executeCampaign(job, archiveRoot, progressPath, execution);
        std::vector<std::filesystem::path> archives;
        for (const auto& point : result.points)
            archives.insert(archives.end(), point.solve.archives.begin(),
                            point.solve.archives.end());
        return indexed({"campaign", result.summary(), std::move(archives),
                        std::move(result.manifest)}, archiveRoot);
    }
    auto result = detail::executeSingleSolve(
        job, archiveRoot, progressPath, execution);
    return indexed({"solve", result.summary(), std::move(result.archives),
                    std::move(result.hierarchyManifest)}, archiveRoot);
}

}  // namespace dynsyslab::workflow
