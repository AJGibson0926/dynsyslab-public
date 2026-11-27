#include "SolveDispatch.hpp"
#include "ExtensionModelSolve.hpp"

#include <dynsyslab/workflow/DelaySolve.hpp>
#include <dynsyslab/workflow/LatticeSolve.hpp>
#include <dynsyslab/workflow/PDESolve.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::detail {

SolveResult executeSingleSolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution) {
    if (job.at("campaign").at("enabled").boolean(false))
        throw std::invalid_argument(
            "single-solve dispatch requires campaign.enabled=false");
    if (isExtensionModelSolve(job))
        return executeExtensionModelSolve(job, archiveRoot, progressPath, execution);
    if (isDelaySolve(job))
        return executeDelaySolve(job, archiveRoot, progressPath, execution);
    if (isLatticeSolve(job))
        return executeLatticeSolve(job, archiveRoot, progressPath, execution);
    if (isPDESolve(job))
        return executePDESolve(job, archiveRoot, progressPath, execution);
    return executeSolve(job, archiveRoot, progressPath, execution);
}

}  // namespace dynsyslab::workflow::detail
