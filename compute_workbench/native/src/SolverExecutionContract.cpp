#include "dynsyslab_compute/SolverExecutionContract.hpp"

#include <dynsyslab/core/Json.hpp>

#include <stdexcept>

namespace dsl::compute {

std::string applySolverExecutionContract(const std::string& jobJson,
                                         const std::string& solverJson,
                                         const std::string& campaignJson) {
    namespace json = dynsyslab::core::json;
    auto root = json::parse(jobJson);
    const auto solver = json::parse(solverJson);
    const auto campaign = json::parse(campaignJson);
    auto apply = [&](json::Value& solve) {
        const auto evaluation = solve.at("solver").at("control_evaluation");
        const auto sampleTime = solve.at("solver").at("control_sample_time");
        const auto legacyPartition = solve.at("solver").at("parallel");
        if (solve.at("system").at("pde_partitioning").isNull()
            && !legacyPartition.isNull())
            solve.at("system").at("pde_partitioning") = legacyPartition;
        solve.at("solver") = solver;
        if (!evaluation.isNull())
            solve.at("solver").at("control_evaluation") = evaluation;
        if (!sampleTime.isNull())
            solve.at("solver").at("control_sample_time") = sampleTime;
        solve.at("campaign") = campaign;
    };
    const std::string kind = root.at("kind").string();
    if (kind == "solve") apply(root);
    else if (kind == "solve_set")
        for (auto& batch : root.at("batches").array())
            for (auto& solve : batch.at("systems").array()) apply(solve);
    else throw std::invalid_argument(
        "Job-level execution settings require a solve or solve_set job");
    return json::dump(root);
}

}  // namespace dsl::compute
