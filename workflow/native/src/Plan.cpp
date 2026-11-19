#include <dynsyslab/workflow/Plan.hpp>

#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/Campaign.hpp>
#include <dynsyslab/workflow/DelaySolve.hpp>
#include <dynsyslab/workflow/LatticeSolve.hpp>
#include <dynsyslab/workflow/KromControlSearch.hpp>
#include <dynsyslab/workflow/PDESolve.hpp>
#include <dynsyslab/workflow/Pipeline.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>

#include "WorkflowExtensions.hpp"
#include "TrainingPreflight.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace dynsyslab::workflow {
namespace {

std::size_t positiveInteger(const json::Value& value, std::size_t fallback,
                            const char* field) {
    const double number = value.number(static_cast<double>(fallback));
    if (number < 1.0 || number != std::floor(number))
        throw std::invalid_argument(std::string(field) + " must be a positive integer");
    return static_cast<std::size_t>(number);
}

dynsyslab::core::SystemConfiguration configuration(const json::Value& job) {
    const auto& source = job.at("system");
    dynsyslab::core::SystemConfiguration result;
    result.identifier = source.at("name").string();
    result.preset = source.at("preset").string();
    if (result.identifier.empty()) throw std::invalid_argument("solve requires system.name");
    for (const auto& parameter : source.at("parameters").object())
        result.parameters.emplace(parameter.first, parameter.second.number());
    return result;
}

struct PlannedSystem {
    dynsyslab::core::SystemDescriptor descriptor;
    dynsyslab::core::CoordinateScaling scaling;
};

PlannedSystem plannedSystem(const json::Value& job,
                            const dynsyslab::core::SystemConfiguration& requested) {
    (void)job;
    const auto resolved=resolveSystemDefinition(requested);
    return {resolved.descriptor,resolved.scaling};
}

void validateController(const json::Value& job) {
    static const std::set<std::string> supported{
        "none", "waveform", "random_pulse", "pid", "lqr", "sdre",
        "mpc", "klqr", "kmpc",
    };
    const std::string kind = job.at("controller").at("kind").string("none");
    if(!job.at("controller").at("settings").at("estimator").isNull()&&kind!="klqr")
        throw std::invalid_argument("Online KLQE composition requires KLQR");
    if (!supported.count(kind))
        throw std::invalid_argument("unknown native controller kind: " + kind);
    if (!job.at("controller").at("settings").at("applied_controller").isNull()
        && (kind != "klqr" || (job.at("controller").at("settings").at("feedback").string() != "estimate"
                               && job.at("controller").at("settings").at("estimator").isNull())))
        throw std::invalid_argument("applied_controller (shadow KLQR) requires KLQR with estimate feedback");
    const auto& estimator = job.at("estimator");
    if (estimator.isNull() && job.at("controller").at("settings").at("feedback").string() == "estimated_state")
        throw std::invalid_argument("estimated_state feedback needs an estimator");
    if (!estimator.isNull()) {
        if (estimator.at("kind").string("klqe") != "klqe")
            throw std::invalid_argument("unknown native estimator kind: " + estimator.at("kind").string());
        if (job.at("controller").at("settings").at("feedback").string() == "estimate" && kind != "klqr")
            throw std::invalid_argument("estimate feedback (KLQE -> controller) requires KLQR");
        if (job.at("controller").at("settings").at("feedback").string() == "estimated_state"
            && estimator.at("output").string().empty())
            throw std::invalid_argument("estimated_state feedback needs an estimator output reconstruction (e.g. state)");
    }
}

bool supportedMethod(const json::Value& stage) {
    static const std::set<std::string> supported{
        "svd", "pca", "pod", "bagged_pca", "robust_pca", "hankel_pod", "hankel_svd", "hankel_matrix",
        "dmd", "dmdc", "edmdc", "hdmd", "hankel_dmd", "sindy", "edmd", "havok",
        "koopman_spectral", "resdmd", "mpedmd", "hermitian_dmd", "tlsdmd", "multdmd", "stochastic_resdmd", "specrkhs",
        "control_regression", "krom_extract", "koopman_validation", "observable_reconstruction",
        "echo_history_koopman", "history_basis", "history_projection", "history_control_regression", "input_map", "control_target", "derivative_benchmark", "mode_selection",
        "coordinate_transform", "sindy_dynamics", "dynamics_validation",
        "boundary_eigenfunctions", "eigenvalue_optimization",
        "observability", "klqe", "linear_analysis", "linear_systems", "controllability", "empirical_observability",
    };
    const std::string method = stage.at("name").string();
    return supported.count(method) || !stage.at("extension").string().empty();
}

json::Value solvePlan(const json::Value& job) {
    const auto requested = configuration(job);
    const auto system = plannedSystem(job, requested);
    const auto& time = job.at("time");
    const double start = time.at("start").number();
    const double stop = time.at("stop").number();
    const double step = time.at("output_step").number();
    if (!(start < stop) || !(step > 0.0))
        throw std::invalid_argument("solve time requires start < stop and positive output_step");
    const double intervals = (stop - start) / step;
    if (std::abs(intervals - std::round(intervals))
        > 1e-9 * std::max(1.0, std::abs(intervals)))
        throw std::invalid_argument("solve interval must be an integer multiple of output_step");
    const std::string coordinates = time.at("coordinate_system").string("automatic");
    if (coordinates != "automatic"
        && coordinates != system.descriptor.preferredCoordinates)
        throw std::invalid_argument("requested coordinate system is unavailable for this model");
    detail::validateTrainingPreflightPlan(
        job, system.descriptor.stateDimension,
        system.descriptor.preferredCoordinates);
    const auto& initial = job.at("initial_conditions");
    const std::size_t trajectories = positiveInteger(initial.at("count"), 1,
        "initial_conditions.count");
    const std::size_t clusters = positiveInteger(initial.at("clusters"), 1,
        "initial_conditions.clusters");
    if (clusters > trajectories)
        throw std::invalid_argument("initial_conditions.clusters exceeds trajectory count");
    const std::string initialMode = initial.at("mode").string("default");
    static const std::set<std::string> continuousModes{
        "default", "explicit", "uniform", "normal"};
    static const std::set<std::string> latticeModes{"default", "explicit", "uniform"};
    const auto& modes = system.descriptor.family == "lattice"
        ? latticeModes : continuousModes;
    if (!modes.count(initialMode))
        throw std::invalid_argument("unknown initial-condition mode");
    validateController(job);
    if (system.descriptor.family != "ode"
        && job.at("controller").at("kind").string("none") != "none")
        throw std::invalid_argument("controllers currently require an ODE system");
    const std::size_t systems = campaignPointCount(job);
    const std::size_t samples = static_cast<std::size_t>(std::llround(intervals)) + 1;
    const std::size_t dimension = system.descriptor.stateDimension;
    const std::string controllerKind =
        job.at("controller").at("kind").string("none");
    const bool controlled = controllerKind != "none";
    const bool hasTarget = controlled && controllerKind != "waveform"
        && controllerKind != "random_pulse";
    const std::size_t scalarsPerSample = 1 + dimension
        + (controlled ? 1 : 0) + (hasTarget ? dimension : 0);
    const double perSystemBytes = static_cast<double>(samples) * trajectories
        * scalarsPerSample * sizeof(double);

    json::Value::Object result;
    result["schema"] = "dynsyslab.workflow-plan/v1";
    result["status"] = "planned";
    result["kind"] = "solve";
    result["system"] = system.descriptor.identifier;
    result["system_family"] = system.descriptor.family;
    result["coordinate_system"] = system.descriptor.preferredCoordinates;
    result["model_time_units"] = system.scaling.modelTimeUnits;
    result["controller"] = controllerKind;
    result["has_control_target"] = hasTarget;
    result["training_preflight"] =
        job.at("preflight").at("enabled").boolean(false);
    result["campaign"] = systems > 1;
    result["system_count"] = static_cast<double>(systems);
    result["trajectory_count"] = static_cast<double>(systems * trajectories);
    result["cluster_count"] = static_cast<double>(systems * clusters);
    result["clusters_per_system"] = static_cast<double>(clusters);
    result["state_dimension"] = static_cast<double>(dimension);
    result["output_samples_per_trajectory"] = static_cast<double>(samples);
    result["estimated_peak_solution_bytes"] = perSystemBytes;
    result["estimated_total_solution_bytes"] = perSystemBytes * systems;
    result["archive_policy"] = "immutable native NPZ plus hierarchy manifest";
    return json::Value(std::move(result));
}

json::Value solveSetPlan(const json::Value& job) {
    const auto& batches = job.at("batches").array();
    if (batches.empty()) throw std::invalid_argument("solve_set requires at least one batch");
    std::size_t configuredSystems = 0, expandedSystems = 0, trajectories = 0, clusters = 0;
    double peakBytes = 0.0, totalBytes = 0.0;
    json::Value::Array plans;
    for (const auto& batch : batches) {
        const auto& systems = batch.at("systems").array();
        if (systems.empty()) throw std::invalid_argument("Every solve_set batch requires a system");
        for (const auto& child : systems) {
            if (child.at("kind").string() != "solve")
                throw std::invalid_argument("solve_set systems must contain solve jobs");
            const auto plan = solvePlan(child);
            ++configuredSystems;
            expandedSystems += static_cast<std::size_t>(plan.at("system_count").number());
            trajectories += static_cast<std::size_t>(plan.at("trajectory_count").number());
            clusters += static_cast<std::size_t>(plan.at("cluster_count").number());
            peakBytes = std::max(peakBytes, plan.at("estimated_peak_solution_bytes").number());
            totalBytes += plan.at("estimated_total_solution_bytes").number();
            plans.push_back(plan);
        }
    }
    json::Value::Object result;
    result["schema"] = "dynsyslab.workflow-plan/v1";
    result["status"] = "planned";
    result["kind"] = "solve_set";
    result["batch_count"] = static_cast<double>(batches.size());
    result["configured_system_count"] = static_cast<double>(configuredSystems);
    result["system_count"] = static_cast<double>(expandedSystems);
    result["trajectory_count"] = static_cast<double>(trajectories);
    result["cluster_count"] = static_cast<double>(clusters);
    result["estimated_peak_solution_bytes"] = peakBytes;
    result["estimated_total_solution_bytes"] = totalBytes;
    result["systems"] = json::Value(std::move(plans));
    result["archive_policy"] = "one hierarchy manifest per proposed batch";
    return json::Value(std::move(result));
}

json::Value postprocessPlan(const json::Value& job) {
    const auto& sources = job.at("sources").array();
    if (sources.empty()) throw std::invalid_argument("post-processing requires a source");
    json::Value::Array methods;
    bool matrixAvailable=false;
    for (const auto& stage : job.at("methods").array()) {
        if (!stage.at("enabled").boolean(true)) continue;
        if (!supportedMethod(stage))
            throw std::invalid_argument("post-processing method is not native: "
                                        + stage.at("name").string());
        const auto name=stage.at("name").string();
        const auto& opts=stage.at("settings");
        if (opts.at("use_previous_hankel").boolean()) {
            if (!matrixAvailable || !opts.at("hankel_archive").isNull() ||
                (name!="hankel_svd" && name!="hankel_pod" && name!="hdmd" && name!="hankel_dmd"
                 && name!="koopman_spectral" && name!="resdmd" && name!="mpedmd" && name!="hermitian_dmd"
                 && name!="tlsdmd" && name!="multdmd" && name!="stochastic_resdmd" && name!="specrkhs"))
                throw std::invalid_argument("use_previous_hankel needs a preceding matrix stage and one Hankel decomposition");
        }
        if (name=="hankel_matrix") matrixAvailable=true;
        methods.emplace_back(name);
    }
    if (methods.empty()) throw std::invalid_argument("post-processing requires an enabled method");
    json::Value::Object result;
    result["schema"] = "dynsyslab.workflow-plan/v1";
    result["status"] = "planned";
    result["kind"] = "postprocess";
    result["source_count"] = static_cast<double>(sources.size());
    result["methods"] = json::Value(methods);
    const auto pooled = static_cast<std::size_t>(std::count_if(
        methods.begin(), methods.end(), [](const json::Value& method) {
            return method.string() == "echo_history_koopman" || method.string() == "history_basis";
        }));
    result["stage_count"] = static_cast<double>(
        sources.size() * (methods.size() - pooled) + pooled);
    result["archive_policy"] = "immutable native derived archive";
    return json::Value(std::move(result));
}

}  // namespace

json::Value planJob(const json::Value& job) {
    validateJob(job);
    const auto kind=job.at("kind").string();
    if(kind=="pipeline")return planPipeline(job);
    detail::loadExtensions(job);
    if(kind=="krom_control_search")return planKromControlSearch(job);
    if(kind=="solve_set")return solveSetPlan(job);
    return kind == "solve" ? solvePlan(job) : postprocessPlan(job);
}

}  // namespace dynsyslab::workflow
