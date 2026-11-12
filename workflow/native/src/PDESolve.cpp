#include <dynsyslab/workflow/PDESolve.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/PDEParallel.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Hierarchy.hpp>

#include "WorkflowSupport.hpp"
#include "AnalyticArchive.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

dynsyslab::core::SystemConfiguration requestedConfiguration(const json::Value& job) {
    dynsyslab::core::SystemConfiguration result;
    result.identifier = job.at("system").at("name").string();
    result.preset = job.at("system").at("preset").string();
    for (const auto& item : job.at("system").at("parameters").object())
        result.parameters[item.first] = item.second.number();
    return result;
}

std::size_t positiveInteger(const json::Value& value, const char* name) {
    const double number = value.number(1.0);
    if (!std::isfinite(number) || number < 1.0 || std::floor(number) != number
        || number > static_cast<double>((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return static_cast<std::size_t>(number);
}

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

std::uint64_t randomSeed(const json::Value& specification) {
    const double value = specification.at("seed").number();
    if (!std::isfinite(value) || value < 0.0 || value > 9007199254740991.0
        || std::floor(value) != value)
        throw std::invalid_argument("PDE seed must be a nonnegative integer");
    return static_cast<std::uint64_t>(value);
}

std::vector<std::vector<double>> initialFields(
    const json::Value& job, const dynsyslab::core::SystemConfiguration& configuration) {
    const auto& spec = job.at("initial_conditions");
    const std::size_t count = positiveInteger(spec.at("count"), "PDE trajectory count");
    const auto base = dynsyslab::core::heatEquationInitialField(configuration);
    const std::string mode = spec.at("mode").string("default");
    if (mode == "default") return std::vector<std::vector<double>>(count, base);
    if (mode == "explicit") {
        const auto& source = spec.at("values").array();
        std::vector<std::vector<double>> result;
        if (count == 1 && source.size() == base.size() && source.front().array().empty())
            result.push_back(numbers(spec.at("values")));
        else
            for (const auto& row : source) result.push_back(numbers(row));
        if (result.size() != count) throw std::invalid_argument("explicit PDE IC count mismatch");
        for (const auto& field : result)
            if (field.size() != base.size())
                throw std::invalid_argument("explicit PDE field size mismatch");
        return result;
    }
    if (mode != "uniform" && mode != "normal")
        throw std::invalid_argument("PDE IC mode must be default, explicit, uniform, or normal");
    const auto& requested = spec.at("ranges").array();
    if (requested.empty()) throw std::invalid_argument("random PDE ICs require ranges");
    std::vector<std::vector<double>> ranges;
    if (requested.size() == 2 && requested.front().array().empty())
        ranges.assign(base.size(), numbers(spec.at("ranges")));
    else {
        for (const auto& row : requested) ranges.push_back(numbers(row));
        if (ranges.size() == 1) ranges.resize(base.size(), ranges.front());
    }
    if (ranges.size() != base.size())
        throw std::invalid_argument("PDE IC ranges require one shared or one per-cell pair");
    for (const auto& bounds : ranges)
        if (bounds.size() != 2 || !std::isfinite(bounds[0]) || !std::isfinite(bounds[1])
            || bounds[0] > bounds[1])
            throw std::invalid_argument("invalid PDE IC range");
    std::mt19937_64 generator(randomSeed(spec));
    std::vector<std::vector<double>> result(count, std::vector<double>(base.size()));
    for (auto& field : result)
        for (std::size_t index = 0; index < field.size(); ++index) {
            const auto bounds = ranges[index];
            if (mode == "uniform") {
                std::uniform_real_distribution<double> distribution(bounds[0], bounds[1]);
                field[index] = distribution(generator);
            } else {
                const double mean = 0.5 * (bounds[0] + bounds[1]);
                const double sigma = (bounds[1] - bounds[0]) / 6.0;
                std::normal_distribution<double> distribution(mean, sigma > 0.0 ? sigma : 1.0);
                field[index] = sigma > 0.0
                    ? std::clamp(distribution(generator), bounds[0], bounds[1]) : mean;
            }
        }
    return result;
}

dynsyslab::core::PDEIntegrationOptions integrationOptions(
    const json::Value& job, const dynsyslab::core::SystemConfiguration& configuration) {
    dynsyslab::core::PDEIntegrationOptions result;
    result.initialTime = job.at("time").at("start").number();
    result.finalTime = job.at("time").at("stop").number();
    result.outputStep = job.at("time").at("output_step").number();
    std::string method = job.at("solver").at("method").string("RK4");
    std::transform(method.begin(), method.end(), method.begin(), [](unsigned char item) {
        return static_cast<char>(std::tolower(item));
    });
    if (method == "euler") result.method = dynsyslab::core::PDEIntegrationMethod::euler;
    else if (method == "rk4") result.method = dynsyslab::core::PDEIntegrationMethod::rk4;
    else throw std::invalid_argument("native PDE workflow supports Euler and RK4");
    const double stable = dynsyslab::core::heatEquationStableStep(configuration);
    result.internalStep = job.at("solver").at("step_size").number(
        std::min(result.outputStep, stable));
    if (result.internalStep > stable * (1.0 + 1.0e-12))
        throw std::invalid_argument(
            "heat-equation step_size exceeds the conservative explicit stability limit");
    return result;
}

std::optional<dynsyslab::core::ParallelMeshPlan> parallelPlan(
    const json::Value& job, const dynsyslab::core::StructuredMesh& mesh) {
    const auto& canonical = job.at("system").at("pde_partitioning");
    const auto& legacy = job.at("solver").at("parallel");
    const auto& source = canonical.isNull() ? legacy : canonical;
    if (source.isNull() || !source.at("enabled").boolean(false)) return {};
    std::vector<std::size_t> shape;
    for (const auto& value : source.at("partition_shape").array())
        shape.push_back(positiveInteger(value, "PDE partition extent"));
    if (shape.empty())
        throw std::invalid_argument("parallel PDE solve requires partition_shape");
    return dynsyslab::core::buildParallelMeshPlan(
        mesh, std::move(shape), positiveInteger(source.at("ghost_width"),
        "PDE ghost width"), positiveInteger(source.at("local_workers"),
        "PDE local worker count"));
}

std::string metadata(const json::Value& job,
                     const dynsyslab::core::SystemConfiguration& configuration,
                     const HierarchyIdentity& hierarchy, std::size_t trajectory,
                     const std::vector<double>& coordinates,
                     const dynsyslab::core::PDEIntegrationOptions& options,
                     const dynsyslab::core::ParallelMeshPlan* parallel,
                     const std::string& created) {
    json::Value::Object root;
    root["format"] = "dynsyslab.solution/v2";
    root["label"] = job.at("label").string("pde-solve");
    root["system"] = configuration.identifier;
    root["system_family"] = "pde";
    root["preset"] = configuration.preset;
    root["system_preset"] = configuration.preset;
    root["coordinate_system"] = "dimensional";
    root["integration_form"] = "dimensional";
    root["execution_backend"] = "native-cpp";
    root["compiled"] = true;
    root["python_required"] = false;
    root["equations_version"] = configuration.equationsVersion;
    root["equations_fingerprint"] = dynsyslab::core::equationsFingerprint(configuration);
    root["hierarchy_manifest"] = hierarchy.manifestPath.generic_u8string();
    root["campaign_id"] = hierarchy.campaignId;
    root["campaign_kind"] = hierarchy.campaignKind;
    root["campaign_point_index"] = job.at("campaign").at("point_index").number(-1.0);
    root["controller_configuration"] = job.at("controller");
    root["batch_id"] = hierarchy.batchId;
    root["system_run_id"] = hierarchy.systemRunId;
    root["parameter_set_key"] = hierarchy.parameterSetKey;
    root["cluster_id"] = hierarchy.clusterIds[
        hierarchy.trajectoryClusters[trajectory]];
    root["trajectory_id"] = hierarchy.trajectoryIds[trajectory];
    root["data_role"] = hierarchy.dataRole;
    root["trajectory_index"] = static_cast<double>(trajectory);
    root["created_utc"] = created;
    root["pde_method"] = options.method == dynsyslab::core::PDEIntegrationMethod::rk4
        ? "RK4" : "Euler";
    root["pde_internal_step"] = options.internalStep;
    root["pde_execution_mode"] = parallel ? "local-partitioned" : "serial";
    if (parallel) {
        json::Value::Array partitionShape;
        for (const auto extent : parallel->partitionShape)
            partitionShape.emplace_back(static_cast<double>(extent));
        root["pde_partition_shape"] = json::Value(std::move(partitionShape));
        root["pde_ghost_width"] = static_cast<double>(parallel->ghostWidth);
        root["pde_local_workers"] = static_cast<double>(parallel->localWorkers);
    }
    json::Value::Object parameters;
    for (const auto& item : configuration.parameters) parameters[item.first] = item.second;
    root["resolved_parameters"] = json::Value(parameters);
    root["resolved_system_parameters"] = json::Value(std::move(parameters));
    json::Value::Object units;
    for (const auto& item : configuration.parameterUnits) units[item.first] = item.second;
    root["parameter_units"] = json::Value(std::move(units));
    json::Value::Array shape;
    shape.emplace_back(static_cast<double>(coordinates.size()));
    root["field_shape"] = json::Value(std::move(shape));
    return json::dump(json::Value(std::move(root)));
}

}  // namespace

bool isPDESolve(const json::Value& job) {
    return job.at("system").at("family").string() == "pde"
        || job.at("system").at("name").string() == "heat_equation_1d";
}

SolveResult executePDESolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution) {
    if (job.at("controller").at("kind").string("none") != "none")
        throw std::invalid_argument("heat-equation control is not defined");
    const auto requested = requestedConfiguration(job);
    const auto resolved = dynsyslab::core::resolveHeatEquation(requested);
    const auto fields = initialFields(job, resolved);
    const auto options = integrationOptions(job, resolved);
    const auto mesh = dynsyslab::core::heatEquationMesh(resolved);
    const auto parallel = parallelPlan(job, mesh);
    if (parallel && execution && execution->limits.maximumThreads
        && parallel->localWorkers > execution->limits.maximumThreads)
        throw dynsyslab::core::ExecutionLimitExceeded(
            "parallel PDE local_workers exceeds the execution thread limit");
    auto rhs = dynsyslab::core::heatEquationRightHandSide(resolved);
    if (parallel) rhs = dynsyslab::core::partitionedPDERightHandSide(
        *parallel, dynsyslab::core::heatEquationSubdomainRightHandSide(resolved));
    const auto coordinates = dynsyslab::core::heatEquationCoordinates(resolved);
    const std::string created = detail::timestamp();
    const auto hierarchy = prepareHierarchy(job, resolved, fields.size(), archiveRoot, created);
    detail::atomicProgress(progressPath, 0.0, 0, fields.size(), "Starting PDE integration");
    SolveResult result;
    result.label = job.at("label").string("pde-solve");
    result.system = resolved.identifier;
    result.trajectories = fields.size();
    const auto directory = archiveRoot / L"solutions" / L"heat_equation_1d";
    std::vector<std::filesystem::path> archives;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        dynsyslab::core::ExecutionContext local;
        if (execution) {
            local.limits = execution->limits;
            local.cancellation = execution->cancellation;
        }
        local.progress = [execution, progressPath, index, count = fields.size()](
                             dynsyslab::core::ExecutionProgress event) {
            const double fraction = (index + std::clamp(event.fraction, 0.0, 1.0)) / count;
            detail::atomicProgress(progressPath, fraction, index, count,
                "PDE trajectory " + std::to_string(index + 1) + " of "
                    + std::to_string(count));
            if (execution) {
                event.fraction = fraction;
                event.completed = index;
                event.planned = count;
                execution->report(std::move(event));
            }
        };
        const auto trajectory = dynsyslab::core::integratePDE(
            mesh, {coordinates.size()}, fields[index], options,
            rhs, &local);
        dynsyslab::core::NpzRecord archive;
        archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
            {trajectory.sampleCount}, trajectory.time);
        archive.arrays["space"] = dynsyslab::core::NpyArray::fromDoubles(
            {coordinates.size()}, coordinates);
        archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
            {trajectory.sampleCount, coordinates.size()}, trajectory.values);
        archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
            metadata(job, resolved, hierarchy, index, coordinates,
                     options, parallel ? &*parallel : nullptr, created));
        detail::appendAnalyticComparison(archive,
            dynsyslab::core::heatEquationDescriptor(coordinates.size()), resolved,
            dynsyslab::core::heatEquationAnalyticSolutions(resolved), trajectory.time,
            trajectory.values, fields[index], options.initialTime, job);
        const std::string label = result.label + (fields.size() == 1 ? ""
            : "-trajectory-" + std::to_string(index));
        const auto path = directory / dynsyslab::core::wide(
            detail::safeName(label, "pde-solve") + "__" + created + ".npz");
        dynsyslab::core::writeNpzAtomicNew(path, archive);
        archives.push_back(path);
        result.samplesPerTrajectory = trajectory.sampleCount;
        result.acceptedSteps += trajectory.acceptedSteps;
    }
    auto descriptor = dynsyslab::core::heatEquationDescriptor(coordinates.size());
    std::vector<double> flattened;
    flattened.reserve(fields.size() * coordinates.size());
    for (const auto& field : fields)
        flattened.insert(flattened.end(), field.begin(), field.end());
    publishHierarchy(hierarchy, job, requested, descriptor, resolved,
                     flattened, archives, false, false, created);
    result.archives = std::move(archives);
    result.hierarchyManifest = hierarchy.manifestPath;
    detail::atomicProgress(progressPath, 1.0, fields.size(), fields.size(),
                           "Archived PDE trajectories");
    return result;
}

}  // namespace dynsyslab::workflow
