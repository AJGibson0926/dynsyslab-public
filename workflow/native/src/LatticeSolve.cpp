#include <dynsyslab/workflow/LatticeSolve.hpp>
#include <dynsyslab/workflow/Hierarchy.hpp>
#include "WorkflowSupport.hpp"
#include "SolveArchive.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <variant>

namespace dynsyslab::workflow {
namespace {

struct LatticeTrajectory {
    std::vector<double> observations;
    std::vector<double> initialObservation;
    std::vector<std::int64_t> initialSites;
    std::vector<std::int64_t> finalSites;
    std::vector<std::size_t> latticeShape;
    std::size_t acceptedMoves{};
};

dynsyslab::core::SystemConfiguration requestedConfiguration(const json::Value& job) {
    dynsyslab::core::SystemConfiguration result;
    const auto& source = job.at("system");
    result.identifier = source.at("name").string();
    result.preset = source.at("preset").string();
    for (const auto& item : source.at("parameters").object())
        result.parameters[item.first] = item.second.number();
    return result;
}

std::size_t integerParameter(const dynsyslab::core::SystemConfiguration& configuration,
                             const char* name, std::size_t minimum) {
    const double value = configuration.parameters.at(name);
    if (!std::isfinite(value) || value < static_cast<double>(minimum)
        || value > static_cast<double>((std::numeric_limits<std::size_t>::max)())
        || std::floor(value) != value)
        throw std::invalid_argument(std::string("lattice parameter must be an integer: ") + name);
    return static_cast<std::size_t>(value);
}

std::uint64_t randomSeed(const json::Value& job) {
    const double value = job.at("initial_conditions").at("seed").number();
    if (!std::isfinite(value) || value < 0.0
        || value > 9007199254740991.0
        || std::floor(value) != value)
        throw std::invalid_argument("lattice seed must be a nonnegative integer");
    return static_cast<std::uint64_t>(value);
}

std::size_t siteCount(std::size_t size, std::size_t dimensions) {
    std::size_t result = 1;
    for (std::size_t index = 0; index < dimensions; ++index) {
        if (result > (std::numeric_limits<std::size_t>::max)() / size)
            throw dynsyslab::core::ExecutionLimitExceeded(
                "lattice site count is not addressable");
        result *= size;
    }
    return result;
}

void flattenSites(const json::Value& value, std::vector<std::int64_t>& result) {
    if (std::holds_alternative<json::Value::Array>(value.storage())) {
        for (const auto& item : value.array()) flattenSites(item, result);
        return;
    }
    if (!std::holds_alternative<double>(value.storage()))
        throw std::invalid_argument("explicit lattice sites must be numeric arrays");
    const double number = value.number();
    if (!std::isfinite(number) || std::floor(number) != number
        || number < static_cast<double>((std::numeric_limits<std::int64_t>::min)())
        || number > static_cast<double>((std::numeric_limits<std::int64_t>::max)()))
        throw std::invalid_argument("explicit lattice sites must contain finite integers");
    result.push_back(static_cast<std::int64_t>(number));
}

std::vector<std::vector<std::int64_t>> initialSites(
    const json::Value& job, std::size_t count, std::size_t sites) {
    const auto& specification = job.at("initial_conditions");
    const std::string mode = specification.at("mode").string("default");
    if (mode == "default" || mode == "uniform")
        return std::vector<std::vector<std::int64_t>>(count);
    if (mode != "explicit")
        throw std::invalid_argument(
            "lattice IC mode must be default, uniform, or explicit");
    const auto& values = specification.at("values");
    std::vector<std::vector<std::int64_t>> result(count);
    if (count == 1) {
        flattenSites(values, result.front());
    } else {
        if (!std::holds_alternative<json::Value::Array>(values.storage())
            || values.array().size() != count)
            throw std::invalid_argument(
                "multiple explicit lattices require one outer row per trajectory");
        for (std::size_t trajectory = 0; trajectory < count; ++trajectory)
            flattenSites(values.array()[trajectory], result[trajectory]);
    }
    for (const auto& item : result)
        if (item.size() != sites)
            throw std::invalid_argument("explicit lattice site count mismatch");
    return result;
}

const dynsyslab::core::SystemDescriptor& descriptor(const std::string& identifier) {
    if (identifier == "ising") return dynsyslab::core::ising2DDescriptor();
    if (identifier == "ising_3d") return dynsyslab::core::ising3DDescriptor();
    return dynsyslab::core::potts2DDescriptor();
}

LatticeTrajectory runTrajectory(
    const dynsyslab::core::SystemConfiguration& configuration,
    std::vector<std::int64_t> sites, std::uint64_t seed, std::size_t sweeps,
    const dynsyslab::core::ExecutionContext* execution,
    const dynsyslab::core::LatticeCheckpoint& checkpoint = {}) {
    const std::size_t size = integerParameter(configuration, "size", 2);
    LatticeTrajectory result;
    if (configuration.identifier == "potts") {
        dynsyslab::core::PottsRequest request;
        request.spins = std::move(sites); request.size = size;
        request.states = integerParameter(configuration, "states", 2);
        request.coupling = configuration.parameters.at("coupling");
        request.field = configuration.parameters.at("field");
        request.temperature = configuration.parameters.at("temperature");
        request.seed = seed; request.sweeps = sweeps; request.checkpoint = checkpoint;
        auto output = dynsyslab::core::runPotts2D(request, execution);
        result.observations = std::move(output.observations);
        result.initialObservation = std::move(output.initialObservation);
        result.initialSites = std::move(output.initialSpins);
        result.finalSites = std::move(output.spins);
        result.acceptedMoves = output.acceptedMoves;
        result.latticeShape = {size, size};
        return result;
    }
    dynsyslab::core::IsingRequest request;
    request.spins = std::move(sites); request.size = size;
    request.dimensions = configuration.identifier == "ising_3d" ? 3 : 2;
    request.coupling = configuration.parameters.at("coupling");
    request.field = configuration.parameters.at("field");
    request.temperature = configuration.parameters.at("temperature");
    request.seed = seed; request.sweeps = sweeps; request.checkpoint = checkpoint;
    auto output = request.dimensions == 3
        ? dynsyslab::core::runIsing3D(request, execution)
        : dynsyslab::core::runIsing2D(request, execution);
    result.observations = std::move(output.observations);
    result.initialObservation = std::move(output.initialObservation);
    result.initialSites = std::move(output.initialSpins);
    result.finalSites = std::move(output.spins);
    result.acceptedMoves = output.acceptedMoves;
    result.latticeShape.assign(request.dimensions, size);
    return result;
}

std::string metadata(
    const json::Value& job, const dynsyslab::core::SystemConfiguration& resolved,
    const dynsyslab::core::SystemDescriptor& description,
    const HierarchyIdentity& hierarchy, const LatticeTrajectory& trajectory,
    std::size_t index, std::uint64_t seed, double start, std::size_t sweeps,
    const std::string& created) {
    std::ostringstream out;
    out << "{\"schema\":\"dynsyslab.solution/v1\",\"label\":\""
        << detail::jsonEscape(job.at("label").string("lattice-solve"))
        << "\",\"system\":\"" << resolved.identifier << "\",\"system_label\":\""
        << detail::jsonEscape(description.label)
        << "\",\"system_family\":\"lattice\",\"system_preset\":\""
        << detail::jsonEscape(resolved.preset) << "\",\"resolved_system_parameters\":{";
    bool first = true;
    for (const auto& item : resolved.parameters) {
        if (!first) out << ',';
        first = false;
        out << '"' << detail::jsonEscape(item.first) << "\":"
            << std::setprecision(17) << item.second;
    }
    const long double attempts = static_cast<long double>(sweeps)
        * static_cast<long double>(trajectory.initialSites.size());
    const double acceptance = attempts > 0.0L
        ? static_cast<double>(trajectory.acceptedMoves / attempts) : 0.0;
    out << "},\"equations_version\":\"" << resolved.equationsVersion
        << "\",\"equations_fingerprint\":\""
        << dynsyslab::core::equationsFingerprint(resolved)
        << "\",\"integration_form\":\"discrete\",\"solver_family\":\"lattice\""
        << ",\"method\":\"metropolis-single-site\",\"sweep_start\":" << start
        << ",\"sweeps\":" << sweeps << ",\"random_engine\":\"std::mt19937_64\""
        << ",\"random_seed\":" << seed << ",\"accepted_moves\":"
        << trajectory.acceptedMoves << ",\"acceptance_fraction\":" << acceptance
        << ",\"trajectory_parity\":\"statistical-not-bitwise-with-numpy-reference\""
        << ",\"execution_backend\":{\"name\":\"native-cpp\",\"compiled\":true,"
           "\"python_required\":false},\"hierarchy_manifest\":\""
        << detail::jsonEscape(hierarchy.manifestPath.generic_u8string())
        << "\",\"campaign_id\":\"" << detail::jsonEscape(hierarchy.campaignId)
        << "\",\"campaign_kind\":\"" << detail::jsonEscape(hierarchy.campaignKind)
        << "\",\"campaign_point_index\":"
        << job.at("campaign").at("point_index").number(-1.0)
        << ",\"controller_configuration\":" << json::dump(job.at("controller"))
        << ",\"batch_id\":\"" << hierarchy.batchId << "\",\"system_run_id\":\""
        << hierarchy.systemRunId << "\",\"parameter_set_key\":\""
        << hierarchy.parameterSetKey << "\",\"cluster_id\":\""
        << hierarchy.clusterIds[hierarchy.trajectoryClusters[index]]
        << "\",\"trajectory_id\":\"" << hierarchy.trajectoryIds[index]
        << "\",\"data_role\":\"" << hierarchy.dataRole
        << "\",\"initial_state\":[" << trajectory.initialObservation[0] << ','
        << trajectory.initialObservation[1] << "],\"trajectory_index\":" << index
        << ",\"created_utc\":\"" << created << "\"}";
    return out.str();
}

}  // namespace

bool isLatticeSolve(const json::Value& job) {
    const std::string family = job.at("system").at("family").string();
    return family == "lattice" || dynsyslab::core::isLatticeSystemIdentifier(
        job.at("system").at("name").string());
}

SolveResult executeLatticeSolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution) {
    if (job.at("controller").at("kind").string("none") != "none")
        throw std::invalid_argument("native lattice control is not defined");
    const auto requested = requestedConfiguration(job);
    const auto resolved = dynsyslab::core::resolveLatticeSystem(requested);
    const auto& description = descriptor(resolved.identifier);
    const double start = job.at("time").at("start").number();
    const double stop = job.at("time").at("stop").number();
    const double duration = stop - start;
    if (!std::isfinite(start) || !std::isfinite(stop) || duration < 1.0
        || std::floor(duration) != duration
        || duration > static_cast<double>((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument(
            "lattice time span must contain a positive integer number of sweeps");
    const std::size_t sweeps = static_cast<std::size_t>(duration);
    const double requestedCount = job.at("initial_conditions").at("count").number(1.0);
    if (!std::isfinite(requestedCount) || requestedCount < 1.0
        || requestedCount > static_cast<double>((std::numeric_limits<std::size_t>::max)())
        || std::floor(requestedCount) != requestedCount)
        throw std::invalid_argument("lattice trajectory count must be a positive integer");
    const std::size_t count = static_cast<std::size_t>(requestedCount);
    const std::size_t size = integerParameter(resolved, "size", 2);
    const std::size_t dimensions = resolved.identifier == "ising_3d" ? 3 : 2;
    auto sites = initialSites(job, count, siteCount(size, dimensions));
    const std::uint64_t seed = randomSeed(job);
    const std::string created = detail::timestamp();
    const auto hierarchy = prepareHierarchy(job, resolved, count, archiveRoot, created);
    detail::atomicProgress(progressPath, 0.0, 0, count, "Starting lattice simulation");
    std::vector<double> initialObservations;
    std::vector<std::filesystem::path> archives;
    const auto directory = archiveRoot / L"solutions" / dynsyslab::core::wide(resolved.identifier);
    for (std::size_t index = 0; index < count; ++index) {
        dynsyslab::core::ExecutionContext local;
        const dynsyslab::core::ExecutionContext* active = execution;
        if (execution) {
            local.limits = execution->limits;
            local.cancellation = execution->cancellation;
            local.progress = [execution, index, count](dynsyslab::core::ExecutionProgress event) {
                event.fraction = (index + event.fraction) / count;
                event.completed = index;
                event.planned = count;
                execution->report(std::move(event));
            };
            active = &local;
        }
        const std::uint64_t trajectorySeed = seed + index;
        const std::string label = job.at("label").string("lattice-solve")
            + (count == 1 ? "" : "-trajectory-" + std::to_string(index));
        const auto path = directory / dynsyslab::core::wide(
            detail::safeName(label, "lattice-solve") + "__" + created + ".npz");
        // history.checkpoint_interval (sweeps): observations so far and the current lattice go to
        // <archive>.partial.npz, removed once the trajectory is archived.
        dynsyslab::core::LatticeCheckpoint checkpoint;
        const auto partial = detail::partialArchivePath(path);
        if (const double interval = detail::checkpointInterval(job); interval > 0.0) {
            checkpoint.interval = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(interval)));
            checkpoint.callback = [&, index](std::size_t done, const std::vector<double>& observations,
                                             const std::vector<std::int64_t>& spins) {
                std::vector<double> times(done);
                for (std::size_t sweep = 0; sweep < done; ++sweep) times[sweep] = start + sweep;
                const std::size_t rank = resolved.identifier == "ising_3d" ? 3 : 2;
                detail::writePartialSeries(partial, times.data(), done, observations.data(), 2,
                    json::dump(json::Value(json::Value::Object{
                        {"schema", std::string("dynsyslab.solution-checkpoint/v1")},
                        {"trajectory_index", static_cast<double>(index)},
                        {"completed_sweeps", static_cast<double>(done)}})),
                    {{"lattice_current", dynsyslab::core::NpyArray::fromInt64s(
                        std::vector<std::size_t>(rank, size), spins)}});
            };
        }
        auto trajectory = runTrajectory(
            resolved, std::move(sites[index]), trajectorySeed, sweeps, active, checkpoint);
        initialObservations.insert(initialObservations.end(),
            trajectory.initialObservation.begin(), trajectory.initialObservation.end());
        std::vector<double> time(sweeps);
        for (std::size_t sweep = 0; sweep < sweeps; ++sweep) time[sweep] = start + sweep;
        dynsyslab::core::NpzRecord archive;
        archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({sweeps}, time);
        archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
            {sweeps, 2}, trajectory.observations);
        archive.arrays["lattice_initial"] = dynsyslab::core::NpyArray::fromInt64s(
            trajectory.latticeShape, trajectory.initialSites);
        archive.arrays["lattice_final"] = dynsyslab::core::NpyArray::fromInt64s(
            trajectory.latticeShape, trajectory.finalSites);
        archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
            metadata(job, resolved, description, hierarchy, trajectory, index,
                     trajectorySeed, start, sweeps, created));
        dynsyslab::core::writeNpzAtomic(path, archive);
        { std::error_code ignored; std::filesystem::remove(dynsyslab::core::extendedLengthPath(partial), ignored); }
        archives.push_back(path);
        detail::atomicProgress(progressPath, static_cast<double>(index + 1) / count,
            index + 1, count, "Archived lattice trajectory");
    }
    publishHierarchy(hierarchy, job, requested, description, resolved,
        initialObservations, archives, false, false, created);
    SolveResult result;
    result.label = job.at("label").string("lattice-solve");
    result.system = resolved.identifier;
    result.trajectories = count;
    result.samplesPerTrajectory = sweeps;
    result.acceptedSteps = sweeps * count;
    result.archives = std::move(archives);
    result.hierarchyManifest = hierarchy.manifestPath;
    return result;
}

}  // namespace dynsyslab::workflow
