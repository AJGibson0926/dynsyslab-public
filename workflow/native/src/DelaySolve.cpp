#include <dynsyslab/workflow/DelaySolve.hpp>
#include <dynsyslab/workflow/Hierarchy.hpp>
#include "WorkflowSupport.hpp"
#include "AnalyticArchive.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

dynsyslab::core::SystemConfiguration requestedConfiguration(const json::Value& job) {
    dynsyslab::core::SystemConfiguration result;
    const auto& source = job.at("system");
    result.identifier = "mackey_glass";
    result.preset = source.at("preset").string();
    for (const auto& item : source.at("parameters").object())
        result.parameters[item.first] = item.second.number();
    return result;
}

std::vector<double> numbers(const json::Value& source) {
    std::vector<double> result;
    for (const auto& value : source.array()) result.push_back(value.number());
    return result;
}

std::vector<double> initialHistories(const json::Value& job,
                                     const dynsyslab::core::SystemConfiguration& resolved) {
    const auto& specification = job.at("initial_conditions");
    const std::size_t count = static_cast<std::size_t>(
        std::max(1.0, specification.at("count").number(1.0)));
    const std::string mode = specification.at("mode").string("default");
    if (mode == "default")
        return std::vector<double>(count, resolved.parameters.at("history_value"));
    if (mode == "explicit") {
        const auto& rows = specification.at("values").array();
        std::vector<double> result;
        if (!rows.empty() && rows.front().array().empty()) {
            result = numbers(specification.at("values"));
        } else {
            result.reserve(rows.size());
            for (const auto& row : rows) {
                const auto values = numbers(row);
                if (values.size() != 1)
                    throw std::invalid_argument("delay explicit history dimension mismatch");
                result.push_back(values.front());
            }
        }
        if (result.size() != count)
            throw std::invalid_argument("delay explicit history count mismatch");
        for (const double value : result)
            if (!std::isfinite(value))
                throw std::invalid_argument("delay history value must be finite");
        return result;
    }
    const auto& ranges = specification.at("ranges").array();
    if (ranges.size() != 1)
        throw std::invalid_argument("delay IC range must contain one state range");
    const auto bounds = numbers(ranges.front());
    if (bounds.size() != 2 || !std::isfinite(bounds[0])
        || !std::isfinite(bounds[1]) || bounds[0] > bounds[1])
        throw std::invalid_argument("invalid delay history range");
    std::mt19937_64 generator(static_cast<std::uint64_t>(
        specification.at("seed").number()));
    std::uniform_real_distribution<double> uniform(bounds[0], bounds[1]);
    std::normal_distribution<double> normal(0.5 * (bounds[0] + bounds[1]),
                                             (bounds[1] - bounds[0]) / 6.0);
    std::vector<double> result(count);
    for (double& value : result)
        value = mode == "uniform" ? uniform(generator)
            : std::clamp(normal(generator), bounds[0], bounds[1]);
    if (mode != "uniform" && mode != "normal")
        throw std::invalid_argument("unknown delay history mode");
    return result;
}

std::string metadata(const json::Value& job,
                     const dynsyslab::core::SystemConfiguration& resolved,
                     const HierarchyIdentity& hierarchy, std::size_t trajectory,
                     double history, double step, const std::string& created) {
    const double delay = resolved.parameters.at("delay");
    const double effectiveDelay = std::llround(delay / step) * step;
    std::ostringstream out;
    out << "{\"format\":\"dynsyslab.solution/v1\",\"label\":\""
        << detail::jsonEscape(job.at("label").string("delay-solve"))
        << "\",\"system\":\"mackey_glass\",\"system_label\":\"Mackey-Glass delay system\""
        << ",\"system_family\":\"dde\",\"system_preset\":\""
        << detail::jsonEscape(resolved.preset) << "\",\"resolved_system_parameters\":{";
    bool first = true;
    for (const auto& item : resolved.parameters) {
        if (!first) out << ',';
        first = false;
        out << '"' << detail::jsonEscape(item.first) << "\":" << std::setprecision(17) << item.second;
    }
    out << "},\"parameter_units\":{";
    first = true;
    for (const auto& item : resolved.parameterUnits) {
        if (!first) out << ',';
        first = false;
        out << '"' << detail::jsonEscape(item.first) << "\":\""
            << detail::jsonEscape(item.second) << '"';
    }
    out << "},\"equations_version\":\"" << resolved.equationsVersion
        << "\",\"equations_fingerprint\":\""
        << dynsyslab::core::equationsFingerprint(resolved)
        << "\",\"integration_form\":\"dimensional\",\"solver_family\":\"dde\""
        << ",\"method\":\"method-of-steps-euler\",\"step_size\":" << step
        << ",\"requested_delay\":" << delay << ",\"effective_grid_delay\":"
        << effectiveDelay << ",\"history_function\":\"constant\""
        << ",\"execution_backend\":{\"name\":\"native-cpp\",\"compiled\":true,"
           "\"python_required\":false},\"hierarchy_manifest\":\""
        << detail::jsonEscape(hierarchy.manifestPath.generic_u8string())
        << "\",\"campaign_id\":\"" << detail::jsonEscape(hierarchy.campaignId)
        << "\",\"campaign_kind\":\"" << detail::jsonEscape(hierarchy.campaignKind)
        << "\",\"campaign_point_index\":"
        << job.at("campaign").at("point_index").number(-1.0)
        << ",\"controller_configuration\":" << json::dump(job.at("controller"))
        << ",\"batch_id\":\""
        << hierarchy.batchId << "\",\"system_run_id\":\"" << hierarchy.systemRunId
        << "\",\"parameter_set_key\":\"" << hierarchy.parameterSetKey
        << "\",\"cluster_id\":\""
        << hierarchy.clusterIds[hierarchy.trajectoryClusters[trajectory]]
        << "\",\"trajectory_id\":\"" << hierarchy.trajectoryIds[trajectory]
        << "\",\"data_role\":\"" << hierarchy.dataRole << "\",\"initial_state\":["
        << history << "],\"trajectory_index\":" << trajectory
        << ",\"created_utc\":\"" << created << "\"}";
    return out.str();
}

}  // namespace

bool isDelaySolve(const json::Value& job) {
    const std::string family = job.at("system").at("family").string();
    const std::string name = job.at("system").at("name").string();
    return family == "dde" || name == "mackey_glass" || name == "mackey-glass";
}

SolveResult executeDelaySolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution) {
    if (job.at("controller").at("kind").string("none") != "none")
        throw std::invalid_argument("Mackey-Glass control is not defined");
    const auto requested = requestedConfiguration(job);
    const auto resolved = dynsyslab::core::resolveMackeyGlass(requested);
    const auto histories = initialHistories(job, resolved);
    dynsyslab::core::DelayIntegrationOptions options;
    options.initialTime = job.at("time").at("start").number();
    options.finalTime = job.at("time").at("stop").number();
    options.step = job.at("solver").at("step_size").number(
        job.at("time").at("output_step").number(0.05));
    const std::string created = detail::timestamp();
    const auto hierarchy = prepareHierarchy(job, resolved, histories.size(), archiveRoot,
                                            created);
    detail::atomicProgress(progressPath, 0.0, 0, histories.size(),
                           "Starting delay integration");
    const auto solution = dynsyslab::core::integrateMackeyGlass(
        resolved, histories, options, execution);
    SolveResult result;
    result.label = job.at("label").string("delay-solve");
    result.system = "mackey_glass";
    result.trajectories = solution.trajectoryCount;
    result.samplesPerTrajectory = solution.sampleCount;
    result.acceptedSteps = solution.acceptedSteps;
    const auto directory = archiveRoot / L"solutions" / L"mackey_glass";
    for (std::size_t trajectory = 0; trajectory < solution.trajectoryCount; ++trajectory) {
        dynsyslab::core::NpzRecord archive;
        archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
            {solution.sampleCount}, solution.time);
        const auto begin = solution.values.begin()
            + static_cast<std::ptrdiff_t>(trajectory * solution.sampleCount);
        archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
            {solution.sampleCount, 1}, std::vector<double>(begin,
                begin + static_cast<std::ptrdiff_t>(solution.sampleCount)));
        archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
            metadata(job, resolved, hierarchy, trajectory, histories[trajectory],
                     options.step, created));
        detail::appendAnalyticComparison(archive,
            dynsyslab::core::mackeyGlassDescriptor(), resolved,
            dynsyslab::core::mackeyGlassAnalyticSolutions(resolved), solution.time,
            std::vector<double>(begin, begin+static_cast<std::ptrdiff_t>(solution.sampleCount)),
            {histories[trajectory]}, options.initialTime, job);
        const std::string label = result.label + (solution.trajectoryCount == 1 ? ""
            : "-trajectory-" + std::to_string(trajectory));
        const auto path = directory / dynsyslab::core::wide(
            detail::safeName(label, "delay-solve") + "__" + created + ".npz");
        dynsyslab::core::writeNpzAtomicNew(path, archive);
        result.archives.push_back(path);
        detail::atomicProgress(progressPath,
            static_cast<double>(trajectory + 1) / solution.trajectoryCount,
            trajectory + 1, solution.trajectoryCount, "Archived delay trajectory");
    }
    publishHierarchy(hierarchy, job, requested,
        dynsyslab::core::mackeyGlassDescriptor(), resolved, histories,
        result.archives, false, false, created);
    result.hierarchyManifest = hierarchy.manifestPath;
    return result;
}

}  // namespace dynsyslab::workflow
