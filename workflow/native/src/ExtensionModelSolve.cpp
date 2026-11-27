#include "ExtensionModelSolve.hpp"

#include <dynsyslab/workflow/Hierarchy.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/DelayRollout.hpp>
#include <dynsyslab/core/Extension.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/Text.hpp>

#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// Solves of extension-supplied delay (method of steps, DelayRollout.hpp), field (method of
// lines, PDE.hpp) and lattice (sweep callback) models. Archives follow the ODE solve layout:
// time, values, metadata_json, one archive per trajectory, and a hierarchy manifest.

namespace dynsyslab::workflow::detail {
namespace {

std::vector<double> numbers(const json::Value& source) {
    std::vector<double> result;
    for (const auto& item : source.array()) result.push_back(item.number());
    return result;
}

std::string stamp() {
    const auto raw = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out; out << std::put_time(&utc, "%Y%m%dT%H%M%S") << 'Z';
    return out.str();
}

std::size_t uniformSamples(double start, double stop, double step) {
    if (!(stop > start) || !(step > 0.0)) throw std::invalid_argument("extension solve needs start < stop and a positive output_step");
    return static_cast<std::size_t>(std::floor((stop - start) / step + 1.0e-9)) + 1;
}

// One trajectory per explicit initial-condition row (default: one zero state).
std::vector<std::vector<double>> initialRows(const json::Value& job, std::size_t dimension) {
    const auto& spec = job.at("initial_conditions");
    std::vector<std::vector<double>> rows;
    for (const auto& row : spec.at("values").array()) {
        auto values = row.array().empty() ? std::vector<double>{row.number()} : numbers(row);
        if (values.size() != dimension) throw std::invalid_argument("extension initial state dimension mismatch");
        rows.push_back(std::move(values));
    }
    if (rows.empty()) rows.assign(static_cast<std::size_t>(std::max(1.0, spec.at("count").number(1.0))),
                                  std::vector<double>(dimension, 0.0));
    return rows;
}

dynsyslab::core::StructuredMesh meshFrom(const json::Value& descriptor) {
    const auto& source = descriptor.at("mesh");
    dynsyslab::core::StructuredMesh mesh;
    for (const auto& extent : source.at("shape").array()) mesh.shape.push_back(static_cast<std::size_t>(extent.number()));
    mesh.spacing = numbers(source.at("spacing"));
    mesh.origin = numbers(source.at("origin"));
    for (const auto& item : source.at("boundary").array()) {
        const auto name = item.string();
        mesh.boundary.push_back(name == "periodic" ? dynsyslab::core::PDEBoundary::periodic
            : name == "edge" ? dynsyslab::core::PDEBoundary::edge : dynsyslab::core::PDEBoundary::zeroValue);
    }
    mesh.validate();
    return mesh;
}

struct Trajectory { std::vector<double> time, values; std::size_t width{}; std::string status{"complete"}; };

}  // namespace

bool isExtensionModelSolve(const json::Value& job) {
    return static_cast<bool>(dynsyslab::core::ExtensionRegistry::instance().model(job.at("system").at("name").string()));
}

SolveResult executeExtensionModelSolve(const json::Value& job, const std::filesystem::path& archiveRoot,
                                       const std::filesystem::path&, const dynsyslab::core::ExecutionContext* execution) {
    const auto model = *dynsyslab::core::ExtensionRegistry::instance().model(job.at("system").at("name").string());
    const auto descriptor = json::parse(model.descriptorJson);
    const double start = job.at("time").at("start").number(), stop = job.at("time").at("stop").number();
    const double step = job.at("time").at("output_step").number();
    std::vector<Trajectory> trajectories;
    std::vector<double> initials;
    std::size_t dimension = 0;
    if (model.family == "dde") {
        dimension = static_cast<std::size_t>(descriptor.at("state_dimension").number());
        dynsyslab::core::DelayRolloutOptions options;
        options.step = step;
        options.lags = numbers(descriptor.at("lags"));
        options.outputSamples = uniformSamples(start, stop, step);
        const auto historyRows = static_cast<std::size_t>(std::ceil(options.lags.back() / step - 1.0e-9)) + 1;
        for (const auto& row : initialRows(job, dimension)) {
            if (execution) execution->checkpoint();
            std::vector<double> history;   // constant prehistory equal to the initial state
            for (std::size_t k = 0; k < std::max<std::size_t>(historyRows, 2); ++k) history.insert(history.end(), row.begin(), row.end());
            const auto rollout = dynsyslab::core::rolloutDelayModel(
                [&](const double* window, double* derivative) {
                    model.delayRhs(window, options.lags.size() * dimension, derivative, dimension);
                }, dimension, history, options);
            Trajectory t; t.width = dimension; t.values = rollout.values;
            for (std::size_t k = 0; k < rollout.samples; ++k) t.time.push_back(start + k * step);
            if (rollout.status != dynsyslab::core::DelayRolloutStatus::complete) t.status = rollout.reason;
            trajectories.push_back(std::move(t));
            initials.insert(initials.end(), row.begin(), row.end());
        }
    } else if (model.family == "pde") {
        const auto mesh = meshFrom(descriptor);
        std::vector<std::size_t> shape;
        for (const auto& extent : descriptor.at("field_shape").array()) shape.push_back(static_cast<std::size_t>(extent.number()));
        std::size_t values = 1;
        for (const auto extent : shape) values *= extent;
        dimension = values;
        dynsyslab::core::PDEIntegrationOptions options;
        options.initialTime = start; options.finalTime = stop; options.outputStep = step;
        const auto method = job.at("solver").at("method").string("RK45");
        if (method == "Euler" || method == "euler" || method == "RK4" || method == "rk4") {
            options.method = method[0] == 'E' || method[0] == 'e' ? dynsyslab::core::PDEIntegrationMethod::euler
                                                                  : dynsyslab::core::PDEIntegrationMethod::rk4;
            options.internalStep = job.at("solver").at("step_size").number(step);
        } else {
            options.method = dynsyslab::core::PDEIntegrationMethod::adaptive;
            options.adaptiveMethod = dynsyslab::core::integrationMethodFromName(method);
            options.relativeTolerance = job.at("solver").at("rtol").number(1.0e-6);
            options.absoluteTolerance = job.at("solver").at("atol").number(1.0e-9);
        }
        const auto initial = model.initialField(json::dump(job.at("system").at("parameters")), values);
        const auto solved = dynsyslab::core::integratePDE(mesh, shape, initial, options,
            [&](double time, const dynsyslab::core::StructuredMesh&, const double* field, std::size_t count, double* derivative) {
                model.fieldRhs(time, field, count, derivative);
            }, execution);
        trajectories.push_back({solved.time, solved.values, values, "complete"});
        initials = initial;
    } else {
        const auto sites = static_cast<std::size_t>(descriptor.at("sites").number());
        const auto observables = descriptor.at("observables").array().size();
        if (!sites || !observables) throw std::invalid_argument("extension lattice needs sites and observables");
        dimension = observables;
        const auto sweeps = uniformSamples(start, stop, step) - 1;
        const double count = job.at("initial_conditions").at("count").number(1.0);
        for (std::size_t member = 0; member < static_cast<std::size_t>(std::max(1.0, count)); ++member) {
            const auto seed = static_cast<std::uint64_t>(job.at("initial_conditions").at("seed").number()) + member;
            std::vector<std::int32_t> state(sites);
            model.latticeInitial(seed, state.data(), sites);
            Trajectory t; t.width = observables;
            std::vector<double> sample(observables);
            for (std::size_t sweep = 0; sweep < sweeps; ++sweep) {
                if (execution) execution->checkpoint();
                model.latticeSweep(seed, sweep, state.data(), sites, sample.data(), observables);
                t.time.push_back(start + (sweep + 1) * step);
                t.values.insert(t.values.end(), sample.begin(), sample.end());
            }
            trajectories.push_back(std::move(t));
            initials.insert(initials.end(), observables, 0.0);
        }
    }

    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = model.identifier;
    configuration.preset = job.at("system").at("preset").string();
    for (const auto& item : job.at("system").at("parameters").object())
        configuration.parameters[item.first] = item.second.number();
    configuration.equationsVersion = "extension." + model.identifier + "/1";
    dynsyslab::core::SystemDescriptor system;
    system.identifier = model.identifier; system.label = model.label;
    system.stateDimension = dimension; system.family = model.family;
    const auto tag = stamp();
    const auto hierarchy = prepareHierarchy(job, configuration, trajectories.size(), archiveRoot, tag);
    SolveResult result;
    result.label = job.at("label").string("extension-solve");
    result.system = model.identifier;
    result.trajectories = trajectories.size();
    const auto directory = archiveRoot / L"solutions" / dynsyslab::core::wide(model.identifier);
    for (std::size_t index = 0; index < trajectories.size(); ++index) {
        const auto& t = trajectories[index];
        dynsyslab::core::NpzRecord archive;
        archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({t.time.size()}, t.time);
        archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles({t.time.size(), t.width}, t.values);
        json::Value::Object meta{{"schema", std::string("dynsyslab.solution/v1")},
            {"label", result.label}, {"system", model.identifier}, {"system_label", model.label},
            {"system_family", model.family}, {"extension_model", true},
            {"integration_status", t.status}, {"time_configuration", job.at("time")},
            {"solver_configuration", job.at("solver")}, {"batch_id", hierarchy.batchId},
            {"trajectory_id", hierarchy.trajectoryIds[index]}, {"trajectory_index", static_cast<double>(index)},
            {"source_job_label", job.at("label").string()}, {"created_utc", tag}};
        archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
        const auto path = directory / dynsyslab::core::wide(result.label + "-trajectory-" + std::to_string(index) + "__" + tag + ".npz");
        dynsyslab::core::writeNpzAtomicNew(path, archive);
        result.archives.push_back(path);
        result.samplesPerTrajectory = std::max(result.samplesPerTrajectory, t.time.size());
    }
    publishHierarchy(hierarchy, job, configuration, system, configuration, initials, result.archives,
                     false, false, tag);
    result.hierarchyManifest = hierarchy.manifestPath;
    return result;
}

}  // namespace dynsyslab::workflow::detail
