#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include "AnalyticArchive.hpp"
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/workflow/EstimatorFactory.hpp>
#include <dynsyslab/workflow/Hierarchy.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/TimeSeriesArchive.hpp>

#include "TrainingPreflight.hpp"
#include "SolveArchive.hpp"
#include "SolveCheckpoint.hpp"
#include "SolveContinuation.hpp"
#include "SolveHistory.hpp"
#include "SolveSettings.hpp"
#include <dynsyslab/core/DimensionalSystem.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <variant>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow {
namespace {

std::string escaped(const std::string& value) {
    std::ostringstream out;
    for (const unsigned char item : value) {
        switch (item) {
        case '"': out << "\\\""; break; case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break; case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (item < 0x20) out << "\\u" << std::hex << std::setw(4)
                                 << std::setfill('0') << static_cast<int>(item) << std::dec;
            else out << static_cast<char>(item);
        }
    }
    return out.str();
}

std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t raw = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out; out << std::put_time(&utc, "%Y%m%dT%H%M%S") << 'Z';
    return out.str();
}

std::string safeName(std::string value) {
    for (char& item : value)
        if (!std::isalnum(static_cast<unsigned char>(item)) && item != '-' && item != '_') item = '-';
    while (value.find("--") != std::string::npos) value.replace(value.find("--"), 2, "-");
    return value.empty() ? "native-solve" : value;
}

void atomicText(const std::filesystem::path& path, const std::string& text) {
    if (path.empty()) return;
    (void)dynsyslab::core::prepareOutputFile(path, "progress file");
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    { std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to create progress record"); stream << text; }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish progress record");
#else
    std::filesystem::rename(temporary, path);
#endif
}

void progress(const std::filesystem::path& path, double fraction, std::size_t completed,
              std::size_t planned, const char* state, const std::string& message) {
    std::ostringstream out;
    out << "{\"schema\":\"dynsyslab.compute-progress/v1\",\"determinate\":true,"
        << "\"percent\":" << std::setprecision(8) << fraction * 100.0
        << ",\"estimated_remaining_seconds\":-1,\"completed_trajectories\":" << completed
        << ",\"planned_trajectories\":" << planned << ",\"state\":\"" << state
        << "\",\"message\":\"" << escaped(message) << "\"}";
    atomicText(path, out.str());
}

dynsyslab::core::SystemConfiguration systemConfiguration(const json::Value& job) {
    const auto& source = job.at("system");
    dynsyslab::core::SystemConfiguration result;
    result.identifier = source.at("name").string();
    result.preset = source.at("preset").string();
    for (const auto& parameter : source.at("parameters").object())
        result.parameters.emplace(parameter.first, parameter.second.number());
    // system.driver: absent/"none" = no external driver; "preset" = the preset's recommended
    // driver; {amplitude, frequency_hz, phase_rad} = an explicit sinusoidal driver.
    const auto& driver = source.at("driver");
    if (driver.string() == "preset") {
        for (const auto& descriptor : dynsyslab::core::systemCatalog())
            if (descriptor.identifier == result.identifier)
                for (const auto& preset : descriptor.presets)
                    if (preset.key == result.preset) result.driver = preset.driver;
    } else if (std::holds_alternative<json::Value::Object>(driver.storage())) {
        result.driver = {driver.at("amplitude").number(), driver.at("frequency_hz").number(),
                         driver.at("phase_rad").number()};
    } else if (!driver.isNull() && driver.string() != "none") {
        throw std::invalid_argument("system.driver must be none, preset or {amplitude, frequency_hz, phase_rad}");
    }
    return result;
}

std::string metadata(const json::Value& job, const dynsyslab::core::ODESystem& system,
                     const std::string& label, const std::vector<double>& initial,
                     std::size_t trajectory, const ControllerSetup& control,
                     const HierarchyIdentity& hierarchy,
                     const detail::TrainingPreflightResult& preflight) {
    const auto& systemSpec = job.at("system");
    const auto coordinates = system.coordinates();
    std::ostringstream out;
    out << "{\"schema\":\"dynsyslab.solution/v1\",\"label\":\"" << escaped(label)
        << "\",\"system\":\"" << escaped(system.descriptor().identifier)
        << "\",\"system_label\":\"" << escaped(system.descriptor().label)
        << "\",\"system_preset\":\"" << escaped(systemSpec.at("preset").string())
        << "\""
        << ",\"system_parameters\":{";
    bool firstParameter = true;
    for (const auto& parameter : systemSpec.at("parameters").object()) {
        if (!firstParameter) out << ',';
        firstParameter = false;
        out << '"' << escaped(parameter.first) << "\":" << std::setprecision(17)
            << parameter.second.number();
    }
    out << "},\"resolved_system_parameters\":{";
    bool firstResolved = true;
    for (const auto& parameter : system.configuration().parameters) {
        if (!firstResolved) out << ',';
        firstResolved = false;
        out << '"' << escaped(parameter.first) << "\":" << std::setprecision(17)
            << parameter.second;
    }
    out << "},\"parameter_units\":{";
    firstResolved = true;
    for (const auto& unit : system.configuration().parameterUnits) {
        if (!firstResolved) out << ',';
        firstResolved = false;
        out << '"' << escaped(unit.first) << "\":\"" << escaped(unit.second) << '"';
    }
    out << "},\"dimensionless_groups\":{";
    firstResolved = true;
    for (const auto& group : system.configuration().dimensionlessGroups) {
        if (!firstResolved) out << ',';
        firstResolved = false;
        out << '"' << escaped(group.first) << "\":" << std::setprecision(17)
            << group.second;
    }
    out << "},\"dimensionless_group_units\":{";
    firstResolved = true;
    for (const auto& unit : system.configuration().dimensionlessGroupUnits) {
        if (!firstResolved) out << ',';
        firstResolved = false;
        out << '"' << escaped(unit.first) << "\":\"" << escaped(unit.second) << '"';
    }
    out << "},\"equations_version\":\"" << escaped(system.configuration().equationsVersion)
        << "\",\"equations_fingerprint\":\""
        << dynsyslab::core::equationsFingerprint(system.configuration())
        << "\",\"integration_form\":\"" << escaped(system.descriptor().preferredCoordinates)
        << "\",\"archive_coordinate_system\":\"" << escaped(system.descriptor().preferredCoordinates)
        << "\",\"forcing_type\":\"" << (control ? "forced" : "free")
        << "\",\"controller_kind\":\"" << escaped(control.kind)
        << "\",\"controller_family\":\"" << escaped(control.family)
        << "\",\"control_evaluation\":\"" << escaped(control.evaluation)
        << "\",\"control_sample_time\":" << std::setprecision(17)
        << control.sampleTime
        << ",\"control_coordinate_system\":\"" << escaped(system.descriptor().preferredCoordinates)
        << "\",\"control_units\":\"" << escaped(coordinates.modelInputUnits)
        << "\",\"control_scale_to_dimensional\":" << std::setprecision(17)
        << coordinates.physicalInputPerModelInput
        << ",\"control_dimensional_units\":\"" << escaped(coordinates.physicalInputUnits)
        << "\",\"control_physical_scale_available\":"
        << (coordinates.physicalInputAvailable ? "true" : "false")
        << ",\"controller_settings_coordinates\":\""
        << escaped(control.settingsCoordinates)
        << "\",\"model_time_units\":\"" << escaped(coordinates.modelTimeUnits)
        << "\",\"time_scale_to_seconds\":" << std::setprecision(17)
        << coordinates.secondsPerModelTime
        << ",\"physical_time_scale_available\":"
        << (coordinates.physicalTimeAvailable ? "true" : "false")
        << ",\"execution_backend\":{\"name\":\"native-cpp\",\"compiled\":true,\"python_required\":false}"
        << ",\"hierarchy_manifest\":\"" << escaped(hierarchy.manifestPath.generic_u8string())
        << "\",\"campaign_id\":\"" << escaped(hierarchy.campaignId)
        << "\",\"campaign_kind\":\"" << escaped(hierarchy.campaignKind)
        << "\",\"campaign_point_index\":"
        << job.at("campaign").at("point_index").number(-1.0)
        << ",\"control_target_space\":\"" << control.targetSpace << "\""
        << ",\"controller_configuration\":"
        << json::dump(job.at("controller"))
        << ",\"model_binding_evidence\":" << json::dump(control.modelBindingEvidence)
        << ",\"solver_configuration\":" << json::dump(job.at("solver"))
        << ",\"time_configuration\":" << json::dump(job.at("time"))
        << ",\"batch_id\":\"" << hierarchy.batchId << "\",\"system_run_id\":\""
        << hierarchy.systemRunId << "\",\"parameter_set_key\":\""
        << hierarchy.parameterSetKey << "\",\"cluster_id\":\""
        << hierarchy.clusterIds[hierarchy.trajectoryClusters[trajectory]]
        << "\",\"trajectory_id\":\"" << hierarchy.trajectoryIds[trajectory]
        << "\",\"data_role\":\"" << escaped(hierarchy.dataRole) << '"'
        << ",\"training_preflight\":" << json::dump(job.at("preflight"))
        << ",\"training_preflight_result\":";
    if (preflight.enabled) {
        json::Value perTrajectory = preflight.report;
        auto& reportObject = perTrajectory.object();
        const auto& entries = preflight.report.at("entries").array();
        reportObject.erase("entries");
        if (trajectory < entries.size()) reportObject["trajectory"] = entries[trajectory];
        out << json::dump(perTrajectory);
    } else out << "{}";
    out << ",\"initial_state\":[";
    for (std::size_t index = 0; index < initial.size(); ++index) {
        if (index) out << ','; out << std::setprecision(17) << initial[index];
    }
    out << "],\"trajectory_index\":" << trajectory << ",\"source_job_label\":\""
        << escaped(job.at("label").string()) << "\",\"created_utc\":\"" << timestamp() << "\"}";
    return out.str();
}

}  // namespace

std::wstring SolveResult::summary() const {
    std::wostringstream out;
    out << L"Native solve completed\r\n\r\nSystem: " << dynsyslab::core::wide(system)
        << L"\r\nTrajectories: " << trajectories << L"\r\nSamples per trajectory: "
        << samplesPerTrajectory << L"\r\nAccepted internal steps: " << acceptedSteps
        << L"\r\nRejected internal steps: " << rejectedSteps << L"\r\n\r\nArchives:";
    for (const auto& path : archives) out << L"\r\n" << path.wstring();
    if (!hierarchyManifest.empty()) out << L"\r\n\r\nHierarchy:\r\n" << hierarchyManifest.wstring();
    return out.str();
}

SolveResult executeSolve(const json::Value& job,
                                     const std::filesystem::path& archiveRoot,
                                     const std::filesystem::path& progressPath,
                                     const dynsyslab::core::ExecutionContext* execution) {
    if (execution) execution->report({"solve", "Preparing native solve", 0.0, 0, 1});
    auto configuration = systemConfiguration(job);
    // time.coordinate_system: "automatic" integrates in the system's preferred coordinates;
    // "dimensional" integrates a dimensionless model in physical units (seconds, SI state)
    // when the system declares its scaling (DimensionalSystem.hpp).
    auto system = dynsyslab::core::systemInCoordinates(
        dynsyslab::core::createSystem(configuration),
        job.at("time").at("coordinate_system").string("automatic"));
    const std::string controlKind = job.at("controller").at("kind").string("none");
    const bool controlled = controlKind != "none";
    const auto initialSet = resolveOdeInitialConditions(job);
    const std::size_t trajectoryCount = initialSet.trajectoryCount;
    const auto& initials = initialSet.values;
    const auto settings = detail::solveSettings(job);
    auto options = settings.integration;
    detail::requireContinuationStart(job, options.initialTime);
    const bool hasPreflight = job.at("preflight").at("enabled").boolean(false);
    const double solveProgressBase = hasPreflight ? 0.1 : 0.0;
    const auto preflight = detail::executeTrainingPreflight(
        job, *system, initials, trajectoryCount,
        [&](double fraction, std::size_t completed, std::size_t planned,
            const std::string& message) {
            progress(progressPath, 0.1 * fraction, completed, planned,
                     "screening", message);
            if (execution) execution->report({"solve/screening", message,
                0.1 * fraction, completed, planned});
        }, execution);
    // Controller and Estimator modules: KLQE -> KLQR output feedback shares one KROM and is
    // composed inside the controller; any other estimator runs beside the input policy.
    auto control = buildController(controllerSpecWithEstimator(job), job.at("solver"),
        *system, trajectoryCount, options.outputStep, initials);
    if (!job.at("estimator").isNull() && !usesEstimatedFeedback(job))
        control = attachEstimator(std::move(control),
            buildEstimator(job.at("estimator"), job.at("solver"), *system, trajectoryCount, options.outputStep),
            trajectoryCount, usesEstimatedState(job) ? system->descriptor().stateDimension : 0);
    const std::string stamp = timestamp();
    const std::string runLabel = job.at("label").string("native-solve");
    const auto directory = archiveRoot / L"solutions" / dynsyslab::core::wide(configuration.identifier);
    const auto archivePath = [&](std::size_t trajectory, const std::string& tag) {
        return directory / dynsyslab::core::wide(safeName(runLabel + (trajectoryCount == 1 ? "" :
            "-trajectory-" + std::to_string(trajectory))) + "__" + tag + ".npz");
    };
    // Run tag: the UTC second, made unique (-2, -3, ...) when an earlier run with this label used
    // it, so runs in the same second never share archives, a manifest or a batch identity.
    std::string runTag = stamp;
    auto hierarchy = prepareHierarchy(job, system->configuration(), trajectoryCount, archiveRoot, runTag);
    const bool namedManifest = !job.at("batch").at("manifest_path").string().empty();
    for (int next = 2; (!namedManifest && std::filesystem::exists(dynsyslab::core::extendedLengthPath(hierarchy.manifestPath)))
             || std::filesystem::exists(dynsyslab::core::extendedLengthPath(archivePath(0, runTag))); ++next) {
        if (next > 10000) throw std::runtime_error("unable to allocate an unused solve archive name");
        runTag = stamp + "-" + std::to_string(next);
        hierarchy = prepareHierarchy(job, system->configuration(), trajectoryCount, archiveRoot, runTag);
    }
    progress(progressPath, solveProgressBase, 0, trajectoryCount,
             "running", "Starting native integration");
    // Checkpoints: the running trajectory's samples so far go to <archive>.partial.npz; an
    // existing history.pause_file pauses the solve after that checkpoint.
    std::vector<std::filesystem::path> partials;
    detail::CheckpointClock checkpointClock(settings.checkpoints);
    if (settings.checkpoints.enabled) {
        options.checkpoint = [&](const dynsyslab::core::IntegrationCheckpoint& point) {
            const auto decision = checkpointClock.at(point.trajectory, point.time);
            if (!decision.write) return true;
            const auto partial = detail::partialArchivePath(archivePath(point.trajectory, runTag));
            json::Value::Object meta{{"schema", std::string("dynsyslab.solution-checkpoint/v1")},
                {"source_job_label", runLabel}, {"trajectory_index", static_cast<double>(point.trajectory)},
                {"time", point.time}, {"completed_samples", static_cast<double>(point.sample + 1)}};
            detail::writePartialArchive(partial, point, json::dump(json::Value(std::move(meta))));
            if (std::find(partials.begin(), partials.end(), partial) == partials.end())
                partials.push_back(partial);
            return !decision.pause;
        };
    }
    dynsyslab::core::ExecutionContext solveExecution;
    const dynsyslab::core::ExecutionContext* activeExecution = execution;
    if (execution && hasPreflight) {
        solveExecution.limits = execution->limits;
        solveExecution.cancellation = execution->cancellation;
        solveExecution.progress = [execution, solveProgressBase](
            dynsyslab::core::ExecutionProgress event) {
            event.fraction = solveProgressBase
                + (1.0 - solveProgressBase)
                    * std::clamp(event.fraction, 0.0, 1.0);
            event.stage = "solve/" + event.stage;
            execution->report(std::move(event));
        };
        activeExecution = &solveExecution;
    }
    dynsyslab::core::ControllerCallLog callLog;
    detail::applyHistoryOptions(job, options, callLog, system->descriptor().stateDimension,
                                control.policy ? control.policy.inputDimension : 1);
    const auto solution = detail::integrateKeepingHistory(runLabel,
        [&](std::size_t trajectory) { return archivePath(trajectory, runTag); }, [&] {
        return dynsyslab::core::integrateControlled(
        *system, initials, trajectoryCount, options, control.policy,
        [&](const dynsyslab::core::IntegrationProgress& value) {
            const std::size_t completed = std::min(value.trajectory, trajectoryCount);
            const double fraction = solveProgressBase
                + (1.0 - solveProgressBase) * value.fraction;
            progress(progressPath, fraction, completed, trajectoryCount, "running",
                     "Integrating trajectory " + std::to_string(value.trajectory + 1));
        }, activeExecution); });

    SolveResult result;
    result.label = job.at("label").string("native-solve"); result.system = configuration.identifier;
    result.trajectories = trajectoryCount;
    result.samplesPerTrajectory = settings.resample.enabled ? settings.resample.times.size() : solution.sampleCount;
    result.acceptedSteps = solution.statistics.acceptedSteps;
    result.rejectedSteps = solution.statistics.rejectedSteps;
    const std::size_t dimension = solution.stateDimension;
    const auto& target = control.target;
    // A failure while publishing removes this run's finished archives: no result without a manifest.
    std::vector<std::size_t> kept;
    try {
    for (std::size_t trajectory = 0; trajectory < trajectoryCount; ++trajectory) {
        const std::size_t computed = trajectory < solution.completedSamples.size()
            ? solution.completedSamples[trajectory] : solution.sampleCount;
        const auto& reason = solution.terminationReasons[trajectory];
        if (computed == 0 || (settings.failure == detail::FailurePolicy::reject
                              && computed < solution.sampleCount && reason != "paused")) {
            continue;  // rejected by a guard, or not started before a pause: see the outcome record
        }
        const std::string itemLabel = result.label + (trajectoryCount == 1 ? "" :
            "-trajectory-" + std::to_string(trajectory));
        dynsyslab::core::NpzRecord archive;
        const auto states = solution.trajectorySeries(trajectory);
        dynsyslab::core::storeTimeSeries(archive, states);
        if(control.appendDiagnostics)control.appendDiagnostics(trajectory,archive);
        if (controlled) {
            dynsyslab::core::storeTimeSeries(
                archive, solution.controlSeries(trajectory), "", "control_values");
            if (control.observableTarget) {
                const auto width=control.observableTarget(solution.time.front()).size();
                std::vector<double> values; values.reserve(solution.sampleCount*width);
                for(const auto time:solution.time) {
                    const auto value=control.observableTarget(time);
                    if(value.size()!=width)throw std::runtime_error("observable target dimension changed");
                    values.insert(values.end(),value.begin(),value.end());
                }
                archive.arrays["control_observable_target_values"]=dynsyslab::core::NpyArray::fromDoubles(
                    {solution.sampleCount,width},values);
                std::string names;
                for(const auto& name:control.targetObservableNames){if(!names.empty())names+='\n';names+=name;}
                archive.arrays["control_observable_target_names"]=dynsyslab::core::NpyArray::fromString(names);
            }
            if (target) {
                std::vector<double> targets(solution.sampleCount * dimension);
                for (std::size_t row = 0; row < solution.sampleCount; ++row) {
                    const auto activeTarget = target(solution.time[row]);
                    std::copy(activeTarget.begin(), activeTarget.end(), targets.begin() + row * dimension);
                }
                archive.arrays["control_target_values"] = dynsyslab::core::NpyArray::fromDoubles(
                    {solution.sampleCount, dimension}, targets);
            }
        }
        const std::vector<double> initial(initials.begin() + static_cast<std::ptrdiff_t>(trajectory * dimension),
            initials.begin() + static_cast<std::ptrdiff_t>((trajectory + 1) * dimension));
        archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
            metadata(job, *system, itemLabel, initial, trajectory, control, hierarchy,
                     preflight));
        detail::appendAnalyticComparison(archive, *system, states, initial, options.initialTime, job);
        auto meta = json::parse(archive.at("metadata_json").string());
        auto& metaObject = meta.object();
        metaObject["integration_method"] = dynsyslab::core::integrationMethodName(options.method);
        metaObject["integration_failure_policy"] = std::string(detail::failurePolicyName(settings.failure));
        metaObject["control_input_dimension"] = static_cast<double>(solution.inputDimension);
        for (const auto& preset : system->descriptor().presets) {
            if (preset.key != configuration.preset) continue;
            if (!preset.reference.empty()) metaObject["preset_reference"] = preset.reference;
            if (preset.processingJson != "{}") metaObject["preset_processing"] = json::parse(preset.processingJson);
        }
        if (configuration.driver.amplitude != 0.0)
            metaObject["system_driver"] = json::Value::Object{{"amplitude", configuration.driver.amplitude},
                {"frequency_hz", configuration.driver.frequencyHz}, {"phase_rad", configuration.driver.phaseRadians}};
        if (auto continued = detail::continuationMetadata(job, trajectory); !continued.object().empty())
            metaObject["continuation"] = std::move(continued);
        if (control.trajectoryMetadata) metaObject["input_design"] = control.trajectoryMetadata(trajectory);
        if (settings.resample.enabled) {
            // Post-solve resampling keeps every computed row as raw_<name>.
            metaObject["interpolation"] = detail::resampleOutputArrays(
                archive, solution.time, computed, settings.resample);
            metaObject["interpolated"] = true;
        } else if (computed < solution.sampleCount) {
            detail::truncateOutputRows(archive, solution.sampleCount, computed);
        }
        if (computed < solution.sampleCount) {
            const json::Value record = json::Value::Object{
                {"completed_samples", static_cast<double>(computed)},
                {"time", solution.time[computed - 1]}, {"reason", reason}};
            metaObject[reason == "paused" ? "integration_paused" : "integration_terminated"] = record;
        }
        detail::appendControllerCalls(archive, metaObject, options, trajectory);
        archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(json::dump(meta));
        const auto path = archivePath(trajectory, runTag);
        dynsyslab::core::writeNpzAtomicNew(path, archive); result.archives.push_back(path);
        kept.push_back(trajectory);
    }
    if (kept.empty()) throw std::runtime_error("every trajectory was rejected by the solve's guards");
    // The manifest lists the archived trajectories; rejected ones are recorded beside it.
    auto published = hierarchy;
    std::vector<double> publishedInitials;
    published.trajectoryIds.clear(); published.trajectoryClusters.clear();
    for (const auto trajectory : kept) {
        published.trajectoryIds.push_back(hierarchy.trajectoryIds[trajectory]);
        published.trajectoryClusters.push_back(hierarchy.trajectoryClusters[trajectory]);
        publishedInitials.insert(publishedInitials.end(),
            initials.begin() + static_cast<std::ptrdiff_t>(trajectory * dimension),
            initials.begin() + static_cast<std::ptrdiff_t>((trajectory + 1) * dimension));
    }
    publishHierarchy(published, job, configuration, system->descriptor(),
                     system->configuration(), publishedInitials, result.archives, controlled,
                     static_cast<bool>(control.target), stamp);
    detail::writeSolveOutcome(hierarchy.manifestPath, job, solution, result.archives, kept, initials);
    } catch (...) {
        for (const auto& written : result.archives) {
            std::error_code ignored; std::filesystem::remove(dynsyslab::core::extendedLengthPath(written), ignored);
        }
        result.archives.clear();
        throw;
    }
    for (const auto& partial : partials) {
        std::error_code ignored;
        std::filesystem::remove(dynsyslab::core::extendedLengthPath(partial), ignored);
    }
    result.hierarchyManifest = hierarchy.manifestPath;
    progress(progressPath, 1.0, trajectoryCount, trajectoryCount, "complete", "Native solve archived");
    if (execution) execution->report({"solve", "Native solve archived", 1.0, 1, 1});
    return result;
}

}  // namespace dynsyslab::workflow
