#include "MetadataPaths.hpp"
#include "PhaseSpaceSerialization.hpp"
#include "PhaseForcingConfiguration.hpp"
#include "WorkflowSupport.hpp"
#include <cmath>
#include <algorithm>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <stdexcept>
namespace dynsyslab::workflow {
json::Value evaluatePhaseSpaceQuery(const json::Value &request,
                                    const std::filesystem::path &archiveRoot,
                                    const dynsyslab::core::ExecutionContext *execution) {
    using namespace dynsyslab::core;
    if (request.at("schema").string() != "dynsyslab.phase-space-query/v1")
        throw std::invalid_argument("Unsupported phase-space query schema");
    const auto mode = request.at("mode").string("lookup");
    if (mode != "lookup" && mode != "compute")
        throw std::invalid_argument("Phase-space mode must be lookup or compute");
    if (archiveRoot.empty())
        throw std::invalid_argument("Phase-space cache requires an archive root");
    const auto &spec = request.at("system");
    SystemConfiguration configuration;
    configuration.identifier = spec.at("name").string(spec.at("identifier").string());
    configuration.preset = spec.at("preset").string();
    for (const auto &p : spec.at("parameters").object()) {
        if (!std::holds_alternative<double>(p.second.storage()) ||
            !std::isfinite(p.second.number()))
            throw std::invalid_argument("System parameters must be finite numbers");
        configuration.parameters[p.first] = p.second.number();
    }
    const auto &driver = spec.at("driver");
    configuration.driver = {driver.at("amplitude").number(), driver.at("frequency_hz").number(),
                            driver.at("phase_radians").number()};
    for (double value : {configuration.driver.amplitude, configuration.driver.frequencyHz, configuration.driver.phaseRadians})
        if (!std::isfinite(value)) throw std::invalid_argument("Driver values must be finite");
    auto system = createSystem(configuration);
    if (request.at("controller").at("kind").string() == "waveform") {
        if (configuration.driver.amplitude != 0)
            throw std::invalid_argument("Specify either a periodic driver or an archived waveform controller");
        configuration = system->configuration();
        configuration.driver = detail::phaseDriverFromController(request.at("controller"),
            request.at("control_evaluation").string("continuous_rhs"), *system);
        system = createSystem(configuration);
    }
    const auto &resolved = system->configuration();
    auto options =
        phaseSpaceOptionsFromJson(request.at("options"), system->descriptor().stateDimension);
    if (options.dynamics == "periodic-forcing") {
        options.forcingPeriod = phaseSpaceForcingPeriod(*system, options.forcingPeriod);
        const auto &requested = request.at("options");
        if (requested.at("maximum_step").isNull()) options.maximumStep = options.forcingPeriod / 32;
        if (requested.at("output_step").isNull()) options.outputStep = options.forcingPeriod / 128;
        if (requested.at("horizon").isNull()) options.horizon = options.forcingPeriod * std::max(8., 2. * options.forcingPeriods);
        if (requested.at("transient_time").isNull()) options.transientTime = 0;
    }
    validatePhaseSpaceOptions(*system, options);
    const auto knowledge = dynamicalKnowledgeJson(loadSystemKnowledge(*system, archiveRoot));
    auto guideIdentity = dynamicalKnowledgeJson(system->dynamicalKnowledge());
    guideIdentity.object().erase("paper_catalogue");
    json::Value::Object parameters, groups, units;
    for (const auto &p : resolved.parameters)
        parameters[p.first] = p.second;
    for (const auto &p : resolved.dimensionlessGroups)
        groups[p.first] = p.second;
    for (const auto &p : resolved.parameterUnits)
        units[p.first] = p.second;
    const auto coordinates = system->coordinates();
    json::Value definition =
        json::Value::Object{{"identifier", resolved.identifier},
                            {"preset", resolved.preset},
                            {"parameters", parameters},
                            {"parameter_units", units},
                            {"dimensionless_groups", groups},
                            {"equations_version", resolved.equationsVersion},
                            {"coordinate_system", system->descriptor().preferredCoordinates},
                            {"time_scale", coordinates.secondsPerModelTime},
                            {"input_scale", coordinates.physicalInputPerModelInput},
                            {"dynamics", options.dynamics},
                            {"driver", json::Value::Object{{"amplitude", resolved.driver.amplitude},
                                {"frequency_hz", resolved.driver.frequencyHz}, {"phase_radians", resolved.driver.phaseRadians}}}};
    const auto systemKey = sha256Text(json::dump(definition));
    json::Value key = json::Value::Object{{"schema", "dynsyslab.phase-cache-key/v1"},
                                          {"method_version", "5"},
                                          {"system", definition},
                                          {"knowledge", guideIdentity},
                                          {"options", phaseSpaceOptionsJson(options)}};
    const auto canonical = json::dump(key), digest = sha256Text(canonical);
    // A distinct extension keeps auxiliary context out of the trajectory data tree.
    const auto path =
        archiveRoot / "system_metadata" / wide(systemKey) / "phase" / wide(digest + ".dslphase");
    const auto ioPath = detail::metadataIoPath(path);
    auto cached = [&]() -> json::Value {
        const auto archive = readNpz(ioPath);
        auto metadata = json::parse(archive.metadataJson());
        if (metadata.at("schema").string() != "dynsyslab.phase-cache/v1" ||
            metadata.at("canonical_key").string() != canonical)
            throw std::runtime_error(
                "Phase-space cache identity mismatch; existing evidence was preserved");
        auto result = metadata.at("result");
        result.object()["dynamical_knowledge"] = knowledge;
        result.object()["cache"] = json::Value::Object{{"status", "hit"},
                                                       {"key", digest},
                                                       {"system_key", systemKey},
                                                       {"path", path.generic_u8string()},
                                                       {"created_utc", metadata.at("created_utc")}};
        return result;
    };
    if (std::filesystem::exists(ioPath))
        return cached();
    if (mode == "lookup")
        return json::Value::Object{
            {"schema", "dynsyslab.phase-space/v1"},
            {"cache", json::Value::Object{{"status", "miss"},
                                          {"key", digest},
                                          {"system_key", systemKey},
                                          {"path", path.generic_u8string()}}}};
    auto result = detail::phaseSpaceResultJson(analyzePhaseSpace(*system, options, execution));
    result.object()["dynamical_knowledge"] = knowledge;
    if (options.dynamics == "periodic-forcing") {
        constexpr double tau = 6.283185307179586476925286766559;
        double sectionPhase = std::fmod(tau * options.forcingSectionTime / options.forcingPeriod, tau);
        if (sectionPhase < 0) sectionPhase += tau;
        result.object()["forcing_clock"] = json::Value::Object{
            {"period", options.forcingPeriod}, {"time_units", coordinates.modelTimeUnits},
            {"time_origin", 0.}, {"phase_origin", 0.}, {"section_time", options.forcingSectionTime},
            {"section_phase", sectionPhase}, {"phase_units", "radians"},
            {"section_convention", "include start; exact returns after integer forcing periods"},
            {"driver", definition.at("driver")},
            {"seconds_per_model_time", coordinates.secondsPerModelTime},
            {"lift", "physical state plus one forcing phase; no symplectic claim"}};
    }
    json::Value::Array sources;
    for (const char *id :
         {"CODE.dynamical-systems.basin-manifold-toolkit",
          "CODE.dynamical-systems.dynamical-systems-model-library",
          "MTH.dynamical-systems.basin-and-separatrix-computation-audit",
          "MTH.dynamical-systems.poincare-return-map-stability-audit",
          "MTH.dynamical-systems.invariant-manifold-parameterization-and-invariance-equation",
          "MTH.dynamical-systems.periodic-forcing-phase-lift-and-viewer-audit"})
        sources.emplace_back(id);
    result.object()["method_provenance"] = json::Value::Object{
        {"peer", "Grothendieck Labs"},
        {"records", sources},
        {"implementation", "Native reformulation with adaptive integration, residual checks and "
                           "explicit finite-search limits"}};
    NpzRecord archive;
    const auto stamp = detail::timestamp();
    archive.arrays["metadata_json"] =
        NpyArray::fromString(json::dump(json::Value::Object{{"schema", "dynsyslab.phase-cache/v1"},
                                                            {"canonical_key", canonical},
                                                            {"created_utc", stamp},
                                                            {"result", result}}));
    try {
        writeNpzAtomicNew(ioPath, archive);
    } catch (const std::exception &) {
        if (std::filesystem::exists(ioPath))
            return cached();
        throw;
    }
    result.object()["cache"] = json::Value::Object{{"status", "computed"},
                                                   {"key", digest},
                                                   {"system_key", systemKey},
                                                   {"path", path.generic_u8string()},
                                                   {"created_utc", stamp}};
    return result;
}
} // namespace dynsyslab::workflow
