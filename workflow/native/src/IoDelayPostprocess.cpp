// Post-processing method "io_delay_model": ridge ARX fit of one named scalar observable from
// forced records on a decimated held-input clock (Ljung, System Identification, 2nd ed., 1999;
// delay-embedded predictors for control: Korda and Mezic, Automatica 93, 2018).
//
// Settings: {"observable":"far_field_acoustic_echo","decimation":1,"ny":10,"nu":10,
//            "relative_ridge":1e-10,"fit_offset":true,"validation_sources":[PATH,...]}
// Sources are the training records. Each record is decimated separately, so no regression row
// crosses a record boundary. See docs/architecture/IO_DELAY_CONTROL.md.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/IoDelayModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;

std::size_t count(const json::Value& value, double fallback, const char* name) {
    const double number = value.number(fallback);
    if (!std::isfinite(number) || number < 1 || std::floor(number) != number || number > 100000)
        throw std::invalid_argument(std::string("io_delay_model ") + name + " must be a positive integer");
    return static_cast<std::size_t>(number);
}

json::Value finiteOrNull(double value) { return std::isfinite(value) ? json::Value(value) : json::Value(); }

struct Loaded { SourceData source; c::IoDelayRecord record; };

Loaded loadRecord(const std::filesystem::path& path, const std::string& observable, std::size_t decimation) {
    Loaded out{loadSource(path), {}};
    const auto& source = out.source;
    if (source.inputs().empty() || source.inputDimension() != 1)
        throw std::invalid_argument("io_delay_model sources need a scalar control history: " + path.u8string());
    const auto system = c::createSystem(systemConfiguration(source));
    const auto measurement = system->observableCatalog().at(observable);
    if (measurement.outputDimension() != 1)
        throw std::invalid_argument("io_delay_model needs a scalar observable");
    std::vector<double> y(source.samples());
    for (std::size_t k = 0; k < source.samples(); ++k)
        y[k] = measurement.evaluateReal(source.values().data() + k * source.dimension()).front();
    out.record = c::decimateHeldRecord(y, source.inputs(), decimation);
    return out;
}

json::Value validationJson(const std::filesystem::path& path, const c::IoDelayValidation& v) {
    return json::Value(json::Value::Object{
        {"source", json::Value(path.generic_u8string())}, {"rows", json::Value(double(v.rows))},
        {"one_step_rmse", finiteOrNull(v.oneStepRmse)}, {"one_step_relative", finiteOrNull(v.oneStepRelative)},
        {"free_run_rmse", finiteOrNull(v.freeRunRmse)}, {"free_run_relative", finiteOrNull(v.freeRunRelative)},
        {"free_run_diverged", json::Value(v.freeRunDiverged)}});
}

}  // namespace

std::filesystem::path runIoDelayModel(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                      const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("io_delay_model needs forced training records");
    const auto observable = settings.at("observable").string("far_field_acoustic_echo");
    const auto decimation = count(settings.at("decimation"), 1, "decimation");
    c::IoDelayFitOptions options;
    options.ny = count(settings.at("ny"), 10, "ny");
    options.nu = count(settings.at("nu"), static_cast<double>(options.ny), "nu");
    options.relativeRidge = settings.at("relative_ridge").number(1e-10);
    options.fitOffset = settings.at("fit_offset").boolean(true);

    std::vector<c::IoDelayRecord> records;
    json::Value::Array trainPaths;
    std::string referenceMetadata;
    c::SystemConfiguration reference;
    double fineStep = 0.0;
    for (const auto& path : sources) {
        auto loaded = loadRecord(path, observable, decimation);
        const double step = uniformDeltaTime(loaded.source);
        const auto configuration = systemConfiguration(loaded.source);
        if (records.empty()) { fineStep = step; reference = configuration; referenceMetadata = loaded.source.metadataJson; }
        else {
            if (std::abs(step - fineStep) > 1e-9 * fineStep)
                throw std::invalid_argument("io_delay_model records must share one sample clock");
            c::requireCompatibleSystems(reference, configuration);
        }
        records.push_back(std::move(loaded.record));
        trainPaths.emplace_back(path.generic_u8string());
    }
    const double sampleTime = fineStep * static_cast<double>(decimation);
    const auto [model, report] = c::fitIoDelayModel(records, options, sampleTime);

    json::Value::Array trainValidation, heldOut, heldPaths;
    for (std::size_t r = 0; r < records.size(); ++r)
        trainValidation.push_back(validationJson(sources[r], c::validateIoDelayModel(model, records[r])));
    for (const auto& item : settings.at("validation_sources").array()) {
        const auto path = std::filesystem::u8path(item.string());
        const auto loaded = loadRecord(path, observable, decimation);
        if (std::abs(uniformDeltaTime(loaded.source) - fineStep) > 1e-9 * fineStep)
            throw std::invalid_argument("io_delay_model validation record uses a different clock");
        c::requireCompatibleSystems(reference, systemConfiguration(loaded.source));
        heldOut.push_back(validationJson(path, c::validateIoDelayModel(model, loaded.record)));
        heldPaths.emplace_back(path.generic_u8string());
    }
    const auto realization = c::ioDelayRealization(model);
    const auto poles = c::systemPoles(realization.system);
    double radius = 0.0;
    for (const auto& pole : poles) radius = std::max(radius, std::abs(pole));

    json::Value::Object meta{
        {"schema", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("io-delay-model")},
        {"label", json::Value(label)}, {"observable", json::Value(observable)},
        {"source_archive", json::Value(sources.front().generic_u8string())},
        {"training_sources", json::Value(std::move(trainPaths))},
        {"validation_sources", json::Value(std::move(heldPaths))},
        {"ny", json::Value(double(model.ny))}, {"nu", json::Value(double(model.nu))},
        {"decimation", json::Value(double(decimation))}, {"source_delta_t", json::Value(fineStep)},
        {"delta_t", json::Value(sampleTime)}, {"relative_ridge", json::Value(options.relativeRidge)},
        {"fit_offset", json::Value(options.fitOffset)},
        {"regression_rows", json::Value(double(report.rows))},
        {"rows_rejected_for_hold", json::Value(double(report.rowsRejectedForHold))},
        {"training_one_step_rmse", json::Value(report.oneStepRmse)},
        {"training_validation", json::Value(std::move(trainValidation))},
        {"held_out_validation", json::Value(std::move(heldOut))},
        {"open_loop_spectral_radius", json::Value(radius)},
        {"state_layout", json::Value("z=(y[k..k-ny+1], u[k-1..k-nu+1])")},
        {"citations", json::Value(json::Value::Array{
            json::Value("Ljung, System Identification: Theory for the User, 2nd ed., 1999"),
            json::Value("Korda and Mezic, Automatica 93 (2018) 149-160")})}};
    c::NpzRecord archive;
    c::storeIoDelayModel(model, archive);
    std::vector<std::complex<double>> poleValues(poles.begin(), poles.end());
    archive.arrays["open_loop_poles"] = c::NpyArray::fromComplex({poleValues.size()}, poleValues);
    archive.arrays["source_metadata_json"] = c::NpyArray::fromString(referenceMetadata);
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const auto directory = kindFolder(root, dynsyslab::core::object_kind::dynamicalModel);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(sources.front().stem().u8string()
        + "__io-delay-model__d" + std::to_string(decimation) + "_ny" + std::to_string(model.ny)
        + "_nu" + std::to_string(model.nu) + "__" + timestamp()) + ".npz");
    stampDataObject(archive, dynsyslab::core::object_kind::dynamicalModel, "io-delay-model", {});
    c::writeNpzAtomicNew(path, archive);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
