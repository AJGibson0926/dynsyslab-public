// Runtime adapter for controller kind "io_delay_lqr": an echo-only discrete LQ tracker on an
// archived ARX model (Anderson and Moore, Optimal Control: Linear Quadratic Methods, 1990;
// Ljung, System Identification, 1999; Korda and Mezic, Automatica 93, 2018).
#include "IoDelayControl.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/IoDelayModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace dynsyslab::workflow::detail {
namespace c = dynsyslab::core;
namespace json = dynsyslab::core::json;

namespace {
std::filesystem::path modelPath(const json::Value& settings) {
    const auto path = settings.at("model").string();
    if (path.empty()) throw std::invalid_argument("io_delay_lqr requires settings.model (an io_delay_model archive)");
    return std::filesystem::u8path(path);
}
}  // namespace

double ioDelayModelSampleTime(const json::Value& settings) {
    return c::loadIoDelayModel(c::readNpz(modelPath(settings))).sampleTime;
}

IoDelayControl buildIoDelayControl(const json::Value& settings, const json::Value& target,
    const c::ODESystem& system, std::size_t streams, double controlSampleTime,
    std::optional<std::pair<double, double>> limits) {
    const auto path = modelPath(settings);
    const auto archive = c::readNpz(path);
    const auto model = c::loadIoDelayModel(archive);
    const auto metadata = json::parse(archive.metadataJson());
    if (metadata.at("method").string() != "io-delay-model")
        throw std::invalid_argument("io_delay_lqr model is not an io_delay_model archive");
    if (std::abs(controlSampleTime - model.sampleTime) > 1e-10 * model.sampleTime)
        throw std::invalid_argument("io_delay_lqr control_sample_time must equal the model clock");
    postprocess_detail::SourceData trained;
    trained.metadataJson = archive.at("source_metadata_json").string();
    c::requireCompatibleSystems(postprocess_detail::systemConfiguration(trained), system.configuration());

    const auto observable = settings.at("measurement_observable").string(metadata.at("observable").string());
    if (observable != metadata.at("observable").string())
        throw std::invalid_argument("io_delay_lqr measurement differs from the fitted observable");
    const auto measurement = system.observableCatalog().at(observable);
    if (measurement.outputDimension() != 1) throw std::invalid_argument("io_delay_lqr needs a scalar measurement");
    if (target.at("space").string() != "observable" || target.at("kind").string() != "from_state")
        throw std::invalid_argument("io_delay_lqr target must be {space: observable, kind: from_state, state_target}");
    const auto states = system.descriptor().stateDimension;
    auto stateTarget = resolveControlTarget(target.at("state_target"), states);

    const double start = settings.at("start_time").number(0.0), stop = settings.at("final_time").number(NAN);
    if (!std::isfinite(start) || !std::isfinite(stop) || !(stop > start))
        throw std::invalid_argument("io_delay_lqr needs final_time > start_time for the known-target recursion");
    const double tail = settings.at("preview_tail").number(5.0);
    if (!std::isfinite(tail) || tail < 0) throw std::invalid_argument("io_delay_lqr preview_tail must be nonnegative");
    const double dt = model.sampleTime;
    const auto ticks = static_cast<std::size_t>(std::ceil((stop + tail - start) / dt - 1e-9)) + 2;
    auto reference = std::make_shared<std::vector<double>>(ticks);
    for (std::size_t k = 0; k < ticks; ++k) {
        const auto x = stateTarget(start + static_cast<double>(k) * dt);
        (*reference)[k] = measurement.evaluateReal(x.data()).front();
        if (!std::isfinite((*reference)[k])) throw std::invalid_argument("io_delay_lqr lifted target is not finite");
    }
    auto tracker = std::make_shared<const c::IoDelayTracker>(model, settings.at("Q").number(1.0),
                                                             settings.at("R").number(1.0), *reference);
    double radius = 0.0;
    for (const auto& pole : tracker->closedLoopPoles()) radius = std::max(radius, std::abs(pole));
    const double warmup = settings.at("warmup_input").number(0.0);
    auto bank = std::make_shared<std::vector<c::IoDelayTrackerStream>>();
    bank->reserve(streams);
    for (std::size_t j = 0; j < streams; ++j) bank->emplace_back(tracker, warmup, limits);

    IoDelayControl out;
    out.observable = observable; out.sampleTime = dt; out.stateTarget = stateTarget;
    out.measuredTarget = [reference, start, dt](double time) {
        const auto k = std::llround((time - start) / dt);
        const auto index = std::clamp<long long>(k, 0, static_cast<long long>(reference->size()) - 1);
        return std::vector<double>{(*reference)[static_cast<std::size_t>(index)]};
    };
    out.law = [bank, measurement, start, dt](std::size_t stream, double time, const double* state) {
        const double y = measurement.evaluateReal(state).front();  // the only use of the state
        const auto tick = std::llround((time - start) / dt);
        if (tick < 0) throw std::runtime_error("io_delay_lqr evaluated before start_time");
        return bank->at(stream).step(static_cast<std::size_t>(tick), y);
    };
    out.evidence = json::Value(json::Value::Object{
        {"controller", json::Value("io_delay_lqr")}, {"model_archive", json::Value(path.generic_u8string())},
        {"model_sha256", json::Value("sha256:" + c::sha256File(path))},
        {"measurement_observable", json::Value(observable)},
        {"online_signals", json::Value("scalar measurement and past applied commands only")},
        {"ny", json::Value(double(model.ny))}, {"nu", json::Value(double(model.nu))},
        {"sample_time", json::Value(dt)}, {"Q", json::Value(settings.at("Q").number(1.0))},
        {"R", json::Value(settings.at("R").number(1.0))}, {"warmup_input", json::Value(warmup)},
        {"final_time", json::Value(stop)}, {"preview_tail", json::Value(tail)},
        {"closed_loop_spectral_radius", json::Value(radius)},
        {"riccati_relative_residual", json::Value(tracker->riccatiResidual())},
        {"citation", json::Value("Anderson and Moore, Optimal Control: Linear Quadratic Methods, 1990")}});
    return out;
}

}  // namespace dynsyslab::workflow::detail
