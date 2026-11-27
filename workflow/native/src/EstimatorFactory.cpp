#include <dynsyslab/workflow/EstimatorFactory.hpp>
#include "ModelBindingOptions.hpp"
#include "BalancedTruncationSettings.hpp"
#include "KlqeProcessCovariance.hpp"
#include "SensorBinding.hpp"
#include "ClassicalEstimator.hpp"
#include "KlqeDesignReuse.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/ReducedOrderEstimation.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include "DerivativeSettings.hpp"
#include <memory>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {
std::vector<double> numbers(const json::Value& v) {
    std::vector<double> result;
    for (const auto& item : v.array()) result.push_back(item.number());
    return result;
}
std::vector<double> square(const json::Value& v, std::size_t n, const char* label) {
    auto values = numbers(v);
    if (values.size() == n * n) return values;
    std::vector<double> result(n * n, 0);
    if (!values.empty() && values.size() != n)
        throw std::invalid_argument(std::string(label) + " must be scalar, diagonal or square");
    for (std::size_t k = 0; k < n; ++k) result[k * n + k] = values.empty() ? v.number(1) : values[k];
    return result;
}
struct Trace {
    std::vector<double> time, lifted, output, measured, predicted, innovation, applied;
    void append(const core::KoopmanEstimate& e, const std::vector<double>& y, double u) {
        time.push_back(e.time);
        lifted.insert(lifted.end(), e.lifted.begin(), e.lifted.end());
        output.insert(output.end(), e.output.begin(), e.output.end());
        measured.insert(measured.end(), y.begin(), y.end());
        predicted.insert(predicted.end(), e.measurement.begin(), e.measurement.end());
        innovation.insert(innovation.end(), e.innovation.begin(), e.innovation.end());
        applied.push_back(u);
    }
    void store(core::NpzRecord& archive) const {
        const auto rows = time.size();
        if (!rows) return;
        archive.arrays["estimator_time"] = core::NpyArray::fromDoubles({rows}, time);
        const auto save = [&](const char* name, const std::vector<double>& values) {
            if (!values.empty()) archive.arrays[name] = core::NpyArray::fromDoubles({rows, values.size() / rows}, values);
        };
        save("lifted_estimate", lifted); save("state_estimate", output); save("measurement_values", measured);
        save("measurement_estimate", predicted); save("innovation", innovation);
        save("estimator_applied_input", applied);
    }
};
}  // namespace

KlqeBinding bindKlqe(const json::Value& spec, const core::KoopmanModel& model, const core::ODESystem& system) {
    KlqeBinding b;
    b.measurement = spec.at("measurement").string();
    if (b.measurement.empty() || !model.hasReconstruction(b.measurement))
        throw std::invalid_argument("KLQE needs the KROM reconstruction named by 'measurement'");
    b.output = spec.at("output").string();
    if (!b.output.empty() && !model.hasReconstruction(b.output))
        throw std::invalid_argument("KLQE output reconstruction is missing: " + b.output);
    const auto n = model.splitDimension(), p = model.reconstruction(b.measurement).outputs;
    auto process = klqeProcessCovariance(spec, model);  // scalar/list, or shaped from the K-ROM residual
    const auto W = std::move(process.matrix);
    b.processCovariance = std::move(process.evidence);
    const auto V = square(spec.at("measurement_covariance"), p, "measurement covariance");
    b.design = core::designKoopmanLQEFromMeasurements(model, b.measurement, W, V, b.output,
        spec.at("require_detectable").boolean(true));
    if (!spec.at("design_archive").isNull())  // a verified, previously published klqe design
        b.processCovariance = reuseKlqeDesign(b.design, spec.at("design_archive").string(), model);
    if (const auto reduction = parseBalancedTruncation(spec.at("balanced_truncation")))
        b.reduction = balancedTruncationEvidence(core::reduceKoopmanLQEDesign(b.design, *reduction, W, V).reduction);
    const auto catalog = system.observableCatalog();
    if (!spec.at("measurement_observables").array().empty())
        for (const auto& name : spec.at("measurement_observables").array()) b.sensors.push_back(catalog.at(name.string()));
    else b.sensors.push_back(catalog.at(spec.at("measurement_observable").string(b.measurement)));
    b.sensor=bindInputSensor(spec,b.sensors,catalog);
    std::size_t width = 0; for (const auto& s : b.sensors) width += s.outputDimension();
    for (const auto& v : spec.at("measurement_derivatives").array()) {
        const double k = v.number(-1);
        if (!(k >= 0) || k != std::floor(k) || k >= double(width))
            throw std::invalid_argument("measurement_derivatives must index sensor outputs");
        b.differentiated.push_back(static_cast<std::size_t>(k));
    }
    width += b.differentiated.size();
    b.derivative.method = spec.at("measurement_derivative_method").string("backward2");
    if (!spec.at("measurement_derivative").isNull()) b.derivative = parseDerivativeSpec(spec.at("measurement_derivative"));
    if (!b.differentiated.empty() && !core::isCausalDerivativeMethod(b.derivative.method))
        throw std::invalid_argument("online measurement derivatives need a causal method (backward, backward2, causal_*, levant)");
    b.noise = parseNoiseSpec(spec.at("measurement_noise"));
    const std::size_t raw = width - b.differentiated.size();
    if (b.noise.kind != "none") {
        const auto& scale = spec.at("measurement_noise").at("scale");
        if (const auto* list = std::get_if<json::Value::Array>(&scale.storage())) for (const auto& v : *list) b.noiseScale.push_back(v.number());
        else b.noiseScale.assign(raw, scale.number(b.noise.relative ? 0.0 : 1.0));
        if (b.noiseScale.size() != raw) throw std::invalid_argument("measurement_noise.scale needs one RMS per sensor output");
        for (double v : b.noiseScale) if (!(v > 0)) throw std::invalid_argument("relative measurement_noise needs a positive scale (sensor RMS)");
    } else b.noiseScale.assign(raw, 1.0);
    const auto transform = spec.at("measurement_transform").string("none");
    if (transform == "eigenfunctions") {
        b.liftSensors = true;
        if (width != model.coordinateDimension())
            throw std::invalid_argument("Eigenfunction measurement needs one sensor output per KROM coordinate");
        if (p != n) throw std::invalid_argument("Eigenfunction measurement needs a split-identity measurement reconstruction");
    } else if (transform != "none") {
        throw std::invalid_argument("measurement_transform must be none or eigenfunctions");
    } else if (width != p) throw std::invalid_argument("Simulated sensor outputs differ from C_y row count");
    b.initial = numbers(spec.at("initial_lifted"));
    if (b.initial.empty()) {
        if (spec.at("initial_estimate").string("zero") != "zero")
            throw std::invalid_argument("estimator initial_estimate must be zero (or give initial_lifted)");
        b.initial.assign(n, 0);
    }
    b.measurementOffset = numbers(spec.at("measurement_offset"));
    b.outputOffset = numbers(spec.at("output_offset"));
    b.measurementFeedthrough = numbers(spec.at("measurement_feedthrough"));
    b.feedthroughScaleOutput = static_cast<int>(spec.at("measurement_feedthrough_scale_output").number(-1));
    return b;
}

EstimatorSetup buildEstimator(const json::Value& spec, const json::Value& solver,
                              const core::ODESystem& system, std::size_t count, double defaultStep) {
    EstimatorSetup result;
    if (spec.isNull()) return result;
    const auto kind = spec.at("kind").string("klqe");
    if (!count) throw std::invalid_argument("estimator needs a trajectory");
    const double dt = spec.at("sample_time").number(solver.at("control_sample_time").number(defaultStep));
    if (!std::isfinite(dt) || dt <= 0) throw std::invalid_argument("estimator sample_time must be positive");
    if (detail::isClassicalEstimator(kind)) return detail::buildClassicalEstimator(spec, dt, system, count);
    if (kind != "klqe") throw std::invalid_argument("Unknown estimator kind: " + kind);
    const auto path = spec.at("krom").string();
    if (path.empty()) throw std::invalid_argument("KLQE estimator requires a KROM archive");
    std::vector<std::size_t> selected;
    for (double v : numbers(spec.at("selected_modes"))) {
        if (!std::isfinite(v) || v < 0 || v != std::floor(v)) throw std::invalid_argument("selected_modes must be nonnegative integers");
        selected.push_back(static_cast<std::size_t>(v));
    }
    auto model = core::loadKoopmanModel(core::wide(path), system.configuration(), selected, modelBindingOptions(spec));
    auto bound = bindKlqe(spec, model, system);
    bound.sampleTime = dt;
    const auto n = model.splitDimension(), p = model.reconstruction(bound.measurement).outputs;
    const auto& design = bound.design;
    const auto& measurementName = bound.measurement;
    const auto& outputName = bound.output;
    auto sensor = bound.sensor;
    const auto representation = spec.at("control_representation").string("none");
    struct Stream { std::unique_ptr<core::KoopmanLQEObserver> observer; Trace trace; SensorDifferentiator derivative; };
    auto streams = std::make_shared<std::vector<Stream>>(count);
    for (auto& s : *streams)
    {
        s.observer = std::make_unique<core::KoopmanLQEObserver>(model, design, representation, bound.initial,
            bound.measurementOffset, bound.outputOffset);
        s.observer->setMeasurementFeedthrough(bound.measurementFeedthrough, bound.feedthroughScaleOutput);
        s.derivative = SensorDifferentiator(bound, std::size_t(&s - &streams->front()));
    }
    result.kind = kind; result.sampleTime = dt;
    const bool lift = bound.liftSensors;
    const auto liftModel = std::make_shared<const core::KoopmanModel>(model);
    result.observe = [streams, sensor, lift, liftModel](std::size_t trajectory, double time, const double* state, double previous) {
        auto& s = streams->at(trajectory);
        auto y = sensor->evaluate(state, previous); // sampled before the next command
        s.derivative.append(y, time);  // numerically differentiated sensor outputs
        if (lift) y = liftModel->evaluateCoordinatesSplit(y.data());  // psi of the sensor values
        const auto estimate = s.observer->observe(time, y, {previous});
        s.trace.append(estimate, y, previous);
        return estimate;
    };
    result.appendDiagnostics = [streams, design, n, p](std::size_t trajectory, core::NpzRecord& archive) {
        streams->at(trajectory).trace.store(archive);
        archive.arrays["estimator_gain"] = core::NpyArray::fromDoubles({n, p}, design.solution.gain);
        archive.arrays["measurement_matrix"] = core::NpyArray::fromDoubles({p, n}, design.measurementMatrix);
        archive.arrays["observer_poles"] = core::NpyArray::fromComplex({n}, design.solution.observerPoles);
        if (design.outputDimension)
            archive.arrays["state_reconstruction_matrix"] = core::NpyArray::fromDoubles({design.outputDimension, n},
                                                                                        design.stateReconstructionMatrix);
    };
    json::Value::Object evidence{{"estimator_kind", json::Value(kind)}, {"estimator_krom", json::Value(path)},
        {"measurement_reconstruction", json::Value(measurementName)}, {"output_reconstruction", json::Value(outputName)},
        {"estimator_control_representation", json::Value(representation)}, {"estimator_sample_time", json::Value(dt)},
        {"estimator_model_identity", json::Value(core::koopmanModelIdentity(model))},
        {"observer_integration", json::Value("exact-zoh-previous-measurement-and-applied-input")}};
    if (!bound.reduction.isNull()) evidence["estimator_balanced_truncation"] = bound.reduction;
    evidence["sensor_feedthrough"] = spec.at("sensor_feedthrough");
    evidence["sensor_input_timing"] = "previous-applied-input-before-new-command";
    if (!bound.processCovariance.isNull()) evidence["estimator_process_covariance"] = bound.processCovariance;
    result.evidence = json::Value(std::move(evidence));
    return result;
}

ControllerSetup attachEstimator(ControllerSetup control, EstimatorSetup estimator, std::size_t count,
                                std::size_t stateDimension) {
    if (!estimator) return control;
    const double dt = estimator.sampleTime;
    if (control.policy && control.policy.evaluation == core::ControlEvaluationMode::sampled
        && std::abs(control.policy.sampleTime - dt) > 1e-9 * dt)
        throw std::invalid_argument("controller and estimator sample clocks differ");
    auto previous = std::make_shared<std::vector<double>>(count, 0.0);
    auto base = control.policy.value;
    auto observe = estimator.observe;
    if (stateDimension) {
        // Separation: observe first (sensors only), then run the controller on x_hat.
        control.policy.value = [base, observe, previous, stateDimension](std::size_t trajectory, double time,
                                                                       const double* state) {
            const auto estimate = observe(trajectory, time, state, previous->at(trajectory));
            if (estimate.output.size() != stateDimension)
                throw std::invalid_argument("estimated_state feedback needs an estimator output with one value per state");
            const double u = base ? base(trajectory, time, estimate.output.data()) : 0.0;
            previous->at(trajectory) = u;
            return u;
        };
    } else control.policy.value = [base, observe, previous](std::size_t trajectory, double time, const double* state) {
        const double u = base ? base(trajectory, time, state) : 0.0;
        observe(trajectory, time, state, previous->at(trajectory));
        previous->at(trajectory) = u;
        return u;
    };
    control.policy.evaluation = core::ControlEvaluationMode::sampled;
    control.policy.sampleTime = dt;
    auto first = control.appendDiagnostics;
    auto second = estimator.appendDiagnostics;
    control.appendDiagnostics = [first, second](std::size_t trajectory, core::NpzRecord& archive) {
        if (first) first(trajectory, archive);
        if (second) second(trajectory, archive);
    };
    if (!std::holds_alternative<json::Value::Object>(control.modelBindingEvidence.storage()))
        control.modelBindingEvidence = json::Value(json::Value::Object{});
    control.modelBindingEvidence.object()["estimator"] = estimator.evidence;
    if (stateDimension) control.modelBindingEvidence.object()["feedback_architecture"] = "klqe-estimated-state/v1";
    return control;
}

bool usesEstimatedState(const json::Value& job) {
    return !job.at("estimator").isNull()
        && job.at("controller").at("settings").at("feedback").string() == "estimated_state";
}

bool usesEstimatedFeedback(const json::Value& job) {
    return !job.at("estimator").isNull() && job.at("controller").at("kind").string() == "klqr"
        && job.at("controller").at("settings").at("feedback").string() == "estimate";
}

json::Value controllerSpecWithEstimator(const json::Value& job) {
    json::Value controller = job.at("controller");
    if (!usesEstimatedFeedback(job)) return controller;
    if (!controller.at("settings").at("estimator").isNull())
        throw std::invalid_argument("give the estimator either at the job level or in controller.settings, not both");
    auto estimator = job.at("estimator");
    const auto kromController = controller.at("settings").at("krom").string();
    const auto kromEstimator = estimator.at("krom").string();
    if (!kromEstimator.empty() && !kromController.empty() && kromEstimator != kromController)
        throw std::invalid_argument("KLQE -> KLQR output feedback needs one shared KROM (estimator.krom == controller krom)");
    controller.object()["settings"].object()["estimator"] = estimator;
    return controller;
}

}  // namespace dynsyslab::workflow
