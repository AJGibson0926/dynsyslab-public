#pragma once

#include <dynsyslab/core/EstimationRuntime.hpp>
#include <dynsyslab/core/NumericalDifferentiation.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <functional>
#include <string>

namespace dynsyslab::workflow {

// The Estimator module, the counterpart of ControlFactory.
//
// A solve job may carry a top-level "estimator" block alongside "controller":
//
//   "estimator": {
//     "kind": "klqe",
//     "krom": "<KROM archive>",              // the estimator's own KROM
//     "selected_modes": [...],               // optional subset of that KROM
//     "measurement": "<reconstruction name>",// C_y (e.g. far_field_acoustic_echo)
//     "measurement_observables": [...],      // simulated sensors (default: the measurement name)
//     "measurement_noise": {noise spec, "scale": [rms per sensor output]},   // optional sensor noise
//     "measurement_derivatives": [k, ...],   // optional: append d/dt of sensor output k, computed
//                                            // online from the (noisy) samples by a causal differentiator
//     "measurement_derivative": {"method": "causal_kalman", "cutoff": 20, ...},  // default backward2
//     "measurement_derivative_method": "backward2" | "backward",   // shorthand for the method only
//     "output": "<reconstruction name>",     // optional decoded output C_x (e.g. state)
//     "process_covariance": 1, "measurement_covariance": 1e-4,   // scalar, diagonal or square
//     "process_covariance": {"from_records": [...], ...},          // or shaped from the K-ROM residual
//     "process_covariance_input_weight": q2,                      // optional + q2 B B^T (KlqeProcessCovariance.hpp)
//     "initial_lifted": [...] | "initial_estimate": "zero",
//     "control_representation": "none" | "<input map name>",     // how the estimator predicts inputs
//     "sample_time": 0.001
//   }
//
// The estimator only reads the plant through its simulated sensors. It never
// evaluates an eigenfunction at the plant state. Its traces are archived with the
// solution as estimator_time, lifted_estimate, state_estimate, measurement_values,
// measurement_estimate, innovation and estimator_applied_input.
//
// Composition with the controller (Solve):
//  - no controller: the estimator runs on the free plant (zero input);
//  - prescribed input (waveform, random_pulse) or a state-feedback controller: the
//    input is held on the estimator clock and the estimator sees the applied input;
//  - KLQR with settings.feedback = "estimate": output feedback. The controller acts
//    on the estimator's lifted estimate (KLQE -> KLQR) and the KROM is shared.
//  - KLQR with feedback = "estimate" and settings.applied_controller (shadow): the plant gets the
//    applied_controller's drive (true state); the KLQE predicts with it and the KLQR force is
//    archived as commanded_input only.
struct EstimatorSetup {
    std::string kind{"none"};
    double sampleTime{};
    // Called once per estimator sample with the plant state (sensor simulation only)
    // and the input that was applied over the interval that just ended.
    std::function<core::KoopmanEstimate(std::size_t trajectory, double time, const double* state,
                                        double previousAppliedInput)> observe;
    std::function<void(std::size_t, core::NpzRecord&)> appendDiagnostics;
    json::Value evidence;

    explicit operator bool() const noexcept { return static_cast<bool>(observe); }
};

// One KLQE bound to a KROM and to the plant's simulated sensors: the design (C_y, C_x,
// gain, intercepts), the sensors, and the initial lifted estimate. Shared by the
// standalone estimator and by KLQE -> KLQR output feedback.
struct KlqeBinding {
    core::KoopmanLQEDesign design;
    std::vector<core::Observable> sensors;
    std::vector<double> initial, measurementOffset, outputOffset;
    std::vector<double> measurementFeedthrough;  // optional D (p x m): y = y0 + C_y z + D u
    int feedthroughScaleOutput = -1;             // optional: scale D by decoded output component
    // measurement_transform "eigenfunctions": the sensors read the KROM's own coordinate
    // observables and the measurement is psi(sensor values) in split real form, so the
    // measurement reconstruction must be the identity on the split lifted state.
    bool liftSensors = false;
    // measurement_derivatives: sensor outputs whose time derivative is appended to y, formed
    // online from the sampled (and possibly noisy) sensor values by a causal differentiator.
    std::vector<std::size_t> differentiated;
    core::DerivativeSpec derivative;
    // measurement_noise: added to every raw sensor output before differentiation.
    core::NoiseSpec noise;
    std::vector<double> noiseScale;
    double sampleTime = 0;  // set by the caller (estimator clock)
    std::string measurement, output;
    // "balanced_truncation": the design's gain is the reduced-order KLQE gain lifted
    // to the full KROM; this records the reduction (null when absent).
    json::Value reduction;
    // How W was formed when it is not a plain scalar/list (records residual, input weight).
    json::Value processCovariance;
};

// Per-trajectory sensor chain: measurement noise on the raw outputs, then the causal derivatives
// of KlqeBinding::differentiated appended in order.
struct SensorDifferentiator {
    std::vector<std::size_t> indices;
    std::vector<core::CausalDifferentiator> differentiators;
    std::vector<core::StreamingNoise> noise;
    SensorDifferentiator() = default;
    explicit SensorDifferentiator(const KlqeBinding& b, std::size_t stream = 0) : indices(b.differentiated) {
        if (!b.noise.kind.empty() && b.noise.kind != "none")
            for (std::size_t k = 0; k < b.noiseScale.size(); ++k) {
                auto spec = b.noise; spec.seed += 1000003u * k + 7919u * stream;
                noise.emplace_back(spec, b.sampleTime, b.noiseScale[k], 0.0);
            }
        for (const auto k : indices)
            differentiators.emplace_back(b.derivative, b.sampleTime, k < b.noiseScale.size() ? b.noiseScale[k] : 1.0);
    }
    void append(std::vector<double>& y, double /*time*/) {
        for (std::size_t k = 0; k < noise.size() && k < y.size(); ++k) y[k] = noise[k].apply(y[k]);
        std::vector<double> d; d.reserve(indices.size());
        for (std::size_t i = 0; i < indices.size(); ++i) d.push_back(differentiators[i].step(y.at(indices[i])));
        y.insert(y.end(), d.begin(), d.end());
    }
};

[[nodiscard]] KlqeBinding bindKlqe(const json::Value& estimatorSpec, const core::KoopmanModel& model,
                                   const core::ODESystem& system);

[[nodiscard]] EstimatorSetup buildEstimator(
    const json::Value& estimatorSpec, const json::Value& solverSpec,
    const core::ODESystem& system, std::size_t trajectoryCount, double defaultSampleTime);

// Runs the estimator alongside an existing input policy (or none), on the estimator clock.
// With stateDimension > 0 (controller.settings.feedback = "estimated_state") the policy is
// evaluated at the estimator's decoded output x_hat instead of the plant state: any
// state-feedback controller (KLQR on its own eigenfunctions, LQR, PID, ...) then runs on a
// state estimated from the estimator's sensors only. Estimation and control may use
// different KROMs and coordinates.
[[nodiscard]] ControllerSetup attachEstimator(ControllerSetup control, EstimatorSetup estimator,
                                              std::size_t trajectoryCount, std::size_t stateDimension = 0);

// True when the job asks the controller to act on the estimated state x_hat.
[[nodiscard]] bool usesEstimatedState(const json::Value& job);

// True when the job asks for KLQE -> KLQR output feedback.
[[nodiscard]] bool usesEstimatedFeedback(const json::Value& job);

// The controller spec that output feedback is built from: the job's estimator block is
// passed to the KLQR controller as settings.estimator (the composed runtime lives in
// EstimatedControl). Other jobs return the controller spec unchanged.
[[nodiscard]] json::Value controllerSpecWithEstimator(const json::Value& job);

}  // namespace dynsyslab::workflow
