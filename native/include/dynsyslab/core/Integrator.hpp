#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace dynsyslab::core {

// ODE methods (docs/architecture/NATIVE_SOLVER.md):
//   euler, rk4: one fixed step per output (or control-clock) interval.
//   bogackiShampine23: embedded RK 3(2), Bogacki & Shampine, Appl. Math. Lett. 2 (1989) 321-325.
//   dormandPrince45: embedded RK 5(4), Dormand & Prince, J. Comput. Appl. Math. 6 (1980) 19-26.
//   dormandPrince853: RK 8(5,3), Hairer, Norsett & Wanner, Solving ODEs I, 2nd ed. (1993), II.10.
//   radau5: Radau IIA order 5, Hairer & Wanner, Solving ODEs II, 2nd ed. (1996), IV.8.
//   bdf, adams: SUNDIALS CVODE, Hindmarsh et al., ACM TOMS 31 (2005) 363-396.
//   lsoda: automatic Adams/BDF switching after Petzold, SIAM J. Sci. Stat. Comput. 4 (1983) 136-148.
enum class IntegrationMethod {
    euler, rk4, bogackiShampine23, dormandPrince45, dormandPrince853, radau5, bdf, adams, lsoda
};

// Accepts "RK45", "dop853", "Radau", "LSODA", ... (case-insensitive); throws on unknown names.
[[nodiscard]] IntegrationMethod integrationMethodFromName(const std::string& name);
[[nodiscard]] std::string integrationMethodName(IntegrationMethod method);

struct IntegrationOptions;
struct ControllerCallLog;   // ControllerCallLog.hpp

// State handed to a checkpoint callback after an output sample.
struct IntegrationCheckpoint {
    std::size_t trajectory{};
    std::size_t sample{};          // index of the sample just stored
    double time{};
    const double* state{};
    std::size_t stateDimension{};
    // This trajectory's samples 0..sample: times[k] and values[k * stateDimension + i];
    // controls[k * inputDimension + j] when the solve is controlled (else null).
    const double* times{};
    const double* values{};
    const double* controls{};
    std::size_t inputDimension{};
};
// Return false to pause: the solve stops after this sample and NativeSolution::paused is set.
using CheckpointCallback = std::function<bool(const IntegrationCheckpoint&)>;

struct IntegrationOptions {
    double initialTime{0.0};
    double finalTime{1.0};
    double outputStep{0.001};
    // Arbitrary, strictly increasing output clock. When set it replaces the uniform
    // outputStep grid: its first entry is the initial time and its last the final time.
    std::vector<double> outputTimes;
    double initialStep{0.0};
    double maximumStep{0.0};
    double relativeTolerance{1.0e-6};
    double absoluteTolerance{1.0e-9};
    // Step budget per trajectory (accepted plus rejected internal steps).
    std::size_t maximumInternalSteps{10000000};
    IntegrationMethod method{IntegrationMethod::dormandPrince45};
    // On an integration failure of one trajectory (left the validity domain, a state guard,
    // step underflow, the step budget, a controller error) keep the samples computed so
    // far, fill the rest with NaN, record the reason and continue with the next trajectory.
    // Cancellation still propagates. Default: the failure aborts the whole solve.
    bool truncateOnFailure{false};
    // Magnitude guards checked after every accepted internal step. Empty vectors (or NaN
    // entries) leave a component unbounded; maximumStateMagnitude <= 0 disables the
    // infinity-norm guard. A violation is an integration failure.
    std::vector<double> stateLowerBounds;
    std::vector<double> stateUpperBounds;
    double maximumStateMagnitude{0.0};
    // Called after an output sample whenever checkpointInterval of model time has passed
    // since the last call (and at every sample when the interval is zero).
    CheckpointCallback checkpoint;
    double checkpointInterval{0.0};
    // Opt-in: every controller evaluation (time, state, input, step fate) is appended here.
    // The caller owns the log and sets its state and input dimensions.
    ControllerCallLog* controllerCallLog{nullptr};
    // With truncateOnFailure false, throw IntegrationFailure (carrying the samples computed
    // so far) instead of the original error.
    bool keepPartialHistoryOnFailure{false};
};

struct IntegrationProgress {
    double fraction{};
    std::size_t trajectory{};
    std::size_t trajectoryCount{};
    std::size_t acceptedSteps{};
    std::size_t rejectedSteps{};
};

struct IntegrationStatistics {
    std::size_t acceptedSteps{};
    std::size_t rejectedSteps{};
    std::size_t rightHandSideCalls{};
    std::size_t controllerCalls{};
    std::size_t methodSwitches{};
};

using ScalarControlLaw = std::function<double(std::size_t trajectory, double time,
                                               const double* state)>;
// Writes inputDimension values into input.
using VectorControlLaw = std::function<void(std::size_t trajectory, double time,
                                            const double* state, double* input)>;
enum class ControlEvaluationMode { continuousRhs, sampled };

// The input applied to the plant. A scalar law drives the system's own input channel
// (ODESystem::rhsWithInput). A vector law or a nonempty additiveMatrix uses additive-state
// actuation: xdot = f(t, x) + B u with B the n x m row-major additiveMatrix.
struct ControlPolicy {
    ScalarControlLaw value;
    VectorControlLaw vector;
    std::size_t inputDimension{1};
    std::vector<double> additiveMatrix;
    ControlEvaluationMode evaluation{ControlEvaluationMode::continuousRhs};
    double sampleTime{};

    explicit operator bool() const noexcept {
        return static_cast<bool>(value) || static_cast<bool>(vector);
    }
    void evaluate(std::size_t trajectory, double time, const double* state,
                  double* input) const;
};

using ProgressCallback = std::function<void(const IntegrationProgress&)>;

struct NativeSolution {
    std::vector<double> time;
    std::vector<double> values;
    std::vector<double> controlValues;   // trajectory x sample x inputDimension
    std::size_t trajectoryCount{};
    std::size_t sampleCount{};
    std::size_t stateDimension{};
    std::size_t inputDimension{1};
    IntegrationStatistics statistics;
    std::string coordinateSystem;
    // Per trajectory: samples actually computed (sampleCount unless truncated or paused)
    // and the failure that ended it (empty when complete).
    std::vector<std::size_t> completedSamples;
    std::vector<std::string> terminationReasons;
    // Set when a checkpoint callback paused the solve: later trajectories have zero
    // completed samples and the paused trajectory stops at its last checkpoint sample.
    bool paused{false};

    [[nodiscard]] const double* state(std::size_t trajectory, std::size_t sample) const;
    [[nodiscard]] double control(std::size_t trajectory, std::size_t sample) const;
    [[nodiscard]] TimeSeries trajectorySeries(std::size_t trajectory) const;
    [[nodiscard]] TimeSeries controlSeries(std::size_t trajectory) const;
};

// The abort-policy failure of an integration with keepPartialHistoryOnFailure: what() is the
// original reason; partial() holds every sample computed before it (the failing trajectory's
// completedSamples and terminationReasons are set; later trajectories have none).
class IntegrationFailure : public std::runtime_error {
public:
    IntegrationFailure(const std::string& reason, NativeSolution partial);
    [[nodiscard]] const NativeSolution& partial() const noexcept { return *partial_; }
private:
    std::shared_ptr<const NativeSolution> partial_;
};

[[nodiscard]] NativeSolution integrate(ODESystem& system,
                                       const std::vector<double>& initialStates,
                                       std::size_t trajectoryCount,
                                       const IntegrationOptions& options,
                                       ScalarControlLaw controller = {},
                                       ProgressCallback progress = {},
                                       const ExecutionContext* execution = nullptr);

[[nodiscard]] NativeSolution integrateControlled(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t trajectoryCount, const IntegrationOptions& options,
    ControlPolicy controller, ProgressCallback progress = {},
    const ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::core
