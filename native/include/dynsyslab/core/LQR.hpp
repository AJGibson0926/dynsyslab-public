#pragma once

#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/PIDCompensator.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ContinuousLQRSolution {
    std::vector<double> riccati; // state x state
    std::vector<double> gain;    // input x state
    std::vector<std::complex<double>> closedLoopPoles;
    double relativeResidual{};
    std::size_t stateDimension{};
    std::size_t inputDimension{};
};

[[nodiscard]] ContinuousLQRSolution solveContinuousLQR(
    const std::vector<double>& generator, const std::vector<double>& input,
    const std::vector<double>& statePenalty, const std::vector<double>& inputPenalty,
    std::size_t stateDimension, std::size_t inputDimension);

// LQR with integral action (LQI). The state is augmented with q, q' = C_i x - y_ref, and
// the gain [K_x K_q] (input x (state + outputs)) minimizes the cost with penalty
// blkdiag(Q, Q_i). Control law: u = -K_x (x - x_ref) - K_q q. (A, B) must be stabilizable
// and (A, B, C_i) must have no transmission zero at s = 0.
[[nodiscard]] ContinuousLQRSolution solveContinuousLQI(
    const std::vector<double>& generator, const std::vector<double>& input,
    const std::vector<double>& outputMatrix, const std::vector<double>& statePenalty,
    const std::vector<double>& integralPenalty, const std::vector<double>& inputPenalty,
    std::size_t stateDimension, std::size_t inputDimension, std::size_t outputDimension);

enum class KoopmanTargetSpace { dynamicalState, observableCoordinates, realLiftedValues };

// A linear output of the real split lifted state, y = offset + C z, compared with
// selected components of a dynamical-state target (e.g. a KROM "state" readout row and
// the target radius). Shared by the KLQR integral-action and PID adapters.
struct LiftedOutputChannel {
    std::vector<double> outputMatrix;          // outputs x split dimension
    std::vector<double> outputOffset;          // outputs
    std::vector<std::size_t> targetComponents; // outputs, indices into the state target
};

struct KoopmanIntegralAction {
    LiftedOutputChannel channel;
    std::vector<double> penalty;          // outputs x outputs (Q_i)
    std::vector<double> integratorLimit;  // per output |q| <= limit (empty = none)
};

struct KoopmanPIDCorrection {
    LiftedOutputChannel channel;
    std::vector<double> proportional, integral, derivative;  // inputs x outputs
    double derivativeFilterTime{0.0};
    std::vector<double> integratorLimit;
    double outputMinimum{-std::numeric_limits<double>::infinity()};
    double outputMaximum{std::numeric_limits<double>::infinity()};
    bool conditionalIntegration{true};
};

using StateTargetFunction = std::function<std::vector<double>(double)>;

class KoopmanLQRController {
public:
    KoopmanLQRController(std::shared_ptr<const KoopmanModel> model,
                         std::vector<double> statePenalty,
                         std::vector<double> inputPenalty,
                         StateTargetFunction target,
                         std::string controlRepresentation = "exact",
                         KoopmanTargetSpace targetSpace = KoopmanTargetSpace::dynamicalState);

    [[nodiscard]] std::vector<double> value(double time, const double* modelState) const;
    [[nodiscard]] std::vector<double> valueFromCoordinates(double time, const double* coordinates) const;
    [[nodiscard]] std::vector<double> valueFromLifted(double time,
        const std::vector<std::complex<double>>& values) const;
    // Feedback with a caller-supplied complex input map (mode x input). Intended
    // for representation "external", whose map is computed per evaluation by the
    // caller (for example from coordinate sensitivities); the gain is re-solved
    // whenever the map changes.
    [[nodiscard]] std::vector<double> valueFromLiftedWithMap(double time,
        const std::vector<std::complex<double>>& values,
        const std::vector<std::complex<double>>& inputMap) const;
    [[nodiscard]] std::size_t gainSolveCount() const noexcept { return gainSolveCount_; }
    [[nodiscard]] const ContinuousLQRSolution& lastSolution() const noexcept;
    [[nodiscard]] const std::string& controlRepresentation() const noexcept;
    // KLQR option `balanced_truncation`: every gain solve balances (A, B, I),
    // designs LQR on the reduced model and applies the lifted gain K_r P to the
    // full lifted state. lastSolution() then holds the lifted gain and the
    // full-KROM closed-loop poles. A solve whose lifted closed loop is not
    // asymptotically stable throws.
    void setBalancedReduction(BalancedReductionOptions options);
    [[nodiscard]] const std::optional<BalancedReductionOptions>& balancedReduction() const noexcept { return reduction_; }
    [[nodiscard]] std::size_t lastReducedOrder() const noexcept { return lastReducedOrder_; }
    [[nodiscard]] double lastReductionErrorBound() const noexcept { return lastReductionErrorBound_; }
    // Integral action (LQI) on a lifted output: q' = y(z) - y_target, gain designed on the
    // augmented KROM. Stateful: evaluate at monotonic sampled times, one controller per
    // trajectory. Requires a dynamical-state target and no balanced truncation.
    void setIntegralAction(KoopmanIntegralAction action);
    // PID correction on a lifted output, added to the KLQR input: u += PID(y_target - y(z)).
    void setPIDCorrection(KoopmanPIDCorrection correction);
    [[nodiscard]] bool stateful() const noexcept { return integral_.has_value() || pid_.has_value(); }
    [[nodiscard]] const std::vector<double>& integralState() const noexcept { return integralState_; }

private:
    std::vector<double> feedback(double time, const std::vector<double>& lifted,
        const std::vector<double>& inputMap) const;
    void solveGain(const std::vector<double>& inputMap) const;
    std::optional<BalancedReductionOptions> reduction_;
    mutable std::size_t lastReducedOrder_{};
    mutable double lastReductionErrorBound_{};
    std::shared_ptr<const KoopmanModel> model_;
    std::vector<double> statePenalty_;
    std::vector<double> inputPenalty_;
    StateTargetFunction target_;
    KoopmanTargetSpace targetSpace_{KoopmanTargetSpace::dynamicalState};
    std::string representation_;
    mutable ContinuousLQRSolution last_;
    mutable std::vector<double> cachedInputMap_;
    mutable std::size_t gainSolveCount_{};
    std::optional<KoopmanIntegralAction> integral_;
    mutable std::vector<double> integralState_;
    mutable double integralTime_{};
    mutable bool integralStarted_{};
    mutable std::optional<PIDCompensator> pid_;
    std::optional<LiftedOutputChannel> pidChannel_;
};

class KoopmanMPCController {
public:
    KoopmanMPCController(std::shared_ptr<const KoopmanModel> model,
                         std::vector<double> statePenalty,
                         std::vector<double> inputPenalty,
                         StateTargetFunction target,
                         std::size_t horizon = 20,
                         double sampleTime = 1.0e-3,
                         std::string controlRepresentation = "exact",
                         std::vector<double> controlLimits = {},
                         KoopmanTargetSpace targetSpace = KoopmanTargetSpace::dynamicalState);

    [[nodiscard]] std::vector<double> value(double time, const double* modelState) const;
    [[nodiscard]] std::vector<double> valueFromCoordinates(double time, const double* coordinates) const;
    [[nodiscard]] std::vector<double> valueFromLifted(double time,
        const std::vector<std::complex<double>>& values) const;
    [[nodiscard]] const std::vector<double>& lastGain() const noexcept;

private:
    std::vector<double> feedback(double time, const std::vector<double>& lifted,
        const std::vector<double>& inputMap) const;
    std::shared_ptr<const KoopmanModel> model_;
    std::vector<double> statePenalty_,inputPenalty_;
    StateTargetFunction target_;
    KoopmanTargetSpace targetSpace_{KoopmanTargetSpace::dynamicalState};
    std::size_t horizon_{};
    double sampleTime_{};
    std::string representation_;
    std::vector<double> controlLimits_;
    mutable std::vector<double> lastGain_;
};

}  // namespace dynsyslab::core
