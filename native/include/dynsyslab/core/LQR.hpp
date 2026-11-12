#pragma once

#include <dynsyslab/core/KoopmanModel.hpp>

#include <cstddef>
#include <functional>
#include <memory>
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

enum class KoopmanTargetSpace { dynamicalState, observableCoordinates, realLiftedValues };

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
    [[nodiscard]] std::size_t gainSolveCount() const noexcept { return gainSolveCount_; }
    [[nodiscard]] const ContinuousLQRSolution& lastSolution() const noexcept;
    [[nodiscard]] const std::string& controlRepresentation() const noexcept;

private:
    std::vector<double> feedback(double time, const std::vector<double>& lifted,
        const std::vector<double>& inputMap) const;
    std::shared_ptr<const KoopmanModel> model_;
    std::vector<double> statePenalty_;
    std::vector<double> inputPenalty_;
    StateTargetFunction target_;
    KoopmanTargetSpace targetSpace_{KoopmanTargetSpace::dynamicalState};
    std::string representation_;
    mutable ContinuousLQRSolution last_;
    mutable std::vector<double> cachedInputMap_;
    mutable std::size_t gainSolveCount_{};
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
