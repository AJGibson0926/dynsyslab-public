#pragma once
// Scalar input-output delay model (ARX) and its discrete LQ tracker.
//
//   y[k+1] = sum_{i<ny} a_i y[k-i] + sum_{j<nu} b_j u[k-j] + c
//
// Non-minimal companion realization, whose state is measured directly:
//   z[k] = (y[k..k-ny+1], u[k-1..k-nu+1]),  z[k+1] = A z[k] + B u[k] + d,  y[k] = C z[k].
//
// References: L. Ljung, System Identification: Theory for the User, 2nd ed., Prentice Hall
// (1999), for ARX least squares. M. Korda and I. Mezic, Automatica 93 (2018) 149-160, for
// linear predictors on delay-embedded outputs in control. B. D. O. Anderson and J. B. Moore,
// Optimal Control: Linear Quadratic Methods, Prentice Hall (1990), for the LQ tracker.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/StateSpace.hpp>

#include <complex>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::core {

// One record on the model clock. u[k] is held over [k, k+1). heldValid[k] is zero when the
// fine-clock input changed inside that interval, so no regression row may use it.
struct IoDelayRecord {
    std::vector<double> y, u;
    std::vector<unsigned char> heldValid;
};

// Keeps every `decimation`-th sample and flags intervals whose input was not held.
[[nodiscard]] IoDelayRecord decimateHeldRecord(const std::vector<double>& y,
    const std::vector<double>& u, std::size_t decimation, double relativeTolerance = 1e-12);

struct IoDelayModel {
    std::size_t ny{}, nu{};
    std::vector<double> a, b;
    double c{};
    double sampleTime{};
    [[nodiscard]] std::size_t stateDimension() const noexcept { return ny + nu - 1; }
};

struct IoDelayFitOptions {
    std::size_t ny{1}, nu{1};
    // Penalty lambda * rows * |theta|^2 on regressors scaled to unit deviation. The offset is
    // not penalized.
    double relativeRidge{1e-10};
    bool fitOffset{true};
};

struct IoDelayFitReport {
    std::size_t rows{}, rowsRejectedForHold{};
    double oneStepRmse{};
};

// Ridge least squares by sequential QR. Regression rows never cross a record boundary.
[[nodiscard]] std::pair<IoDelayModel, IoDelayFitReport> fitIoDelayModel(
    const std::vector<IoDelayRecord>& records, const IoDelayFitOptions& options, double sampleTime);

struct IoDelayRealization {
    StateSpaceModel system;  // discrete (A, B, C), non-minimal
    std::vector<double> offset;  // d = c e_0
};
[[nodiscard]] IoDelayRealization ioDelayRealization(const IoDelayModel& model);

struct IoDelayValidation {
    std::size_t rows{};
    double oneStepRmse{}, oneStepRelative{};   // relative to the deviation of y about its mean
    double freeRunRmse{}, freeRunRelative{};   // simulated from recorded u and the first window
    bool freeRunDiverged{};
};
[[nodiscard]] IoDelayValidation validateIoDelayModel(const IoDelayModel& model,
                                                     const IoDelayRecord& record);

void storeIoDelayModel(const IoDelayModel& model, NpzRecord& archive);
[[nodiscard]] IoDelayModel loadIoDelayModel(const NpzRecord& archive);

// Discrete LQ tracker for cost sum q (y[k] - r[k])^2 + rho u[k]^2 over a known reference.
// Gain: steady Riccati solution from designLQR. Feedforward: backward adjoint recursion
//   w = v[k+1] - P d,  f[k] = S^-1 B^T w,  v[k] = C^T q r[k] + (A - B K)^T w,
// with S = rho + B^T P B, giving u[k] = -K z[k] + f[k] (Anderson and Moore 1990, ch. 4).
class IoDelayTracker {
public:
    IoDelayTracker(const IoDelayModel& model, double stateWeight, double inputWeight,
                   const std::vector<double>& reference);
    [[nodiscard]] const IoDelayModel& model() const noexcept { return model_; }
    [[nodiscard]] const std::vector<double>& gain() const noexcept { return gain_; }
    [[nodiscard]] const std::vector<double>& feedforward() const noexcept { return feedforward_; }
    [[nodiscard]] const std::vector<std::complex<double>>& closedLoopPoles() const noexcept { return poles_; }
    [[nodiscard]] double riccatiResidual() const noexcept { return residual_; }
private:
    IoDelayModel model_;
    std::vector<double> gain_, feedforward_;
    std::vector<std::complex<double>> poles_;
    double residual_{};
};

// Online state of one trajectory: reads only y[k] and its own past commands.
class IoDelayTrackerStream {
public:
    IoDelayTrackerStream(std::shared_ptr<const IoDelayTracker> tracker, double warmupInput,
                         std::optional<std::pair<double, double>> limits = std::nullopt);
    // Repeating a tick returns the held command without changing the history.
    [[nodiscard]] double step(std::size_t tick, double measurement);
private:
    std::shared_ptr<const IoDelayTracker> tracker_;
    std::vector<double> z_;
    std::optional<std::pair<double, double>> limits_;
    double warmup_{}, last_{};
    std::size_t samples_{}, required_{}, lastTick_{};
    bool started_{};
};

}  // namespace dynsyslab::core
