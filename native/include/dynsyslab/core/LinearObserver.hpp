#pragma once
// Classical state estimators for the solve loop, on a linear(ized) plant model
//   dx/dt = A dx + B du,   dy = C dx,   dx = x - x_e, dy = y - y_e.
// Luenberger observer (Luenberger 1964, "Observing the state of a linear system", IEEE Trans.
// Mil. Electron. 8:74-80) with a given gain L, or the steady-state continuous Kalman-Bucy /
// linear-quadratic estimator (Kalman & Bucy 1961, "New results in linear filtering and
// prediction theory", J. Basic Eng. 83:95-108) whose gain L = P C^T V^-1 solves the filter
// Riccati equation for process covariance W and measurement covariance V. The observer
//   dxhat/dt = (A - L C) xhat + B du + L dy
// is advanced exactly over each sample interval with u and y held (zero-order hold).
// A static map estimates x_hat = M y + c directly.

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct LinearObserverModel {
    std::vector<double> A;  // n x n, row-major
    std::vector<double> B;  // n x m (may be empty: no input)
    std::vector<double> C;  // p x n
    std::vector<double> L;  // n x p
    std::size_t states{}, inputs{}, outputs{};
    std::vector<double> stateEquilibrium;        // x_e (n), default zeros
    std::vector<double> measurementEquilibrium;  // y_e (p), default C x_e
    std::vector<std::complex<double>> observerPoles;  // eig(A - L C)
};

// The continuous LQE gain for (A, C, W, V); throws when (A, C) is not detectable.
[[nodiscard]] LinearObserverModel designLinearQuadraticEstimator(
    std::vector<double> A, std::vector<double> B, std::vector<double> C,
    const std::vector<double>& processCovariance, const std::vector<double>& measurementCovariance,
    std::size_t states, std::size_t inputs, std::size_t outputs);

// A Luenberger observer with a given gain L; validates the shapes and computes the poles.
[[nodiscard]] LinearObserverModel luenbergerObserver(
    std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> L,
    std::size_t states, std::size_t inputs, std::size_t outputs);

class LinearObserver {
public:
    LinearObserver(LinearObserverModel model, double sampleTime, std::vector<double> initialState);
    // y measured at `time`; `input` is the input applied over the interval that just ended.
    // The first call only records the measurement. Returns x_hat = x_e + dx_hat.
    std::vector<double> observe(double time, const std::vector<double>& y,
                                const std::vector<double>& input);
    [[nodiscard]] std::vector<double> predictedMeasurement() const;  // y_e + C dx_hat
    [[nodiscard]] const LinearObserverModel& model() const noexcept { return model_; }

private:
    LinearObserverModel model_;
    double sampleTime_{};
    std::vector<double> transition_, response_;  // ZOH of (A - L C, [B L])
    std::vector<double> deviation_, lastInput_, lastMeasurement_;
    double lastTime_{};
    bool started_{};
};

// x_hat = M y + c (n x p map and n offset).
struct StaticMapEstimator {
    std::vector<double> map;
    std::vector<double> offset;
    std::size_t states{}, outputs{};
    [[nodiscard]] std::vector<double> estimate(const std::vector<double>& y) const;
};

}  // namespace dynsyslab::core
