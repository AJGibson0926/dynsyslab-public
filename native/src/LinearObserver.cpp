#include <dynsyslab/core/LinearObserver.hpp>

#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {
namespace {

void shape(const std::vector<double>& matrix, std::size_t rows, std::size_t columns, const char* name) {
    if (matrix.size() != rows * columns)
        throw std::invalid_argument(std::string(name) + " has the wrong shape for the observer");
    for (double value : matrix)
        if (!std::isfinite(value)) throw std::invalid_argument(std::string(name) + " must be finite");
}

std::vector<double> closedLoop(const LinearObserverModel& m) {
    auto result = m.A;
    for (std::size_t i = 0; i < m.states; ++i)
        for (std::size_t j = 0; j < m.states; ++j)
            for (std::size_t k = 0; k < m.outputs; ++k)
                result[i * m.states + j] -= m.L[i * m.outputs + k] * m.C[k * m.states + j];
    return result;
}

void finish(LinearObserverModel& m) {
    shape(m.A, m.states, m.states, "A");
    shape(m.B, m.states, m.inputs, "B");
    shape(m.C, m.outputs, m.states, "C");
    shape(m.L, m.states, m.outputs, "L");
    if (m.stateEquilibrium.empty()) m.stateEquilibrium.assign(m.states, 0.0);
    if (m.stateEquilibrium.size() != m.states) throw std::invalid_argument("observer x_e has the wrong size");
    if (m.measurementEquilibrium.empty()) {
        m.measurementEquilibrium.assign(m.outputs, 0.0);
        for (std::size_t k = 0; k < m.outputs; ++k)
            for (std::size_t j = 0; j < m.states; ++j)
                m.measurementEquilibrium[k] += m.C[k * m.states + j] * m.stateEquilibrium[j];
    }
    if (m.measurementEquilibrium.size() != m.outputs) throw std::invalid_argument("observer y_e has the wrong size");
    m.observerPoles = numerics::eigensystem(closedLoop(m), m.states, false).eigenvalues;
}

}  // namespace

LinearObserverModel designLinearQuadraticEstimator(
    std::vector<double> A, std::vector<double> B, std::vector<double> C,
    const std::vector<double>& W, const std::vector<double>& V,
    std::size_t states, std::size_t inputs, std::size_t outputs) {
    const auto solution = solveContinuousLQE(A, C, W, V, states, outputs, true);
    LinearObserverModel model;
    model.A = std::move(A); model.B = std::move(B); model.C = std::move(C); model.L = solution.gain;
    model.states = states; model.inputs = inputs; model.outputs = outputs;
    finish(model);
    return model;
}

LinearObserverModel luenbergerObserver(
    std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> L,
    std::size_t states, std::size_t inputs, std::size_t outputs) {
    LinearObserverModel model;
    model.A = std::move(A); model.B = std::move(B); model.C = std::move(C); model.L = std::move(L);
    model.states = states; model.inputs = inputs; model.outputs = outputs;
    finish(model);
    return model;
}

LinearObserver::LinearObserver(LinearObserverModel model, double sampleTime, std::vector<double> initial)
    : model_(std::move(model)), sampleTime_(sampleTime) {
    if (!(sampleTime > 0.0)) throw std::invalid_argument("observer sample time must be positive");
    if (initial.empty()) initial = model_.stateEquilibrium;
    if (initial.size() != model_.states) throw std::invalid_argument("observer initial state has the wrong size");
    deviation_.resize(model_.states);
    for (std::size_t i = 0; i < model_.states; ++i) deviation_[i] = initial[i] - model_.stateEquilibrium[i];
    // Driving inputs of the observer: [du; dy], held over each interval.
    const std::size_t drive = model_.inputs + model_.outputs;
    std::vector<double> input(model_.states * drive, 0.0);
    for (std::size_t i = 0; i < model_.states; ++i) {
        for (std::size_t j = 0; j < model_.inputs; ++j) input[i * drive + j] = model_.B[i * model_.inputs + j];
        for (std::size_t k = 0; k < model_.outputs; ++k)
            input[i * drive + model_.inputs + k] = model_.L[i * model_.outputs + k];
    }
    auto hold = numerics::zeroOrderHold(closedLoop(model_), input, model_.states, drive, sampleTime);
    transition_ = std::move(hold.transition);
    response_ = std::move(hold.response);
}

std::vector<double> LinearObserver::observe(double time, const std::vector<double>& y,
                                            const std::vector<double>& input) {
    if (y.size() != model_.outputs) throw std::invalid_argument("observer measurement has the wrong size");
    if (started_) {
        const double elapsed = time - lastTime_;
        if (std::abs(elapsed - sampleTime_) > 1e-6 * sampleTime_ + 1e-12)
            throw std::invalid_argument("observer measurements must arrive on its sample clock");
        const std::size_t drive = model_.inputs + model_.outputs;
        std::vector<double> held(drive, 0.0);
        for (std::size_t j = 0; j < model_.inputs && j < input.size(); ++j) held[j] = input[j];
        for (std::size_t k = 0; k < model_.outputs; ++k)
            held[model_.inputs + k] = lastMeasurement_[k] - model_.measurementEquilibrium[k];
        std::vector<double> next(model_.states, 0.0);
        for (std::size_t i = 0; i < model_.states; ++i) {
            for (std::size_t j = 0; j < model_.states; ++j) next[i] += transition_[i * model_.states + j] * deviation_[j];
            for (std::size_t j = 0; j < drive; ++j) next[i] += response_[i * drive + j] * held[j];
        }
        deviation_ = std::move(next);
    }
    started_ = true; lastTime_ = time; lastMeasurement_ = y; lastInput_ = input;
    std::vector<double> estimate(model_.states);
    for (std::size_t i = 0; i < model_.states; ++i) estimate[i] = model_.stateEquilibrium[i] + deviation_[i];
    return estimate;
}

std::vector<double> LinearObserver::predictedMeasurement() const {
    auto result = model_.measurementEquilibrium;
    for (std::size_t k = 0; k < model_.outputs; ++k)
        for (std::size_t j = 0; j < model_.states; ++j) result[k] += model_.C[k * model_.states + j] * deviation_[j];
    return result;
}

std::vector<double> StaticMapEstimator::estimate(const std::vector<double>& y) const {
    if (y.size() != outputs || map.size() != states * outputs)
        throw std::invalid_argument("static map estimator shape mismatch");
    std::vector<double> result = offset.empty() ? std::vector<double>(states, 0.0) : offset;
    for (std::size_t i = 0; i < states; ++i)
        for (std::size_t k = 0; k < outputs; ++k) result[i] += map[i * outputs + k] * y[k];
    return result;
}

}  // namespace dynsyslab::core
