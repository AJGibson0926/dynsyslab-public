#include <dynsyslab/core/ODE.hpp>

#include "Models.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

const double* TrajectoryBatch::state(std::size_t trajectory, std::size_t sample) const {
    if (trajectory >= batch || sample >= samples) throw std::out_of_range("trajectory state index");
    return values.data() + (trajectory * samples + sample) * stateDimension;
}

TimeSeries TrajectoryBatch::trajectorySeries(std::size_t trajectory) const {
    if (trajectory >= batch) throw std::out_of_range("trajectory time series index");
    return TimeSeries::view(time,
        values.data() + trajectory * samples * stateDimension,
        samples, stateDimension);
}

std::size_t modelStateDimension(const std::string& model) { return detail::stateDimension(model); }
std::size_t modelParameterCount(const std::string& model) { return detail::parameterCount(model); }
std::vector<std::string> fixedStepModels() {
    return {"duffing", "van_der_pol", "rpe_dimensionless"};
}

TrajectoryBatch integrateFixed(const FixedIntegrationRequest& request) {
    if (!std::isfinite(request.initialTime) || !std::isfinite(request.finalTime) ||
        request.finalTime <= request.initialTime)
        throw std::invalid_argument("native integration requires finite initialTime < finalTime");
    if (!std::isfinite(request.stepSize) || request.stepSize <= 0.0)
        throw std::invalid_argument("native integration stepSize must be positive and finite");
    const std::size_t dimension = detail::stateDimension(request.model);
    if (request.parameters.size() != detail::parameterCount(request.model))
        throw std::invalid_argument("native model parameter array has the wrong length");
    if (!request.batch || request.initialStates.size() != request.batch * dimension)
        throw std::invalid_argument("native initial states have inconsistent dimensions");

    TrajectoryBatch result;
    result.batch = request.batch;
    result.stateDimension = dimension;
    result.samples = static_cast<std::size_t>(
                         std::ceil((request.finalTime - request.initialTime) / request.stepSize)) + 1;
    result.time.resize(result.samples);
    result.values.resize(result.batch * result.samples * dimension);
    result.time[0] = request.initialTime;
    for (std::size_t item = 0; item < result.batch; ++item)
        std::copy_n(request.initialStates.data() + item * dimension, dimension,
                    result.values.data() + item * result.samples * dimension);

    double current = request.initialTime;
    for (std::size_t step = 0; step + 1 < result.samples; ++step) {
        const double dt = std::min(request.stepSize, request.finalTime - current);
        const double nextTime = current + dt;
        result.time[step + 1] = nextTime;
        for (std::size_t item = 0; item < result.batch; ++item) {
            const auto inputOffset = (item * result.samples + step) * dimension;
            const auto outputOffset = inputOffset + dimension;
            const double* y = result.values.data() + inputOffset;
            double* next = result.values.data() + outputOffset;
            detail::validateState(request.model, y);
            double k1[2]{};
            detail::evaluateRHS(request.model, request.parameters.data(), current, y, k1);
            if (request.method == FixedMethod::euler) {
                next[0] = y[0] + dt * k1[0];
                next[1] = y[1] + dt * k1[1];
            } else {
                const double y2[2]{y[0] + 0.5 * dt * k1[0], y[1] + 0.5 * dt * k1[1]};
                double k2[2]{};
                detail::evaluateRHS(request.model, request.parameters.data(), current + 0.5 * dt, y2, k2);
                const double y3[2]{y[0] + 0.5 * dt * k2[0], y[1] + 0.5 * dt * k2[1]};
                double k3[2]{};
                detail::evaluateRHS(request.model, request.parameters.data(), current + 0.5 * dt, y3, k3);
                const double y4[2]{y[0] + dt * k3[0], y[1] + dt * k3[1]};
                double k4[2]{};
                detail::evaluateRHS(request.model, request.parameters.data(), nextTime, y4, k4);
                next[0] = y[0] + (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
                next[1] = y[1] + (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
            }
            detail::validateState(request.model, next);
        }
        current = nextTime;
    }
    return result;
}

}  // namespace dynsyslab::core
