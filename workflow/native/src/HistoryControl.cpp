#include "HistoryControl.hpp"
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>
namespace dynsyslab::workflow::detail {
namespace c = dynsyslab::core;

dynsyslab::core::ScalarControlLaw historyControlLaw(
    dynsyslab::core::ScalarControlLaw coordinateLaw,
    dynsyslab::core::FrozenHistoryTransform transform,
    dynsyslab::core::Observable measurement,std::size_t count,double warmupInput) {
    if(!count||!std::isfinite(warmupInput)||measurement.outputDimension()!=1)
        throw std::invalid_argument("History control requires scalar measurements and a finite warmup input");
    auto streams=std::make_shared<std::vector<dynsyslab::core::CausalHistoryProjector>>();
    streams->reserve(count);for(std::size_t j=0;j<count;++j)streams->emplace_back(transform);
    return [streams,measurement=std::move(measurement),coordinateLaw=std::move(coordinateLaw),warmupInput](
        std::size_t stream,double time,const double* state) {
        const auto measured=measurement.evaluateReal(state);
        const auto q=streams->at(stream).push(time,measured.front());
        return q?coordinateLaw(stream,time,q->data()):warmupInput;
    };
}

std::vector<double> newestSampleWeights(const c::FrozenHistoryTransform& transform) {
    const auto& basis = transform.basis();
    if (!basis.delays || basis.leftSingularVectors.size() != basis.delays * basis.rank)
        throw std::invalid_argument("History transform basis is malformed");
    std::vector<double> weights;
    for (const auto column : transform.retainedCoordinates()) {
        if (column >= basis.rank || !(basis.singularValues.at(column) > 0))
            throw std::invalid_argument("History transform retains an invalid coordinate");
        weights.push_back(basis.leftSingularVectors[(basis.delays - 1) * basis.rank + column]
                          / basis.singularValues[column]);
    }
    return weights;
}

dynsyslab::core::ScalarControlLaw historyExactControlLaw(
    std::shared_ptr<const c::KoopmanLQRController> controller,
    std::shared_ptr<const c::KoopmanModel> model,
    c::FrozenHistoryTransform transform, c::Observable measurement,
    std::shared_ptr<const c::ODESystem> system, std::size_t count, double warmupInput) {
    if (!count || !std::isfinite(warmupInput) || measurement.outputDimension() != 1)
        throw std::invalid_argument("History control requires scalar measurements and a finite warmup input");
    if (model->inputDimension() != 1) throw std::invalid_argument("history_exact supports one scalar input");
    auto weights = newestSampleWeights(transform);
    if (weights.size() != model->coordinateDimension())
        throw std::invalid_argument("History transform coordinates differ from the KROM");
    auto streams = std::make_shared<std::vector<c::CausalHistoryProjector>>();
    streams->reserve(count); for (std::size_t j = 0; j < count; ++j) streams->emplace_back(transform);
    const auto states = system->descriptor().stateDimension;
    return [=, measurement = std::move(measurement), weights = std::move(weights)](
        std::size_t stream, double time, const double* state) {
        const auto measured = measurement.evaluateReal(state);
        const auto q = streams->at(stream).push(time, measured.front());
        if (!q) return warmupInput;
        const auto jacobian = measurement.jacobian(state);          // 1 x states
        const auto input = system->inputMatrix(state);              // states x 1
        double lgy = 0.0;
        for (std::size_t axis = 0; axis < states; ++axis) lgy += jacobian[axis].real() * input[axis];
        std::vector<double> sensitivity(weights.size());
        for (std::size_t k = 0; k < weights.size(); ++k) sensitivity[k] = weights[k] * lgy;
        const auto map = model->controlMapFromCoordinateSensitivity(q->data(), sensitivity);
        return controller->valueFromLiftedWithMap(time, model->evaluateCoordinates(q->data()), map).front();
    };
}

dynsyslab::core::StateTargetFunction historyCoordinateTarget(
    c::StateTargetFunction stateTarget, c::FrozenHistoryTransform transform, c::Observable measurement) {
    if (!stateTarget) throw std::invalid_argument("History-coordinate target needs a state target");
    if (measurement.outputDimension() != 1) throw std::invalid_argument("History target needs a scalar measurement");
    struct Cache { std::mutex lock; std::map<long long, double> y; std::map<long long, std::vector<double>> q; };
    auto cache = std::make_shared<Cache>();
    const double dt = transform.sampleInterval();
    if (!(dt > 0)) throw std::invalid_argument("History-coordinate target needs a clocked transform");
    return [=](double time) {
        const auto k = std::llround(time / dt);
        std::lock_guard<std::mutex> guard(cache->lock);
        if (const auto found = cache->q.find(k); found != cache->q.end()) return found->second;
        const auto& basis = transform.basis();
        auto sample = [&](long long index) {
            if (const auto found = cache->y.find(index); found != cache->y.end()) return found->second;
            const auto x = stateTarget(std::max(0.0, double(index) * dt));
            return cache->y[index] = measurement.evaluateReal(x.data()).front();
        };
        std::vector<double> q;
        for (const auto column : transform.retainedCoordinates()) {
            double sum = 0.0;
            for (std::size_t j = 0; j < basis.delays; ++j)
                sum += basis.leftSingularVectors[j * basis.rank + column]
                     * sample(k - static_cast<long long>(basis.delays - 1 - j));
            q.push_back(sum / basis.singularValues[column]);
        }
        return cache->q[k] = q;
    };
}
}
