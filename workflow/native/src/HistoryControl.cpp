#include "HistoryControl.hpp"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <stdexcept>
#include <vector>
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
    if (!basis.rows || basis.leftSingularVectors.size() != basis.rows * basis.rank)
        throw std::invalid_argument("History transform basis is malformed");
    std::vector<double> weights;
    for (const auto column : transform.retainedCoordinates()) {
        if (column >= basis.rank || !(basis.singularValues.at(column) > 0))
            throw std::invalid_argument("History transform retains an invalid coordinate");
        weights.push_back(basis.leftSingularVectors[(basis.rows - 1) * basis.rank + column]
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
    // y[i] is the measured target at sample i >= 0; earlier samples equal y[0]. The last
    // projection is cached because the integrator can ask for one sample more than once.
    struct Cache { std::mutex lock; std::vector<double> y; bool hasLast{}; long long lastK{}; std::vector<double> last; };
    auto cache = std::make_shared<Cache>();
    const double dt = transform.sampleInterval();
    if (!(dt > 0)) throw std::invalid_argument("History-coordinate target needs a clocked transform");
    const auto& fitted = transform.basis();
    const auto columns = transform.retainedCoordinates();
    auto weights = std::make_shared<std::vector<double>>();
    weights->reserve(columns.size() * fitted.rows);
    for (const auto column : columns) for (std::size_t j = 0; j < fitted.rows; ++j)
        weights->push_back(fitted.leftSingularVectors[j * fitted.rank + column]);
    return [=](double time) {
        const auto k = std::llround(time / dt);
        std::lock_guard<std::mutex> guard(cache->lock);
        if (cache->hasLast && cache->lastK == k) return cache->last;
        const auto& basis = transform.basis();
        const auto newest = std::max(k, 0LL);
        for (auto i = static_cast<long long>(cache->y.size()); i <= newest; ++i) {
            const auto x = stateTarget(std::max(0.0, double(i) * dt));
            cache->y.push_back(measurement.evaluateReal(x.data()).front());
        }
        const auto delays = static_cast<long long>(basis.rows);
        const long long oldest = k - (delays - 1);
        const auto repeated = static_cast<std::size_t>(std::min(delays, std::max(-oldest, 0LL)));
        const double* history = repeated < basis.rows
            ? cache->y.data() + (oldest + static_cast<long long>(repeated)) : nullptr;
        std::vector<double> q;
        q.reserve(columns.size());
        for (std::size_t axis = 0; axis < columns.size(); ++axis) {
            const double* weight = weights->data() + axis * basis.rows;
            double sum = 0.0;
            for (std::size_t j = 0; j < repeated; ++j) sum += weight[j] * cache->y[0];
            for (std::size_t j = repeated; j < basis.rows; ++j) sum += weight[j] * history[j - repeated];
            q.push_back(sum / basis.singularValues[columns[axis]]);
        }
        cache->hasLast = true; cache->lastK = k; cache->last = q;
        return q;
    };
}
}
