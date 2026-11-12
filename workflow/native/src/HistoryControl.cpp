#include "HistoryControl.hpp"
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow::detail {
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
}
