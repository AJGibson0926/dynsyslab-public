#pragma once
#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <memory>
namespace dynsyslab::workflow::detail {
// Simulator measurement adapter only: the controller receives coordinates, not x.
dynsyslab::core::ScalarControlLaw historyControlLaw(
    dynsyslab::core::ScalarControlLaw coordinateLaw,
    dynsyslab::core::FrozenHistoryTransform transform,
    dynsyslab::core::Observable measurement,std::size_t streams,double warmupInput);
template<class Controller>
dynsyslab::core::ScalarControlLaw coordinateControllerBank(
    const Controller& prototype,std::size_t count) {
    auto controllers=std::make_shared<std::vector<Controller>>(count,prototype);
    return [controllers](std::size_t stream,double time,const double* q) {
        return controllers->at(stream).valueFromCoordinates(time,q).front();
    };
}
}
