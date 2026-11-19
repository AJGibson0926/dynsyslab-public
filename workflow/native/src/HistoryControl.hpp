#pragma once
#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/System.hpp>
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

// d(xi)/d(y_newest) for xi = (U^T h)/sigma with h ordered oldest to newest.
std::vector<double> newestSampleWeights(const dynsyslab::core::FrozenHistoryTransform& transform);

// control_representation "history_exact": state-informed input map
//   B_psi = grad_xi psi(xi) * m_newest * (dy/dx g(x)),
// the history-coordinate analogue of the exact state-Jacobian map. The law
// still forms xi causally from the measured observable; it additionally reads
// the simulator state for L_g y, so it is not a measurement-only controller.
dynsyslab::core::ScalarControlLaw historyExactControlLaw(
    std::shared_ptr<const dynsyslab::core::KoopmanLQRController> controller,
    std::shared_ptr<const dynsyslab::core::KoopmanModel> model,
    dynsyslab::core::FrozenHistoryTransform transform,
    dynsyslab::core::Observable measurement,
    std::shared_ptr<const dynsyslab::core::ODESystem> system,
    std::size_t streams,double warmupInput);

// History-coordinate reference: xi_ref(t_k) is the frozen projection of the
// measured observable along the state target, sampled on the transform clock.
// Before t = 0 the state target is evaluated at t = 0 (its rest/initial state).
dynsyslab::core::StateTargetFunction historyCoordinateTarget(
    dynsyslab::core::StateTargetFunction stateTarget,
    dynsyslab::core::FrozenHistoryTransform transform,
    dynsyslab::core::Observable measurement);
}
