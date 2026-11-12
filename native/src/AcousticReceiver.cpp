#include <dynsyslab/core/AcousticReceiver.hpp>
#include <dynsyslab/core/System.hpp>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {

double AcousticReceiver::rangeMetres() const {
    const auto& p = positionMetres;
    const double range = std::hypot(p[0], p[1], p[2]);
    if (!(range > 0) || !std::isfinite(range))
        throw std::invalid_argument("Receiver position must have a positive finite range");
    return range;
}

double AcousticReceiver::travelTimeSeconds() const {
    if (!(soundSpeedMetresPerSecond > 0) || !std::isfinite(soundSpeedMetresPerSecond))
        throw std::invalid_argument("Receiver sound speed must be positive and finite");
    return rangeMetres() / soundSpeedMetresPerSecond;
}

std::vector<Observable> compactForceDipolePressures(const AcousticReceiver& receiver,
    const std::array<AnalyticExpression,3>& force, const std::array<AnalyticExpression,3>& forceRate,
    const std::string& stateCoordinates) {
    const double range = receiver.rangeMetres();
    const double travel = receiver.travelTimeSeconds();
    const auto& p = receiver.positionMetres;
    const auto project = [&](const std::array<AnalyticExpression,3>& components) {
        return components[0]*(p[0]/range) + components[1]*(p[1]/range) + components[2]*(p[2]/range);
    };
    if (force[0].stateDimension() != forceRate[0].stateDimension())
        throw std::invalid_argument("Force and force rate must use the same state dimension");
    const double inductionScale = 1.0 / (4.0 * std::acos(-1.0) * range * range);
    return {
        Observable::fromExpressions({"dipole_induction_pressure","Dipole induction pressure","Pa",stateCoordinates,{},false},
            {project(force)*inductionScale}),
        Observable::fromExpressions({"dipole_radiation_pressure","Dipole radiation pressure","Pa",stateCoordinates,{},false},
            {project(forceRate)*(inductionScale*travel)})
    };
}

ObservableCatalog rpeReceiverPressures(const ODESystem& system, const AcousticReceiver& receiver) {
    if (system.descriptor().identifier != "rpe")
        throw std::invalid_argument("The spherical receiver catalog requires the RPE physical state");
    const double range = receiver.rangeMetres();
    const auto& config = system.configuration();
    const double radiusRatio = config.parameters.at("equilibrium_radius_m") / range;
    const double pressureScale = config.dimensionlessGroups.at("pressure_scale_pa");
    const auto r = AnalyticExpression::state(2,0), v = AnalyticExpression::state(2,1);
    const auto zero = AnalyticExpression::constant(0.,2);
    // A centered spherical pulsation exerts zero net directional force.
    const std::array<AnalyticExpression,3> zeroForce{zero,zero,zero};
    auto terms = compactForceDipolePressures(receiver,zeroForce,zeroForce,"dimensionless");
    const auto echo = system.observableCatalog().at("far_field_acoustic_echo").expressions().at(0);
    terms.push_back(Observable::fromExpressions(
        {"monopole_pressure","Leading monopole pressure","Pa","dimensionless",{},false},{echo/range}));
    terms.push_back(Observable::fromExpressions(
        {"spherical_near_field_pressure","Spherical near-field pressure correction","Pa","dimensionless",{},false},
        {-0.5*pressureScale*std::pow(radiusRatio,4)*r.power(4)*v.power(2)}));
    return ObservableCatalog(std::move(terms));
}
}
