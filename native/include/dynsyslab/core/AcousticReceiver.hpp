#pragma once

#include <dynsyslab/core/Observable.hpp>
#include <array>

namespace dynsyslab::core {
class ODESystem;

struct AcousticReceiver {
    std::array<double,3> positionMetres; // Measured from the stationary source center.
    double soundSpeedMetresPerSecond;
    [[nodiscard]] double rangeMetres() const;
    [[nodiscard]] double travelTimeSeconds() const;
};

// Compact stationary point force on the fluid, in N, and its physical-time
// derivative, in N/s. Returns separate induction and radiation pressures in Pa.
std::vector<Observable> compactForceDipolePressures(const AcousticReceiver& receiver,
    const std::array<AnalyticExpression,3>& force, const std::array<AnalyticExpression,3>& forceRate,
    const std::string& stateCoordinates);

// Explicit receiver catalog, evaluated at emission time with D > R(t).
// The monopole uses the RPE's intrinsic zero-input acceleration.
ObservableCatalog rpeReceiverPressures(const ODESystem& system, const AcousticReceiver& receiver);
}
