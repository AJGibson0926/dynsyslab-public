#pragma once
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/NumericalDifferentiation.hpp>
#include <dynsyslab/core/System.hpp>
#include <memory>
namespace dynsyslab::workflow::detail {
// KLQR setting "online_measurement": the controller receives one scalar signal y and forms the
// KROM coordinates (y, y') with a causal differentiator. The state argument is used only to
// simulate the sensor. The input map must be a stored constant or bilinear map, since an exact
// map would need the state and the input field online.
dynsyslab::core::ScalarControlLaw measuredJetControlLaw(
    const dynsyslab::core::KoopmanLQRController& prototype,
    dynsyslab::core::Observable measurement,
    dynsyslab::core::DerivativeSpec derivative,
    double sampleTime,std::size_t streams,std::size_t warmupSamples,double warmupInput);
}
