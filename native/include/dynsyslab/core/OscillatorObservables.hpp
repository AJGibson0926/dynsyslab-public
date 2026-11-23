#pragma once
#include <dynsyslab/core/Observable.hpp>

namespace dynsyslab::core {
// Colbrook, Rigged DMD (2024), PENDULUM_example.m; mpEDMD (2023),
// mpEDMD_pendulum.m. Real components retain exact expression derivatives.
// Phase is in radians; velocityScale converts native velocity to the
// nondimensional argument and participates in all analytic derivatives.
enum class OscillatorEnvelope { sech, gaussianVelocity };
std::vector<Observable> phaseEnvelopeObservables(
    const std::string& name, std::size_t harmonic, OscillatorEnvelope envelope,
    const std::string& coordinates = "dimensionless", double velocityScale = 1.0);
// Colbrook spectralcomputations, chapter_11/Duffing_oscillator_convergence.m.
Observable gaussianSumObservable(const std::string& name = "gaussian_sum",
                                const std::string& coordinates = "dimensionless");
} // namespace dynsyslab::core
