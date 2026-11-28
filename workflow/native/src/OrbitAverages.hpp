#pragma once
// Orbit records against a system's known spectrum (OrbitAverages.cpp): the Fourier coefficients of
// catalogue observables in a spectral band's exact phase along each record, and transverse weights
// for a family of records.
//
// On a periodic orbit of a shell band, a scalar observable is g = sum_n a_n e^{i n phi} with e^{i phi} the
// band's (normalized) profile; a_n are g's continuous Koopman modes there (Mezic, J. Nonlinear Sci. 30,
// 2020, section 11). They are computed as weighted Birkhoff averages
//   a_n = sum_k w(k/K) g(x_k) e^{-i n phi(x_k)} / sum_k w(k/K),  w(s) = exp(-1/(s(1 - s))),
// which converge faster than any power of the record length on periodic and quasiperiodic orbits (Das &
// Yorke, Nonlinearity 31, 2018), so a record needs no whole number of periods and no resampling.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/math/HarmonicAnalysis.hpp>

#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

struct OrbitSeries {
    std::vector<math::FourierSeries> channels;  // one per observable, angular frequency Omega
    double frequency{};                         // Omega: the band's level at the record's first state
    double invariant{std::numeric_limits<double>::quiet_NaN()};  // a conserved quantity there, when asked
};

// The library band `band` of the record's system (core::libraryBands).
[[nodiscard]] core::SpectralBand recordBand(const SourceData& record, const std::string& band);
// Fourier coefficients |n| <= harmonics of each catalogue observable in the band's phase along the record.
[[nodiscard]] OrbitSeries orbitSeries(const SourceData& record, const core::SpectralBand& band,
                                      const std::vector<std::string>& observables, int harmonics,
                                      const std::string& invariant = {});
// Transverse weights of a family of records, normalized to sum one: "uniform"; a list; or
// {"canonical_area": {"invariant": id, "band": id}}, the trapezoidal dE/Omega of the canonical area
// dE dphi/Omega (Arnold, Mathematical Methods of Classical Mechanics, section 50), with E the system's
// catalogue invariant and Omega its band's level, each at the record's first state.
[[nodiscard]] std::vector<double> recordWeights(const json::Value& rule, const std::vector<SourceData>& records);

}  // namespace dynsyslab::workflow::postprocess_detail
