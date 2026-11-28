#pragma once
// The parts of a Koopman spectrum beyond its eigenvalues (Spectrum.hpp). By the Lebesgue decomposition
// of the spectral measure of an observable g, the spectrum has an atomic part, its eigenpairs, and a
// continuous part (Halmos, Lectures on Ergodic Theory, 1956; Mezic, Nonlinear Dyn. 41, 2005). The
// continuous part further splits into absolutely and singularly continuous parts; a band records
// which when its method distinguishes them, and is simply continuous otherwise. This header holds
//   the eigenfunctions of the atomic part: named representations (one output per pair) and their
//     values at the training samples;
//   the spectral bands: connected pieces of the continuous part, subsets of C parametrized by a real
//     parameter p (a support and, optionally, a grid of p), with their generalized eigenfunctions
//     psi_p, the operator-level counterpart of the eigenfunctions (not square integrable; they live
//     in a rigged Hilbert space: Gelfand & Vilenkin, Generalized Functions 4, 1964);
//   the decomposition of observables over both: an observable's Koopman modes, its coefficients c_j on
//     the eigenfunctions (g ~ offset + sum_j c_j psi_j), and its spectral measure on each band, with
//     its wave packets there (the smoothed projection of g onto the band's generalized eigenfunctions).
// The parallel is exact: an eigenpair is to its Koopman mode what a band's generalized eigenfunctions
// are to an observable's spectral measure and wave packets.
// Each fact is stored once. The atoms of g's measure are |c_j|^2 ||psi_j||^2 at the eigenvalues, so they
// are derived from the Koopman modes (Spectrum::atoms), not stored with the measure. A band lists
// no observables and an observable lists no measures: both directions are views over one table
// (Spectrum::decomposition, Spectrum::bandMeasures, Eigenpair::modes).
//
// Data methods estimate a smoothed measure: Rigged DMD's resolvent-smoothed measure of a unitary compression (Colbrook, Drysdale & Horning,
// arXiv:2405.00782, 2024; Colbrook, SIAM J. Numer. Anal. 61, 2023) and the autocorrelation Fejer
// estimate (Korda, Putinar & Mezic, Appl. Comput. Harmon. Anal. 48, 2020).
#include <dynsyslab/core/FunctionRepresentation.hpp>
#include <dynsyslab/core/Reconstruction.hpp>

#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace dynsyslab::core {

class FunctionExpansion;

// Values of every eigenfunction at the training samples the method produced them on.
struct EigenfunctionSamples {
    std::vector<double> time;                  // label time of each sample
    std::vector<std::size_t> sampleIndices;    // row of the input record each sample is labelled by
    std::vector<std::complex<double>> values;  // samples x modes, row-major
    std::size_t modes{};
    [[nodiscard]] std::size_t samples() const noexcept { return time.size(); }
};

// One representation of every eigenfunction (one output per pair), named by representationName.
struct EigenfunctionRepresentation {
    std::string name;
    std::shared_ptr<const FunctionRepresentation> function;
};

// Values of SpectralBand::continuity, when a method distinguishes them.
namespace continuity {
inline constexpr const char* absolute = "absolutely-continuous";
inline constexpr const char* singular = "singularly-continuous";
}  // namespace continuity

// How a band's parameter p maps to C: "unit-circle", p an angle and the multiplier e^{ip} of a sampled
// operator; "imaginary-axis", p a frequency and the rate ip of a generator.
namespace band_support {
inline constexpr const char* unitCircle = "unit-circle";
inline constexpr const char* imaginaryAxis = "imaginary-axis";
}  // namespace band_support

// The generalized eigenfunctions psi_p of a band, one per value of its parameter, with eigenvalue
// bandPoint(band, p). Two forms, by what is set:
//   sampled    `functions` only, one output per grid point: psi_p at each grid p;
//   level set  `level` and `functions`: psi_p(x) = delta(level(x) - p) functions(x), with `functions`
//              one output (the profile) and `level` the scalar map whose level sets carry them, so that
//              grad(profile).f = bandPoint(band, level(x)) profile. An integrable system's phase
//              e^{i n phi} with frequency n Omega(H) is of this form, level(x) = n Omega(H(x)) (Mezic,
//              J. Nonlinear Sci. 30, 2020, section 11).
//   explicit   `family`, when known: psi_p(x) as one expression in the state and p (the last
//              coordinate of its domain), with distributions such as delta(level(x) - p) written out
//              (it is 0 off its support and has no pointwise value on it).
struct GeneralizedEigenfunctions {
    std::shared_ptr<const FunctionExpansion> functions, level, family;
    [[nodiscard]] bool empty() const noexcept { return !functions; }
    [[nodiscard]] bool levelSet() const noexcept { return bool(level); }
};

// A band of the continuous part: a subset of C given by a parameter on a support, optionally sampled
// on a grid, with its generalized eigenfunctions when they are known, and how it was obtained.
struct SpectralBand {
    std::string name, support{band_support::unitCircle}, method, assumptions;
    std::string continuity;    // optional: absolutely or singularly continuous, when known; empty otherwise
    std::vector<double> grid;  // values of the parameter (angles, frequencies); optional, needed for measures
    double bandwidth{};        // smoothing width, when the estimate is smoothed
    std::size_t order{};       // kernel order (Rigged DMD) or number of lags (Fejer)
    GeneralizedEigenfunctions eigenfunctions;
};

// The point of C at parameter p of the band (a multiplier on the unit circle, a rate on the imaginary axis).
[[nodiscard]] inline std::complex<double> bandPoint(const SpectralBand& band, double p) {
    if (band.support == band_support::unitCircle) return std::polar(1.0, p);
    if (band.support == band_support::imaginaryAxis) return {0.0, p};
    throw std::invalid_argument("unknown spectral band support " + band.support);
}

// The spectral measure of output `output` of observable `observable` on band `band`:
// mu_g(A) = <E(A) g, g>, with E the operator's projection-valued spectral resolution. A spectral
// measure is always an observable's; the operator-level objects are the eigenpairs and bands.
struct SpectralMeasure {
    std::string observable;
    std::size_t output{};
    std::string band;
    std::vector<double> density;                // on the band's grid
    std::vector<std::complex<double>> moments;  // <K^k g, g>, when estimated
    double totalMass{};                         // ||g||^2
    std::string assumptions;                    // about this observable (e.g. a projected readout)
    // The wave packets of g at each grid point (the smoothed projection of g onto the band's generalized
    // eigenfunctions there), one output per grid point, on the operator's frozen dictionary; optional.
    std::shared_ptr<const FunctionExpansion> packets;
};

// An atom of an observable's measure: the pair (by stable mode id) and its mass |c_j|^2 ||psi_j||^2.
struct SpectralAtom {
    std::size_t mode{};
    double mass{};
};

// An observable's decomposition over a spectrum: pointers into it, valid while it lives unchanged.
struct ObservableDecomposition {
    const ReconstructionResult* modes{};             // Koopman modes: outputs x pairs, with offset and errors
    std::vector<const SpectralMeasure*> measures;  // one per (band, output)
};

// A measure estimated from data alone (no krom): the band it lies on and the observable's measure.
struct SpectralMeasureEstimate {
    SpectralBand band;
    SpectralMeasure measure;
};

// Keys of the functions on spectral sets, as archives name them.
[[nodiscard]] inline std::string measureKey(const std::string& band, const std::string& observable, std::size_t output) {
    return band + "/" + observable + "/" + std::to_string(output);
}
[[nodiscard]] inline std::string measureFunctionKey(const SpectralMeasure& measure) {
    return "measure/" + measureKey(measure.band, measure.observable, measure.output);
}
[[nodiscard]] inline const char* pseudospectrumFunctionKey() { return "pseudospectrum"; }
[[nodiscard]] inline std::string bandFunctionKey(const std::string& band) { return "band/" + band; }
[[nodiscard]] inline std::string bandLevelKey(const std::string& band) { return "band/" + band + "/level"; }
[[nodiscard]] inline std::string bandFamilyKey(const std::string& band) { return "band/" + band + "/family"; }

}  // namespace dynsyslab::core
