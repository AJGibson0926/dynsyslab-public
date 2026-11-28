#pragma once
// The compressed Koopman operator of a spectral fit and its frozen dictionary, plus the diagnostics
// a spectral fit adds. A KoopmanModel (krom) holds one as compression(): the fit's eigenvalues are its
// Spectrum, its eigenfunctions an expansion on this dictionary, its residuals on the eigenpairs. The
// compression keeps only the operator they were derived from, for what needs the operator itself:
// Rigged DMD (Colbrook, Drysdale & Horning 2024), residuals on new data (Colbrook & Townsend 2024),
// rollout through the operator, and readouts of dictionary features. EDMD-family compressions
// (Williams, Kevrekidis & Rowley 2015) need not be diagonalizable; the operator is kept as fitted.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Pseudospectrum.hpp>
#include <dynsyslab/core/SpectralDecomposition.hpp>
#include <complex>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace dynsyslab::core {
// A frozen dictionary is part of the fitted model, never rebuilt on evaluation data.
struct SpectralDictionary {
    std::string kind{"linear"}; // see OBSERVABLE_DICTIONARIES.md for all families
    std::size_t dimension{};
    bool constant{true};
    std::vector<double> center, scale, nodes; // nodes: centers x dimension
    std::vector<std::size_t> exponents;       // features x dimension
    double bandwidth{1};
    double smoothness{2}, amplitude{1}; // Matern nu; radial kernel amplitude
    void validate() const;
    std::size_t size() const;
    std::vector<double> evaluate(const std::vector<double>& points) const;
};
struct SpectralInvariance {
    std::vector<double> principalAngles;
    double projectionDefect{};
    std::size_t sourceRank{}, successorRank{};
};
struct StochasticSpectralResiduals {
    std::vector<double> sampleMean, innovation, total;
    std::size_t replicates{};
    std::string interpretation;
};
// Eigenvectors of K* acting on kernel sections are not Koopman eigenfunctions.
struct RkhsAdjointSpectrum {
    std::size_t dimension{}, nodesCount{}, rank{};
    std::vector<double> nodes, successors;
    double bandwidth{};
    std::vector<std::complex<double>> eigenvalues, sectionCoefficients;
    std::vector<double> residuals;                               // one per adjoint eigenpair
    Pseudospectrum pseudospectrum;                               // on the requested grid
    std::vector<std::complex<double>> pseudospectrumMinimizers;  // grid x nodes
    std::string assumptions;
};
// What a spectral fit computes beyond its spectrum and functions, kept on a krom under this name:
// invariance of the trial space, stochastic residuals and RKHS adjoint spectra.
struct SpectralDiagnostics {
    StochasticSpectralResiduals stochastic;
    std::vector<RkhsAdjointSpectrum> rkhsAdjoints;
    SpectralInvariance invariance;
};
// A least-squares readout of dictionary features: outputs = features^T coefficients.
struct FeatureReadout {
    std::size_t outputs{};
    std::vector<double> coefficients; // features x outputs
    double relativeTrainingError{};
};
struct KoopmanCompression {
    SpectralDictionary dictionary;
    std::string method, assumptions, coordinateContract{"measured"};
    double sampleInterval{};
    std::size_t rank{}, samples{};
    std::vector<double> basis;   // dictionary features x retained rank; the X basis is W-orthonormal
    std::vector<double> action;  // rank x rank, the coefficient action of the Koopman operator
    std::map<std::string, FeatureReadout> readouts;
    void validate() const;
    std::vector<double> reconstruct(const std::string& readout, const std::vector<double>& points) const;
};
json::Value spectralDictionaryJson(const SpectralDictionary&);
SpectralDictionary spectralDictionaryFromJson(const json::Value&);
json::Value spectralDiagnosticsJson(const SpectralDiagnostics&);
SpectralDiagnostics spectralDiagnosticsFromJson(const json::Value&);
// A continuous component, and an observable's measure on one (values; its wave packets are a spectral
// function, archived apart).
json::Value spectralBandJson(const SpectralBand&);
SpectralBand spectralBandFromJson(const json::Value&);
json::Value spectralMeasureJson(const SpectralMeasure&);
SpectralMeasure spectralMeasureFromJson(const json::Value&);
json::Value koopmanCompressionJson(const KoopmanCompression&);
KoopmanCompression koopmanCompressionFromJson(const json::Value&);
// The operator's identity: its dictionary, basis, action and readouts. References to a krom's
// operator-derived observables (feature readouts, generalized eigenfunctions) are checked against it.
std::string koopmanCompressionIdentity(const KoopmanCompression&);
}  // namespace dynsyslab::core
