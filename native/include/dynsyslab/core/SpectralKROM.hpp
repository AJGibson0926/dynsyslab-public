#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Observable.hpp>
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
struct SpectralMeasure {
    std::string name, domain{"unit-circle"}, method, assumptions;
    std::vector<double> grid, density, atomLocations, atomMasses;
    std::vector<std::complex<double>> moments;
    double totalMass{}, bandwidth{};
    std::size_t order{};
};
struct SpectralResiduals {
    std::string scope{"empirical-training"}, definition;
    std::vector<double> eigenpair, minimum;
    std::vector<std::complex<double>> grid, minimizingCoefficients;
    std::size_t features{};
};
struct SpectralInvariance {
    std::vector<double> principalAngles;
    double projectionDefect{};
    std::size_t sourceRank{}, successorRank{};
};
struct SpectralWavePackets {
    std::string observable, dictionaryIdentity, assumptions;
    std::vector<double> grid;
    std::vector<std::complex<double>> coefficients; // grid x features
    std::size_t features{}, order{};
    double bandwidth{};
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
    SpectralResiduals residuals;
    std::string assumptions;
};
struct SpectralData {
    std::vector<SpectralMeasure> measures;
    StochasticSpectralResiduals stochastic;
    std::vector<RkhsAdjointSpectrum> rkhsAdjoints;
    SpectralResiduals residuals;
    SpectralInvariance invariance;
    std::vector<SpectralWavePackets> wavePackets;
};
struct SpectralReconstruction {
    std::size_t outputs{};
    std::vector<double> coefficients; // features x outputs (dictionary readout)
    double relativeTrainingError{};
};
// A spectral KROM stores an operator in a frozen dictionary and separate spectral
// data. It need not be diagonalizable. Continuous spectral densities are not modes.
struct SpectralKROM {
    SpectralDictionary dictionary;
    std::string method, assumptions, coordinateContract{"measured"};
    double sampleInterval{};
    std::size_t rank{}, samples{};
    std::vector<double> basis;       // dictionary features x retained rank; X basis is W-orthonormal
    std::vector<double> action;      // rank x rank, coefficient action of Koopman
    std::vector<std::complex<double>> multipliers, eigenfunctions; // features x modes
    SpectralData spectral;
    std::map<std::string, SpectralReconstruction> reconstructions;
    void validate() const;
    std::vector<std::complex<double>> evaluateEigenfunctions(const std::vector<double>& points) const;
    std::vector<double> reconstruct(const std::string& name,const std::vector<double>& points) const;
};
json::Value spectralDictionaryJson(const SpectralDictionary&);
SpectralDictionary spectralDictionaryFromJson(const json::Value&);
json::Value spectralDataJson(const SpectralData&);
SpectralData spectralDataFromJson(const json::Value&);
json::Value spectralKromJson(const SpectralKROM&);
SpectralKROM spectralKromFromJson(const json::Value&);
} // namespace dynsyslab::core
