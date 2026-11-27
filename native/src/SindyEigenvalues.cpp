// Koopman eigenvalues from data in the SINDy polynomial library: Kaiser, Kutz & Brunton,
// "Data-driven discovery of Koopman eigenfunctions for control", Mach. Learn.: Sci. Technol. 2:035023
// (2021), Sec. 5.1, Algorithm 1. Initial guesses are the eigenvalues of the library generator
// K = Theta^+ Gamma (generator EDMD: Klus, Nüske, Peitz, Niemann, Clementi & Schütte, Physica D
// 406:132416, 2020). See SINDy.hpp for the steps and the adaptations.
#include <dynsyslab/core/SINDy.hpp>
#include "SindyLibrary.hpp"
#include <dynsyslab/numerics/Eigensystem.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
using namespace sindy_detail;
namespace {

using Complex = std::complex<double>;

// K = A^-1 Bᵀ with A = ThetaᵀTheta, B = GammaᵀTheta (so Bᵀ = ThetaᵀGamma): Theta K ≈ Gamma.
std::vector<double> libraryGenerator(const GramBlocks& g) {
    const std::size_t n = g.n;
    std::vector<double> K(n * n);
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<Complex> a(g.A.begin(), g.A.end()), b(n);
        for (std::size_t row = 0; row < n; ++row) b[row] = g.B[column * n + row];  // column of Bᵀ
        const auto x = solveDense(std::move(a), std::move(b), n);
        for (std::size_t row = 0; row < n; ++row) K[row * n + column] = x[row].real();
    }
    return K;
}

Complex rayleigh(const std::vector<double>& K, const std::vector<Complex>& xi) {
    const std::size_t n = xi.size();
    Complex numerator{}; double denominator = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        Complex row{};
        for (std::size_t j = 0; j < n; ++j) row += K[i * n + j] * xi[j];
        numerator += std::conj(xi[i]) * row;
        denominator += std::norm(xi[i]);
    }
    return numerator / denominator;
}

// ||phi(t_k) - c e^{lambda (t_k - t_0)}|| / ||phi||, c by least squares.
double validationError(const std::vector<double>& theta, std::size_t samples, std::size_t features,
                       const std::vector<Complex>& xi, Complex lambda, const std::vector<double>& time) {
    std::vector<Complex> phi(samples), model(samples);
    Complex cross{}; double modelNorm = 0.0, phiNorm = 0.0;
    for (std::size_t k = 0; k < samples; ++k) {
        for (std::size_t f = 0; f < features; ++f) phi[k] += theta[k * features + f] * xi[f];
        model[k] = std::exp(lambda * (time[k] - time.front()));
        cross += std::conj(model[k]) * phi[k];
        modelNorm += std::norm(model[k]);
        phiNorm += std::norm(phi[k]);
    }
    if (!(phiNorm > 0) || !(modelNorm > 0) || !std::isfinite(modelNorm)) return INFINITY;
    const Complex c = cross / modelNorm;
    double error = 0.0;
    for (std::size_t k = 0; k < samples; ++k) error += std::norm(phi[k] - c * model[k]);
    return std::sqrt(error / phiNorm);
}

}  // namespace

KoopmanSINDyEigenvalueEstimate estimateKoopmanEigenvaluesInCoordinates(
    const std::vector<double>& dynamicalStates, std::size_t samples, std::size_t dynamicalStateDimension,
    const std::vector<double>& time, const Observable& coordinates, const std::vector<double>& dynamicalVelocities,
    const KoopmanSINDyOptions& options, const KoopmanSINDyEigenvalueOptions& estimate) {
    if (coordinates.stateDimension() != dynamicalStateDimension || samples < 3 || time.size() != samples)
        throw std::invalid_argument("SINDy eigenvalue estimate: inputs do not match");
    const auto values = coordinates.evaluateRealBatch(dynamicalStates, samples);
    const std::size_t dimension = coordinates.outputDimension();
    const auto velocities = dynamicalVelocities.empty() ? finiteDifference(values, samples, dimension, time)
        : coordinates.pushForwardRealBatch(dynamicalStates, dynamicalVelocities, samples);
    std::vector<double> center, scale;
    normalization(values, samples, dimension, options.normalizeStates, center, scale);
    const auto exponents = polynomialExponentTable(dimension, options.polynomialDegree, options.includeConstant);
    const std::size_t features = exponents.size() / dimension;
    const auto theta = featureMatrix(values, samples, dimension, exponents, center, scale);
    const auto gamma = featureDerivative(values, velocities, samples, dimension, exponents, center, scale);
    const auto gram = gramBlocks(theta, gamma, samples, features);
    const auto K = libraryGenerator(gram);
    const auto guesses = numerics::eigensystem(K, features, false).eigenvalues;

    KoopmanSINDyEigenvalueEstimate out;
    out.candidates = guesses.size();
    double radius = 1.0;
    for (const auto value : guesses) radius = std::max(radius, std::abs(value));
    for (const auto guess : guesses) {
        Complex lambda = guess;
        std::vector<Complex> xi;
        std::size_t steps = 0;
        for (; steps < estimate.refinementIterations; ++steps) {
            auto fitted = sparseNullVectorGram(gram, lambda, options.threshold, options.maximumIterations);
            if (!fitted) {
                std::vector<Complex> generator(samples * features);
                for (std::size_t i = 0; i < generator.size(); ++i) generator[i] = gamma[i] - lambda * theta[i];
                fitted = sparseNullVector(generator, samples, features, options.threshold, options.maximumIterations);
            }
            xi = std::move(*fitted);
            const Complex next = rayleigh(K, xi);
            const bool settled = std::abs(next - lambda) <= estimate.refinementTolerance * std::max(1.0, std::abs(lambda));
            lambda = next;
            if (settled) { ++steps; break; }
        }
        const bool trivial = std::abs(lambda) <= 1e-8 * radius;
        const bool duplicate = std::any_of(out.eigenvalues.begin(), out.eigenvalues.end(), [&](Complex kept) {
            return std::abs(kept - lambda) <= 1e-6 * std::max(1.0, std::abs(lambda));
        });
        const double error = validationError(theta, samples, features, xi, lambda, time);
        if ((trivial && !estimate.keepTrivial) || duplicate || !(error <= estimate.validationTolerance)) {
            ++out.rejected;
            continue;
        }
        out.eigenvalues.push_back(lambda);
        out.initial.push_back(guess);
        out.validationError.push_back(error);
        out.iterations.push_back(steps);
    }
    return out;
}

}  // namespace dynsyslab::core
