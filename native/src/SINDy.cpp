#include <dynsyslab/core/SINDy.hpp>
#include "PolynomialEvaluation.hpp"
#include "SindyLibrary.hpp"
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Nullspace.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace dynsyslab::core {
using namespace sindy_detail;


namespace {

// The library the generator equation is written in: the normalized polynomial table, or a given
// ExpansionDictionary on the raw coordinates.
struct Library {
    std::vector<std::size_t> exponents;
    std::vector<double> center, scale;
    std::shared_ptr<const ExpansionDictionary> dictionary;
    std::size_t dimension{}, features{};

    std::vector<double> theta(const std::vector<double>& q, std::size_t rows) const {
        if (!dictionary) return featureMatrix(q, rows, dimension, exponents, center, scale);
        std::vector<double> out(rows * features);
        for (std::size_t k = 0; k < rows; ++k) {
            const auto row = dictionary->evaluate(q.data() + k * dimension);
            std::copy(row.begin(), row.end(), out.begin() + std::ptrdiff_t(k * features));
        }
        return out;
    }
    std::vector<double> gamma(const std::vector<double>& q, const std::vector<double>& v, std::size_t rows,
                              const KoopmanSINDyOptions& options) const {
        if (options.directionalDerivative == "finite-difference") {
            // Directional central difference along the vector field (step relative to the state size).
            std::vector<double> plus(q.size()), minus(q.size()), steps(rows);
            for (std::size_t k = 0; k < rows; ++k) {
                double largest = 1.0;
                for (std::size_t i = 0; i < dimension; ++i) largest = std::max(largest, std::abs(q[k*dimension+i]));
                steps[k] = options.directionalStep * largest;
                for (std::size_t i = 0; i < dimension; ++i) {
                    plus[k*dimension+i] = q[k*dimension+i] + steps[k] * v[k*dimension+i];
                    minus[k*dimension+i] = q[k*dimension+i] - steps[k] * v[k*dimension+i];
                }
            }
            auto out = theta(plus, rows);
            const auto low = theta(minus, rows);
            for (std::size_t k = 0; k < rows; ++k)
                for (std::size_t f = 0; f < features; ++f)
                    out[k*features+f] = (out[k*features+f] - low[k*features+f]) / (2.0 * steps[k]);
            return out;
        }
        if (options.directionalDerivative != "analytic")
            throw std::invalid_argument("SINDy directional derivative must be analytic or finite-difference");
        if (!dictionary) return featureDerivative(q, v, rows, dimension, exponents, center, scale);
        std::vector<double> out(rows * features);
        for (std::size_t k = 0; k < rows; ++k) {
            const auto g = dictionary->gradient(q.data() + k * dimension);  // features x dimension
            for (std::size_t f = 0; f < features; ++f)
                for (std::size_t i = 0; i < dimension; ++i) out[k*features+f] += g[f*dimension+i] * v[k*dimension+i];
        }
        return out;
    }
};

// Trajectory lengths covering the rows (one trajectory when none are given).
std::vector<std::size_t> segmentLengths(const KoopmanSINDyOptions& options, std::size_t samples) {
    if (options.segments.empty()) return {samples};
    std::size_t total = 0;
    for (const auto length : options.segments) {
        if (!length) throw std::invalid_argument("SINDy trajectory segments must not be empty");
        total += length;
    }
    if (total != samples) throw std::invalid_argument("SINDy trajectory segments do not cover the samples");
    return options.segments;
}

// The smallest right singular vector of the unnormalized generator, entries below the absolute
// threshold zeroed ("single-absolute": one SVD, no refit).
std::vector<std::complex<double>> singleAbsoluteVector(const std::vector<double>& theta, const std::vector<double>& gamma,
                                                       std::size_t rows, std::size_t features,
                                                       std::complex<double> lambda, double threshold) {
    std::vector<std::complex<double>> generator(rows * features);
    for (std::size_t index = 0; index < generator.size(); ++index) generator[index] = gamma[index] - lambda * theta[index];
    auto vector = numerics::smallestRightSingularVector(generator, rows, features);
    for (auto& value : vector) if (std::abs(value) < threshold) value = 0.0;
    if (std::all_of(vector.begin(), vector.end(), [](auto z) { return z == std::complex<double>{}; }))
        throw std::invalid_argument("single-absolute SINDy threshold removed every coefficient");
    return vector;
}

}  // namespace

EigenfunctionFitResult fitKoopmanEigenfunctions(
    const std::vector<double>& stateValues, std::size_t samples, std::size_t stateDimension,
    const std::vector<double>& time, const std::vector<std::complex<double>>& eigenvalues,
    const std::vector<double>& velocityValues, const KoopmanSINDyOptions& options) {
    if (samples < 3 || !stateDimension || stateValues.size() != samples * stateDimension
        || time.size() != samples || eigenvalues.empty() || (!options.polynomialDegree && !options.dictionary)
        || !options.maximumIterations || options.threshold < 0.0)
        throw std::invalid_argument("invalid Koopman-with-SINDy inputs");
    if (options.sparsification != "iterative-relative" && options.sparsification != "single-absolute")
        throw std::invalid_argument("SINDy sparsification must be iterative-relative or single-absolute");
    const auto segments = segmentLengths(options, samples);
    // Time increases within each trajectory; the sampling step is the common one when every step agrees.
    double dt = 0.0;
    bool uniform = true;
    for (std::size_t start = 0, s = 0; s < segments.size(); start += segments[s++])
        for (std::size_t k = start + 1; k < start + segments[s]; ++k) {
            const double step = time[k] - time[k-1];
            if (!(step > 0)) throw std::invalid_argument("SINDy time must increase");
            if (dt == 0.0) dt = step;
            else if (std::abs(step - dt) > 1e-7 * dt) uniform = false;
        }
    Library library;
    library.dimension = stateDimension;
    library.dictionary = options.dictionary;
    if (library.dictionary) {
        if (library.dictionary->dimension() != stateDimension)
            throw std::invalid_argument("SINDy dictionary dimension differs from the coordinates");
        library.features = library.dictionary->size();
    } else {
        normalization(stateValues, samples, stateDimension, options.normalizeStates, library.center, library.scale);
        if (!options.rescale.empty()) {  // q' = (q - c) / s * m, stored as the scale s / m
            if (options.rescale.size() != 1 && options.rescale.size() != stateDimension)
                throw std::invalid_argument("SINDy rescale needs one value or one per component");
            for (std::size_t i = 0; i < stateDimension; ++i) {
                const double m = options.rescale[options.rescale.size() == 1 ? 0 : i];
                if (!(m > 0) || !std::isfinite(m)) throw std::invalid_argument("SINDy rescale must be positive and finite");
                library.scale[i] /= m;
            }
        }
        library.exponents = polynomialExponentTable(stateDimension, options.polynomialDegree, options.includeConstant);
        library.features = library.exponents.size() / stateDimension;
    }
    const auto theta = library.theta(stateValues, samples);
    std::vector<double> velocities;
    if (velocityValues.empty()) {
        velocities.resize(stateValues.size());
        for (std::size_t start = 0, s = 0; s < segments.size(); start += segments[s++]) {
            if (segments[s] < 3) throw std::invalid_argument("SINDy finite differences need three samples per trajectory");
            const auto first = stateValues.begin() + std::ptrdiff_t(start * stateDimension);
            const std::vector<double> part(first, first + std::ptrdiff_t(segments[s] * stateDimension));
            const std::vector<double> clock(time.begin() + std::ptrdiff_t(start), time.begin() + std::ptrdiff_t(start + segments[s]));
            const auto rates = finiteDifference(part, segments[s], stateDimension, clock);
            std::copy(rates.begin(), rates.end(), velocities.begin() + std::ptrdiff_t(start * stateDimension));
        }
    } else {
        if (velocityValues.size() != stateValues.size()) throw std::invalid_argument("SINDy velocity shape mismatch");
        velocities = velocityValues;
    }
    const auto gamma = library.gamma(stateValues, velocities, samples, options);
    const auto& exponents = library.exponents;
    const auto& center = library.center;
    const auto& scale = library.scale;
    struct {
        std::size_t featureCount{}, directFits{};
        std::vector<std::complex<double>> coefficients;  // eigenvalue x feature
        std::vector<double> residualRms, residualMaximum;
    } result;
    result.featureCount = library.features;
    result.coefficients.resize(eigenvalues.size() * result.featureCount);
    result.residualRms.resize(eigenvalues.size()); result.residualMaximum.resize(eigenvalues.size());
    std::vector<bool> assigned(eigenvalues.size(), false);
    const bool anchored = options.anchorStates && options.anchorValues;
    std::optional<GramBlocks> gram;
    const bool single = options.sparsification == "single-absolute";
    if (single && anchored) throw std::invalid_argument("single-absolute sparsification has no anchored form");
    if ((options.gramNullspace && !single) || anchored) gram = gramBlocks(theta, gamma, samples, result.featureCount);
    // Anchored fit: Theta at the anchor states, its Gram matrix, and the reference values (see SINDy.hpp).
    const std::size_t F = result.featureCount;
    std::vector<double> anchorTheta, anchorGram;
    std::size_t anchors = 0;
    if (anchored) {
        anchors = options.anchorStates->size() / stateDimension;
        if (!anchors || options.anchorStates->size() != anchors * stateDimension
            || options.anchorValues->size() != anchors * eigenvalues.size())
            throw std::invalid_argument("anchored SINDy: anchor states and values do not match the eigenvalues");
        anchorTheta = library.theta(*options.anchorStates, anchors);
        anchorGram.assign(F * F, 0.0);
        for (std::size_t k = 0; k < anchors; ++k)
            for (std::size_t i = 0; i < F; ++i)
                for (std::size_t j = 0; j < F; ++j) anchorGram[i*F+j] += anchorTheta[k*F+i] * anchorTheta[k*F+j];
    }
    for (std::size_t mode = 0; mode < eigenvalues.size(); ++mode) {
        if (assigned[mode]) continue;
        std::optional<std::vector<std::complex<double>>> fitted;
        if (anchored) {  // each mode against its own reference values; no conjugate reuse
            std::vector<std::complex<double>> rhs(F);
            for (std::size_t k = 0; k < anchors; ++k)
                for (std::size_t i = 0; i < F; ++i)
                    rhs[i] += anchorTheta[k*F+i] * (*options.anchorValues)[k * eigenvalues.size() + mode];
            fitted = anchoredVector(anchorGram, *gram, eigenvalues[mode], rhs, options.anchorWeight,
                                    options.threshold, options.maximumIterations);
        } else if (single) {
            fitted = singleAbsoluteVector(theta, gamma, samples, F, eigenvalues[mode], options.threshold);
        } else if (options.gramNullspace) fitted = sparseNullVectorGram(*gram, eigenvalues[mode],
                                                options.threshold, options.maximumIterations);
        if (!fitted) {  // tall-SVD reference path (also the fallback for unresolved Gram modes)
            std::vector<std::complex<double>> generator(samples*result.featureCount);
            for(std::size_t index=0;index<generator.size();++index)
                generator[index]=gamma[index]-eigenvalues[mode]*theta[index];
            fitted = sparseNullVector(generator,samples,
                result.featureCount,options.threshold,options.maximumIterations);
        }
        const auto& coefficient = *fitted;
        std::copy(coefficient.begin(), coefficient.end(),
                  result.coefficients.begin() + static_cast<std::ptrdiff_t>(mode * result.featureCount));
        assigned[mode] = true; ++result.directFits;
        if (!anchored && options.reuseConjugateFits && std::abs(eigenvalues[mode].imag()) > 1e-12)
            for (std::size_t partner = mode + 1; partner < eigenvalues.size(); ++partner)
                if (!assigned[partner] && conjugatePartner(eigenvalues[mode], eigenvalues[partner])) {
                    for (std::size_t feature = 0; feature < result.featureCount; ++feature)
                        result.coefficients[partner * result.featureCount + feature]
                            = std::conj(coefficient[feature]);
                    assigned[partner] = true; break;
                }
    }
    for (std::size_t mode = 0; mode < eigenvalues.size(); ++mode) {
        double residualSquared=0.0;
        for (std::size_t sample = 0; sample < samples; ++sample) {
            std::complex<double> residual{};
            for(std::size_t feature=0;feature<result.featureCount;++feature) {
                const auto coefficient=
                    result.coefficients[mode*result.featureCount+feature];
                residual+=(gamma[sample*result.featureCount+feature]
                    -eigenvalues[mode]*theta[sample*result.featureCount+feature])
                    *coefficient;
            }
            const double magnitude=std::abs(residual);
            residualSquared+=magnitude*magnitude;
            result.residualMaximum[mode]=std::max(
                result.residualMaximum[mode],magnitude);
        }
        result.residualRms[mode] = std::sqrt(residualSquared / samples);
    }
    // Each pair's generator residual RMS_k |Dpsi(q_k) qdot_k - lambda psi(q_k)| (Kaiser, Kutz & Brunton 2021);
    // the largest pointwise value stays with the fit.
    EigenfunctionFitDiagnostics diagnostics{"sindy-generator", samples, result.directFits, std::move(result.residualMaximum)};
    auto functions = library.dictionary
        ? FunctionExpansion(library.dictionary, std::move(result.coefficients), FunctionDomain::physicalState(stateDimension))
        : FunctionExpansion({stateDimension, exponents, center, scale}, std::move(result.coefficients),
                            FunctionDomain::physicalState(stateDimension));
    auto spectrum = Spectrum::fromContinuousRates(eigenvalues, uniform ? dt : 0.0);
    spectrum.setResiduals(result.residualRms, residual_definition::generator, residual_scope::empiricalTraining);
    return {std::move(functions), std::move(spectrum), std::move(diagnostics)};
}

EigenfunctionFitResult fitKoopmanEigenfunctionsInCoordinates(
    const std::vector<double>& dynamicalStates, std::size_t samples,
    std::size_t dynamicalStateDimension, const std::vector<double>& time,
    const std::vector<std::complex<double>>& eigenvalues,
    const Observable& coordinates, const std::vector<double>& dynamicalVelocities,
    const KoopmanSINDyOptions& options) {
    if (coordinates.stateDimension()!=dynamicalStateDimension)
        throw std::invalid_argument("SINDy coordinate observable is not state aligned");
    auto values=coordinates.evaluateRealBatch(dynamicalStates,samples);
    const auto velocities=dynamicalVelocities.empty()?std::vector<double>{}
        :coordinates.pushForwardRealBatch(dynamicalStates,dynamicalVelocities,samples);
    auto mapped=options;  // anchor states are dynamical states: map them into the fit's coordinates
    std::vector<double> anchorCoordinates;
    if (options.anchorStates) {
        anchorCoordinates=coordinates.evaluateRealBatch(*options.anchorStates,
            options.anchorStates->size()/dynamicalStateDimension);
        mapped.anchorStates=&anchorCoordinates;
    }
    auto result=fitKoopmanEigenfunctions(values,samples,coordinates.outputDimension(),time,
                                         eigenvalues,velocities,mapped);
    // The domain is y = h(x): the observable itself when it has closed forms, else its catalog names.
    auto names=coordinates.descriptor().outputNames;
    if (names.empty() && coordinates.outputDimension()==1) names.push_back(coordinates.descriptor().name);
    auto domain=!coordinates.expressions().empty()
        ? FunctionDomain::fromStateMap(Observable::fromExpressionDefinition(coordinates.expressionDefinition()))
        : FunctionDomain::fromCatalog(names, coordinates.outputDimension());
    if (!domain.stateMap) domain.physicalStateDimension=dynamicalStateDimension;
    if (!names.empty()) domain.names=names;
    result.functions=result.functions.withCoordinates(std::move(domain));
    return result;
}


}  // namespace dynsyslab::core
