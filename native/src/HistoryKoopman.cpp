#include <dynsyslab/core/HistoryKoopman.hpp>
#include <dynsyslab/core/HankelOperator.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace dynsyslab::core {
namespace {

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ComplexMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
    Eigen::Dynamic, Eigen::RowMajor>;
using Vector = Eigen::VectorXd;


double power(double value, std::size_t exponent) {
    double result = 1.0;
    for (std::size_t count = 0; count < exponent; ++count) result *= value;
    return result;
}

std::vector<double> featureRow(const double* coordinates,
                               const PolynomialDictionary& basis) {
    const auto dimension = basis.coordinateDimension, features = basis.featureCount();
    std::vector<double> result(features, 1.0);
    for (std::size_t feature = 0; feature < features; ++feature)
        for (std::size_t axis = 0; axis < dimension; ++axis)
            result[feature] *= power(
                (coordinates[axis] - basis.center[axis]) / basis.scale[axis],
                basis.exponents[feature * dimension + axis]);
    return result;
}

std::vector<std::complex<double>> evaluateRun(
    const PolynomialDictionary& contract,
    const std::vector<std::complex<double>>& coefficient,
    const std::vector<double>& run) {
    const std::size_t samples = run.size() / contract.coordinateDimension;
    std::vector<std::complex<double>> result(samples);
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const auto features = featureRow(
            run.data() + sample * contract.coordinateDimension, contract);
        for (std::size_t feature = 0; feature < features.size(); ++feature)
            result[sample] += features[feature] * coefficient[feature];
    }
    return result;
}

double phaseRSquared(const std::vector<std::complex<double>>& values,
                     double sampleTime) {
    if (values.size() < 2) return -std::numeric_limits<double>::infinity();
    std::vector<double> phase(values.size());
    phase[0] = std::arg(values[0]);
    constexpr double pi = 3.14159265358979323846;
    for (std::size_t index = 1; index < values.size(); ++index) {
        double delta = std::arg(values[index]) - std::arg(values[index - 1]);
        while (delta > pi) delta -= 2.0 * pi;
        while (delta < -pi) delta += 2.0 * pi;
        phase[index] = phase[index - 1] + delta;
    }
    const double meanTime = sampleTime * static_cast<double>(values.size() - 1) / 2.0;
    const double meanPhase = std::accumulate(phase.begin(), phase.end(), 0.0)
        / static_cast<double>(phase.size());
    double timeSquared = 0.0, covariance = 0.0;
    for (std::size_t index = 0; index < phase.size(); ++index) {
        const double centeredTime = sampleTime * static_cast<double>(index) - meanTime;
        timeSquared += centeredTime * centeredTime;
        covariance += centeredTime * (phase[index] - meanPhase);
    }
    if (timeSquared <= 0.0) return -std::numeric_limits<double>::infinity();
    const double slope = covariance / timeSquared;
    double residual = 0.0, total = 0.0;
    for (std::size_t index = 0; index < phase.size(); ++index) {
        const double fitted = meanPhase + slope
            * (sampleTime * static_cast<double>(index) - meanTime);
        const double difference = phase[index] - fitted;
        residual += difference * difference;
        const double centered = phase[index] - meanPhase;
        total += centered * centered;
    }
    return 1.0 - residual / std::max(total, std::numeric_limits<double>::min());
}

}  // namespace

HankelLeftBasis fitPooledHankelBasis(
    const std::vector<std::vector<double>>& observations,
    const PooledHankelOptions& options) {
    if (observations.empty() || options.delays < 2 || !options.rank
        || !options.subspaceIterations)
        throw std::invalid_argument("invalid pooled Hankel options");
    std::size_t observationCount = 0;
    double squared = 0.0;
    for (const auto& signal : observations) {
        if (signal.size() < options.delays)
            throw std::invalid_argument("pooled Hankel trajectory is shorter than its delay window");
        for (const double value : signal) {
            if (!std::isfinite(value))
                throw std::invalid_argument("pooled Hankel observations must be finite");
            squared += value * value;
        }
        observationCount += signal.size();
    }
    const double scale = options.normalizeByRootMeanSquare
        ? std::sqrt(squared / static_cast<double>(observationCount)) : 1.0;
    if (!(scale > std::numeric_limits<double>::epsilon()))
        throw std::invalid_argument("pooled Hankel observations have zero scale");
    // The shared Hankel module: one scalar channel scaled by the pooled RMS.
    std::vector<HankelSegment> segments;
    segments.reserve(observations.size());
    for (const auto& source : observations) segments.push_back({source, source.size(), 1});
    const HankelOperator hankel(segments, {{0, 0, scale}}, options.delays);
    const auto fitted = fitHankelLeftBasis(hankel, {options.rank, options.oversampling,
        options.subspaceIterations, options.seed});
    HankelLeftBasis result;
    result.rows = options.delays;
    result.rank = options.rank;
    result.trajectories = observations.size();
    result.columns = hankel.columns();
    result.observationScale = scale;
    result.totalEnergy = fitted.totalEnergy;
    result.leftSingularVectors = fitted.leftSingularVectors;
    result.singularValues = fitted.singularValues;
    result.cumulativeEnergy = fitted.cumulativeEnergy;
    return result;
}


PooledPolynomialEDMD fitPooledPolynomialEDMD(
    const std::vector<std::vector<double>>& coordinateRuns,
    std::size_t coordinateDimension,
    const PooledPolynomialEDMDOptions& options) {
    if (coordinateRuns.empty() || !coordinateDimension || !options.degree
        || !options.pairStride || !(options.sampleTime > 0.0)
        || options.relativeRidge < 0.0 || !std::isfinite(options.relativeRidge))
        throw std::invalid_argument("invalid pooled polynomial EDMD options");
    std::size_t samples = 0, pairs = 0;
    for (const auto& run : coordinateRuns) {
        if (run.size() % coordinateDimension)
            throw std::invalid_argument("pooled EDMD run shape is inconsistent");
        const std::size_t count = run.size() / coordinateDimension;
        if (count <= options.pairStride)
            throw std::invalid_argument("pooled EDMD run is shorter than its pair stride");
        if (!std::all_of(run.begin(), run.end(), [](double value) { return std::isfinite(value); }))
            throw std::invalid_argument("pooled EDMD coordinates must be finite");
        samples += count;
        pairs += count - options.pairStride;
    }
    PolynomialDictionary contract;
    contract.coordinateDimension = coordinateDimension;
    contract.center.assign(coordinateDimension, 0.0);
    contract.scale.assign(coordinateDimension, 0.0);
    for (const auto& run : coordinateRuns)
        for (std::size_t index = 0; index < run.size(); ++index)
            contract.center[index % coordinateDimension] += run[index];
    for (double& value : contract.center) value /= static_cast<double>(samples);
    for (const auto& run : coordinateRuns)
        for (std::size_t index = 0; index < run.size(); ++index) {
            const std::size_t axis = index % coordinateDimension;
            const double difference = run[index] - contract.center[axis];
            contract.scale[axis] += difference * difference;
        }
    for (double& value : contract.scale) {
        value = std::sqrt(value / static_cast<double>(samples));
        if (!(value > std::numeric_limits<double>::epsilon()))
            throw std::invalid_argument("pooled EDMD coordinate normalization is singular");
    }
    contract.exponents = polynomialExponentTable(
        coordinateDimension, options.degree, options.includeConstant);
    const std::size_t featureCount = contract.featureCount();

    std::vector<std::vector<double>> rawFeatures;
    rawFeatures.reserve(coordinateRuns.size());
    std::vector<double> featureScaleSquared(featureCount, 0.0);
    for (const auto& run : coordinateRuns) {
        const std::size_t count = run.size() / coordinateDimension;
        std::vector<double> features(count * featureCount);
        for (std::size_t sample = 0; sample < count; ++sample) {
            const auto row = featureRow(
                run.data() + sample * coordinateDimension, contract);
            for (std::size_t feature = 0; feature < featureCount; ++feature) {
                features[sample * featureCount + feature] = row[feature];
                featureScaleSquared[feature] += row[feature] * row[feature];
            }
        }
        rawFeatures.push_back(std::move(features));
    }
    std::vector<double> featureScale(featureCount);
    for (std::size_t feature = 0; feature < featureCount; ++feature) {
        featureScale[feature] = std::sqrt(
            featureScaleSquared[feature] / static_cast<double>(samples));
        if (!(featureScale[feature] > std::numeric_limits<double>::epsilon()))
            featureScale[feature] = 1.0;
    }
    Matrix minus(static_cast<Eigen::Index>(pairs),
                 static_cast<Eigen::Index>(featureCount));
    Matrix plus = minus;
    std::size_t outputRow = 0;
    for (std::size_t runIndex = 0; runIndex < coordinateRuns.size(); ++runIndex) {
        const std::size_t count = coordinateRuns[runIndex].size() / coordinateDimension;
        const auto& features = rawFeatures[runIndex];
        for (std::size_t sample = 0; sample + options.pairStride < count; ++sample) {
            for (std::size_t feature = 0; feature < featureCount; ++feature) {
                minus(static_cast<Eigen::Index>(outputRow), static_cast<Eigen::Index>(feature))
                    = features[sample * featureCount + feature] / featureScale[feature];
                plus(static_cast<Eigen::Index>(outputRow), static_cast<Eigen::Index>(feature))
                    = features[(sample + options.pairStride) * featureCount + feature]
                    / featureScale[feature];
            }
            ++outputRow;
        }
    }
    const Matrix gram = minus.transpose() * minus / static_cast<double>(pairs);
    const Matrix cross = minus.transpose() * plus / static_cast<double>(pairs);
    const double ridge = options.relativeRidge * gram.trace()
        / static_cast<double>(featureCount);
    Matrix regularized = gram;
    regularized.diagonal().array() += ridge;
    const Matrix operation = regularized.completeOrthogonalDecomposition().solve(cross);
    Eigen::ComplexEigenSolver<Matrix> eigensolver(operation);
    if (eigensolver.info() != Eigen::Success)
        throw std::runtime_error("pooled polynomial EDMD eigensystem did not converge");
    const ComplexMatrix eigenvectors = eigensolver.eigenvectors();

    PooledPolynomialEDMD result;
    result.trajectories = coordinateRuns.size();
    result.pairs = pairs;
    result.pairStride = options.pairStride;
    result.pairTime = options.sampleTime * static_cast<double>(options.pairStride);
    result.ridge = ridge;
    result.featureScale = featureScale;
    result.operatorMatrix.assign(operation.data(), operation.data() + operation.size());
    result.discreteEigenvalues.assign(eigensolver.eigenvalues().data(),
        eigensolver.eigenvalues().data() + eigensolver.eigenvalues().size());
    result.samples = samples;
    result.coefficients.resize(featureCount * featureCount);
    for (std::size_t mode = 0; mode < featureCount; ++mode)
        for (std::size_t feature = 0; feature < featureCount; ++feature)
            result.coefficients[mode * featureCount + feature]
                = eigenvectors(static_cast<Eigen::Index>(feature),
                    static_cast<Eigen::Index>(mode)) / featureScale[feature];
    result.basis = std::move(contract);
    return result;
}

std::vector<EDMDModeCandidate> selectPolynomialEDMDModes(
    const PooledPolynomialEDMD& model,
    const std::vector<std::vector<double>>& coordinateRuns,
    const EDMDModeSearchOptions& options) {
    const auto& contract = model.basis;
    const std::size_t featureCount = contract.featureCount();
    if (coordinateRuns.empty() || model.discreteEigenvalues.size() * featureCount
        != model.coefficients.size())
        throw std::invalid_argument("invalid pooled EDMD mode search inputs");
    std::vector<EDMDModeCandidate> result;
    for (std::size_t mode = 0; mode < model.discreteEigenvalues.size(); ++mode) {
        const auto discrete = model.discreteEigenvalues[mode];
        const auto continuous = std::log(discrete) / model.pairTime;
        const double magnitude = std::abs(discrete);
        if (magnitude < options.minimumDiscreteMagnitude
            || magnitude > options.maximumDiscreteMagnitude
            || continuous.real() < options.minimumContinuousReal
            || continuous.real() > options.maximumContinuousReal
            || continuous.imag() < options.minimumContinuousImaginary
            || continuous.imag() > options.maximumContinuousImaginary
            || (options.positiveImaginaryOnly && discrete.imag() <= 0.0))
            continue;
        std::vector<std::complex<double>> coefficient(
            model.coefficients.begin() + std::ptrdiff_t(mode * featureCount),
            model.coefficients.begin() + std::ptrdiff_t((mode + 1) * featureCount));
        std::vector<std::vector<std::complex<double>>> valuesByRun;
        double normSquared = 0.0;
        std::size_t valueCount = 0;
        for (const auto& run : coordinateRuns) {
            valuesByRun.push_back(evaluateRun(contract, coefficient, run));
            for (const auto value : valuesByRun.back()) normSquared += std::norm(value);
            valueCount += valuesByRun.back().size();
        }
        const double rms = std::sqrt(normSquared / static_cast<double>(valueCount));
        if (!(rms > std::numeric_limits<double>::min())) continue;
        for (auto& value : coefficient) value /= rms;
        for (auto& run : valuesByRun) for (auto& value : run) value /= rms;
        const auto& reference = valuesByRun.front();
        const auto pivot = *std::max_element(reference.begin(), reference.end(),
            [](auto left, auto right) { return std::abs(left) < std::abs(right); });
        if (std::abs(pivot) > std::numeric_limits<double>::min()) {
            const auto rotation = std::exp(std::complex<double>(0.0, -std::arg(pivot)));
            for (auto& value : coefficient) value *= rotation;
            for (auto& run : valuesByRun) for (auto& value : run) value *= rotation;
        }
        double residualSum = 0.0;
        double minimumPhase = std::numeric_limits<double>::infinity();
        for (const auto& values : valuesByRun) {
            double residualSquared = 0.0, targetSquared = 0.0;
            for (std::size_t index = 0; index + model.pairStride < values.size(); ++index) {
                const auto residual = values[index + model.pairStride]
                    - discrete * values[index];
                residualSquared += std::norm(residual);
                targetSquared += std::norm(values[index + model.pairStride]);
            }
            residualSum += std::sqrt(residualSquared
                / std::max(targetSquared, std::numeric_limits<double>::min()));
            minimumPhase = std::min(minimumPhase,
                phaseRSquared(values, model.pairTime
                    / static_cast<double>(model.pairStride)));
        }
        result.push_back({mode, discrete, continuous, std::move(coefficient),
            residualSum / static_cast<double>(coordinateRuns.size()), minimumPhase});
    }
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        if (left.meanPairResidual != right.meanPairResidual)
            return left.meanPairResidual < right.meanPairResidual;
        return left.minimumPhaseRSquared > right.minimumPhaseRSquared;
    });
    return result;
}

EigenfunctionFitResult selectedPolynomialEigenfunction(
    const PooledPolynomialEDMD& model,
    const EDMDModeCandidate& candidate) {
    const auto& basis = model.basis;
    if (candidate.coefficients.size() != basis.featureCount())
        throw std::invalid_argument("selected EDMD coefficient shape is inconsistent");
    // The independent coordinates xi are produced by a stateful history projector, not by an
    // instantaneous Observable: the domain is supplied values, tagged with that provider.
    auto domain = FunctionDomain::supplied(basis.coordinateDimension);
    for (std::size_t axis = 0; axis < basis.coordinateDimension; ++axis) domain.names.push_back("xi" + std::to_string(axis + 1));
    domain.providerIdentity = "causal-history";
    // The pair's residual is the one-step residual normalized by the successor and averaged over runs
    // (selectPolynomialEDMDModes), a definition of its own; it is not converted to ResDMD's.
    EigenfunctionFitDiagnostics diagnostics{"pooled-polynomial-edmd", model.samples, 1, {}};
    auto spectrum = Spectrum::fromStoredRates({candidate.continuousEigenvalue}, model.pairTime, EigenvalueDomain::discreteTime);
    spectrum.setResidual(0, EigenpairResidual{candidate.meanPairResidual, residual_definition::oneStepSuccessorMean, residual_scope::empiricalTraining});
    return {FunctionExpansion(basis, candidate.coefficients, std::move(domain)), std::move(spectrum), std::move(diagnostics)};
}

}  // namespace dynsyslab::core
