#pragma once
#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <dynsyslab/numerics/NonlinearLeastSquares.hpp>
namespace dynsyslab::core {
// A real linear transform of real past measurements. Complex features are
// represented by adjacent real/imaginary columns. No information is added.
struct PeriodicHistoryFeatureMap {
    std::vector<std::size_t> channels;
    std::vector<double> timeOffsets; // nonpositive, at least one zero
    std::size_t outputs{};
    std::vector<double> coefficients; // (offset/channel) x output
    std::vector<double> noiseStandardDeviations; // independent input noise
    bool whiten{true};
    double relativeTolerance{1e-10};
};
struct PeriodicFeatureIdentificationOptions {
    std::size_t tableSamples{8192},phaseSeeds{256},starts{8},parameterSubdivisions{1};
    numerics::NonlinearLeastSquaresOptions optimization;
};
struct PeriodicFeatureIdentification {
    numerics::NonlinearLeastSquaresResult fit;
    double frequency{},historyRelativeError{};
    std::vector<double> current,forecast; // forecast time/channel
    std::vector<double> jacobianSingularValues; // parameter/unit-range, phase/2pi
};
struct PeriodicFeatureIdentificationBatch {
    std::size_t inputs{},outputs{},rank{};
    std::vector<double> singularValues,effectiveCoefficients;
    std::vector<PeriodicFeatureIdentification> estimates;
};
// Fixed atlas and features; coarse seeds use training knots and uniform phases,
// followed by local bounded least squares. No true parameter/phase is accepted.
// Whitening uses the full possibly singular transformed noise covariance, via
// SVD, not independent errors assigned to correlated output features.
PeriodicFeatureIdentificationBatch identifyPeriodicHistoryFeatures(
    const math::PeriodicAtlas&,const PeriodicHistoryFeatureMap&,
    const std::vector<double>& histories,const std::vector<double>& forecastTimes={},
    const PeriodicFeatureIdentificationOptions& = {});
}
