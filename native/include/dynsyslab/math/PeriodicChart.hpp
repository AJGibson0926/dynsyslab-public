#pragma once
#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <array>
namespace dynsyslab::math {
struct PeriodicChartValue {
    std::vector<double> value, parameterDerivative, phaseDerivative;
};
// Real periodic chart with tensor cubic Hermite interpolation. Phase slopes come
// from the supplied Fourier series; parameter slopes use three-point differences
// (one-sided at endpoints, secants for two knots). Derivatives are of this interpolant, not of the unknown
// exact chart. No extrapolation. No global injectivity certificate.
class PeriodicChart {
  public:
    PeriodicChart(const PeriodicAtlas&,std::size_t phaseSamples=8192);
    PeriodicChartValue evaluate(double parameter,double phase) const;
    std::pair<double,double> frequency(double parameter) const; // value, derivative
    const std::vector<double>& parameters() const { return parameters_; }
    std::size_t channels() const { return channels_; }
    double maximumRealityProjectionL1() const { return realityProjection_; }
  private:
    struct Entry { double value{},phase{},parameter{},mixed{}; };
    std::vector<double> parameters_,frequencies_,frequencySlopes_;
    std::vector<Entry> table_;
    std::size_t samples_,channels_;
    double realityProjection_{};
    std::pair<std::size_t,double> bracket(double) const;
};
struct PlanarPullback {
    std::array<double,2> coordinates;
    double determinant{}, normalizedDeterminant{};
};
// Solve [X_parameter X_phase] b = B in the selected two real state channels.
// Reject singular/nearly parallel chart columns. Column normalization removes
// parameter units, but the diagnostic still depends on physical state scaling.
PlanarPullback pullbackPlanarField(const PeriodicChartValue&,
    std::array<std::size_t,2> stateChannels,std::array<double,2> field,
    double minimumNormalizedDeterminant=1e-10);
}
