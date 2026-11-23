#pragma once
#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <dynsyslab/numerics/NonlinearLeastSquares.hpp>
namespace dynsyslab::core {
struct PeriodicFilterOptions {
    std::vector<double> parameters; // quadrature/search grid inside the trained atlas
    std::size_t phaseCount{256}, tableSamples{4096};
    std::vector<std::size_t> measurementChannels, outputChannels;
    std::vector<double> measurementVariances, parameterPrior;
    double frozenFrequency{}; // zero uses the atlas; positive is a diagnostic baseline
};
struct PeriodicFilterResult {
    std::vector<double> outputs; // sample x requested output channels
    std::vector<double> parameterMean, parameterVariance, effectiveSampleSize, logEvidence;
    std::vector<double> finalWeights; // parameter-major, then initial phase
};
// Causal, unforced periodic-family filter. Fixed parameter/initial-phase particles
// advance with their own frequencies. No resampling/process noise. Readouts are
// interpolated on a uniform phase table. Not a controlled estimator and not KLQE.
PeriodicFilterResult filterPeriodicAtlas(const math::PeriodicAtlas&,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const PeriodicFilterOptions&);
// Refine a constant orbit parameter and initial phase using a finite past window.
// Starts are [parameter,phase]. Nonlinear local fits; caller supplies global seeds.
// Forecasting after the window requires no additional measurements.
numerics::NonlinearLeastSquaresResult identifyPeriodicOrbit(const math::PeriodicAtlas&,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const std::vector<std::vector<double>>& starts,const PeriodicFilterOptions&,
    const numerics::NonlinearLeastSquaresOptions& = {});
}
