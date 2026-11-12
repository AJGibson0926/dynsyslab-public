#pragma once
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/core/Integrator.hpp>

namespace dynsyslab::core {
// Original SINDy: qdot = Theta(q) C. Separate from Koopman eigenfunction fitting.
struct SindyDynamicsOptions {
    double threshold{0.001}; // applied after feature and response RMS normalization
    double singularValueTolerance{1e-10};
    std::size_t maximumIterations{20};
};
struct SindyDynamicsFitDiagnostics {
    // Per output: least-squares error and support size BEFORE each pruning step.
    std::vector<std::vector<double>> relativeResidual;
    std::vector<std::vector<std::size_t>> supportSize;
    json::Value toJson() const;
};
class SindyDynamicsModel {
public:
    SindyDynamicsModel(Observable library,std::vector<double> coefficients);
    const Observable& library() const noexcept {return library_;}
    const std::vector<double>& coefficients() const noexcept {return coefficients_;}
    std::size_t dimension() const noexcept {return library_.stateDimension();}
    std::vector<double> evaluate(const double* q) const;
    TimeSeries evaluate(const TimeSeries& q) const;
    Observable vectorField() const;
    json::Value definition() const;
    static SindyDynamicsModel fromDefinition(const json::Value& definition);
private:
    Observable library_;
    std::vector<double> coefficients_; // feature x output, row-major, original library units
};
// Normalization is encoded in the returned expression library, making the fitted
// model self-contained. Empty center/scale mean raw coordinates.
Observable polynomialDynamicsLibrary(std::size_t dimension,std::size_t degree=4,
    const std::vector<double>& center={},const std::vector<double>& scale={});
SindyDynamicsModel fitSindyDynamics(const TimeSeries& coordinates,const TimeSeries& derivatives,
    const Observable& library,const SindyDynamicsOptions& options={},
    SindyDynamicsFitDiagnostics* diagnostics=nullptr);
struct DynamicsError {
    std::vector<double> referenceRms,rmse,relativeRms,maximumAbsolute;
    json::Value toJson() const;
};
DynamicsError dynamicsError(const TimeSeries& expected,const TimeSeries& predicted);
NativeSolution integrateSindyDynamics(const SindyDynamicsModel& model,
    const std::vector<double>& initialState,const IntegrationOptions& options);
}
