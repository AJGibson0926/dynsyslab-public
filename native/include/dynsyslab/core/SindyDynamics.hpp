#pragma once
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/core/Integrator.hpp>

namespace dynsyslab::core {
// Original SINDy: qdot = Theta(q) C. Separate from Koopman eigenfunction fitting.
// Sequentially thresholded least squares (SINDy: Brunton, Proctor & Kutz, PNAS 113:3932, 2016) on
// rows that are either the instantaneous rates or, in weak (integral) form, state increments over
// physical-time windows, q(t + D) - q(t) = int_t^{t+D} Theta(q) dt xi (integral SINDy: Schaeffer &
// McCalla, Phys. Rev. E 96:023302, 2017; weak SINDy: Messenger & Bortz, Multiscale Model. Simul.
// 19:1474, 2021), with Theta piecewise linear between samples and the window end interpolated.
struct SindyDynamicsOptions {
    double threshold{0.001}; // applied after feature and response RMS normalization
    double singularValueTolerance{1e-10};
    std::size_t maximumIterations{20};
    // Tikhonov ridge on the normalized problem, min ||A xi - b||^2 + ridge ||xi||^2 (ridge
    // regression: Hoerl & Kennard, Technometrics 12:55, 1970). Zero: plain least squares.
    double ridge{0.0};
    // Per-sample weights on the instantaneous rows (empty: equal); a row's squared residual is
    // multiplied by its weight.
    std::vector<double> sampleWeights;
    // Weak form: physical window lengths D (empty: instantaneous rows), windows starting every
    // `weakStride` samples. Only for a vector field (the response is the predictors' own rate).
    std::vector<double> weakDurations;
    std::size_t weakStride{1};
};
// Nodal weights of the trapezoid rule on an increasing clock (the integral of the piecewise-linear
// interpolant): a sample's share of the time it represents.
std::vector<double> trapezoidWeights(const TimeSeries& series);
struct SindyDynamicsFitDiagnostics {
    // Per output: least-squares error and support size BEFORE each pruning step.
    std::vector<std::vector<double>> relativeResidual;
    std::vector<std::vector<std::size_t>> supportSize;
    json::Value toJson() const;
};
class SindyDynamicsModel {
public:
    // outputs == 0 means one output per input coordinate (a vector field z'=f(z)).
    // A different count fits another coordinate set's rates, w'=g(z).
    SindyDynamicsModel(Observable library,std::vector<double> coefficients,std::size_t outputs=0);
    const Observable& library() const noexcept {return library_;}
    const std::vector<double>& coefficients() const noexcept {return coefficients_;}
    std::size_t dimension() const noexcept {return library_.stateDimension();}
    std::size_t outputs() const noexcept {return outputs_;}
    bool isVectorField() const noexcept {return outputs_==dimension();}
    std::vector<double> evaluate(const double* q) const;
    TimeSeries evaluate(const TimeSeries& q) const;
    Observable vectorField() const;
    json::Value definition() const;
    static SindyDynamicsModel fromDefinition(const json::Value& definition);
private:
    Observable library_;
    std::vector<double> coefficients_; // feature x output, row-major, original library units
    std::size_t outputs_{};
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
