#pragma once
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/math/PeriodicChart.hpp>
#include <functional>
namespace dynsyslab::core {
using PlanarInputField=std::function<std::array<double,2>(std::array<double,2>)>;
struct PeriodicTransportOptions {
    std::array<std::size_t,2> stateChannels{0,1};
    std::size_t phaseSamples{8192};
    int inputMaximumHarmonic{2047};
    double minimumNormalizedDeterminant{1e-10};
    double frozenFrequency{}; // positive replaces only the unforced phase rate
};
struct PeriodicTransportDiagnostics {
    double minimumAbsoluteNormalizedDeterminant{1};
    double maximumInputProjectionError{}; // sampled physical B reconstruction
    double inputProjectionRelativeL2{};
    std::array<double,2> worstInputLocation{}; // parameter, phase
    double minimumNormalizedDeterminant{1},maximumNormalizedDeterminant{-1};
    double chartRealityProjectionL1{},inputRealityProjectionL1{};
};
// Precompute b = (DX)^(-1)B, then advance only [orbit parameter,phase].
// The physical field callback is discarded after construction. No physical
// system RHS is used in the rollout. This is nonlinear, not a finite LTI KROM.
class PeriodicTransportSystem final:public ODESystem {
  public:
    PeriodicTransportSystem(const math::PeriodicAtlas&,PlanarInputField,
                            PeriodicTransportOptions = {});
    const SystemDescriptor& descriptor()const override {return descriptor_;}
    const SystemConfiguration& configuration()const override {return configuration_;}
    std::vector<double> initialState()const override;
    void rhs(double,const double*,double*)const override;
    void rhsWithInput(double,const double*,double,double*)const override;
    std::vector<double> jacobian(double,const double*)const override;
    std::vector<double> inputMatrix(const double*)const override;
    bool invalid(double,const double*)const override;
    math::PeriodicChartValue readout(const double* state)const;
    std::array<double,2> physicalState(const double* state)const;
    const PeriodicTransportDiagnostics& diagnostics()const {return diagnostics_;}
  private:
    math::PeriodicChart chart_;
    std::unique_ptr<math::PeriodicChart> input_;
    PeriodicTransportOptions options_;
    PeriodicTransportDiagnostics diagnostics_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
};
}
