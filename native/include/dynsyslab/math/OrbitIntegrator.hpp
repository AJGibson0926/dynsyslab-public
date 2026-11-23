#pragma once
#include <cstddef>
#include <functional>
#include <vector>
namespace dynsyslab::math {
// Autonomous first-order field. The output vector is presized by the caller.
using FlowField=std::function<void(const std::vector<double>&,std::vector<double>&)>;
struct FlowOptions {
    double relativeTolerance{1e-11},absoluteTolerance{1e-13},initialStep{1e-2},maximumStep{0};
    std::size_t maximumSteps{2000000};
};
// One Gragg-Bulirsch-Stoer step: modified midpoint with step counts 2,4,...,16 and
// Aitken-Neville extrapolation in h^2. Gragg, SIAM J. Numer. Anal. 2 (1965);
// Bulirsch & Stoer, Numer. Math. 8 (1966); Hairer, Norsett & Wanner, Solving
// Ordinary Differential Equations I, section II.9. The core integrator is not reused
// because it binds ODESystem and output grids; event refinement needs local steps.
void extrapolatedStep(const FlowField&,const std::vector<double>& x,double h,
    std::vector<double>& out,std::vector<double>* lowerOrder=nullptr);
// Advances x by one accepted step of at most limit; updates the suggested h.
double adaptiveFlowStep(const FlowField&,std::vector<double>& x,double& h,
    const FlowOptions&,double limit);
void integrateFlow(const FlowField&,std::vector<double>& x,double duration,const FlowOptions&);
// Tables have strictly increasing abscissae. Local four-point Lagrange interpolation;
// integrals use two-point Gauss per panel, exact for the local cubic.
double interpolateTable(const std::vector<double>& x,const std::vector<double>& f,double at);
std::vector<double> cumulativeTableIntegral(const std::vector<double>& x,const std::vector<double>& f);
std::vector<double> tableDerivative(const std::vector<double>& x,const std::vector<double>& f);
void requireIncreasing(const std::vector<double>& x,const char* name);
}
