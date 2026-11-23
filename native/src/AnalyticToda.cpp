#include "AnalyticSystems.hpp"
#include <cmath>
#include <stdexcept>
// Two-site periodic Toda reduces to the cosh oscillator; the open chain
// has a hyperbolic scattering solution. General periodic Toda integrability:
// Henrici/Kappeler https://doi.org/10.1088/0951-7715/21/12/001.
namespace dynsyslab::core::detail {
std::vector<AnalyticSolutionFamily> todaAnalyticSolutions(const SystemConfiguration& c) {
    if(c.identifier!="toda_lattice")return {};
    AnalyticSolutionDescriptor d{"toda_two_site","Two-site Toda solution",
        "R=q_1-q_2,\\quad C=(q_1+q_2)/2",
        "Two sites, zero damping and forcing; periodic boundary has a torus only at zero total momentum",
        "cosh/elliptic periodic relative motion or elementary open-chain scattering",
        {"https://www.mat.univie.ac.at/~gerald/ftp/book-jac/","https://arxiv.org/abs/0802.4252"}};
    return {{d,[=](const AnalyticState& x,double t0)->AnalyticTrajectory {
        if(x.size()!=4)return {{},"General Toda theta-function evaluator is not implemented; this evaluator requires two sites"};
        if(c.driver.amplitude!=0||c.parameters.at("damping")!=0)return {{},"Analytic Toda evaluator requires zero forcing and damping"};
        const double alpha=c.parameters.at("alpha"),R=x[0]-x[1],v=x[2]-x[3];
        const double C=(x[0]+x[1])/2,P=x[2]+x[3];
        const auto assemble=[=](double t,const AnalyticState& relative) {
            const double center=C+P*(t-t0)/2;
            return AnalyticState{center+relative[0]/2,center-relative[0]/2,(P+relative[1])/2,(P-relative[1])/2};
        };
        if(c.preset=="open") {
            const double w=std::hypot(v/2,alpha*std::exp(R/2));
            if(!std::isfinite(w)||w==0)return {{},"Open Toda relative energy outside representable range"};
            const double u0=std::asinh(-v/(2*alpha*std::exp(R/2)));
            const auto logcosh=[](double u){const double a=std::abs(u);return a+std::log1p(std::exp(-2*a))-std::log(2.);};
            return {[=](double t) {
                const double u=u0+w*(t-t0);
                return assemble(t,{R-2*(logcosh(u)-logcosh(u0)),-2*w*std::tanh(u)});
            },{}};
        }
        if(R==0&&v==0)return {[=](double t){return assemble(t,{0,0});},{}};
        SystemConfiguration relative;relative.identifier="cosh_oscillator";
        relative.parameters["strength"]=2*alpha*alpha;
        const auto s=createSystem(relative);
        const auto orbit=s->analyticSolutions().at(0).instantiate({R,v},t0);
        if(!orbit)return orbit;
        AnalyticTrajectory result{[=](double t){return assemble(t,orbit.evaluate(t));},{}};
        if(P==0) {
            result.period=orbit.period;result.angularFrequencies=orbit.angularFrequencies;
            result.torus=[=](const std::vector<double>& theta){return assemble(t0,orbit.torus(theta));};
        }
        return result;
    }}};
}
}
