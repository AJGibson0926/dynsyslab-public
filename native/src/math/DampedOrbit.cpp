#include <dynsyslab/math/DampedOrbit.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
void finite(double v,const char* m){if(!std::isfinite(v))throw std::invalid_argument(m);}
}
DampedAverage dampedAverage(const std::vector<double>& q,double frequency){
    if(q.size()<3)throw std::invalid_argument("Damped average needs at least three phase samples");
    if(!(frequency>0)||!std::isfinite(frequency))throw std::invalid_argument("Orbit frequency must be positive");
    std::vector<C> samples(q.begin(),q.end());const int maximum=int((q.size()-1)/2);
    const auto series=periodicFourierSeries(samples,frequency,maximum);
    DampedAverage out;out.mean=series.coefficients.at(0).real();out.chiSeries.angularFrequency=frequency;
    // Harmonic primitive of the zero-mean part; the n=0 coefficient is dropped.
    for(const auto& [n,c]:series.coefficients)if(n!=0)out.chiSeries.coefficients[n]=c/(C(0,double(n))*frequency);
    const auto chi=uniformFourierSamples(out.chiSeries,q.size());out.chi.reserve(q.size());
    for(const auto& z:chi)out.chi.push_back(z.real());
    return out;
}
DampedPhase dampedPhase(const DampedTable& t,double omega0,double realLambda,double d){
    requireIncreasing(t.energy,"Damped table energy");const auto n=t.energy.size();
    if(t.frequency.size()!=n||t.damping.size()!=n)throw std::invalid_argument("Damped table sizes differ");
    finite(omega0,"Equilibrium frequency must be finite");finite(realLambda,"Eigenvalue must be finite");
    if(!(d>0)||!std::isfinite(d))throw std::invalid_argument("Damping scale must be positive");
    const bool origin=t.damping[0]==0;
    for(std::size_t i=origin?1:0;i<n;++i)if(!(t.damping[i]>0)||!std::isfinite(t.damping[i]))
        throw std::invalid_argument("Averaged damping must be positive on the table");
    if(origin&&n<3)throw std::invalid_argument("A table starting at the equilibrium needs three nodes");
    std::vector<double> slope(n),rate(n);
    for(std::size_t i=origin?1:0;i<n;++i){slope[i]=(t.frequency[i]-omega0)/t.damping[i];rate[i]=-realLambda/(d*t.damping[i]);}
    if(origin)slope[0]=slope[1]+(slope[2]-slope[1])*(t.energy[0]-t.energy[1])/(t.energy[2]-t.energy[1]);
    DampedPhase out;out.energy=t.energy;out.phaseSlope=slope;out.phaseFunction=cumulativeTableIntegral(t.energy,slope);
    out.equilibriumFrequency=omega0;out.realEigenvalue=realLambda;out.dampingScale=d;out.amplitude.assign(n,0);
    const std::size_t first=origin?1:0;
    const std::vector<double> e(t.energy.begin()+first,t.energy.end()),g(rate.begin()+first,rate.end());
    const auto integral=cumulativeTableIntegral(e,g);
    for(std::size_t i=0;i<e.size();++i)out.amplitude[i+first]=integral[i]-integral.back();
    if(origin)out.amplitude[0]=realLambda<0?-std::numeric_limits<double>::infinity():std::numeric_limits<double>::infinity();
    return out;
}
C wkbEigenfunction(const DampedPhase& p,double energy,double phase,double chi){
    if(energy<p.energy.front()||energy>p.energy.back())throw std::invalid_argument("Energy outside the damped table");
    if(!std::isfinite(phase)||!std::isfinite(chi))throw std::invalid_argument("Phase and chi must be finite");
    const bool origin=!std::isfinite(p.amplitude.front());
    if(origin&&energy<=p.energy.front())return 0;
    const std::vector<double> e(p.energy.begin()+(origin?1:0),p.energy.end()),a(p.amplitude.begin()+(origin?1:0),p.amplitude.end());
    if(energy<e.front())throw std::invalid_argument("Energy lies in the unresolved first panel at the equilibrium");
    const double amplitude=interpolateTable(e,a,energy),S=interpolateTable(p.energy,p.phaseFunction,energy),
        slope=interpolateTable(p.energy,p.phaseSlope,energy);
    return std::exp(C(amplitude,phase+S/p.dampingScale+slope*chi));
}
LinearMode linearMode(double a,double b,double c,double d){
    for(double v:{a,b,c,d})finite(v,"Jacobian entries must be finite");
    const double half=.5*(a+d),gap=.5*(a-d),disc=gap*gap+b*c;
    const C lambda=disc<0?C(half,std::sqrt(-disc)):C(half+std::sqrt(disc));
    std::array<C,2> w{C(c),lambda-a};if(std::abs(w[0])+std::abs(w[1])==0)w={lambda-d,C(b)};
    std::array<C,2> v{C(b),lambda-a};if(std::abs(v[0])+std::abs(v[1])==0)v={lambda-d,C(c)};
    if(std::abs(w[0])+std::abs(w[1])==0){w={C(1),C(0)};v={C(1),C(0)};}
    const C pairing=w[0]*v[0]+w[1]*v[1];
    if(std::abs(pairing)<1e-14*(std::abs(w[0])+std::abs(w[1]))*(std::abs(v[0])+std::abs(v[1])))
        throw std::invalid_argument("Defective linear part has no normalized left eigenvector");
    return {lambda,{w[0]/pairing,w[1]/pairing}};
}
LaplaceValue laplaceEigenfunction(const PlanarField& drift,const PlanarPoint& star,const LinearMode& mode,
    const PlanarPoint& x0,const LaplaceOptions& o){
    if(!drift)throw std::invalid_argument("Laplace limit needs a drift");
    if(!(mode.eigenvalue.real()<0))throw std::invalid_argument("Laplace limit needs a stable eigenvalue");
    if(!(o.settleRadius>0)||!(o.maximumTime>0))throw std::invalid_argument("Laplace options must be positive");
    // Integrate the deviation so relative tolerance remains meaningful near x*.
    FlowOptions flow=o.flow;flow.absoluteTolerance=std::min(flow.absoluteTolerance,o.settleRadius*1e-9);
    const FlowField field=[&](const std::vector<double>& u,std::vector<double>& du){
        const auto f=drift({star[0]+u[0],star[1]+u[1]});du[0]=f[0];du[1]=f[1];};
    std::vector<double> u{x0[0]-star[0],x0[1]-star[1]};double t=0,h=flow.initialStep>0?flow.initialStep:1e-2;
    for(std::size_t step=0;std::hypot(u[0],u[1])>o.settleRadius;++step){
        if(step>=flow.maximumSteps||t>=o.maximumTime)throw std::runtime_error("Trajectory did not reach the equilibrium");
        t+=adaptiveFlowStep(field,u,h,flow,o.maximumTime-t);
    }
    const C projected=mode.left[0]*u[0]+mode.left[1]*u[1];
    return {std::exp(-mode.eigenvalue*t)*projected,t,std::hypot(u[0],u[1])};
}
C annularEigenfunction(C psi,C lambda,int n,double w){
    if(psi==C(0))throw std::invalid_argument("Annular eigenfunction is undefined at psi=0");
    if(lambda.real()==0||!std::isfinite(w))throw std::invalid_argument("Annular eigenfunction needs Re(lambda)!=0 and finite w");
    const C s=(C(0,w)-double(n)*lambda)/lambda.real();
    return std::exp(double(n)*std::log(psi)+s*std::log(std::abs(psi)));
}
double shellWeight(double d,double D,int n,double slope){
    if(!(d>0)||!(D>0)||n==0||slope==0||!std::isfinite(slope))throw std::invalid_argument("Shell weight needs d,D>0 and a nonzero twist");
    return std::sqrt(2*std::acos(-1.)*d*D/std::abs(n*slope));
}
}
