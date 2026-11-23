#include <dynsyslab/math/MelnikovTools.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;using V=std::vector<double>;
struct Leg {C core,tail;double time{},radius{},decay{};};
void validate(const ConnectionProblem& p){
    if(!p.drift)throw std::invalid_argument("Connection needs a drift");
    if(!(p.forwardRate>0)||!(p.backwardRate>0)||!std::isfinite(p.forwardRate)||!std::isfinite(p.backwardRate))
        throw std::invalid_argument("Connection rates must be positive");
    if(!(p.tailRadius>0)||!(p.maximumTime>0))throw std::invalid_argument("Connection tail radius and time must be positive");
}
// Integrates x, signed time and the complex quadrature along one direction.
Leg integrateLeg(const ConnectionProblem& p,const ConnectionIntegrand& g,const PlanarPoint& saddle,
    double rate,double sign,double oscillation){
    const FlowField field=[&](const V& x,V& dx){
        const auto f=p.drift({x[0],x[1]});const auto value=g({x[0],x[1]},x[2]);
        dx[0]=sign*f[0];dx[1]=sign*f[1];dx[2]=sign;dx[3]=value.real();dx[4]=value.imag();};
    auto flow=p.flow;const double cap=.5/rate;flow.maximumStep=flow.maximumStep>0?std::min(flow.maximumStep,cap):cap;
    V x{p.midpoint[0],p.midpoint[1],0,0,0};double s=0,h=std::min(flow.initialStep>0?flow.initialStep:1e-2,cap);
    std::vector<std::pair<double,double>> history{{0.,std::abs(g(p.midpoint,0))}};
    const auto radius=[&]{return std::hypot(x[0]-saddle[0],x[1]-saddle[1]);};
    for(std::size_t step=0;radius()>p.tailRadius;++step){
        if(step>=flow.maximumSteps||s>=p.maximumTime)throw std::runtime_error("Connection did not reach its saddle; check the midpoint");
        s+=adaptiveFlowStep(field,x,h,flow,p.maximumTime-s);
        history.emplace_back(s,std::abs(g({x[0],x[1]},x[2])));
    }
    Leg leg;leg.core={x[3],x[4]};leg.time=s;leg.radius=radius();
    const C end=g({x[0],x[1]},x[2]);const double magnitude=std::abs(end);double decay=rate;
    const double target=s-1/rate;
    auto it=std::upper_bound(history.begin(),history.end(),target,[](double a,const auto& b){return a<b.first;});
    if(it!=history.begin()){--it;
        if(it->second>0&&magnitude>0&&s>it->first)decay=std::clamp(std::log(it->second/magnitude)/(s-it->first),rate,4*rate);}
    leg.decay=decay;
    // Forward tail g/(r-iw); backward tail g/(r+iw). sign carries the conjugate direction.
    leg.tail=magnitude>0?end/C(decay,-sign*oscillation):C(0);
    return leg;
}
}
ConnectionIntegral connectionIntegral(const ConnectionProblem& p,const ConnectionIntegrand& g,double oscillation){
    validate(p);if(!g)throw std::invalid_argument("Connection needs an integrand");
    if(!std::isfinite(oscillation))throw std::invalid_argument("Oscillation frequency must be finite");
    const auto f=integrateLeg(p,g,p.forwardSaddle,p.forwardRate,1,oscillation);
    const auto b=integrateLeg(p,g,p.backwardSaddle,p.backwardRate,-1,oscillation);
    ConnectionIntegral r;r.core=f.core+b.core;r.forwardTail=f.tail;r.backwardTail=b.tail;
    r.value=r.core+f.tail+b.tail;r.forwardTime=f.time;r.backwardTime=b.time;r.forwardRadius=f.radius;
    r.backwardRadius=b.radius;r.forwardDecay=f.decay;r.backwardDecay=b.decay;return r;
}
ConnectionIntegral melnikovWedge(const ConnectionProblem& p,const PlanarField& Y){
    if(!Y)throw std::invalid_argument("Melnikov wedge needs a perturbation");
    return connectionIntegral(p,[&](const PlanarPoint& x,double){const auto f=p.drift(x);const auto y=Y(x);return C(f[0]*y[1]-f[1]*y[0]);});
}
ConnectionIntegral dampingIntegral(const ConnectionProblem& p,const PlanarScalar& Q){
    if(!Q)throw std::invalid_argument("Damping integral needs a dissipation");
    return connectionIntegral(p,[&](const PlanarPoint& x,double){return C(Q(x));});
}
ConnectionIntegral forcingIntegral(const ConnectionProblem& p,const PlanarField& Y,double w){
    if(!Y)throw std::invalid_argument("Forcing integral needs a forcing direction");
    if(!std::isfinite(w))throw std::invalid_argument("Forcing frequency must be finite");
    return connectionIntegral(p,[&](const PlanarPoint& x,double t){const auto f=p.drift(x);const auto y=Y(x);
        return C(f[0]*y[1]-f[1]*y[0])*std::exp(C(0,w*t));},w);
}
double melnikovThreshold(double D,double A){
    if(!(A>0)||!std::isfinite(D)||!std::isfinite(A))throw std::invalid_argument("Threshold needs a positive forcing amplitude");
    return std::abs(D)/A;
}
LobeFlux lobeFlux(double d,double D,double F,double A,double w){
    if(!(d>=0)||!(F>0)||!(A>0)||!(w>0)||!std::isfinite(D))throw std::invalid_argument("Lobe flux needs d>=0 and F,A,w>0");
    LobeFlux r;r.ratio=d*std::abs(D)/(F*A);r.transverse=r.ratio<1;
    if(!r.transverse)return r;
    const double angle=std::acos(r.ratio);
    r.area=2*(F*A*std::sin(angle)-d*std::abs(D)*angle)/w;
    r.fluxPerTime=r.area*w/(2*std::acos(-1.));return r;
}
std::vector<double> captureProbabilities(const std::vector<double>& M){
    if(M.empty())throw std::invalid_argument("Capture needs at least one lobe");
    double total=0;for(double v:M){if(!(v>=0)||!std::isfinite(v))throw std::invalid_argument("Lobe integrals must be finite and nonnegative");total+=v;}
    if(!(total>0))throw std::invalid_argument("Capture needs a positive total lobe integral");
    std::vector<double> out;out.reserve(M.size());for(double v:M)out.push_back(v/total);return out;
}
std::vector<ShellMelnikovRow> shellMelnikovComparison(const std::vector<int>& n,const V& e,const V& T,const V& y,double A){
    if(n.size()!=e.size()||n.size()!=T.size()||n.size()!=y.size())throw std::invalid_argument("Shell comparison sizes differ");
    if(!(A>0)||!std::isfinite(A))throw std::invalid_argument("Shell comparison needs a positive amplitude");
    std::vector<ShellMelnikovRow> out;out.reserve(n.size());
    for(std::size_t i=0;i<n.size();++i){
        if(!(T[i]>0)||!(y[i]>=0)||!std::isfinite(T[i])||!std::isfinite(y[i]))throw std::invalid_argument("Shell periods and magnitudes must be valid");
        const double product=T[i]*y[i];out.push_back({n[i],e[i],T[i],y[i],product,product/A});
    }
    return out;
}
}
