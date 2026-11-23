#include "AnalyticSystems.hpp"
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Euler solution: Jacobi (1849); modern explicit review:
// https://doi.org/10.3390/dynamics3030032 (The Free Euler Rigid Body Revisited).
// Kepler: physical mean anomaly, Kepler inverse rather than eccentric-anomaly time.
// Fourier-Bessel inversion: https://www.math.nthu.edu.tw/~kchen/teaching/CelMech2024-4.pdf
namespace dynsyslab::core::detail {
namespace {
constexpr double pi=3.14159265358979323846;
double p(const SystemConfiguration&c,const char*k){return c.parameters.at(k);}
AnalyticTrajectory no(const char*r){return {{},r};}
AnalyticTrajectory periodic(std::function<AnalyticState(double)> f,double w,double t0){
    AnalyticTrajectory r{[=](double t){return f(w*(t-t0));},{},2*pi/w};r.angularFrequencies={w};
    r.torus=[=](const std::vector<double>&x){if(x.size()!=1||!std::isfinite(x[0]))throw std::invalid_argument("Expected one finite phase");return f(x[0]);};return r;
}
double keplerInverse(double M,double e){
    M=std::remainder(M,2*pi);double lo=-pi,hi=pi,x=M;
    for(int j=0;j<80;++j){const double f=x-e*std::sin(x)-M;
        if(std::abs(f)<2e-15)return x;if(f>0)hi=x;else lo=x;
        const double candidate=x-f/(1-e*std::cos(x));x=candidate>lo&&candidate<hi?candidate:(lo+hi)/2;}
    return x;
}
}
std::vector<AnalyticSolutionFamily> celestialAnalyticSolutions(const SystemConfiguration&c){
    const auto& id=c.identifier;if(id!="kepler_planar"&&id!="euler_top")return {};
    AnalyticSolutionDescriptor d{id+"_regular",id+" regular periodic solution",
        id=="euler_top"?"M_j=A_j\\operatorname{sn/cn/dn}(\\nu t+u_0,m)":"M=E_a-e\\sin E_a,\\quad M'=\\sqrt{\\mu/a^3}",
        "Unforced and undamped; bound noncollision Kepler orbit or regular Euler momentum orbit",
        id=="euler_top"?"Jacobi elliptic evaluation":"analytic parametric orbit with monotone physical-clock inversion",
        {id=="euler_top"?"https://doi.org/10.3390/dynamics3030032":"https://www.math.nthu.edu.tw/~kchen/teaching/CelMech2024-4.pdf"}};
    return {{d,[=](const AnalyticState&x,double t0)->AnalyticTrajectory{
        if(c.driver.amplitude!=0||p(c,"damping")!=0)return no("Periodic celestial evaluator requires zero forcing and damping");
        if(id=="kepler_planar"){
            const double mu=p(c,"mu"),R=std::hypot(x[0],x[1]),L=x[0]*x[3]-x[1]*x[2];
            if(R==0||L==0)return no("Collision/radial Kepler orbit is outside the regular domain");
            const double E=(x[2]*x[2]+x[3]*x[3])/2-mu/R;
            if(E>=0)return no("Kepler Fourier circle requires negative energy");
            double ex=x[3]*L/mu-x[0]/R,ey=-x[2]*L/mu-x[1]/R,e=std::hypot(ex,ey);
            if(e>=1)return no("Kepler eccentricity must be below one");
            const double a=-mu/(2*E),b=a*std::sqrt(1-e*e),w=std::sqrt(mu/(a*a*a)),sign=L>0?1.:-1.;
            double ux,uy,mean;
            if(e<1e-12){e=0;ux=x[0]/R;uy=x[1]/R;mean=0;}
            else{ux=ex/e;uy=ey/e;const double cosE=(x[0]*ux+x[1]*uy)/a+e;
                const double sinE=sign*(-x[0]*uy+x[1]*ux)/b,anomaly=std::atan2(sinE,cosE);mean=anomaly-e*std::sin(anomaly);}
            const double vx=-sign*uy,vy=sign*ux;
            return periodic([=](double phase){const double u=keplerInverse(mean+phase,e),co=std::cos(u),si=std::sin(u),speed=w/(1-e*co);
                return AnalyticState{a*(co-e)*ux+b*si*vx,a*(co-e)*uy+b*si*vy,
                    speed*(-a*si*ux+b*co*vx),speed*(-a*si*uy+b*co*vy)};},w,t0);
        }
        const double I1=p(c,"inertia_1"),I2=p(c,"inertia_2"),I3=p(c,"inertia_3");
        const double C=x[0]*x[0]+x[1]*x[1]+x[2]*x[2],E=(x[0]*x[0]/I1+x[1]*x[1]/I2+x[2]*x[2]/I3)/2;
        if(C==0||std::abs(2*E*I2-C)<1e-12*C)return no("Euler equilibrium/separatrix has no regular periodic phase");
        const bool low=2*E*I2<C;
        const double a=std::sqrt(std::max(0.,I1*(2*E*I3-C)/(I3-I1)));
        const double z=std::sqrt(std::max(0.,I3*(C-2*E*I1)/(I3-I1)));
        const double b=std::sqrt(low?I2*(2*E*I3-C)/(I3-I2):I2*(C-2*E*I1)/(I2-I1));
        if(a<1e-14||b<1e-14||z<1e-14)return no("Principal-axis equilibrium has no regular phase circle");
        const double nu=std::sqrt(low?(I3-I2)*(C-2*E*I1)/(I1*I2*I3):(I2-I1)*(2*E*I3-C)/(I1*I2*I3));
        const double m=low?(I2-I1)*(2*E*I3-C)/((I3-I2)*(C-2*E*I1)):(I3-I2)*(C-2*E*I1)/((I2-I1)*(2*E*I3-C));
        const double sign=(low?x[2]:x[0])<0?-1.:1.,K=special::ellipticK(m),w=pi*nu/(2*K);
        const double u0=special::ellipticF(std::atan2(x[1]/b,low?x[0]/a:x[2]/z),m);
        return periodic([=](double phase){const auto j=special::jacobi(u0+sign*2*K*phase/pi,m);
            return low?AnalyticState{a*j.cn,b*j.sn,sign*z*j.dn}:AnalyticState{sign*a*j.dn,b*j.sn,z*j.cn};},w,t0);
    }}};
}
} // namespace
