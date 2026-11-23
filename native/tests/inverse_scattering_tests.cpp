#include <dynsyslab/math/InverseScattering.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::math;
void check(bool b,const char*m){if(!b)throw std::runtime_error(m);}
int main(){try{
    auto free=rightScattering(std::vector<double>(100,0),-5,.1,1.3);
    check(std::abs(free.reflection)<1e-13&&std::abs(free.transmission-1.)<1e-13,"Free scattering phase");
    auto hill=hillMonodromy(std::vector<double>(100,0),.1,1.69);
    check(std::abs(hill.discriminant-2*std::cos(13.))<1e-12,"Free Hill discriminant");
    const ReflectionlessData data{{1},{std::log(2.)}};
    for(double x:{-5.,-.3,0.,.7,3.})check(std::abs(reflectionlessPotential(evolveKdV(data,.2),x)+2/std::pow(std::cosh(x-.8),2))<1e-12,"One-soliton inverse and evolution");
    double previous=1;for(int n:{1000,2000,4000}){const double h=24./n;std::vector<double>v;
        for(int j=0;j<n;++j)v.push_back(-2/std::pow(std::cosh(-12+(j+.5)*h),2));
        const auto s=rightScattering(v,-12,h,1.1);check(s.fluxDefect<1e-10,"Scattering flux");
        check(std::abs(s.reflection)<previous*.4,"Reflectionless potential grid convergence");previous=std::abs(s.reflection);}
    std::vector<double>v;const int n=1600;for(int j=1;j<=n;++j)v.push_back(-2/std::pow(std::cosh(-12+j*24./(n+1)),2));
    const auto eig=dirichletNegativeSpectrum(v,24);check(!eig.empty()&&std::abs(eig[0]+1)<3e-5,"Dirichlet bound-state convergence");
    ReflectionlessData two{{.7,1.2},{std::log(1.4),std::log(2.4)}};
    for(double x:{-1.,0.,1.}){const double h=.002,dt=1e-5,t=.1;
        auto at=[&](double a,double b){return reflectionlessPotential(evolveKdV(two,b),a);};
        const double u=at(x,t),ux=(at(x+h,t)-at(x-h,t))/(2*h);
        const double uxxx=(at(x+2*h,t)-2*at(x+h,t)+2*at(x-h,t)-at(x-2*h,t))/(2*h*h*h);
        const double ut=(at(x,t+dt)-at(x,t-dt))/(2*dt);check(std::abs(ut-6*u*ux+uxxx)<.002,"Two-soliton KdV residual");}
    bool rejected=false;try{rightScattering({0},0,1,0);}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Zero-energy scattering rejected");
    std::cout<<"Inverse scattering tests passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
