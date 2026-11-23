#include <dynsyslab/math/PeriodicClock.hpp>
#include <dynsyslab/math/AnalyticFourier.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::math;using C=std::complex<double>;
void check(bool x){if(!x)throw std::runtime_error("Periodic clock test failed");}
void near(C a,C b,double e=1e-12){check(std::abs(a-b)<e*(1+std::abs(b)));}
template<class F>void rejects(F f){bool failed=false;try{f();}catch(const std::exception&){failed=true;}check(failed);}
int main(){try{
    FourierSeries cosine{2,{{-1,.5},{1,.5}}},sine{2,{{-1,C(0,.5)},{1,C(0,-.5)}}};
    auto constant=reparameterizePeriodicFourier({cosine,sine},{2,{{0,3.}}},128,12);
    near(constant.channels[0].angularFrequency,2./3);near(constant.physicalPeriod,3*std::acos(-1.));
    near(constant.channels[0].coefficients.at(1),.5);near(constant.channels[1].coefficients.at(1),C(0,-.5));
    // Independent reference: Kepler's equation M=alpha-e*sin(alpha). Known
    // Fourier-Bessel coefficients test signs, normalization and clock Jacobian.
    for(double e:{.1,.6,.9}){
        FourierSeries w{2,{{0,1.},{-1,-e/2},{1,-e/2}}};
        auto r=reparameterizePeriodicFourier({cosine,sine},w,2048,32);
        auto c=keplerAnomalyFourier(false,e,32,2),s=keplerAnomalyFourier(true,e,32,2);
        for(int n=-32;n<=32;++n){near(r.channels[0].coefficients.at(n),c.coefficients.at(n));
            near(r.channels[1].coefficients.at(n),s.coefficients.count(n)?s.coefficients.at(n):C{});}
        for(int n:{-3,0,1,7}){
            const auto t=periodicClockCoefficientExpansion(cosine,w,n,36);
            check(std::abs(t.coefficient-c.coefficients.at(n))<=t.exactArithmeticRemainderBound+5e-12);
        }
        check(r.clockPolynomialLowerBound>0);
    }
    // General non-even clock: constant physical observable stays constant,
    // including a nonzero phase-correction constant chosen so theta(0)=0.
    FourierSeries w{2,{{0,1.},{-1,C(0,.2)},{1,C(0,-.2)}}};
    auto r=reparameterizePeriodicFourier({{2,{{0,1.}}}},w,256,20);
    near(evaluate(r.phaseCorrection,0),0);
    for(int n=-20;n<=20;++n)near(r.channels[0].coefficients.at(n),n==0?1.:0.);
    rejects([&]{reparameterizePeriodicFourier({cosine},{2,{{0,1.},{1,.1}}},128,12);});
    rejects([&]{reparameterizePeriodicFourier({cosine},{2,{{0,1.},{-1,1.},{1,1.}}},128,12);});
    rejects([&]{reparameterizePeriodicFourier({cosine},{1,{{0,1.}}},128,12);});
    rejects([&]{reparameterizePeriodicFourier({cosine},{2,{{0,1.}}},20,12);});
    std::cout<<"Periodic clock, Kepler coefficients and series remainder checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
