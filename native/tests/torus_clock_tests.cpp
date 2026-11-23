#include <dynsyslab/math/TorusClock.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::math;
namespace {
using C=std::complex<double>;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void near(C a,C b,double tol=1e-10){check(std::abs(a-b)<tol,"Torus clock comparison failed");}
template<class F> void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,"Invalid clock accepted");}
void clockChecks() {
    TorusFourierSeries clock{{2.},{{{0},1.3},{{-1},.2},{{1},.2}}};
    auto analysis=analyzeTorusClock(clock);
    near(analysis.lowerBound,.9);near(analysis.meanRate,1.3);
    const std::vector<double> times={0.,3.,-2.,.1,5.};
    auto run=evolveTorusClock(clock,{.7},times,{{1},{-3}},1e-12,1e-13);
    for(const auto& row:run.samples) {
        near(row.straightenedPhases[0],run.samples[0].straightenedPhases[0]+2*row.time/1.3,2e-12);
        for(std::size_t i=0;i<2;++i) {
            const double k=i==0?1:-3;
            near(row.correctedCharacters[i],run.samples[0].correctedCharacters[i]*std::exp(C(0,k*2*row.time/1.3)),5e-12);
        }
    }
    // Independent integration of the time-changed ODE; no primitive or inverse.
    double angle=.7;const double dt=1e-4;
    auto drift=[](double a){return 2/(1.3+.4*std::cos(a));};
    for(int j=0;j<30000;++j) {
        const double a=drift(angle),b=drift(angle+dt*a/2),c=drift(angle+dt*b/2),d=drift(angle+dt*c);
        angle+=dt*(a+2*b+2*c+d)/6;
    }
    near(run.samples[1].phases[0],angle,2e-11);
    // Group property and negative times.
    auto restarted=evolveTorusClock(clock,run.samples[1].phases,{-3.,2.});
    near(restarted.samples[0].phases[0],.7,2e-11);
    near(restarted.samples[1].phases[0],run.samples.back().phases[0],2e-11);
    // Corrected eigenfunctions are orthonormal in rho/mean(rho) Haar measure.
    C weighted{},unweighted{};const double pi=std::acos(-1.);
    for(int j=0;j<4096;++j) {
        const double a=2*pi*j/4096,eta=a+2*evaluate(analysis.primitive,{a}).real()/1.3;
        const C psi=std::exp(C(0,eta));
        weighted+=psi*(1.3+.4*std::cos(a))/1.3/4096.;
        unweighted+=psi/4096.;
    }
    near(weighted,0,1e-12);check(std::abs(unweighted)>.1,"Haar change is visible");
    // Cover refinement changes the factor coordinate, not the physical result.
    auto lifted=clock;lifted.angularFrequencies[0]/=2;lifted.coefficients.clear();
    for(const auto& [k,c]:clock.coefficients)lifted.coefficients[{2*k[0]}]=c;
    auto other=evolveTorusClock(lifted,{.35},times,{{2},{-6}});
    for(std::size_t j=0;j<times.size();++j) {
        near(other.samples[j].baseElapsed,run.samples[j].baseElapsed);
        near(other.samples[j].correctedCharacters[0],run.samples[j].correctedCharacters[0]);
    }
    // Two-angle clocks share one elapsed-time variable.
    TorusFourierSeries two{{1.,std::sqrt(2.)},{{{0,0},1.},{{1,-1},.1},{{-1,1},.1}}};
    auto twoRun=evolveTorusClock(two,{.1,.2},{0,7},{{2,3}});
    const double frequency=2+3*std::sqrt(2.);
    near(twoRun.samples[1].correctedCharacters[0],twoRun.samples[0].correctedCharacters[0]*std::exp(C(0,7*frequency)));
    auto invalid=clock;invalid.coefficients[{1}]=C(.2,.1);rejects([&]{analyzeTorusClock(invalid);});
    invalid=clock;invalid.coefficients[{0}]=.3;rejects([&]{analyzeTorusClock(invalid);});
    two.angularFrequencies={1,1};rejects([&]{analyzeTorusClock(two);});
    rejects([&]{evolveTorusClock(clock,{.1,.2},{0});});
    rejects([&]{evolveTorusClock(clock,{.1},{0},{{1,1}});});
    rejects([&]{evolveTorusClock(clock,{.1},{0},{},1e-12,0);});
    // Finite approximations can have increasing primitive norms despite rapidly
    // converging positive clocks. This check does not certify an infinite limit.
    for(int depth:{1,4,12}) {
        TorusFourierSeries hierarchy{{std::ldexp(1.,-depth)},{{{0},1.}}};
        for(int j=1;j<=depth;++j) {
            const int k=1<<(depth-j);const double c=.2*std::ldexp(1.,-j)/2;
            hierarchy.coefficients[{k}]=c;hierarchy.coefficients[{-k}]=c;
        }
        auto a=analyzeTorusClock(hierarchy);
        near(meanSquare(a.primitive),depth*.2*.2/2);
        check(a.lowerBound>=.8,"Uniform positive clock bound");
    }
}
}
int main(){try{clockChecks();std::cout<<"Torus clock conjugacy, independent ODE and invariant-measure checks passed\n";return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
