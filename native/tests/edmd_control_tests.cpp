#include <dynsyslab/core/EdmdControl.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;
void require(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
template<class F> void rejects(F f) {
    bool failed=false; try {f();} catch(const std::invalid_argument&) {failed=true;}
    require(failed,"invalid EDMDc input accepted");
}
void defaultsAndSharedRepresentation() {
    constexpr std::size_t count=240;
    std::vector<double> x(count*2),y(count*2),u(count);
    // Four-degree lifting of a forced linear system is NOT additive in u.
    // Nevertheless its degree-one rows have an exact additive representation.
    for(std::size_t k=0;k<count;++k) {
        x[2*k]=10+2*std::sin(.37*k); x[2*k+1]=-4+3*std::cos(.23*k);
        u[k]=std::cos(.71*k)+.3*std::sin(.57*k);
        y[2*k]=.8*x[2*k]+.1*x[2*k+1]+.2*u[k];
        y[2*k+1]=-.1*x[2*k]+.9*x[2*k+1]-.15*u[k];
    }
    const auto fit=core::extendedDynamicModeDecompositionWithControl(x,y,u,count,2,1,.1);
    require(fit.dictionary.degree==4&&fit.dictionary.includeConstant,"wrong defaults");
    require(fit.lifted.dimension==15&&fit.lifted.fullRowRank,"default library lacks mixed monomials");
    require(fit.dictionary.exponents==core::polynomialExponentTable(2,4,true),"not SINDy order");
    require(fit.snapshotOneStepRelativeError<1e-10,"degree-one snapshot decoder failed");
    require(fit.lifted.oneStepRelativeError>1e-5,"input-dependent lifting error was hidden");
    require(fit.lifted.hasEigenfunctions,"expected nondefective eigenbasis");
    const auto polynomial=fit.expansion();
    const auto values=core::evaluateEdmdControlEigenfunctions(fit,x,count);
    require(polynomial.samples==0&&polynomial.residualRms.empty(),"fabricated SINDy fit diagnostics");
    for(std::size_t mode=0;mode<fit.lifted.rank;++mode) {
        require(std::abs(polynomial.discreteEigenvalues()[mode]-fit.lifted.discreteEigenvalues[mode])<1e-10,
            "shared spectrum conversion failed");
        for(std::size_t k=0;k<10;++k)
            require(std::abs(core::evaluatePolynomialEigenfunction(polynomial,mode,x.data()+2*k)
                -values[k*fit.lifted.rank+mode])<1e-10,"shared polynomial evaluator disagrees");
    }
    // Out-of-sample evaluation must use frozen training normalization.
    const auto theta=core::liftEdmdControlSnapshots(fit.dictionary,{50,-30},1);
    const double a=(50-fit.dictionary.center[0])/fit.dictionary.scale[0];
    const double b=(-30-fit.dictionary.center[1])/fit.dictionary.scale[1];
    require(std::abs(theta[4]-a*b)<1e-12,"mixed term / frozen scaling incorrect");
    rejects([&]{(void)core::liftEdmdControlSnapshots(fit.dictionary,{1},1);});
    auto bad=fit.dictionary; bad.scale[0]=0;
    rejects([&]{(void)core::liftEdmdControlSnapshots(bad,{1,2},1);});
    bad=fit.dictionary; bad.exponents[1]=100;
    rejects([&]{(void)core::liftEdmdControlSnapshots(bad,{1,2},1);});
    core::EdmdControlOptions options; options.polynomialDegree=0;
    rejects([&]{(void)core::extendedDynamicModeDecompositionWithControl(x,y,u,count,2,1,.1,options);});
    options.polynomialDegree=33;
    rejects([&]{(void)core::extendedDynamicModeDecompositionWithControl(x,y,u,count,2,1,.1,options);});
    options.polynomialDegree=4; options.regression.rank=3;
    const auto reduced=core::extendedDynamicModeDecompositionWithControl(x,y,u,count,2,1,.1,options);
    require(reduced.lifted.rank==3&&reduced.lifted.dimension==15,"rank confused with dictionary size");
    options.polynomialDegree=1; options.normalize=false; options.includeConstant=false;
    const auto linear=core::dynamicModeDecompositionWithControl(x,y,u,count,2,1,.1,options.regression);
    const auto equivalent=core::extendedDynamicModeDecompositionWithControl(x,y,u,count,2,1,.1,options);
    const auto p=core::predictDmdControl(linear,x,u,count);
    const auto q=core::predictEdmdControl(equivalent,x,u,count);
    for(std::size_t k=0;k<p.size();++k) require(std::abs(p[k]-q[k])<1e-12,"degree-one DMDc parity failed");
}
void nonlinearAdvantage() {
    constexpr std::size_t count=200;
    std::vector<double> x(count),y(count),u(count);
    for(std::size_t k=0;k<count;++k) {
        x[k]=std::sin(.39*k); u[k]=std::cos(.73*k);
        y[k]=.55*x[k]+.12*x[k]*x[k]+.05*u[k];
    }
    const auto lifted=core::extendedDynamicModeDecompositionWithControl(x,y,u,count,1,1,.1);
    const auto linear=core::dynamicModeDecompositionWithControl(x,y,u,count,1,1,.1);
    require(lifted.snapshotOneStepRelativeError<1e-10&&linear.oneStepRelativeError>.05,
        "polynomial EDMDc did not capture nonlinear snapshot mapping");
    std::vector<double> testX{-.87,.33,.77},testU{.17,-.43,.29};
    const auto predicted=core::predictEdmdControl(lifted,testX,testU,3);
    for(std::size_t i=0;i<3;++i)
        require(std::abs(predicted[i]-(.55*testX[i]+.12*testX[i]*testX[i]+.05*testU[i]))<1e-9,
            "held-out nonlinear prediction failed");
}
void freePolynomialSemigroup() {
    std::vector<double> x(200),y(200),u(200);
    // Known degree-four invariant dictionary. Independent recorded inputs have
    // zero actual effect, so the exact identified B must vanish.
    for(std::size_t k=0;k<x.size();++k) {
        x[k]=std::sin(.37*k); y[k]=.8*x[k]; u[k]=std::cos(.63*k);
    }
    const auto fit=core::extendedDynamicModeDecompositionWithControl(x,y,u,200,1,1,.1);
    require(fit.lifted.oneStepRelativeError<1e-12,"closed polynomial dictionary was not recovered");
    for(const double b:fit.lifted.controlMatrix) require(std::abs(b)<1e-12,"spurious input coupling");
    const auto initial=core::evaluateEdmdControlEigenfunctions(fit,{.7},1);
    double state=.7;
    for(int k=0;k<30;++k) {
        const auto observed=core::evaluateEdmdControlEigenfunctions(fit,{state},1);
        for(std::size_t j=0;j<fit.lifted.rank;++j)
            require(std::abs(observed[j]-initial[j]*std::pow(fit.lifted.discreteEigenvalues[j],k))<1e-10,
                "held-out free eigenfunction did not follow its complex exponential");
        state*=.8;
    }
}
} // namespace
int main() {
    try {defaultsAndSharedRepresentation(); nonlinearAdvantage(); freePolynomialSemigroup();
        std::cout<<"EDMDc: PASS\n"; return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
