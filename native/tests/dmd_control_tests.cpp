#include <dynsyslab/core/DmdControl.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class F> void rejects(F f) {
    bool threw=false; try { f(); } catch(const std::invalid_argument&) { threw=true; }
    require(threw,"invalid data were accepted");
}

// Independent pairs are also a pooled multi-trajectory dataset: no endpoints
// are connected. A rotated contraction has genuinely complex eigenfunctions.
void testRecoveryAndEigenfunctions() {
    constexpr std::size_t count=151,n=2,m=2;
    const double a[]{.91,-.18,.18,.91},b[]{.4,-.1,.2,.7};
    std::vector<double> x(count*n),y(count*n),u(count*m);
    for(std::size_t k=0;k<count;++k) {
        x[2*k]=std::sin(.37*k); x[2*k+1]=std::cos(.23*k);
        u[2*k]=std::cos(.71*k); u[2*k+1]=std::sin(.59*k);
        for(std::size_t i=0;i<n;++i)
            for(std::size_t j=0;j<n;++j) y[2*k+i]+=a[2*i+j]*x[2*k+j]+b[2*i+j]*u[2*k+j];
    }
    const auto fit=core::dynamicModeDecompositionWithControl(x,y,u,count,n,m,.02);
    require(fit.rank==2&&fit.augmentedRank==4&&fit.fullRowRank,"wrong DMDc ranks");
    require(fit.oneStepRelativeError<1e-12&&fit.projectionRelativeError<1e-12,"linear recovery failed");
    require(fit.hasEigenfunctions,"complex linear eigenfunctions missing");
    const auto predicted=core::predictDmdControl(fit,x,u,count);
    const auto values=core::evaluateDmdControlEigenfunctions(fit,x,count);
    const auto polynomial=fit.expansion();
    for(std::size_t mode=0;mode<fit.rank;++mode)
        require(std::abs(core::evaluatePolynomialEigenfunction(polynomial,mode,x.data())
            -values[mode])<1e-12,"DMDc shared polynomial export failed");
    require(values.size()==count*fit.rank,"wrong eigenfunction evaluation shape");
    for(std::size_t i=0;i<y.size();++i) require(std::abs(predicted[i]-y[i])<1e-12,"prediction API failed");
    for(std::size_t mode=0;mode<n;++mode) {
        require(std::abs(fit.discreteEigenvalues[mode].real()-.91)<1e-12
            &&std::abs(std::abs(fit.discreteEigenvalues[mode].imag())-.18)<1e-12,"wrong spectrum");
        for(std::size_t j=0;j<n;++j) {
            std::complex<double> la=0,lb=0;
            for(std::size_t i=0;i<n;++i) {
                const auto l=fit.eigenfunctionCoefficients[mode*n+i];
                la+=l*a[i*n+j]; lb+=l*b[i*m+j];
            }
            require(std::abs(la-fit.discreteEigenvalues[mode]*fit.eigenfunctionCoefficients[mode*n+j])<1e-12,
                "left eigenfunction coefficients violate L A = Lambda L");
            require(std::abs(lb-fit.modalInputMatrix[mode*m+j])<1e-12,"modal B is not L B");
        }
        for(std::size_t k=0;k<count;++k) {
            std::complex<double> psi=0,next=0,drive=0;
            for(std::size_t j=0;j<n;++j) {
                psi+=fit.eigenfunctionCoefficients[mode*n+j]*x[k*n+j];
                next+=fit.eigenfunctionCoefficients[mode*n+j]*y[k*n+j];
                drive+=fit.modalInputMatrix[mode*m+j]*u[k*m+j];
            }
            require(std::abs(next-fit.discreteEigenvalues[mode]*psi-drive)<1e-12,
                "forced eigenfunction evolution failed");
        }
    }
    // Reconstruct physical A/B from reduced matrices without storing n-by-n A.
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<n;++j) {
        double aa=0,bb=0;
        for(std::size_t p=0;p<n;++p) {
            bb+=fit.basis[i*n+p]*fit.controlMatrix[p*m+j];
            for(std::size_t q=0;q<n;++q)
                aa+=fit.basis[i*n+p]*fit.operatorMatrix[p*n+q]*fit.basis[j*n+q];
        }
        require(std::abs(aa-a[i*n+j])<1e-12&&std::abs(bb-b[i*n+j])<1e-12,"A/B recovery failed");
    }
    core::DmdControlOptions reduced; reduced.rank=1; reduced.augmentedRank=3;
    const auto shortFit=core::dynamicModeDecompositionWithControl(x,y,u,count,n,m,.02,reduced);
    require(shortFit.rank==1&&shortFit.augmentedRank==3,"rank truncation ignored");
    require(shortFit.projectionRelativeError>.1,"rank reduction lost no energy unexpectedly");
    // Held-out, zero-input trajectory: ideal modal behavior is mu^k psi(q0).
    std::vector<double> state{.43,-.69};
    auto ideal=core::evaluateDmdControlEigenfunctions(fit,state,1);
    for(int k=0;k<100;++k) {
        state=core::predictDmdControl(fit,state,{0,0},1);
        const auto actual=core::evaluateDmdControlEigenfunctions(fit,state,1);
        for(std::size_t j=0;j<n;++j) {
            ideal[j]*=fit.discreteEigenvalues[j];
            require(std::abs(actual[j]-ideal[j])<1e-12,"unforced semigroup law failed");
        }
    }
    // A Jordan block may be identified accurately without providing a safe
    // eigenfunction coordinate chart. A/B must remain available in that case.
    for(std::size_t k=0;k<count;++k) {
        y[2*k]=.8*x[2*k]+x[2*k+1]+u[2*k];
        y[2*k+1]=.8*x[2*k+1]+u[2*k+1];
    }
    const auto jordan=core::dynamicModeDecompositionWithControl(x,y,u,count,n,m,.02);
    require(jordan.oneStepRelativeError<1e-12&&!jordan.hasEigenfunctions,
        "defective eigenbasis was advertised as usable");
    rejects([&]{(void)core::evaluateDmdControlEigenfunctions(jordan,x,count);});
}

void testIdentifiabilityAndValidation() {
    const std::vector<double> x{1,2,3,4,5},u{.2,-.4,.7,-.1,.3};
    std::vector<double> y;
    for(std::size_t k=0;k<x.size();++k) y.push_back(.8*x[k]+.3*u[k]);
    const auto fit=core::dynamicModeDecompositionWithControl(x,y,u,5,1,1,1);
    require(fit.fullRowRank&&fit.rank==1&&std::abs(fit.operatorMatrix[0]-.8)<1e-12,"scalar DMDc failed");
    const auto ambiguous=core::dynamicModeDecompositionWithControl(x,y,x,5,1,1,1);
    require(!ambiguous.fullRowRank&&ambiguous.augmentedNumericalRank==1,"feedback collinearity not detected");
    rejects([&]{(void)core::dynamicModeDecompositionWithControl(x,y,std::vector<double>(5),5,1,1,1);});
    rejects([&]{(void)core::dynamicModeDecompositionWithControl(x,y,u,4,1,1,1);});
    rejects([&]{(void)core::dynamicModeDecompositionWithControl(x,y,u,5,1,1,0);});
    auto bad=x; bad[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{(void)core::dynamicModeDecompositionWithControl(bad,y,u,5,1,1,1);});
    core::DmdControlOptions opts; opts.singularValueTolerance=-1;
    rejects([&]{(void)core::dynamicModeDecompositionWithControl(x,y,u,5,1,1,1,opts);});
    // Zero eigenvalue is a valid discrete model; no log(0) in the core.
    const auto zero=core::dynamicModeDecompositionWithControl(x,u,u,5,1,1,1);
    require(std::abs(zero.discreteEigenvalues[0])<1e-14,"zero discrete multiplier mishandled");
}
} // namespace
int main() {
    try {
        testRecoveryAndEigenfunctions(); testIdentifiabilityAndValidation();
        std::cout<<"dynsyslab_dmd_control_tests: PASS\n"; return 0;
    } catch(const std::exception& e) {
        std::cerr<<"dynsyslab_dmd_control_tests: FAIL: "<<e.what()<<'\n'; return 1;
    }
}
