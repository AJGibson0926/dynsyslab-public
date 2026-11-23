#include <dynsyslab/math/TorusFourier.hpp>
#include <dynsyslab/math/Perturbation.hpp>
#include <dynsyslab/math/HamiltonianPerturbation.hpp>
#include <dynsyslab/core/TorusAnalysis.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
void check(bool ok,const char*m){if(!ok)throw std::runtime_error(m);}
int main(){try{
    const double pi=std::acos(-1.);std::vector<std::complex<double>> samples;
    for(int j=0;j<16;++j)for(int k=0;k<18;++k)samples.push_back(std::cos(2*pi*j/16)+2*std::sin(4*pi*k/18));
    auto f=math::torusFourierSeries(samples,{16,18},{1,std::sqrt(2.)},{5,5});
    check(std::abs(math::evaluate(f,{.3,.7})-(std::cos(.3)+2*std::sin(1.4)))<1e-12,"2D Fourier evaluation");
    check(std::abs(math::meanSquare(f)-2.5)<1e-12,"Torus Parseval");
    auto df=math::torusTimeDerivative(f);
    check(std::abs(math::evaluate(df,{.3,.7})-(-std::sin(.3)+4*std::sqrt(2.)*std::cos(1.4)))<1e-11,"Torus generator");
    auto solved=math::solveCohomological(f,1e-12);
    auto restored=math::torusTimeDerivative(solved.solution);
    check(std::abs(math::evaluate(restored,{.8,.2})+math::evaluate(solved.retained,{.8,.2})-math::evaluate(f,{.8,.2}))<1e-12,"Cohomological identity");
    math::TorusFourierSeries resonant{{1,2},{{{2,-1},3},{{1,0},2}}};
    auto split=math::solveCohomological(resonant,1e-10);check(split.retained.coefficients.size()==1,"Resonant divisor retained");
    const auto records=math::smallDivisors({1,2},4,2);check(!records.empty()&&records[0].divisor==0,"Small divisor finds resonance");
    const auto normal=math::homogeneousNormalForm({-1.,-2.},{{1,{2,0},3.},{0,{2,0},2.}},1e-12);
    check(normal.retained.size()==1&&normal.change.size()==1&&std::abs(normal.change[0].coefficient+2.)<1e-12,"Normal form resonance and sign");
    check(math::conformalSymplecticDefect({0,1,-2,-.1},{0,1,-1,0},2,.1)<1e-15,"Conformal symplectic defect");
    check(math::conformalSymplecticDefect({0,1,-2,-.1},{0,1,-1,0},2)>0.1,"Damping detected");
    check(std::abs(math::sampledMelnikov({0,.5,1},{1,0,1,0,1,0},{0,0,0,.5,0,1})-.5)<1e-15,"Melnikov determinant and quadrature");
    // Independently test a corrected torus against Hamilton equations of
    // H=I^2/2+epsilon*cos(theta), not just against the coefficient formula.
    const double action=1.3;
    math::TorusFourierSeries perturbation{{action},{{{-1},.5},{{1},.5}}}, gradient{{action},{}};
    auto hstep=math::firstOrderHamiltonian(perturbation,{gradient},{1},1e-12);
    auto residual=[&](double epsilon) {
        const auto dI=math::torusPartialDerivative(hstep.actionCorrection[0],0);
        const auto dt=math::torusPartialDerivative(hstep.angleCorrection[0],0);
        double maximum=0;
        for(int j=0;j<101;++j) {
            const std::vector<double> theta{2*pi*j/101};
            const double mappedI=action+epsilon*math::evaluate(hstep.actionCorrection[0],theta).real();
            const double mappedTheta=theta[0]+epsilon*math::evaluate(hstep.angleCorrection[0],theta).real();
            const double a=epsilon*action*math::evaluate(dI,theta).real()-epsilon*std::sin(mappedTheta);
            const double b=action*(1+epsilon*math::evaluate(dt,theta).real())-mappedI;
            maximum=std::max(maximum,std::hypot(a,b));
        }return maximum;
    };
    const auto ratio=residual(.01)/residual(.005);
    check(ratio>3.9&&ratio<4.1,"Hamiltonian torus residual has second-order scaling");
    auto resonantStep=math::firstOrderHamiltonian(
        {{1,2},{{{2,-1},.5},{{-2,1},.5}}},
        {{{1,2},{}},{{1,2},{}}},{1,0,0,1},1e-12);
    check(resonantStep.generator.coefficients.empty()&&resonantStep.retained.coefficients.size()==2,
        "Hamiltonian resonance retained");
    auto vanishingAmplitude=math::firstOrderHamiltonian({{2},{}},
        {{{2},{{{1},.5},{{-1},.5},{{0},3}}}},{1},1e-12);
    check(vanishingAmplitude.generator.coefficients.size()==2&&
        std::abs(math::evaluate(vanishingAmplitude.angleCorrection[0],{.7}).real()-std::sin(.7)/2)<1e-12&&
        vanishingAmplitude.meanFrequencyShift[0]==3,"Zero amplitude still has nonzero action derivative");
    if(math::twistDiagnosticsAvailable()) {
        auto twist=math::hamiltonianTwist({1,1.4142},{1,0,0,2});
        check(twist.kolmogorovNumericallyNondegenerate&&twist.isoenergeticNumericallyNondegenerate,"Nondegenerate twist");
        auto isochronous=math::hamiltonianTwist({1},{0});
        check(!isochronous.kolmogorovNumericallyNondegenerate&&isochronous.isoenergeticNumericallyNondegenerate,"Distinct twist conditions");
        auto degenerate=math::hamiltonianTwist({1,2},{0,0,0,0});
        check(!degenerate.kolmogorovNumericallyNondegenerate&&!degenerate.isoenergeticNumericallyNondegenerate,"Degenerate two-torus");
    }
    core::SystemConfiguration c;c.identifier="coupled_duffing_hamiltonian";auto system=core::createSystem(c);
    auto orbit=system->analyticSolutions()[0].instantiate(system->initialState(),0);
    auto data=core::sampleInvariantTorus(*system,orbit,{system->observableCatalog().at("q_1")},{8,9});
    check(data.maximumInvarianceResidual<1e-7&&data.maximumWrapResidual<1e-9,"Product torus invariance");
    c.parameters["coupling"]=.1;auto perturbed=core::createSystem(c);
    auto mismatch=core::sampleInvariantTorus(*perturbed,orbit,{perturbed->observableCatalog().at("q_1")},{8,9});
    check(mismatch.maximumInvarianceResidual>1e-3,"Unperturbed torus defect detects Hamiltonian perturbation");
    bool rejected=false;try{math::torusFourierSeries(samples,{16,18},{1,2},{8,5});}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Nyquist rejected");
    std::cout<<"Torus analysis tests passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
