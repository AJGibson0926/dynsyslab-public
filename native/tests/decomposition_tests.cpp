#include <dynsyslab/core/Decomposition.hpp>
#include <dynsyslab/core/PCA.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        constexpr std::size_t samples=240;
        std::vector<double> time(samples),states(samples*2);
        for (std::size_t i=0;i<samples;++i) {
            time[i]=0.01*static_cast<double>(i);
            states[2*i]=std::sin(time[i]); states[2*i+1]=std::cos(time[i]);
        }
        dynsyslab::core::EDMDOptions edmdOptions;
        edmdOptions.dictionary.polynomialDegree=2; edmdOptions.regularization=1.0e-10;
        const auto edmd=dynsyslab::core::extendedDynamicModeDecomposition(
            states,samples,2,0.01,edmdOptions);
        require(edmd.featureCount==6,"quadratic two-state dictionary must have six features");
        require(edmd.spectrum.values().size()==6,"EDMD spectrum shape is wrong");
        require(edmd.eigenfunctionValues.size()==samples*6,"EDMD values shape is wrong");
        double predictionError=0.0,referenceEnergy=0.0;
        auto features=[](double x,double y) {
            return std::vector<double>{1.0,y,x,y*y,x*y,x*x};
        };
        for(std::size_t sample=0;sample+1<samples;++sample) {
            const auto source=features(states[2*sample],states[2*sample+1]);
            const auto target=features(states[2*(sample+1)],states[2*(sample+1)+1]);
            for(std::size_t output=0;output<6;++output) {
                double predicted=0.0;
                for(std::size_t input=0;input<6;++input)
                    predicted+=source[input]*edmd.operatorMatrix[input*6+output];
                const double residual=predicted-target[output];
                predictionError+=residual*residual;
                referenceEnergy+=target[output]*target[output];
            }
        }
        require(std::sqrt(predictionError/referenceEnergy)<1.0e-10,
                "EDMD one-step dictionary reconstruction changed");
        const auto havok=dynsyslab::core::havokDecomposition(
            states,time,samples,2,dynsyslab::core::HAVOKOptions{30,8,0});
        require(havok.rank==8 && havok.embeddedSamples==211,"HAVOK geometry is wrong");
        require(havok.generator.size()==49 && havok.input.size()==7,"HAVOK model shape is wrong");
        require(havok.coordinates.size()==211*8,"HAVOK coordinate shape is wrong");
        for (const double value:havok.generator) require(std::isfinite(value),"HAVOK generator is non-finite");
        const std::vector<double> collinear{
            1.0,2.0, 2.0,4.0, 3.0,6.0, 4.0,8.0};
        const auto pca=dynsyslab::core::principalComponentAnalysis(
            collinear,4,2,{1,true});
        require(pca.rank==1 && pca.components.size()==2
            && pca.scores.size()==4,"PCA result geometry is wrong");
        require(std::abs(pca.explainedVarianceRatio[0]-1.0)<1.0e-12,
                "PCA must explain the collinear sample exactly");
        dynsyslab::core::BaggedPCAOptions bagOptions;
        bagOptions.rank=1; bagOptions.bags=8; bagOptions.seed=73;
        const auto bagged=dynsyslab::core::baggedPrincipalComponentAnalysis(
            collinear,4,2,bagOptions);
        require(bagged.seed==73 && bagged.bags==8,
                "bagged PCA did not retain reproducibility metadata");
        require(std::abs(bagged.decomposition.explainedVarianceRatio[0]-1.0)<1.0e-12,
                "bagged PCA must preserve a collinear direction");
        std::vector<double> corrupted{
            1.0,.5,-1.0, 2.0,11.0,-2.0, 3.0,1.5,-3.0,
            4.0,2.0,-4.0, 5.0,2.5,-12.0};
        const auto robust=dynsyslab::core::robustPrincipalComponentAnalysis(
            corrupted,5,3,{});
        require(robust.converged && robust.relativeResidual<1.0e-7,
                "robust PCA did not converge");
        require(std::abs(robust.sparse[4]-10.0)<1.0e-4
            && std::abs(robust.sparse[14]+7.0)<1.0e-4,
                "robust PCA did not isolate the reference outliers");
        require(std::abs(robust.lowRank[13]-2.5)<1.0e-4,
                "robust PCA low-rank reference parity changed");
        const auto hankel=dynsyslab::core::hankelProperOrthogonalDecomposition(
            states,samples,2,{10,2,4,true});
        require(hankel.embeddedSamples==222 && hankel.embeddedFeatures==20
            && hankel.decomposition.rank==4,
                "Hankel POD geometry is wrong");
        std::cout<<"native EDMD and HAVOK tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"decomposition test failure: "<<error.what()<<'\n'; return 1;
    }
}
