#include <dynsyslab/core/Controllability.hpp>
#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/core/EstimationSelection.hpp>
#include <dynsyslab/core/EstimationRealization.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
int main(){try{
    const std::vector<double> A{-1,0,0,-2},B{1,2},C{1,2};
    for(auto method:{ObservabilityMethod::classicalMatrix,ObservabilityMethod::discreteGramian,
        ObservabilityMethod::continuousFiniteGramian,ObservabilityMethod::continuousInfiniteGramian,
        ObservabilityMethod::discreteInfiniteGramian}) {
        ObservabilityOptions o;o.method=method;o.blockCount=7;o.horizon=0.8;
        const bool discrete=method==ObservabilityMethod::discreteGramian||method==ObservabilityMethod::discreteInfiniteGramian;
        const auto M=discrete?std::vector<double>{0.5,0,0,0.25}:A;
        const auto obs=analyzeObservability(M,C,2,1,o);
        const auto ctr=analyzeControllability(M,B,2,1,o);
        check(obs.rank==2&&ctr.rank==2,"full-rank dual pair lost rank");
        for(std::size_t k=0;k<4;++k)check(std::abs(obs.gramian[k]-ctr.gramian[k])<1e-12,"Gramian duality");
        if(method==ObservabilityMethod::continuousInfiniteGramian) {
            check(std::abs(ctr.gramian[0]-0.5)<1e-12,"continuous diagonal");
            check(std::abs(ctr.gramian[1]-2.0/3)<1e-12,"continuous cross term");
        }
        if(method==ObservabilityMethod::discreteInfiniteGramian) {
            check(std::abs(ctr.gramian[0]-1.0/0.75)<1e-12,"discrete diagonal");
            check(std::abs(ctr.gramian[1]-2.0/0.875)<1e-12,"discrete cross term");
        }
    }
    ObservabilityOptions longWindow;longWindow.method=ObservabilityMethod::continuousFiniteGramian;longWindow.horizon=100;
    check(std::abs(analyzeObservability({-100.},{1.},1,1,longWindow).gramian[0]-0.005)<1e-12,
        "stiff stable finite Gramian overflowed");
    const auto pbh=pbhObservability({0.5,0,0,1.1},{1,0},2,1,true);
    check(!pbh.fullRank&&!pbh.stableHiddenModes,"discrete PBH missed hidden unstable mode");
    const auto stableHidden=pbhControllability({-1,0,0,-2},{1,0},2,1);
    check(!stableHidden.fullRank&&stableHidden.stableHiddenModes,"PBH stabilizability");
    const auto hidden=analyzeControllability(A,{1,0},2,1);
    check(hidden.rank==1&&std::isinf(hidden.conditionNumber),"rank deficient conditioning");
    check(hidden.leftDirections.size()==4,"controllable directions missing");
    bool rejected=false;ObservabilityOptions o;o.method=ObservabilityMethod::discreteInfiniteGramian;
    try{(void)analyzeObservability({1,0,0,0.5},C,2,1,o);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"marginal infinite Gramian accepted");
    SystemConfiguration config;config.identifier="duffing";config.preset="double_well_reference";
    auto model=KoopmanModel::fromPolynomialCoefficients(config,2,{1,0,0,1},
        {{1,0},{0,0},{0,0},{1,0}},{{-1,0},{-2,0}});
    const std::vector<std::complex<double>> psi{1.,0.,0.,1.,2.,3.},y{4.,-2.,2.};
    const auto fit=fitAndSetObservableReconstruction(model,"sensor",psi,3,y,1);
    check(fit.residualRms<1e-12&&model.hasReconstruction("sensor"),"arbitrary reconstruction fit");
    check(std::abs(fit.matrix[0]-4.)<1e-12&&std::abs(fit.matrix[1]+2.)<1e-12,"reconstruction coefficients");
    const auto obs=analyzeObservability(model,"sensor");check(obs.rank==2,"KROM observability");
    ObservabilityOptions modalOptions;modalOptions.method=ObservabilityMethod::continuousFiniteGramian;
    const auto modal=analyzeModalObservability(model,"sensor",modalOptions);
    const auto direct=analyzeObservability(model,"sensor",modalOptions);
    for(std::size_t i=0;i<4;++i)check(std::abs(modal.gramian[i]-direct.gramian[i])<1e-10,"modal Gramian formula");
    model.setControlMap("constant",ConstantControlMap{{1.,2.},2,1,"test"});
    check(analyzeControllability(model,"constant").rank==2,"KROM controllability");
    model.setReconstruction("state",{{1.,0.,0.,1.},2,2});
    EstimationSelectionOptions selection;selection.measurement="sensor";
    selection.outputScales={1,1};selection.horizon=0.2;selection.sampleTime=0.01;
    const std::vector<double> t{0.,0.1,0.2};
    std::vector<std::complex<double>> exactPsi;std::vector<double> exactX,exactY;
    for(double time:t){double a=std::exp(-time),b=std::exp(-2*time);
        exactPsi.insert(exactPsi.end(),{a,b});exactX.insert(exactX.end(),{a,b});exactY.push_back(4*a-2*b);}
    const auto assessment=assessEstimationSubset(model,t,exactPsi,exactX,exactY,selection);
    check(assessment.stateReconstructionMse<1e-24&&assessment.autonomousPropagationMse<1e-24,"exact subset score");
    check(assessment.posteriorStateNoiseMse>0&&std::isfinite(assessment.score),"subset noise score");
    auto rescaled=model;
    rescaled.setReconstruction("state",{{.1,0.,0.,10.},2,2});
    rescaled.setReconstruction("sensor",{{.4,-20.},1,2});
    auto scaledPsi=exactPsi;
    for(std::size_t k=0;k<t.size();++k){scaledPsi[2*k]*=10.;scaledPsi[2*k+1]*=.1;}
    const auto scaledAssessment=assessEstimationSubset(rescaled,t,scaledPsi,exactX,exactY,selection);
    check(std::abs(scaledAssessment.score-assessment.score)<1e-10,"selection changed under modal rescaling");
    auto hiddenModel=model;hiddenModel.setReconstruction("sensor",{{1.,0.},1,2});
    std::vector<double> firstChannel;for(std::size_t k=0;k<t.size();++k)firstChannel.push_back(exactX[2*k]);
    const auto hiddenAssessment=assessEstimationSubset(hiddenModel,t,exactPsi,exactX,firstChannel,selection);
    const double hiddenPrior=.5*(1+std::exp(-.4)+std::exp(-.8))/3*std::exp(-.8);
    check(hiddenAssessment.posteriorStateNoiseMse>=hiddenPrior*(1-1e-10),"hidden state uncertainty disappeared through Gramian pseudoinverse");
    ScatteredEigenfunctionData paired;paired.samples=3;paired.stateDimension=2;
    paired.trainingStates={1,0,0,1,-1,0};paired.continuousEigenvalues={{-1,2},{-1,-2}};
    paired.trainingValues={{1,0},{0,2},{0,1},{2,0},{-1,0},{0,-2}};
    KoopmanModel pairModel(paired,config);
    pairModel.setReconstruction("y",{{{2,3},{-1,4}},1,2});
    const auto realModel=makeRealModalRealization(pairModel,{0});
    check(realModel.model.splitDimension()==2,"conjugate pair double counted");
    const auto before=reconstructObservable(pairModel.reconstruction("y"),paired.trainingValues,3);
    const auto after=reconstructObservable(realModel.model.reconstruction("y"),
        realModel.model.scatteredEigenfunctions()->trainingValues,3);
    for(std::size_t k=0;k<3;++k)check(std::abs(before[k].real()-after[k].real())<1e-12,"real modal reconstruction changed");
    ScatteredEigenfunctionData nominal;nominal.samples=3;nominal.stateDimension=2;
    nominal.trainingStates=paired.trainingStates;nominal.continuousEigenvalues={{-0.01,-1e-9}};
    nominal.trainingValues={{0,1},{0,2},{0,3}};
    KoopmanModel nominalModel(nominal,config);nominalModel.setReconstruction("y",{{{0,-2}},1,1});
    const auto snapped=makeRealModalRealization(nominalModel,{0},1e-6,1e-8);
    check(snapped.model.splitDimension()==1,"nominal real mode not corrected within tolerance");
    const auto snappedValues=reconstructObservable(snapped.model.reconstruction("y"),snapped.model.scatteredEigenfunctions()->trainingValues,3);
    check(std::abs(snappedValues[2].real()-6)<1e-12,"real-phase correction changed reconstruction");
    std::cout<<"dynsyslab_linear_analysis_tests: PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
