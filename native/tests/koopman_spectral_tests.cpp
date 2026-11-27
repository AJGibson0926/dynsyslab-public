#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <chrono>
namespace {
namespace c=dynsyslab::core;constexpr double pi=3.14159265358979323846;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"Expected validation rejection");}
c::KoopmanSnapshotPairs rotation(double factor=1) {
    c::KoopmanSnapshotPairs data;data.samples=128;data.dimension=2;
    for(std::size_t i=0;i<data.samples;++i) {
        const double t=2*pi*i/data.samples;
        data.x.insert(data.x.end(),{std::cos(t),std::sin(t)});
        data.y.insert(data.y.end(),{factor*std::cos(t+.31),factor*std::sin(t+.31)});
    }
    return data;
}
c::SpectralDictionary linear(){c::SpectralDictionary d;d.dimension=2;d.constant=false;return d;}
void rotationAndRollout() {
    for(const std::string method:{"resdmd","edmd","mpedmd","tlsdmd"}) {
        auto data=rotation();c::KoopmanSpectralOptions options;options.method=method;
        options.acknowledgeMeasurePreserving=true;options.rank=2;options.sampleInterval=.1;
        options.residualGrid={std::exp(std::complex<double>(0,.31)),{0,0}};
        auto model=c::fitKoopmanSpectral(data,linear(),options);
        require(model.rank==2,"Rotation rank");
        for(auto z:model.multipliers)require(std::abs(std::abs(z)-1)<1e-10,"Rotation eigenvalue radius");
        for(auto res:model.spectral.residuals.eigenpair)require(res<1e-10,"Analytic rotation eigenpair residual");
        require(model.spectral.residuals.pseudospectrum.minimum[0]<1e-7&&std::abs(model.spectral.residuals.pseudospectrum.minimum[1]-1)<1e-10,"Pseudospectral distances for normal rotation");
        require(model.spectral.invariance.projectionDefect<1e-10,"Invariant coordinate span");
        c::fitSpectralReconstruction(model,"position",data.x,data.x,2);
        auto rollout=c::rolloutSpectralKrom(model,"position",{1,0},40);
        for(std::size_t k=0;k<=40;++k)require(std::abs(rollout[2*k]-std::cos(.31*k))<1e-9&&std::abs(rollout[2*k+1]-std::sin(.31*k))<1e-9,"Autonomous rotation forecast");
        auto restored=c::readSpectralKromArchive(c::makeSpectralKromArchive(model));
        require(c::spectralKromIdentity(model)==c::spectralKromIdentity(restored),"Archive identity round trip");
        c::ObservableDecompositionReference ref;ref.modelIdentity=c::spectralKromIdentity(model);ref.component="position";
        auto observable=c::spectralKromObservable(model,ref,"position");
        double point[]{.3,.4};auto values=observable.evaluate(point);
        require(std::abs(values[0]-.3)<1e-10&&std::abs(values[1]-.4)<1e-10,"Bound observable readout");
        ref.modelIdentity="stale";rejects([&]{c::spectralKromObservable(model,ref,"bad");});
        auto held=rotation(.8);auto residual=c::koopmanResiduals(model,held);
        for(auto value:residual.eigenpair)require(std::abs(value-.2)<1e-10,"Held-out residual detects changed dynamics");
    }
}
void measures() {
    // Nonnormal repeated-root Schur factorization: reconstruct the original
    // matrix, so discarding triangular couplings cannot pass this test.
    const std::vector<double> nonnormal{1,2,-.3,0,1,4,0,0,.5};
    const auto schur=dynsyslab::numerics::complexSchur(nonnormal,3);
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j){std::complex<double> a{},identity{};
        for(std::size_t k=0;k<3;++k){identity+=schur.unitary[i*3+k]*std::conj(schur.unitary[j*3+k]);
            for(std::size_t l=0;l<3;++l)a+=schur.unitary[i*3+k]*schur.triangular[k*3+l]*std::conj(schur.unitary[j*3+l]);}
        require(std::abs(a-nonnormal[i*3+j])<1e-12,"Full Schur reconstruction");
        require(std::abs(identity-double(i==j))<1e-12,"Unitary Schur basis");}
    c::KoopmanSpectralOptions options;options.method="mpedmd";options.acknowledgeMeasurePreserving=true;
    auto model=c::fitKoopmanSpectral(rotation(),linear(),options);
    std::vector<double> angles;for(std::size_t i=0;i<2048;++i)angles.push_back(-pi+2*pi*i/2048);
    c::riggedSpectralMeasure(model,{1,0},"cos",angles,.05,1);
    const auto& measure=model.spectral.measures.front();double mass=0;
    for(auto p:measure.density){require(p>=0,"First-order Poisson measure positivity");mass+=p*2*pi/angles.size();}
    require(std::abs(mass-.5)<1e-10,"Integrated spectral mass equals observable norm");
    const auto packets=model.spectral.wavePackets.front();
    auto restored=c::readSpectralKromArchive(c::makeSpectralKromArchive(model));
    require(restored.spectral.wavePackets.front().coefficients==packets.coefficients,"Wave packet archive preservation");
    c::riggedSpectralMeasure(model,{1,0},"cos-high",angles,.05,4);
    mass=0;for(auto p:model.spectral.measures.back().density)mass+=p*2*pi/angles.size();
    require(std::abs(mass-.5)<1e-8,"Higher-order resolvent mass");
    const dynsyslab::math::RationalCircleKernel kernel(.05,4);
    for(std::size_t j=0;j<angles.size();j+=47){
        require(std::abs(model.spectral.measures.back().density[j]-.25*(kernel(angles[j]-.31)+kernel(angles[j]+.31)))<2e-9,"Analytic two-atom density");
        const auto& coeff=model.spectral.wavePackets.back().coefficients;
        const double phase=.72;const auto value=std::cos(phase)*coeff[2*j]+std::sin(phase)*coeff[2*j+1];
        const auto exact=.5*(std::exp(std::complex<double>(0,phase))*kernel(angles[j]-.31)+std::exp(std::complex<double>(0,-phase))*kernel(angles[j]+.31));
        require(std::abs(value-exact)<2e-9,"Analytic rotation packet");
    }
    auto fixed=rotation();fixed.y=fixed.x;auto repeated=c::fitKoopmanSpectral(fixed,linear(),options);
    c::riggedSpectralMeasure(repeated,{1,0},"invariant",{0.,.2},.05,4);
    require(std::abs(repeated.spectral.measures.back().density[0]-.5*kernel(0))<2e-9,"Repeated invariant eigenvalue");
    std::vector<std::complex<double>> signal;for(std::size_t i=0;i<400;++i)signal.push_back(std::exp(std::complex<double>(0,.4*i)));
    auto correlation=c::correlationSpectralMeasure(signal,100,angles,true);mass=0;
    for(auto p:correlation.density){require(p>-1e-10,"Fejer positivity");mass+=p*2*pi/angles.size();}
    require(std::abs(mass-1)<1e-10,"Fejer mass");rejects([&]{c::correlationSpectralMeasure(signal,10,angles,false);});
}
void periodicComparison() {
    namespace m=dynsyslab::math;
    m::FourierSeries cosine{1,{{-1,.5},{1,.5}}},sine{1,{{-1,{0,.5}},{1,{0,-.5}}}};
    std::vector<m::PeriodicSnapshotOrbit> training{{{cosine,sine},2}},validation=training;
    for(auto&channel:validation[0].channels)for(auto&[n,a]:channel.coefficients)a*=2;
    c::KoopmanSpectralOptions o;o.method="mpedmd";o.acknowledgeMeasurePreserving=true;o.sampleInterval=.2;
    const std::vector<double> offsets{-.2,0,.2};const std::vector<std::size_t> channels{1,0};
    auto model=c::fitPeriodicDelaySpectrum(training,channels,offsets,32,o);
    std::vector<double> angles;for(unsigned j=0;j<128;++j)angles.push_back(-pi+2*pi*j/128);
    const auto reference=c::periodicSpectralReference(training,validation,0,.2,angles,{.2,.5},.08,4);
    const auto result=c::comparePeriodicDelaySpectrum(model,channels,offsets,validation,reference);
    require(result.rank==2&&result.features==6,"Redundant delay dictionary rank");
    require(result.densityErrorNorm<1e-9,"Exact invariant periodic density");
    for(auto error:result.complexErrorNorm)require(error<1e-9,"Held-amplitude complex packet error");
    const m::FourierSeries fast{1,{{-19,.5},{19,.5}}};
    const std::vector<m::PeriodicSnapshotOrbit> aliased{{{fast},1}};
    const auto audit=m::auditPeriodicQuadrature(aliased,{19,40});
    require(audit.aliasFreePhaseCount==39&&std::abs(audit.exactMass[0]-.5)<1e-12,"Finite polynomial resolution rule");
    require(std::abs(audit.sampledMass[0]-1)<1e-12&&std::abs(audit.sampledMass[1]-.5)<1e-12,"Quadrature aliasing detected");
    rejects([&]{c::fitPeriodicDelaySpectrum(aliased,{0},offsets,32,o);});
    rejects([&]{c::comparePeriodicDelaySpectrum(model,channels,{-1,-.2,.2},validation,reference);});
}
void existingModelAttachments() {
    c::PolynomialDictionary dictionary;dictionary.coordinateDimension=1;dictionary.exponents={1};
    c::FunctionDomain coordinates;coordinates.dimension=1;coordinates.names={"x"};
    c::KoopmanModel model(c::EigenfunctionFitResult{c::FunctionExpansion(dictionary, {{1,0}}, coordinates), c::Spectrum::fromContinuousRates({{-1,0}}), {}});
    // Eigenpair, Spectrum and Pseudospectrum round trip through the krom archive (Spectrum.hpp).
    c::SpectralData data;data.residuals.eigenpair={.01};data.residuals.definition=c::residual_definition::oneStep;
    data.residuals.scope=c::residual_scope::empiricalTraining;
    data.residuals.pseudospectrum={data.residuals.definition,data.residuals.scope,{{0.5,0.1},{-1,0}},{0.3,0.02},{{1,0},{0,1}},1};
    c::SpectralMeasure measure;measure.name="x";measure.method="test-atom";measure.atomLocations={0};measure.atomMasses={2};measure.totalMass=2;
    data.measures.push_back(measure);data.invariance.projectionDefect=0.25;
    auto evidence=model.spectrum();c::attachSpectralData(evidence,data);
    model.setSpectralEvidence(evidence);model.setSpectralDiagnostics(data);
    const auto path=std::filesystem::temp_directory_path()/("dynsyslab-spectral-attachment-"+
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    c::saveKoopmanModel(path,model);
    const auto loaded=c::loadKoopmanModel(path);
    std::filesystem::remove(path);
    const auto pair=loaded.eigenpair(0);
    require(pair.residual&&pair.residual->value==.01&&pair.residual->definition=="one-step"&&pair.residual->scope=="empirical-training"
            &&pair.eigenfunction&&pair.eigenvalue.continuous==std::complex<double>(-1,0),"Eigenpair round trip");
    const auto& pseudo=loaded.spectrum().pseudospectrum();
    require(pseudo&&pseudo->grid.size()==2&&pseudo->minimum[1]==0.02&&pseudo->minimizingCoefficients.size()==2&&pseudo->levelSet(0.1).size()==1,
            "Pseudospectrum round trip");
    require(loaded.spectrum().continuous().front().atomMasses.front()==2,"Continuous spectrum round trip");
    require(loaded.spectralDiagnostics().invariance.projectionDefect==0.25,"Spectral diagnostics round trip");
    const auto subset=loaded.extract({0});
    require(subset.spectrum().residual(0)->value==.01&&subset.spectrum().pseudospectrum()&&!subset.spectrum().continuous().empty(),
            "Extraction keeps the selected pairs' residuals and the operator's evidence");
    require(subset.spectralDiagnostics().invariance.projectionDefect==0,"Extraction clears diagnostics of the parent trial space");
    {   // A pair with no residual stays empty.
        c::KoopmanModel plain(c::EigenfunctionFitResult{c::FunctionExpansion(dictionary, {{1,0}}, coordinates), c::Spectrum::fromContinuousRates({{-2,0}}), {}});
        const auto other=std::filesystem::temp_directory_path()/("dynsyslab-plain-krom-"+
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
        c::saveKoopmanModel(other,plain);
        const auto archive=c::readNpz(other);
        const auto reread=c::loadKoopmanModel(other);std::filesystem::remove(other);
        require(!archive.contains("eigenpair_residual")&&!reread.eigenpair(0).residual&&!reread.spectrum().pseudospectrum(),
                "A method with no residual publishes pairs with empty residuals");
    }
    auto observable=c::Observable::fromExpressions({"x","x","","measured",{},false},
        {c::AnalyticExpression::state(1,0)});
    c::ObservableDecompositionReference ref{"test.npz","identity","x","spectral-measure",0};
    auto associated=observable.withDecomposition(ref);
    require(c::Observable::fromExpressionDefinition(associated.expressionDefinition()).descriptor().decompositions.size()==1,"Observable association serialization");
    require(associated.partialDerivative(0).descriptor().decompositions.empty(),"Derivative cannot inherit a different function's decomposition");
}
void stochasticAndRkhs() {
    c::KoopmanSnapshotPairs data;data.samples=4;data.dimension=1;data.x={1,2,3,4};data.y={.8,1.6,2.4,3.2};
    data.replicateCount=2;data.successorReplicates={.9,1.8,2.7,3.6,.7,1.4,2.1,2.8};
    c::SpectralDictionary d;d.dimension=1;d.constant=false;c::KoopmanSpectralOptions o;o.method="stochastic_resdmd";
    auto model=c::fitKoopmanSpectral(data,d,o);
    require(std::abs(model.multipliers[0]-.8)<1e-10,"Stochastic conditional mean multiplier");
    const auto& s=model.spectral.stochastic;
    require(s.sampleMean[0]<1e-10&&std::abs(s.total[0]-.1)<1e-10&&std::abs(s.innovation[0]-std::sqrt(.02))<1e-10,"Innovation versus conditional mean residual");
    data.y=data.x;data.weights.clear();auto rkhs=c::rkhsAdjointSpectrum(data,1,0,1e-10,{{1,0}});
    for(auto z:rkhs.eigenvalues)require(std::abs(z-std::complex<double>(1))<1e-9,"RKHS adjoint identity spectrum");
    for(double r:rkhs.residuals.eigenpair)require(r<1e-6,"RKHS exact identity residual");
    c::SpectralData associated;associated.rkhsAdjoints.push_back(rkhs);
    require(c::spectralDataFromJson(c::spectralDataJson(associated)).rkhsAdjoints[0].nodes==data.x,"RKHS typed serialization");
    c::SpectralDictionaryOptions polynomial;polynomial.kind="polynomial";polynomial.degree=3;
    auto dictionary=c::buildSpectralDictionary(rotation().x,2,polynomial);require(dictionary.size()==10,"Automatic total-degree dictionary");
}
void dictionariesAndStructure() {
    auto data=rotation(.9);c::KoopmanSpectralOptions options;
    auto model=c::fitKoopmanSpectral(data,linear(),options);
    for(auto z:model.multipliers)require(std::abs(std::abs(z)-.9)<1e-10,"Damped spectrum");
    options.method="mpedmd";rejects([&]{c::fitKoopmanSpectral(data,linear(),options);});
    options.method="hermitian_dmd";options.acknowledgeSelfAdjoint=true;
    for(std::size_t i=0;i<data.samples;++i){data.y[2*i]=.9*data.x[2*i];data.y[2*i+1]=.7*data.x[2*i+1];}
    auto hermitian=c::fitKoopmanSpectral(data,linear(),options);
    std::vector<double> eig;for(auto z:hermitian.multipliers)eig.push_back(z.real());std::sort(eig.begin(),eig.end());
    require(std::abs(eig[0]-.7)<1e-10&&std::abs(eig[1]-.9)<1e-10,"Self-adjoint spectrum");
    c::SpectralDictionary d;d.dimension=1;d.kind="partition";d.nodes={0,1,2};
    c::KoopmanSnapshotPairs cycle;cycle.dimension=1;cycle.samples=6;cycle.x={0,1,2,0,1,2};cycle.y={1,2,0,1,2,0};
    options.method="multdmd";auto mult=c::fitKoopmanSpectral(cycle,d,options);
    for(auto z:mult.multipliers)require(std::abs(std::pow(z,3)-std::complex<double>(1))<1e-9,"Multiplicative 3-cycle spectrum");
    c::fitSpectralReconstruction(mult,"index",cycle.x,cycle.x,1);
    auto rollout=c::rolloutSpectralKrom(mult,"index",{0},6);
    for(std::size_t i=0;i<rollout.size();++i)require(std::abs(rollout[i]-double(i%3))<1e-9,"Partition cycle forecast");
    d.kind="gaussian";d.bandwidth=.5;auto basis=d.evaluate({0,1});
    require(std::abs(basis[0]-1)<1e-14&&std::abs(basis[1]-std::exp(-2.))<1e-14,"Frozen Gaussian sections");
    d.kind="polynomial";d.exponents={0,1,2};auto p=d.evaluate({2});require(p==std::vector<double>({1,2,4}),"Polynomial dictionary");
    auto page=c::koopmanPageRank({0,1,1,0},2,{1,0},.8);
    require(std::abs(page[0]-1/1.8)<1e-10,"Personalized PageRank orientation");
}
}
int main(){try{rotationAndRollout();measures();periodicComparison();existingModelAttachments();stochasticAndRkhs();dictionariesAndStructure();std::cout<<"spectral Koopman tests passed\n";return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
