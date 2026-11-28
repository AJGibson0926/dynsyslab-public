#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/KromSpectralObservables.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <filesystem>
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
std::vector<double> residualsOf(const c::Spectrum& s){std::vector<double> out;for(std::size_t j=0;j<s.size();++j)out.push_back(s.residual(j)->value);return out;}
const c::SpectralMeasure& measureOf(const c::KoopmanModel& m,const std::string& band,const std::string& observable="cos"){
    const auto* found=m.spectrum().measure(band,observable);require(found!=nullptr,"measure present");return *found;}
const std::vector<std::complex<double>>& packetsOf(const c::KoopmanModel& m,const std::string& band,const std::string& observable="cos"){
    const auto& f=measureOf(m,band,observable).packets;require(f!=nullptr,"generalized eigenfunctions present");return f->coefficients();}
std::filesystem::path scratch(const std::string& stem){return std::filesystem::temp_directory_path()/(stem+"-"+
    std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");}
void rotationAndRollout() {
    for(const std::string method:{"resdmd","edmd","mpedmd","tlsdmd"}) {
        auto data=rotation();c::KoopmanSpectralOptions options;options.method=method;
        options.acknowledgeMeasurePreserving=true;options.rank=2;options.sampleInterval=.1;
        options.residualGrid={std::exp(std::complex<double>(0,.31)),{0,0}};
        auto model=c::fitKoopmanSpectral(data,linear(),options);
        require(model.compression()->rank==2&&model.liftedDimension()==2,"Rotation rank");
        for(auto z:model.spectrum().discreteMultipliers())require(std::abs(std::abs(z)-1)<1e-10,"Rotation eigenvalue radius");
        for(auto res:residualsOf(model.spectrum()))require(res<1e-10,"Analytic rotation eigenpair residual");
        const auto& pseudo=*model.spectrum().pseudospectrum();
        require(pseudo.minimum[0]<1e-7&&std::abs(pseudo.minimum[1]-1)<1e-10,"Pseudospectral distances for normal rotation");
        require(pseudo.minimizers&&pseudo.minimizers->outputCount()==2
                &&model.spectralFunction(c::pseudospectrumFunctionKey())==pseudo.minimizers.get(),"The pseudospectrum holds its minimizers");
        require(model.spectralDiagnostics().invariance.projectionDefect<1e-10,"Invariant coordinate span");
        c::fitFeatureReadout(model,"position",data.x,data.x,2);
        auto rollout=c::rolloutCompression(*model.compression(),"position",{1,0},40);
        for(std::size_t k=0;k<=40;++k)require(std::abs(rollout[2*k]-std::cos(.31*k))<1e-9&&std::abs(rollout[2*k+1]-std::sin(.31*k))<1e-9,"Autonomous rotation forecast");
        const auto path=scratch("dynsyslab-spectral-krom");c::saveKoopmanModel(path,model);
        const auto restored=c::loadKoopmanModel(path);std::filesystem::remove(path);
        require(c::kromObservableIdentity(model)==c::kromObservableIdentity(restored)
                &&restored.compression()->readouts.count("position"),"The krom keeps its compression through the archive");
        c::ObservableDecompositionReference ref;ref.modelIdentity=c::kromObservableIdentity(model);ref.component="position";
        auto observable=c::kromSpectralObservable(model,ref,"position");
        double point[]{.3,.4};auto values=observable.evaluate(point);
        require(std::abs(values[0]-.3)<1e-10&&std::abs(values[1]-.4)<1e-10,"Bound observable readout");
        ref.modelIdentity="stale";rejects([&]{c::kromSpectralObservable(model,ref,"bad");});
        auto held=rotation(.8);auto residual=c::koopmanResiduals(model,held);
        for(auto value:residualsOf(residual))require(std::abs(value-.2)<1e-10,"Held-out residual detects changed dynamics");
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
    c::riggedSpectralMeasure(model,{1,0},"cos",0,"kernel-1",angles,.05,1);
    const auto measure=measureOf(model,"kernel-1");double mass=0;
    require(model.spectrum().bands().size()==1&&model.spectrum().band("kernel-1")->grid==angles
            &&model.spectrum().bands().front().continuity.empty(),"One spectral band on the grid");
    for(auto p:measure.density){require(p>=0,"First-order Poisson measure positivity");mass+=p*2*pi/angles.size();}
    require(std::abs(mass-.5)<1e-10,"Integrated spectral mass equals observable norm");
    const auto packets=packetsOf(model,"kernel-1");
    require(packets.size()==angles.size()*2 && measure.packets->outputCount()==angles.size()
            && model.spectralFunction(c::measureFunctionKey(measure))==measure.packets.get(),
            "The measure holds its generalized eigenfunctions, one per angle");
    {
        const auto path=scratch("dynsyslab-rigged-krom");c::saveKoopmanModel(path,model);
        const auto restored=c::loadKoopmanModel(path);std::filesystem::remove(path);
        require(packetsOf(restored,"kernel-1")==packets&&measureOf(restored,"kernel-1").density==measure.density
                &&restored.spectrum().koopmanModes("cos"),
                "Measures and their generalized eigenfunctions survive the archive");
    }
    auto shifted=model;{auto c=*shifted.compression();c.action[1]*=1.01;shifted.setCompression(c);}  // no longer orthogonal
    rejects([&]{c::riggedSpectralMeasure(shifted,{1,0},"cos",0,"not-unitary",angles,.05,1);});
    rejects([&]{c::riggedSpectralMeasure(model,{1,0},"cos",0,"kernel-1",angles,.05,1);});  // one measure per observable and band
    rejects([&]{c::riggedSpectralMeasure(model,{1,0},"sin",0,"kernel-1",{0.},.05,1);});    // the band's grid is fixed
    // The compression's measure is atomic: one atom per eigenpair, masses |c_k|^2 summing to ||g||^2,
    // derived from the observable's Koopman modes (the rigged stage stores those, not the atoms).
    const auto atoms=model.spectrum().atoms("cos");
    require(atoms.size()==2 && std::abs(atoms[0].mass+atoms[1].mass-measure.totalMass)<1e-12
            && std::abs(atoms[0].mass-.25)<1e-10 && std::abs(atoms[1].mass-.25)<1e-10 && atoms[0].mode!=atoms[1].mode,
            "Atoms reference the two eigenpairs with the observable's masses");
    {   // Views: a pair lists the observables projected on it; an observable its modes and measures.
        const auto pair=model.eigenpair(0);
        require(pair.modes.size()==1&&pair.modes[0].observable=="cos"&&std::norm(pair.modes[0].value)==atoms[0].mass,
                "The pair lists the observable's Koopman mode on it");
        const auto view=model.spectrum().decomposition("cos");
        require(view.modes&&view.measures.size()==1&&view.measures[0]->band=="kernel-1"
                &&model.spectrum().bandMeasures("kernel-1").size()==1,"The observable's decomposition");
    }
    c::riggedSpectralMeasure(model,{1,0},"cos",0,"kernel-4",angles,.05,4);
    require(model.spectrum().decomposition("cos").measures.size()==2,"One observable, two bands");
    mass=0;for(auto p:measureOf(model,"kernel-4").density)mass+=p*2*pi/angles.size();
    require(std::abs(mass-.5)<1e-8,"Higher-order resolvent mass");
    const dynsyslab::math::RationalCircleKernel kernel(.05,4);
    for(std::size_t j=0;j<angles.size();j+=47){
        require(std::abs(measureOf(model,"kernel-4").density[j]-.25*(kernel(angles[j]-.31)+kernel(angles[j]+.31)))<2e-9,"Analytic two-atom density");
        const auto& coeff=packetsOf(model,"kernel-4");
        const double phase=.72;const auto value=std::cos(phase)*coeff[2*j]+std::sin(phase)*coeff[2*j+1];
        const auto exact=.5*(std::exp(std::complex<double>(0,phase))*kernel(angles[j]-.31)+std::exp(std::complex<double>(0,-phase))*kernel(angles[j]+.31));
        require(std::abs(value-exact)<2e-9,"Analytic rotation packet");
    }
    auto fixed=rotation();fixed.y=fixed.x;auto repeated=c::fitKoopmanSpectral(fixed,linear(),options);
    c::riggedSpectralMeasure(repeated,{1,0},"invariant",0,"kernel-4",{0.,.2},.05,4);
    require(std::abs(measureOf(repeated,"kernel-4","invariant").density[0]-.5*kernel(0))<2e-9,"Repeated invariant eigenvalue");
    std::vector<std::complex<double>> signal;for(std::size_t i=0;i<400;++i)signal.push_back(std::exp(std::complex<double>(0,.4*i)));
    auto correlation=c::correlationSpectralMeasure(signal,100,angles,true);mass=0;
    for(auto p:correlation.measure.density){require(p>-1e-10,"Fejer positivity");mass+=p*2*pi/angles.size();}
    require(correlation.band.grid==angles&&correlation.measure.band==correlation.band.name,"Fejer band");
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
    auto evidence=model.spectrum();
    evidence.setResiduals({.01},c::residual_definition::oneStep,c::residual_scope::empiricalTraining);
    evidence.setPseudospectrum(c::Pseudospectrum{c::residual_definition::oneStep,c::residual_scope::empiricalTraining,{{0.5,0.1},{-1,0}},{0.3,0.02}});
    c::SpectralBand band;band.name="arc";band.method="test";band.grid={0.,.5};
    evidence.setBands({band});
    c::SpectralMeasure measure;measure.observable="x";measure.band="arc";measure.density={.1,.2};measure.totalMass=2;
    evidence.addMeasure(measure);
    model.setSpectralEvidence(evidence);
    c::ReconstructionResult modes;modes.outputs=1;modes.modes=1;modes.matrix={std::sqrt(2.)};
    model.setReconstruction("x",modes);  // the observable's Koopman mode: its atom has mass 2
    c::SpectralDiagnostics diagnostics;diagnostics.invariance.projectionDefect=0.25;model.setSpectralDiagnostics(diagnostics);
    const auto path=std::filesystem::temp_directory_path()/("dynsyslab-spectral-attachment-"+
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    c::saveKoopmanModel(path,model);
    const auto loaded=c::loadKoopmanModel(path);
    std::filesystem::remove(path);
    const auto pair=loaded.eigenpair(0);
    require(pair.residual&&pair.residual->value==.01&&pair.residual->definition=="one-step"&&pair.residual->scope=="empirical-training"
            &&pair.eigenfunction&&pair.eigenvalue.continuous==std::complex<double>(-1,0),"Eigenpair round trip");
    const auto& pseudo=loaded.spectrum().pseudospectrum();
    require(pseudo&&pseudo->grid.size()==2&&pseudo->minimum[1]==0.02&&pseudo->levelSet(0.1).size()==1,
            "Pseudospectrum round trip");
    require(loaded.spectrum().band("arc")&&loaded.spectrum().measure("arc","x")->density[1]==.2
            &&std::abs(loaded.spectrum().atoms("x").front().mass-2)<1e-12,"Continuous spectrum round trip");
    require(loaded.spectralDiagnostics().invariance.projectionDefect==0.25,"Spectral diagnostics round trip");
    const auto subset=loaded.extract({0});
    require(subset.spectrum().residual(0)->value==.01&&subset.spectrum().pseudospectrum()&&!subset.spectrum().bands().empty(),
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
    require(std::abs(model.spectrum().discreteMultipliers()[0]-.8)<1e-10,"Stochastic conditional mean multiplier");
    const auto& s=model.spectralDiagnostics().stochastic;
    require(s.sampleMean[0]<1e-10&&std::abs(s.total[0]-.1)<1e-10&&std::abs(s.innovation[0]-std::sqrt(.02))<1e-10,"Innovation versus conditional mean residual");
    data.y=data.x;data.weights.clear();auto rkhs=c::rkhsAdjointSpectrum(data,1,0,1e-10,{{1,0}});
    for(auto z:rkhs.eigenvalues)require(std::abs(z-std::complex<double>(1))<1e-9,"RKHS adjoint identity spectrum");
    for(double r:rkhs.residuals)require(r<1e-6,"RKHS exact identity residual");
    c::SpectralDiagnostics associated;associated.rkhsAdjoints.push_back(rkhs);
    require(c::spectralDiagnosticsFromJson(c::spectralDiagnosticsJson(associated)).rkhsAdjoints[0].nodes==data.x,"RKHS typed serialization");
    c::SpectralDictionaryOptions polynomial;polynomial.kind="polynomial";polynomial.degree=3;
    auto dictionary=c::buildSpectralDictionary(rotation().x,2,polynomial);require(dictionary.size()==10,"Automatic total-degree dictionary");
}
void dictionariesAndStructure() {
    auto data=rotation(.9);c::KoopmanSpectralOptions options;
    auto model=c::fitKoopmanSpectral(data,linear(),options);
    for(auto z:model.spectrum().discreteMultipliers())require(std::abs(std::abs(z)-.9)<1e-10,"Damped spectrum");
    options.method="mpedmd";rejects([&]{c::fitKoopmanSpectral(data,linear(),options);});
    options.method="hermitian_dmd";options.acknowledgeSelfAdjoint=true;
    for(std::size_t i=0;i<data.samples;++i){data.y[2*i]=.9*data.x[2*i];data.y[2*i+1]=.7*data.x[2*i+1];}
    auto hermitian=c::fitKoopmanSpectral(data,linear(),options);
    std::vector<double> eig;for(auto z:hermitian.spectrum().discreteMultipliers())eig.push_back(z.real());std::sort(eig.begin(),eig.end());
    require(std::abs(eig[0]-.7)<1e-10&&std::abs(eig[1]-.9)<1e-10,"Self-adjoint spectrum");
    c::SpectralDictionary d;d.dimension=1;d.kind="partition";d.nodes={0,1,2};
    c::KoopmanSnapshotPairs cycle;cycle.dimension=1;cycle.samples=6;cycle.x={0,1,2,0,1,2};cycle.y={1,2,0,1,2,0};
    options.method="multdmd";auto mult=c::fitKoopmanSpectral(cycle,d,options);
    for(auto z:mult.spectrum().discreteMultipliers())require(std::abs(std::pow(z,3)-std::complex<double>(1))<1e-9,"Multiplicative 3-cycle spectrum");
    c::fitFeatureReadout(mult,"index",cycle.x,cycle.x,1);
    auto rollout=c::rolloutCompression(*mult.compression(),"index",{0},6);
    for(std::size_t i=0;i<rollout.size();++i)require(std::abs(rollout[i]-double(i%3))<1e-9,"Partition cycle forecast");
    d.kind="gaussian";d.bandwidth=.5;auto basis=d.evaluate({0,1});
    require(std::abs(basis[0]-1)<1e-14&&std::abs(basis[1]-std::exp(-2.))<1e-14,"Frozen Gaussian sections");
    d.kind="polynomial";d.exponents={0,1,2};auto p=d.evaluate({2});require(p==std::vector<double>({1,2,4}),"Polynomial dictionary");
    auto page=c::koopmanPageRank({0,1,1,0},2,{1,0},.8);
    require(std::abs(page[0]-1/1.8)<1e-10,"Personalized PageRank orientation");
}
}
// A library krom: the pendulum's energy as an eigenpair and its libration phase as a band whose
// generalized eigenfunctions are on level sets (profile e^{i phi}, level map Omega(H)); both functions
// survive the archive and a mode selection.
void libraryBands() {
    c::SystemConfiguration config;config.identifier="nonlinear_pendulum";config.parameters={{"damping_n_m_s",0.0}};
    const auto system=c::createSystem(config);
    auto bank=c::libraryEigenfunctionBank(*system,{"energy","libration_phase"});
    const c::KoopmanModel model(bank,config);
    require(model.liftedDimension()==1&&model.spectrum().bands().size()==1,"one pair and one band");
    const auto path=scratch("dynsyslab-library-bands");c::saveKoopmanModel(path,model);
    const auto loaded=c::loadKoopmanModel(path);std::filesystem::remove(path);
    const auto* band=loaded.spectrum().band("libration_phase");
    require(band&&band->eigenfunctions.levelSet()&&band->support==c::band_support::imaginaryAxis,"the band and its level-set form");
    const double x[]{.8,.1};
    const auto* original=model.spectrum().band("libration_phase");
    require(std::abs(band->eigenfunctions.functions->evaluate(x).front()-original->eigenfunctions.functions->evaluate(x).front())<1e-12
            &&std::abs(band->eigenfunctions.level->evaluate(x).front()-original->eigenfunctions.level->evaluate(x).front())<1e-12,
            "profile and level map round trip");
    require(loaded.extract({0}).spectrum().band("libration_phase")->eigenfunctions.levelSet(),"selection keeps the band");
    {   // The band's generalized eigenfunctions written out: delta(Omega(H(x)) - p) e^{i phi(x)}, 0 off the
        // level set p = Omega(H(x)) and without a value on it; it survives the archive.
        const auto& family=band->eigenfunctions.family;
        const double level=band->eigenfunctions.level->evaluate(x).front().real();
        const double off[]{x[0],x[1],level+0.1},on[]{x[0],x[1],level};
        bool refused=false;  // a distribution has no pointwise value on its support: the evaluation refuses
        try{(void)family->evaluate(on);}catch(const std::invalid_argument&){refused=true;}
        require(family&&family->coordinateDimension()==3&&family->evaluate(off).front()==std::complex<double>(0,0)&&refused,
                "the explicit family is a distribution on the level sets");
    }
    // The inviscid RPE's known spectrum from its catalogue: the energy as an eigenpair (eigenvalue 0)
    // and the phase as a band whose expressions (shared definitions, elliptic integrals) round trip.
    c::SystemConfiguration rpe;rpe.identifier="rpe";rpe.parameters={{"mu_pa_s",0.0}};
    const auto bubble=c::createSystem(rpe);
    const c::KoopmanModel known(c::libraryEigenfunctionBank(*bubble),bubble->configuration());
    const auto rpePath=scratch("dynsyslab-rpe-known-spectrum");c::saveKoopmanModel(rpePath,known);
    const auto back=c::loadKoopmanModel(rpePath);std::filesystem::remove(rpePath);
    const auto* phase=back.spectrum().band("phase");
    const double y[]{1.6,-0.7};
    require(back.liftedDimension()==1&&phase&&phase->eigenfunctions.levelSet()
            &&std::abs(phase->eigenfunctions.functions->evaluate(y).front()-known.spectrum().band("phase")->eigenfunctions.functions->evaluate(y).front())<1e-13
            &&std::abs(back.evaluate(y).front()-known.evaluate(y).front())<1e-12,"the RPE's known spectrum as a krom");
}
int main(){try{rotationAndRollout();measures();periodicComparison();existingModelAttachments();libraryBands();stochasticAndRkhs();dictionariesAndStructure();std::cout<<"spectral Koopman tests passed\n";return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
