// Case-only data selection and native API composition. No fitting algorithms,
// ODE derivatives, unforced trajectories, or previously fitted eigenfunctions.
#include <dynsyslab/core/EdmdControlArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
namespace c=dynsyslab::core;
using J=c::json::Value;using Path=std::filesystem::path;
constexpr double dt=.01,pi=3.14159265358979323846;
J zjson(std::complex<double> z){return J::Array{z.real(),z.imag()};}
void writeJson(const Path& p,const J& value){
    if(std::filesystem::exists(p))throw std::runtime_error("Refusing overwrite");
    std::ofstream f(p);f<<c::json::dump(value)<<'\n';if(!f)throw std::runtime_error("write failed");
}
struct Pairs{std::vector<double>x,y,u,time,offsets;J::Array sources;};
Pairs readBank(const Path& folder,const std::string& prefix,int count){
    Pairs out;
    for(int i=1;i<=count;++i){
        const auto stem=prefix+std::to_string(i);
        const auto measurement=folder/(stem+".measurements.npz"),truth=folder/(stem+".truth.npz");
        const auto a=c::readNpz(measurement),b=c::readNpz(truth);const auto meta=c::json::parse(a.metadataJson());
        const auto t=a.at("time").doubles(),s=b.at("states").doubles(),u=a.at("input").doubles();
        if(t!=b.at("time").doubles()||s.size()!=2*t.size()||u.size()!=t.size())throw std::runtime_error("Unaligned source records");
        const auto stride=std::size_t(std::llround(dt/meta.at("sensor_step").number()));
        if(!stride||std::abs(meta.at("control_step").number()-dt)>1e-12)throw std::runtime_error("Input tick mismatch");
        const double start=meta.at("drive_start").number(),stop=meta.at("drive_stop").number();
        const auto initial=out.u.size();out.offsets.push_back(double(initial));
        for(std::size_t k=0;k+stride<t.size();k+=stride){
            if(t[k]<start-1e-10||t[k+stride]>stop+1e-10)continue;
            if(std::abs(t[k+stride]-t[k]-dt)>1e-10)throw std::runtime_error("Nonuniform pairs");
            // Verify each actual transition sees precisely the recorded held input.
            for(std::size_t h=k;h<k+stride;++h)
                if(std::abs(u[h]-u[k])>1e-12)throw std::runtime_error("Input changes inside fitted interval");
            out.x.insert(out.x.end(),{s[2*k],s[2*k+1]});out.y.insert(out.y.end(),{s[2*(k+stride)],s[2*(k+stride)+1]});
            out.u.push_back(u[k]);out.time.push_back(t[k]-start);
        }
        if(out.u.size()==initial)throw std::runtime_error("No active forcing interval");
        out.sources.push_back(J::Object{{"id",stem},{"states",truth.generic_string()},
            {"inputs",measurement.generic_string()},{"drive_start",start},{"drive_stop",stop},
            {"pairs",double(out.u.size()-initial)},{"amplitude",meta.at("amplitude_model")},
            {"frequency",meta.at("frequency_cycles_per_tau")},{"phase",meta.at("phase")}});
        out.sources.back().object()["initial_radius"]=meta.at("initial_radius");
    }return out;
}
// Publish the existing EDMDc/v1 contract, using the same public archive types
// as workflow DmdControlPostprocess. Raw datasets stay at their source paths.
void saveFit(const Path& path,const c::EdmdControlResult& model,J meta){
    const auto& f=model.lifted;const auto& d=model.dictionary;c::NpzRecord a;
    auto real=[&](const char* k,std::vector<std::size_t> s,const std::vector<double>& v){a.arrays[k]=c::NpyArray::fromDoubles(s,v);};
    auto complex=[&](const char* k,std::vector<std::size_t> s,const std::vector<std::complex<double>>& v){a.arrays[k]=c::NpyArray::fromComplex(s,v);};
    real("operator",{f.rank,f.rank},f.operatorMatrix);real("control_matrix",{f.rank,f.inputDimension},f.controlMatrix);
    real("projection_basis",{f.dimension,f.rank},f.basis);
    real("augmented_singular_values",{f.augmentedSingularValues.size()},f.augmentedSingularValues);
    real("output_singular_values",{f.outputSingularValues.size()},f.outputSingularValues);
    complex("eigenvalues",{f.rank},f.discreteEigenvalues);complex("modes",{f.dimension,f.rank},f.modes);
    complex("exact_modes",{f.dimension,f.rank},f.exactModes);
    complex("eigenfunction_coefficients",{f.rank,f.dimension},f.eigenfunctionCoefficients);
    complex("modal_input_matrix",{f.rank,f.inputDimension},f.modalInputMatrix);
    a.arrays["polynomial_exponents"]=c::NpyArray::fromInt64s({f.dimension,d.dimension},{d.exponents.begin(),d.exponents.end()});
    real("state_center",{d.dimension},d.center);real("state_scale",{d.dimension},d.scale);
    auto& m=meta.object();m["model_schema"]="dynsyslab.edmdc/v1";m["format"]="dynsyslab.derived/v1";
    m["rank"]=double(f.rank);m["feature_count"]=double(f.dimension);m["snapshot_dimension"]=double(d.dimension);
    m["input_dimension"]=double(f.inputDimension);m["pairs"]=double(f.pairs);m["delta_t"]=f.deltaTime;
    m["polynomial_degree"]=double(d.degree);m["include_constant"]=d.includeConstant;m["normalize"]=true;
    m["eigenfunctions_available"]=f.hasEigenfunctions;m["full_row_rank"]=f.fullRowRank;
    m["augmented_rank"]=double(f.augmentedRank);m["augmented_numerical_rank"]=double(f.augmentedNumericalRank);
    m["one_step_relative_error"]=f.oneStepRelativeError;m["projection_relative_error"]=f.projectionRelativeError;
    m["retained_condition_number"]=f.retainedConditionNumber;m["eigenvector_condition_number"]=f.eigenvectorConditionNumber;
    m["snapshot_one_step_relative_error"]=model.snapshotOneStepRelativeError;
    a.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(meta));c::writeNpzAtomicNew(path,a);
}
J score(const c::EdmdControlResult& fit,const Pairs& p,std::size_t mode){
    const auto predicted=c::predictEdmdControl(fit,p.x,p.u,p.u.size());
    const auto a=c::evaluateEdmdControlEigenfunctions(fit,p.x,p.u.size());
    const auto b=c::evaluateEdmdControlEigenfunctions(fit,p.y,p.u.size());
    double error=0,energy=0,increment=0,pe=0,py=0,pdelta=0;
    for(std::size_t j=0;j<p.y.size();++j){error+=std::pow(predicted[j]-p.y[j],2);energy+=p.y[j]*p.y[j];increment+=std::pow(p.y[j]-p.x[j],2);}
    for(std::size_t k=0;k<p.u.size();++k){const auto x=a[k*fit.lifted.rank+mode],y=b[k*fit.lifted.rank+mode];
        const auto next=fit.lifted.discreteEigenvalues[mode]*x+fit.lifted.modalInputMatrix[mode]*p.u[k];
        pe+=std::norm(next-y);py+=std::norm(y);pdelta+=std::norm(y-x);}
    return J::Object{{"pairs",double(p.u.size())},{"state_one_step_relative_l2",std::sqrt(error/energy)},
        {"state_error_over_increment_l2",std::sqrt(error/increment)},
        {"selected_mode_forced_one_step_relative_l2",std::sqrt(pe/py)},
        {"selected_mode_error_over_increment_l2",std::sqrt(pe/pdelta)}};
}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::runtime_error("SINUSOIDAL_BANK OUTPUT_DIRECTORY");const Path folder=std::filesystem::absolute(argv[1]),out=std::filesystem::absolute(argv[2]);
    if(std::filesystem::exists(out/"results.json"))throw std::runtime_error("Completed experiment is immutable");std::filesystem::create_directories(out);
    const auto train=readBank(folder,"train-",60),test=readBank(folder,"heldout-",9);
    std::cout<<"Forced training pairs: "<<train.u.size()<<"; held-out forced pairs: "<<test.u.size()<<std::endl;
    c::EdmdControlOptions options;options.polynomialDegree=4;options.includeConstant=true;options.normalize=true;
    const auto fit=c::extendedDynamicModeDecompositionWithControl(train.x,train.y,train.u,train.u.size(),2,1,dt,options);
    if(!fit.lifted.hasEigenfunctions)throw std::runtime_error("Eigenbasis is not usable");
    // Select the positive-frequency oscillatory mode closest to one cycle/tau.
    // No unforced data, derivative, prior spectrum or control result enters this rule.
    std::size_t selected=fit.lifted.rank;double best=std::numeric_limits<double>::infinity();J::Array spectrum;
    for(std::size_t j=0;j<fit.lifted.rank;++j){const auto mu=fit.lifted.discreteEigenvalues[j];const auto rate=std::log(mu)/dt;
        spectrum.push_back(J::Object{{"index",double(j)},{"mu",zjson(mu)},{"rate_for_interpretation",zjson(rate)}});
        if(mu.imag()>1e-8&&std::abs(rate.imag()-2*pi)<best){best=std::abs(rate.imag()-2*pi);selected=j;}}
    if(selected==fit.lifted.rank)throw std::runtime_error("No complex oscillatory mode");
    J meta=J::Object{{"label","RPE forced-only sine-bank EDMDc degree 4"},{"method","edmdc"},
        {"eigenvalue_domain","discrete-time"},{"training_policy","Only active sinusoidal-forcing intervals; no warmup, ringdown, or unforced records"},
        {"training_sources",train.sources},{"validation_sources",test.sources},{"input_alignment","input[k] held throughout X[k] to Y[k]"},
        {"integration_form","dimensionless"},{"delta_t",dt},{"control_ready",false}};
    saveFit(out/"edmdc.npz",fit,meta);
    const auto reloaded=c::loadEdmdControlModel(out/"edmdc.npz");const auto poly=reloaded.polynomialEigenfunctions();
    // Physics metadata ONLY; it is not used in the EDMDc regression above.
    const auto provenance=c::json::load((folder/"provenance.json").wstring());
    auto config=provenance.at("source_metadata").at("system_configuration");
    // The provenance file originated in a different study: copy its physics,
    // not that study's unforced initial condition. Each sine-training run starts at r=1.
    config.at("initial_state")=J::Array{1.0,0.0};
    c::SystemConfiguration cfg;cfg.identifier=config.at("identifier").string();cfg.preset=config.at("preset").string();
    for(const auto& [k,v]:config.at("parameters").object())cfg.parameters[k]=v.number();
    const auto selectedModel=c::KoopmanModel(poly,cfg).extract({selected});
    auto selectedMeta=meta;selectedMeta.at("label")="RPE forced-only EDMDc polynomial eigenfunction";
    selectedMeta.object()["selected_mode"]=double(selected);
    if(std::filesystem::exists(out/"selected-polynomial.npz"))throw std::runtime_error("Polynomial archive exists");
    c::savePolynomialKoopmanModel(out/"selected-polynomial.npz",selectedModel,c::json::dump(selectedMeta),{
        {"training_states",c::NpyArray::fromDoubles({train.u.size(),2},train.x)},
        {"training_time",c::NpyArray::fromDoubles({train.time.size()},train.time)},
        {"training_trajectory_offsets",c::NpyArray::fromDoubles({train.offsets.size()},train.offsets)},
        {"source_mode_indices",c::NpyArray::fromDoubles({1},{double(selected)})}});
    double parity=0;const auto direct=c::evaluateEdmdControlEigenfunctions(fit,train.x,train.u.size());
    const auto loaded=c::loadKoopmanModel(out/"selected-polynomial.npz");
    for(std::size_t k=0;k<train.u.size();++k)parity=std::max(parity,std::abs(direct[k*fit.lifted.rank+selected]-loaded.evaluate(&train.x[2*k])[0]));
    if(parity>1e-10)throw std::runtime_error("Shared polynomial export differs");
    J::Array powers,coefficients,center,scale;
    for(std::size_t k=0;k<selectedModel.eigenfunctions().featureCount;++k){powers.push_back(J::Array{double(poly.exponents[2*k]),double(poly.exponents[2*k+1])});coefficients.push_back(zjson(selectedModel.eigenfunctions().coefficients[k]));}
    for(auto v:poly.center)center.emplace_back(v);for(auto v:poly.scale)scale.emplace_back(v);
    const J result=J::Object{{"schema","dynsyslab.forced-rpe-edmdc-case/v1"},{"configuration",config},
        {"metadata",meta},{"spectrum",spectrum},{"selected_index",double(selected)},
        {"mu",zjson(fit.lifted.discreteEigenvalues[selected])},{"modal_input",zjson(fit.lifted.modalInputMatrix[selected])},
        {"selection","Positive-frequency mode nearest one cycle per nondimensional time unit"},
        {"rank",double(fit.lifted.rank)},{"augmented_rank",double(fit.lifted.augmentedRank)},
        {"augmented_condition",fit.lifted.retainedConditionNumber},{"eigenvector_condition",fit.lifted.eigenvectorConditionNumber},
        {"full_row_rank",fit.lifted.fullRowRank},{"lifted_training_one_step_error",fit.lifted.oneStepRelativeError},
        {"training_score",score(fit,train,selected)},{"heldout_forced_score",score(fit,test,selected)},
        {"export_max_absolute_error",parity},{"exponents",powers},{"coefficients_real_imag",coefficients},
        {"center",center},{"scale",scale},{"normalization","Fit on forced X only; frozen for all evaluation"},
        {"limitations","Eigenfunctions of fitted A, not certified unforced Koopman eigenfunctions. Forced evolution includes modal_input*u. No control attempt or continuous B conversion."}};
    writeJson(out/"results.json",result);std::cout<<c::json::dump(result)<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
