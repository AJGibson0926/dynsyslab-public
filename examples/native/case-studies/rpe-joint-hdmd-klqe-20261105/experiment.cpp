// Experiment orchestration and candidate proposals only. Fitting, reconstruction,
// Gramians, subset assessments and KLQE all use the standard native APIs.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/EstimationSelection.hpp>
#include <dynsyslab/core/EstimationRealization.hpp>
#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/core/Json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
namespace c=dynsyslab::core;
using J=c::json::Value;using Path=std::filesystem::path;
void jsonFile(const Path& p,const J& j){std::ofstream f(p);f<<c::json::dump(j);if(!f)throw std::runtime_error("write failed");}
J finite(double x){return std::isfinite(x)?J(x):J{};}
J indices(const std::vector<std::size_t>& x){J::Array a;for(auto k:x)a.emplace_back(double(k));return a;}
std::vector<double> scales(const std::vector<double>& x,std::size_t p) {
    const auto N=x.size()/p;std::vector<double> result(p);
    for(std::size_t j=0;j<p;++j){double mean=0;for(std::size_t k=0;k<N;++k)mean+=x[k*p+j]/double(N);
        for(std::size_t k=0;k<N;++k)result[j]+=std::pow(x[k*p+j]-mean,2)/double(N);
        result[j]=std::max(1e-8,std::sqrt(result[j]));}return result;
}
J metrics(const std::vector<double>& actual,const std::vector<double>& predicted,std::size_t p) {
    J::Array rmse;for(std::size_t j=0;j<p;++j){double e=0;
        for(std::size_t k=0;k<actual.size()/p;++k)e+=std::pow(actual[k*p+j]-predicted[k*p+j],2);
        rmse.emplace_back(std::sqrt(e/double(actual.size()/p)));}return rmse;
}
J intervalMetrics(const std::vector<double>& actual,const std::vector<double>& predicted,
    std::size_t p,const std::vector<double>& time,double begin,double end) {
    std::vector<double> a,b;
    for(std::size_t k=0;k<time.size();++k)if(time[k]>=begin&&time[k]<=end) {
        a.insert(a.end(),actual.begin()+k*p,actual.begin()+(k+1)*p);
        b.insert(b.end(),predicted.begin()+k*p,predicted.begin()+(k+1)*p);
    }
    return metrics(a,b,p);
}
void inspect(const Path& path,const Path& out) {
    auto model=c::loadKoopmanModel(path);const auto record=c::readNpz(path);
    const auto* data=model.scatteredEigenfunctions();const auto N=data->samples;
    const auto system=c::createSystem(model.systemConfiguration());
    const auto echo=system->observableCatalog().at("far_field_acoustic_echo").evaluateRealBatch(data->trainingStates,N);
    const auto x=c::reconstructObservable(model.reconstruction("state"),data->trainingValues,N);
    const auto y=c::reconstructObservable(model.reconstruction("far_field_acoustic_echo"),data->trainingValues,N);
    std::vector<double> xr,yr;for(auto z:x)xr.push_back(z.real());for(auto z:y)yr.push_back(z.real());
    const auto singular=record.at("hankel_singular_values").doubles();
    std::size_t growing=0;double maxGrowth=-1e300;
    for(auto l:model.continuousEigenvalues()){if(l.real()>1e-6)++growing;maxGrowth=std::max(maxGrowth,l.real());}
    jsonFile(out/"inspection.json",J::Object{{"source_krom",path.generic_u8string()},
        {"mode_count",double(model.eigenfunctionCount())},{"growing_modes",double(growing)},
        {"maximum_growth_rate",maxGrowth},{"singular_first",singular.front()},
        {"singular_last",singular.back()},{"state_reconstruction_rmse",metrics(data->trainingStates,xr,2)},
        {"echo_reconstruction_rmse",metrics(echo,yr,1)}});
    c::NpzRecord view;
    view.arrays["time"]=record.at("training_time");view.arrays["true_states"]=c::NpyArray::fromDoubles({N,2},data->trainingStates);
    view.arrays["state_reconstruction"]=c::NpyArray::fromDoubles({N,2},xr);
    view.arrays["echo"]=c::NpyArray::fromDoubles({N,1},echo);
    view.arrays["echo_reconstruction"]=c::NpyArray::fromDoubles({N,1},yr);
    view.arrays["hankel_singular_values"]=record.at("hankel_singular_values");
    view.arrays["eigenvalues"]=record.at("continuous_eigenvalues");
    c::writeNpzAtomicNew(out/"full-reference.npz",view);
}
void modal(const Path& path,const Path& out) {
    const auto model=c::loadKoopmanModel(path);const auto* data=model.scatteredEigenfunctions();
    const auto n=model.eigenfunctionCount();std::vector<double> rms(n);
    for(std::size_t j=0;j<n;++j){for(std::size_t k=0;k<data->samples;++k)rms[j]+=std::norm(data->trainingValues[k*n+j]);
        rms[j]=std::max(1e-30,std::sqrt(rms[j]/double(data->samples)));}
    c::ObservabilityOptions options;options.method=c::ObservabilityMethod::continuousFiniteGramian;
    options.horizon=0.001; // whole bank includes severely unstable spurious modes
    const auto result=c::analyzeModalObservability(model,"far_field_acoustic_echo",options,rms);
    c::NpzRecord record;
    record.arrays["gramian"]=c::NpyArray::fromComplex({n,n},result.gramian);
    record.arrays["directions"]=c::NpyArray::fromComplex({n,n},result.directions);
    record.arrays["gramian_eigenvalues"]=c::NpyArray::fromDoubles({n},result.gramianEigenvalues);
    record.arrays["singular_values"]=c::NpyArray::fromDoubles({n},result.singularValues);
    record.arrays["coordinate_scales"]=c::NpyArray::fromDoubles({n},rms);
    record.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(J::Object{
        {"source_krom",path.generic_u8string()},{"horizon",options.horizon},{"rank",double(result.rank)},
        {"rank_threshold",result.rankTolerance},{"coordinate_space","complex modal; conjugate constraints not yet imposed"}}));
    c::writeNpzAtomicNew(out/"full-modal-gramian.npz",record);
    std::cout<<"Full modal Gramian rank "<<result.rank<<" / "<<n<<std::endl;
}
struct Proposal{std::size_t source;double importance,visibility,defect,joint;};
void fullReal(const Path& path,const Path& out) {
    const auto parent=c::loadKoopmanModel(path);std::vector<std::size_t> selected;
    for(std::size_t j=0;j<parent.eigenfunctionCount();++j) {
        const auto l=parent.continuousEigenvalues()[j];
        if(l.imag()>0||std::abs(l.imag())<=1e-7*std::max(1.0,std::abs(l)))selected.push_back(j);
    }
    const auto r=c::makeRealModalRealization(parent,selected,1e-3,1e-7);
    c::saveKoopmanModel(out/"candidate.npz",r.model,c::json::dump(J::Object{
        {"format","dynsyslab.derived/v1"},{"method","klqe-full-real"},
        {"source_krom",std::filesystem::absolute(path).generic_u8string()},
        {"source_modes",indices(selected)},{"coefficients_refitted",false},
        {"real_rate_tolerance",1e-7},{"sampled_conjugacy_tolerance",1e-3}}));
}
void assess(const Path& path,const Path& out) {
    const auto model=c::loadKoopmanModel(path);const auto meta=c::json::parse(c::readNpz(path).metadataJson());
    const auto parent=c::readNpz(Path(meta.at("source_krom").string()));
    const auto time=parent.at("training_time").doubles();const auto* data=model.scatteredEigenfunctions();
    const auto system=c::createSystem(model.systemConfiguration());
    const auto y=system->observableCatalog().at("far_field_acoustic_echo").evaluateRealBatch(data->trainingStates,time.size());
    c::EstimationSelectionOptions o;o.measurement="far_field_acoustic_echo";o.outputScales=scales(data->trainingStates,2);
    o.measurementStandardDeviation=0.01*scales(y,1)[0];o.horizon=2;o.sampleTime=0.0001;
    const auto a=c::assessEstimationSubset(model,time,data->trainingValues,data->trainingStates,y,o);
    jsonFile(out/"assessment.json",J::Object{{"state_mse",a.stateReconstructionMse},
        {"echo_mse",a.measurementReconstructionMse},{"dynamics_mse",a.autonomousPropagationMse},
        {"noise_mse",a.posteriorStateNoiseMse},{"score",a.score},{"rank",double(a.observability.rank)}});
    c::NpzRecord g;const auto n=model.splitDimension();
    g.arrays["gramian_eigenvalues"]=c::NpyArray::fromDoubles({n},a.observability.gramianEigenvalues);
    g.arrays["directions"]=c::NpyArray::fromDoubles({n,n},a.observability.rightDirections);
    c::writeNpzAtomicNew(out/"gramian.npz",g);
}
void screen(const Path& path,const Path& out) {
    auto parent=c::loadKoopmanModel(path);const auto record=c::readNpz(path);
    const auto* data=parent.scatteredEigenfunctions();const auto time=record.at("training_time").doubles();
    const auto x=data->trainingStates;const auto N=time.size();
    const auto system=c::createSystem(parent.systemConfiguration());
    const auto y=system->observableCatalog().at("far_field_acoustic_echo").evaluateRealBatch(x,N);
    c::EstimationSelectionOptions options;options.measurement="far_field_acoustic_echo";
    options.outputScales=scales(x,2);options.measurementStandardDeviation=0.01*scales(y,1)[0];
    options.horizon=2;options.sampleTime=0.0001;options.dynamicsWeight=1;
    std::vector<Proposal> pool;J::Array profiles;std::size_t rejected=0;
    for(std::size_t j=0;j<parent.eigenfunctionCount();++j) {
        const auto l=parent.continuousEigenvalues()[j];
        const bool nominalReal=std::abs(l.imag())<=1e-7*std::max(1.0,std::abs(l));
        if((l.imag()<0&&!nominalReal)||l.real()>0.01||l.real()<-50||std::abs(l.imag())>1000)continue;
        try{
            const auto r=c::makeRealModalRealization(parent,{j},1e-3,1e-7);const auto& model=r.model;
            const auto& psi=model.scatteredEigenfunctions()->trainingValues;
            double defect=0;
            for(std::size_t k=1;k<N;++k)defect+=std::norm(psi[k]-std::exp(l*(time[k]-time[k-1]))*psi[k-1]);
            defect=std::sqrt(defect/double(N-1));double importance=0;
            for(std::size_t row=0;row<2;++row)importance+=std::norm(model.reconstruction("state").matrix[row]/options.outputScales[row]);
            importance=std::sqrt(importance);
            const double visibility=std::abs(model.reconstruction(options.measurement).matrix[0])/options.measurementStandardDeviation;
            const bool eligible=defect<0.5;
            profiles.emplace_back(J::Object{{"source",double(j)},{"real",l.real()},{"imag",l.imag()},
                {"importance",importance},{"visibility",visibility},{"propagation_defect",defect},{"eligible",eligible}});
            if(eligible)pool.push_back({j,importance,visibility,defect,importance*visibility/(1+visibility)/(1+defect)});
        }catch(const std::exception&){++rejected;}
    }
    jsonFile(out/"profiles.json",J::Object{{"modes",profiles},{"failed_realization_checks",double(rejected)},
        {"eligible_representatives",double(pool.size())}});
    if(pool.empty())throw std::runtime_error("No eligible mode groups");
    J::Array results;std::set<std::vector<std::size_t>> seen;std::size_t id=0;
    for(const std::string order:{"reconstruction","observability","joint"}) {
        auto sorted=pool;std::sort(sorted.begin(),sorted.end(),[&](const auto& a,const auto& b){
            return order=="reconstruction"?a.importance>b.importance:order=="observability"?a.visibility>b.visibility:a.joint>b.joint;});
        for(std::size_t count:{1,2,4,8,16,32,64}) {
            if(count>sorted.size())continue;std::vector<std::size_t> selected;
            for(std::size_t k=0;k<count;++k)selected.push_back(sorted[k].source);
            std::sort(selected.begin(),selected.end());if(!seen.insert(selected).second)continue;
            J::Object row{{"id",double(id)},{"proposal",order},{"representatives",indices(selected)}};
            try{
                auto r=c::makeRealModalRealization(parent,selected,1e-3,1e-7);
                const auto score=c::assessEstimationSubset(r.model,time,
                    r.model.scatteredEigenfunctions()->trainingValues,x,y,options);
                row["state_mse"]=finite(score.stateReconstructionMse);
                row["echo_mse"]=finite(score.measurementReconstructionMse);
                row["noise_mse"]=finite(score.posteriorStateNoiseMse);
                row["dynamics_mse"]=finite(score.autonomousPropagationMse);row["score"]=finite(score.score);
                row["observability_rank"]=double(score.observability.rank);
                row["real_dimension"]=double(r.model.splitDimension());
                const auto modelPath=out/("candidate-"+std::to_string(id)+".npz");
                J::Object meta{{"format","dynsyslab.derived/v1"},{"method","klqe-subset"},
                    {"source_krom",path.generic_u8string()},{"source_modes",indices(selected)},
                    {"conjugate_modes",indices(r.conjugateModes)},{"coefficients_refitted",false},
                    {"real_rate_tolerance",1e-7},{"sampled_conjugacy_tolerance",1e-3}};
                c::saveKoopmanModel(modelPath,r.model,c::json::dump(meta));row["krom"]=modelPath.generic_u8string();
                c::NpzRecord diag;const auto n=r.model.splitDimension();
                diag.arrays["observability_gramian"]=c::NpyArray::fromDoubles({n,n},score.observability.gramian);
                diag.arrays["gramian_eigenvalues"]=c::NpyArray::fromDoubles({n},score.observability.gramianEigenvalues);
                diag.arrays["directions"]=c::NpyArray::fromDoubles({n,n},score.observability.rightDirections);
                c::writeNpzAtomicNew(out/("candidate-"+std::to_string(id)+"-diagnostics.npz"),diag);
                std::cout<<"Candidate "<<id<<" "<<order<<" "<<count<<" score "<<score.score<<std::endl;
            }catch(const std::exception& e){row["error"]=e.what();std::cout<<"Candidate failed: "<<e.what()<<std::endl;}
            results.emplace_back(std::move(row));++id;
        }
    }
    J::Array sx;for(double s:options.outputScales)sx.emplace_back(s);
    jsonFile(out/"screening.json",J::Object{{"candidates",results},{"state_scales",sx},
        {"noise_sd",options.measurementStandardDeviation},{"information_horizon",options.horizon},
        {"information_sample_time",options.sampleTime},{"noise_prior","unit isotropic RMS-normalized modal prior"}});
}
void run(const Path& modelPath,const Path& sourcePath,const Path& out,double process,double noise) {
    const auto model=c::loadKoopmanModel(modelPath);const auto source=c::readNpz(sourcePath);
    const auto time=source.at("time").doubles();const auto x=source.at(source.contains("states")?"states":"values").doubles();
    const auto N=time.size();if(x.size()!=2*N)throw std::runtime_error("state shape");
    const auto system=c::createSystem(model.systemConfiguration());
    auto y=system->observableCatalog().at("far_field_acoustic_echo").evaluateRealBatch(x,N);
    const auto clean=y;std::mt19937_64 rng(20261105);std::normal_distribution<double> normal(0,noise>0?noise:1.0);
    if(noise>0)for(double& value:y)value+=normal(rng);
    const auto n=model.splitDimension();std::vector<double> Q(n*n),initial(n);
    for(std::size_t j=0;j<n;++j)Q[j*n+j]=process;
    const auto design=c::designKoopmanLQEFromMeasurements(model,"far_field_acoustic_echo",Q,{std::max(1e-10,noise*noise*0.0001)},"state");
    const auto result=c::runKoopmanLQEFromLifted(model,design,time,y,initial,x);
    c::NpzRecord record;
    record.arrays["time"]=c::NpyArray::fromDoubles({N},time);
    record.arrays["values"]=c::NpyArray::fromDoubles({N,2},result.stateEstimate);
    record.arrays["true_states"]=c::NpyArray::fromDoubles({N,2},x);
    record.arrays["measurement_values"]=c::NpyArray::fromDoubles({N,1},y);
    record.arrays["clean_echo"]=c::NpyArray::fromDoubles({N,1},clean);
    record.arrays["measurement_estimate"]=c::NpyArray::fromDoubles({N,1},result.measurementEstimate);
    record.arrays["observer_poles"]=c::NpyArray::fromComplex({n},design.solution.observerPoles);
    record.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(J::Object{
        {"format","dynsyslab.derived/v1"},{"method","klqe"},{"system","rpe"},
        {"source_krom",modelPath.generic_u8string()},{"source_archive",sourcePath.generic_u8string()},
        {"process_variance",process},{"measurement_noise_sd",noise},{"initial_estimate","zero lifted coordinates"}}));
    // Retain Psi comparison values at the fitted sample times, not a second
    // massive copy of the full-rate observer state.
    const auto* data=model.scatteredEigenfunctions();
    const auto nativeMeta=c::json::parse(c::readNpz(modelPath).metadataJson());
    const auto parentPath=nativeMeta.at("source_krom").string();
    if(!parentPath.empty()) {
        const auto parent=c::readNpz(Path(parentPath));
        const auto origin=c::json::parse(parent.metadataJson()).at("source_archive").string();
        if(std::filesystem::absolute(sourcePath).lexically_normal()!=std::filesystem::absolute(Path(origin)).lexically_normal()) {
            c::writeNpzAtomicNew(out/"run.npz",record);
            jsonFile(out/"metrics.json",J::Object{{"state_rmse",metrics(x,result.stateEstimate,2)},
                {"echo_rmse",metrics(clean,result.measurementEstimate,1)},
                {"state_rmse_early_0_20",intervalMetrics(x,result.stateEstimate,2,time,0,20)},
                {"state_rmse_late_60_80",intervalMetrics(x,result.stateEstimate,2,time,60,80)},
                {"echo_rmse_early_0_20",intervalMetrics(clean,result.measurementEstimate,1,time,0,20)},
                {"echo_rmse_late_60_80",intervalMetrics(clean,result.measurementEstimate,1,time,60,80)},
                {"care_residual",design.solution.relativeResidual}});
            return;
        }
        const auto comparison=parent.at("training_time").doubles();
        std::vector<double> estimated;estimated.reserve(comparison.size()*n);
        for(double t:comparison){const auto k=std::min<std::size_t>(N-1,std::size_t(std::llround((t-time[0])/(time[1]-time[0]))));
            estimated.insert(estimated.end(),result.liftedEstimate.begin()+k*n,result.liftedEstimate.begin()+(k+1)*n);}
        record.arrays["comparison_time"]=parent.at("training_time");
        record.arrays["estimated_psi_real"]=c::NpyArray::fromDoubles({comparison.size(),n},estimated);
        record.arrays["reference_psi"]=c::NpyArray::fromComplex({data->samples,model.eigenfunctionCount()},data->trainingValues);
    }
    c::writeNpzAtomicNew(out/"run.npz",record);
    jsonFile(out/"metrics.json",J::Object{{"state_rmse",metrics(x,result.stateEstimate,2)},
        {"echo_rmse",metrics(clean,result.measurementEstimate,1)},
                {"state_rmse_early_0_20",intervalMetrics(x,result.stateEstimate,2,time,0,20)},
                {"state_rmse_late_60_80",intervalMetrics(x,result.stateEstimate,2,time,60,80)},
                {"echo_rmse_early_0_20",intervalMetrics(clean,result.measurementEstimate,1,time,0,20)},
                {"echo_rmse_late_60_80",intervalMetrics(clean,result.measurementEstimate,1,time,60,80)},
                {"care_residual",design.solution.relativeResidual}});
}
int main(int argc,char** argv){try{
    if(argc<4)throw std::runtime_error("usage: inspect|screen KROM OUTPUT; run KROM OUTPUT TRAJECTORY Q NOISE_SD");
    const std::string action=argv[1];const Path source=argv[2],out=argv[3];std::filesystem::create_directories(out);
    if(action=="inspect")inspect(source,out);else if(action=="screen")screen(source,out);else if(action=="modal")modal(source,out);
    else if(action=="full-real")fullReal(source,out);
    else if(action=="assess")assess(source,out);
    else if(action=="run"&&argc==7)run(source,Path(argv[4]),out,std::stod(argv[5]),std::stod(argv[6]));
    else throw std::runtime_error("invalid action");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
