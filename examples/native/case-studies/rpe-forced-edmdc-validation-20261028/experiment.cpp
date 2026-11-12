// Frozen-fit validation only. Numerical metrics come from the shared core;
// data verification/publication uses the standard post-processing workflow.
#include <dynsyslab/core/EdmdControlArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/workflow/KromValidationPlan.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
namespace c=dynsyslab::core;
namespace w=dynsyslab::workflow;
using J=c::json::Value;using Path=std::filesystem::path;
J zjson(std::complex<double> z){return J::Array{z.real(),z.imag()};}
void writeJson(const Path& path,const J& value){
    if(std::filesystem::exists(path))throw std::runtime_error("Refusing overwrite");
    std::ofstream f(path);f<<c::json::dump(value)<<'\n';if(!f)throw std::runtime_error("Write failed");
}
J metrics(const c::KoopmanValidationResult& result){
    J::Array rows;
    for(std::size_t j=0;j<result.modes.size();++j){const auto& m=result.modes[j];J::Array horizons;
        for(const auto& h:m.horizons)horizons.emplace_back(J::Object{{"horizon",h.horizon},{"relative_l2",h.statistics.relativeL2}});
        rows.emplace_back(J::Object{{"index",double(j)},{"one_step",m.oneStep.relativeL2},
            {"change_relative",m.oneStepChangeDefined?J(m.oneStepChangeRelativeL2):J{}},
            {"anchored",m.anchoredLongTerm.relativeL2},{"generator",m.generator.relativeL2},
            {"one_step_absolute_rms",m.oneStep.absoluteRms},{"horizons",std::move(horizons)}});
    }
    return J::Object{{"samples",double(result.samples)},{"outside_training_fraction",result.outsideTrainingFraction},
        {"maximum_relative_excursion",result.maximumRelativeExcursion},{"modes",std::move(rows)}};
}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::runtime_error("SOURCE_ROOT OUTPUT_DIRECTORY");
    const Path root=std::filesystem::absolute(argv[1]),out=std::filesystem::absolute(argv[2]);
    if(std::filesystem::exists(out/"all-modes.npz"))throw std::runtime_error("Experiment exists; select a new output directory");
    std::filesystem::create_directories(out);
    const auto original=root/"archives/experiments/rpe-forced-edmdc-20261022-v2";
    const auto old=c::loadKoopmanModel(original/"selected-polynomial.npz");
    const auto edmdc=c::loadEdmdControlModel(original/"edmdc.npz");
    auto model=c::KoopmanModel(edmdc.polynomialEigenfunctions(),old.systemConfiguration());
    const auto prior=c::json::load((original/"results.json").wstring());
    for(const auto& row:prior.at("metadata").at("training_sources").array())
        model.addDatasetReference({row.at("states").string(),"training","forced",row.at("inputs").string()});
    for(const auto& row:prior.at("metadata").at("validation_sources").array())
        model.addDatasetReference({row.at("states").string(),"test","forced",row.at("inputs").string()});
    std::vector<Path> sources;
    for(int i=0;i<5;++i){
        sources.push_back(root/("archives/solutions/rpe/rpe-viscous-echo-long-validation-five-ics-trajectory-"
            +std::to_string(i)+"__20261002T001403Z.npz"));
        if(!std::filesystem::is_regular_file(sources.back()))throw std::runtime_error("Missing validation archive");
        model.addDatasetReference({sources.back().generic_u8string(),"validation","unforced",{}});
    }
    std::vector<double> low(2,std::numeric_limits<double>::infinity()),high(2,-std::numeric_limits<double>::infinity());
    {const auto a=c::readNpz(original/"selected-polynomial.npz");const auto x=a.at("training_states").doubles();
        for(std::size_t k=0;k<x.size();++k){low[k%2]=std::min(low[k%2],x[k]);high[k%2]=std::max(high[k%2],x[k]);}}
    const auto modelPath=out/"all-modes.npz";
    const auto plan=w::planAutonomousKoopmanValidation(model,modelPath,root/"archives");
    writeJson(out/"validation-plan.json",plan);
    J meta=J::Object{{"method","edmdc"},{"extraction_algorithm","edmdc"},{"label","Frozen forced EDMDc: free validation"},
        {"source_edmdc",(original/"edmdc.npz").generic_u8string()},
        {"source_selected_krom",(original/"selected-polynomial.npz").generic_u8string()},
        {"delta_t",.01},{"integration_form","dimensionless"},{"control_ready",false},
        {"validation_policy","unforced held-out trajectories; no coefficient, eigenvalue, normalization, or B refit"},
        {"autonomous_validation_plan",plan}};
    const std::map<std::string,c::NpyArray> evidence{
        {"training_minimum",c::NpyArray::fromDoubles({2},low)},
        {"training_maximum",c::NpyArray::fromDoubles({2},high)}};
    c::saveKoopmanModel(modelPath,model,c::json::dump(meta),evidence);
    c::saveKoopmanModel(out/"selected-mode-10.npz",model.extract({10}),c::json::dump(meta),evidence);
    if(c::loadKoopmanModel(modelPath).eigenfunctions().coefficients!=model.eigenfunctions().coefficients)
        throw std::runtime_error("Frozen coefficients changed on round-trip");
    J::Array reports,spectrum;const double equilibrium[]{1.,0.};const auto eq=model.evaluate(equilibrium);
    for(std::size_t j=0;j<model.liftedDimension();++j)spectrum.emplace_back(J::Object{
        {"index",double(j)},{"mu",zjson(edmdc.lifted.discreteEigenvalues[j])},
        {"principal_rate",zjson(model.continuousEigenvalues()[j])},{"psi_at_equilibrium",zjson(eq[j])}});
    for(std::size_t i=0;i<sources.size();++i){
        const J settings=J::Object{{"krom",modelPath.generic_u8string()},
            {"dataset_role","validation"},{"retain_validation_pointer",false},
            {"validation_delta_t",.01},{"horizons",J::Array{.1,1.,5.,10.}}};
        J job=J::Object{{"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
            {"label","Frozen EDMDc free validation IC "+std::to_string(i+1)},
            {"sources",J::Array{sources[i].generic_u8string()}},
            {"methods",J::Array{J::Object{{"name","koopman_validation"},{"settings",settings}}}},
            {"archive",J::Object{{"enabled",true}}}};
        const auto published=w::executeJob(job,root/"archives");
        if(published.archives.size()!=1)throw std::runtime_error("Validation publication failed");
        const auto a=c::readNpz(sources[i]);const auto rawT=a.at("time").doubles(),rawX=a.at("values").doubles();
        const auto sourceMeta=c::json::parse(a.metadataJson());
        std::vector<double> t,x;
        for(std::size_t k=0;k<rawT.size();k+=100){t.push_back(rawT[k]);x.insert(x.end(),{rawX[2*k],rawX[2*k+1]});}
        if(std::abs(t[1]-t[0]-.01)>1e-10)throw std::runtime_error("Validation sampling mismatch");
        J::Array windows;
        for(const auto& interval:std::vector<std::pair<double,double>>{{0,80},{0,5},{5,20},{20,80}}){
            std::vector<double> wt,wx;
            for(std::size_t k=0;k<t.size();++k)if(t[k]>=interval.first-1e-10&&t[k]<=interval.second+1e-10){
                wt.push_back(t[k]);wx.insert(wx.end(),{x[2*k],x[2*k+1]});}
            auto record=metrics(c::validateKoopmanSemigroup(model,wt,wx,wt.size(),2,{.1,1.,5.,10.},low,high));
            record.object()["start"]=interval.first;record.object()["stop"]=interval.second;windows.push_back(record);
        }
        // Display evidence for the previously selected mode only; source states
        // and all training trajectories remain pointers, not duplicated arrays.
        std::vector<std::complex<double>> observed,predicted;const auto p0=model.evaluate(x.data())[10];
        for(std::size_t k=0;k<t.size();++k){observed.push_back(model.evaluate(&x[2*k])[10]);
            predicted.push_back(std::pow(edmdc.lifted.discreteEigenvalues[10],double(k))*p0);}
        c::NpzRecord display;display.arrays["time"]=c::NpyArray::fromDoubles({t.size()},t);
        display.arrays["observed_psi"]=c::NpyArray::fromComplex({t.size()},observed);
        display.arrays["ideal_psi"]=c::NpyArray::fromComplex({t.size()},predicted);
        display.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(J::Object{
            {"format","dynsyslab.derived/v1"},{"method","koopman-validation-display"},
            {"source_archive",sources[i].generic_u8string()},{"source_krom",modelPath.generic_u8string()},{"source_mode",10.}}));
        c::writeNpzAtomicNew(out/("ic-"+std::to_string(i+1)+"-display.npz"),display);
        reports.emplace_back(J::Object{{"source",sources[i].generic_u8string()},
            {"initial_state",sourceMeta.at("initial_state")},{"validation_archive",published.archives.front().generic_u8string()},
            {"windows",windows}});
        std::cout<<"Validated IC "<<i+1<<": "<<c::json::dump(sourceMeta.at("initial_state"))<<std::endl;
    }
    writeJson(out/"results.json",J::Object{{"schema","dynsyslab.rpe-edmdc-autonomous-validation/v1"},
        {"krom",modelPath.generic_u8string()},{"selected_mode",10.},{"delta_t",.01},{"spectrum",spectrum},
        {"training_minimum",J::Array{low[0],low[1]}},{"training_maximum",J::Array{high[0],high[1]}},
        {"training_refit",false},{"generator_method","sampled finite differences; principal logarithm branch diagnostic"},
        {"validation",reports}});
    std::cout<<"Saved "<<(out/"results.json").generic_string()<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
