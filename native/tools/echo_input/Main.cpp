#include "Experiment.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace echo_experiment;
namespace {
constexpr double pi=3.14159265358979323846;
void sineRecord(const core::SystemConfiguration& config,const std::filesystem::path& folder,
    const Frozen& f,const std::string& id,double amplitude,double frequency,double phase,
    double duration,double stop,const std::string& role,double radius=1.0) {
    if(std::filesystem::exists(folder/(id+".measurements.npz"))) return;
    std::cout<<"generate "<<id<<" A="<<amplitude<<" f="<<frequency<<std::endl;
    double held=0; std::size_t tick=0;
    simulate(config,folder,id,f,duration,radius,[&](double t,double) {
        if(tick++%f.stride==0) held=(t>=6.0-1e-9 && t<stop-1e-9)
            ? amplitude*std::sin(2*pi*frequency*(t-6.0)+phase):0.0;
        return held;
    },Json::Object{{"role",role},{"amplitude_model",amplitude},{"frequency_cycles_per_tau",frequency},
        {"phase",phase},{"drive_start",6.0},{"drive_stop",stop},{"initial_radius",radius}});
}
Json score(const core::HistoryInputMap& map,const std::vector<core::HistoryInputSample>& rows) {
    const auto m=core::scoreHistoryInputMap(map,rows);
    return Json::Object{{"relative_residual",m.relativeResidual},{"rmse",m.rmse},
        {"zero_input_baseline_rmse",m.unforcedBaselineRmse},{"pairs",double(rows.size())}};
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=3 && argc!=4) throw std::invalid_argument("usage: dynsyslab_echo_input_experiment FROZEN_KROM OUTPUT_DIRECTORY [--refine]");
        const std::filesystem::path source=argv[1],folder=argv[2];
        std::filesystem::create_directories(folder);
        if(argc==3 && std::filesystem::exists(folder/"results.json"))
            throw std::invalid_argument("completed experiment is immutable; choose a new output directory");
        const auto model=core::loadKoopmanModel(source);
        const auto archive=core::readNpz(source);
        const auto f=frozen(archive,model.eigenfunctions());
        const auto config=model.systemConfiguration();
        if(argc==4) {
            if(std::string(argv[3])=="--short-history") {
                shortHistoryStudy(config,folder,f); return 0;
            }
            if(std::string(argv[3])=="--eigenfunctions") {
                eigenfunctionSweep(config,folder,f,archive); return 0;
            }
            if(std::string(argv[3])!="--refine") throw std::invalid_argument("unknown experiment option");
            refineExperiment(config,folder,f,archive); return 0;
        }
        if(!std::filesystem::exists(folder/"provenance.json")) saveJson(folder/"provenance.json",Json::Object{
            {"source_krom",source.generic_string()},{"source_metadata",core::json::parse(archive.metadataJson())},
            {"representation","frozen two-latent degree-five polynomial; no forced eigenfunction training"},
            {"identification_contract","measurement-only files; echo history + input; no states or EOM"},
            {"driver","sampled sine, zero-order held for 0.01 tau"}});
        std::vector<core::HistoryInputSample> train,small,validation;
        Json::Array datasets;
        int run=0;
        for(double amplitude:{0.003,0.01,0.03,0.08,0.2,0.4})
            for(double frequency:{0.65,0.85,1.0,1.15,1.35})
                for(double phase:{0.0,pi/2}) {
                    const std::string id="train-"+std::to_string(++run);
                    const double stop=6+16/frequency,duration=stop+6;
                    sineRecord(config,folder,f,id,amplitude,frequency,phase,duration,stop,
                        "forced input-response identification, not eigenfunction training",1.0);
                    const auto measured=readMeasurements(folder/(id+".measurements.npz"));
                    const auto lifted=lift(measured,f);
                    train.insert(train.end(),lifted.pairs.begin(),lifted.pairs.end());
                    if(amplitude<=0.03) small.insert(small.end(),lifted.pairs.begin(),lifted.pairs.end());
                    datasets.emplace_back(Json::Object{{"id",id},{"split","fit"},
                        {"amplitude",amplitude},{"frequency",frequency},{"phase",phase}});
                }
        run=0;
        for(double amplitude:{0.015,0.12,0.3}) for(double frequency:{0.75,0.93,1.25}) {
            const std::string id="heldout-"+std::to_string(++run);
            const double stop=6+16/frequency;
            sineRecord(config,folder,f,id,amplitude,frequency,0.37,stop+6,stop,"held-out waveform",1.005);
            const auto lifted=lift(readMeasurements(folder/(id+".measurements.npz")),f);
            validation.insert(validation.end(),lifted.pairs.begin(),lifted.pairs.end());
            datasets.emplace_back(Json::Object{{"id",id},{"split","heldout"},
                {"amplitude",amplitude},{"frequency",frequency},{"phase",0.37}});
        }
        std::vector<std::pair<std::string,core::HistoryInputMap>> maps;
        const std::array<double,2> scale{f.polynomial.scale[0],f.polynomial.scale[1]};
        maps.emplace_back("constant-small",core::fitHistoryInputMap(small,f.multiplier,f.controlStep,0,scale));
        maps.emplace_back("bilinear-small",core::fitHistoryInputMap(small,f.multiplier,f.controlStep,1,scale));
        for(std::size_t degree=0;degree<=1;++degree)
            maps.emplace_back("pooled-degree-"+std::to_string(degree),
                core::fitHistoryInputMap(train,f.multiplier,f.controlStep,degree,scale));
        Json::Array fits,results;
        for(const auto& [name,map]:maps) {
            std::cout<<"fit "<<name<<" heldout="<<core::json::dump(score(map,validation))<<std::endl;
            fits.emplace_back(Json::Object{{"name",name},{"degree",double(map.degree)},
                {"design_rank",double(map.rank)},{"design_condition",map.conditionNumber},
                {"fit_small",score(map,small)},{"fit_pooled",score(map,train)},
                {"heldout",score(map,validation)}});
            core::NpzRecord saved;
            saved.arrays["coefficients"]=core::NpyArray::fromComplex({map.coefficients.size()},map.coefficients);
            saved.arrays["coordinate_scale"]=core::NpyArray::fromDoubles({2},{scale[0],scale[1]});
            const auto standard=map.canonicalRepresentation();
            if(const auto* c=std::get_if<core::ConstantControlMap>(&standard)) {
                saved.arrays["matrix"]=core::NpyArray::fromComplex({1,1},c->matrix);
                saved.arrays["discrete_matrix"]=core::NpyArray::fromComplex({1,1},c->discreteMatrix);
                saved.arrays["sample_time"]=core::NpyArray::fromDoubles({1},{c->sampleTime});
            } else {
                const auto& b=std::get<core::BilinearControlMap>(standard);
                saved.arrays["operators"]=core::NpyArray::fromComplex({1,1,1},b.operators);
                saved.arrays["conjugate_operators"]=core::NpyArray::fromComplex({1,1,1},b.conjugateOperators);
                saved.arrays["offsets"]=core::NpyArray::fromComplex({1,1},b.offsets);
            }
            saved.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(Json::Object{
                {"schema","dynsyslab.history-input-map/v1"},{"degree",double(map.degree)},
                {"sample_time",map.sampleTime},{"meaning","discrete response B_d; not an exact state Jacobian"},
                {"source_krom",source.generic_string()}}));
            if(!std::filesystem::exists(folder/(name+".map.npz")))
                core::writeNpzAtomicNew(folder/(name+".map.npz"),saved);
        }
        for(const auto& targetSpec:std::vector<std::pair<std::string,std::array<double,2>>>{
            {"mild",{0.012,0.93}},{"nonlinear",{0.3,0.83}}}) {
            const auto& name=targetSpec.first; const auto spec=targetSpec.second;
            sineRecord(config,folder,f,"target-"+name,spec[0],spec[1],0,60,60,"echo-only control reference");
            const auto target=readMeasurements(folder/("target-"+name+".measurements.npz"));
            const auto reference=lift(target,f);
            for(const auto& [mapName,map]:maps) for(double penalty:{1.0,0.01,0.0001}) {
                const auto id=name+"-"+mapName+"-klqr-r"+std::to_string(penalty);
                std::cout<<"control "<<id<<std::endl;
                try {
                    if(std::filesystem::exists(folder/(id+".measurements.npz"))) {
                        const auto actual=readMeasurements(folder/(id+".measurements.npz"));
                        results.emplace_back(Json::Object{{"id",id},{"late",metrics(actual,target,30.0)},
                            {"whole_driven",metrics(actual,target,6.0)},{"input_penalty_multiplier",penalty}});
                    } else results.push_back(controlTrial(config,folder,id,f,map,target,reference,penalty));
                    std::cout<<core::json::dump(results.back())<<std::endl;
                } catch(const std::exception& e) {
                    results.emplace_back(Json::Object{{"id",id},{"failure",e.what()}});
                    std::cout<<"failed: "<<e.what()<<std::endl;
                }
            }
        }
        saveJson(folder/"results.json",Json::Object{{"datasets",datasets},{"maps",fits},{"controls",results}});
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}
