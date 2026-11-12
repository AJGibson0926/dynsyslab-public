#include "Experiment.hpp"
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace echo_experiment {
namespace {
struct FreeRun { std::vector<double> time,coordinates; };
FreeRun projectFree(const Measurements& m,const Frozen& f) {
    if(std::any_of(m.input.begin(),m.input.end(),[](double u){return u!=0.0;}))
        throw std::invalid_argument("eigenfunction training must be unforced");
    const auto all=core::projectHistory(m.echo,f.basis);
    FreeRun run;
    for(std::size_t k=0;k<all.windows;k+=f.stride) {
        run.time.push_back(m.time[k+f.basis.delays-1]);
        run.coordinates.insert(run.coordinates.end(),all.coordinates.begin()+k*f.basis.rank,
            all.coordinates.begin()+(k+1)*f.basis.rank);
    }
    return run;
}
FreeRun subset(const FreeRun& in,std::size_t rank,double start) {
    FreeRun out;
    for(std::size_t k=0;k<in.time.size();++k) if(in.time[k]>=start) {
        out.time.push_back(in.time[k]);
        out.coordinates.insert(out.coordinates.end(),in.coordinates.begin()+k*4,
            in.coordinates.begin()+k*4+rank);
    }
    return out;
}
Json validate(const core::KoopmanSINDyResult& p,const core::SystemConfiguration& config,
    const std::vector<FreeRun>& runs) {
    const core::KoopmanModel model(p,config);
    double one=0,four=0,pair=0;
    for(const auto& run:runs) {
        const auto x=subset(run,p.stateDimension,5.0);
        const auto v=core::validateKoopmanSemigroup(model,x.time,x.coordinates,
            x.time.size(),p.stateDimension,{1.0,4.0});
        pair+=v.modes[0].oneStep.relativeL2;
        one+=v.modes[0].horizons[0].statistics.relativeL2;
        four+=v.modes[0].horizons[1].statistics.relativeL2;
    }
    return Json::Object{{"pair",pair/runs.size()},{"one_tau",one/runs.size()},
        {"four_tau",four/runs.size()},{"validation_interval","5 to 80 tau, every candidate"}};
}
void retest(const std::string& modelName,const core::SystemConfiguration& config,
    const std::filesystem::path& folder,const std::filesystem::path& inputs,
    const Frozen& f,Json::Array& fits,Json::Array& controls) {
    std::vector<core::HistoryInputSample> small,pooled;
    for(int k=1;k<=60;++k) {
        const auto x=lift(readMeasurements(inputs/("train-"+std::to_string(k)+".measurements.npz")),f);
        pooled.insert(pooled.end(),x.pairs.begin(),x.pairs.end());
        if(k<=30) small.insert(small.end(),x.pairs.begin(),x.pairs.end());
    }
    for(const auto& name:std::vector<std::string>{"mild","nonlinear"}) {
        const auto target=readMeasurements(inputs/("target-"+name+".measurements.npz"));
        const auto targetLift=lift(target,f);
        std::vector<core::HistoryInputSample> local;
        for(int k=1;k<=8;++k) {
            const auto measured=readMeasurements(inputs/("local-"+name+"-"+std::to_string(k)+".measurements.npz"));
            const auto x=lift(measured,f);
            for(std::size_t j=0;j<x.pairs.size();++j)
                if(x.time[j]>=12 && x.time[j]<measured.time.back()-6.01) local.push_back(x.pairs[j]);
        }
        for(const auto& set:std::vector<std::string>{"broad","local"}) {
            const auto& rows=set=="local"?local:(name=="mild"?small:pooled);
            // Same canonical conjugate-complete bilinear map family as the existing framework.
            const auto map=core::fitHistoryInputMap(rows,f.multiplier,f.controlStep,1,{1,1});
            const auto fit=core::scoreHistoryInputMap(map,rows);
            const auto check=core::scoreHistoryInputMap(map,targetLift.pairs);
            const auto id=modelName+"-"+name+"-"+set;
            fits.emplace_back(Json::Object{{"id",id},{"fit_relative_residual",fit.relativeResidual},
                {"target_relative_residual",check.relativeResidual},{"rank",double(map.rank)},
                {"condition",map.conditionNumber}});
            core::NpzRecord saved;
            const auto b=std::get<core::BilinearControlMap>(map.canonicalRepresentation());
            saved.arrays["operators"]=core::NpyArray::fromComplex({1,1,1},b.operators);
            saved.arrays["conjugate_operators"]=core::NpyArray::fromComplex({1,1,1},b.conjugateOperators);
            saved.arrays["offsets"]=core::NpyArray::fromComplex({1,1},b.offsets);
            saved.arrays["coefficients"]=core::NpyArray::fromComplex({map.coefficients.size()},map.coefficients);
            saved.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(fits.back()));
            if(!std::filesystem::exists(folder/(id+".map.npz"))) core::writeNpzAtomicNew(folder/(id+".map.npz"),saved);
            for(double penalty:{100.0,10.0,1.0,0.01}) {
                const auto trial=id+"-klqr-r"+std::to_string(penalty);
                std::cout<<"control "<<trial<<std::endl;
                try {
                    if(std::filesystem::exists(folder/(trial+".measurements.npz"))) {
                        const auto actual=readMeasurements(folder/(trial+".measurements.npz"));
                        controls.emplace_back(Json::Object{{"id",trial},{"late",metrics(actual,target,30)}});
                    } else controls.push_back(controlTrial(config,folder,trial,f,map,target,targetLift,penalty));
                    std::cout<<core::json::dump(controls.back())<<std::endl;
                } catch(const std::exception& e) { controls.emplace_back(Json::Object{{"id",trial},{"failure",e.what()}}); }
            }
        }
    }
}
}
void eigenfunctionSweep(const core::SystemConfiguration& config,const std::filesystem::path& folder,
    const Frozen& original,const core::NpzRecord& archive,const std::filesystem::path& inputFolder) {
    if(std::filesystem::exists(folder/"eigenfunction-sweep.json"))
        throw std::invalid_argument("completed eigenfunction sweep is immutable");
    std::vector<FreeRun> train,validation;
    const auto inputs=inputFolder.empty()?folder:inputFolder;
    Json::Array freeRecords,candidates,selected,fits,controls;
    for(const auto& split:std::vector<std::string>{"train","validation"}) {
        int k=0;
        for(double radius:split=="train"?std::vector<double>{0.65,0.8,0.9,1.05,1.2,1.4,1.7,2.0}
            :std::vector<double>{0.7,0.85,1.1,1.3,1.55,1.85}) {
            const auto id="free-"+split+"-"+std::to_string(++k);
            std::cout<<"unforced "<<id<<std::endl;
            if(!std::filesystem::exists(inputs/(id+".measurements.npz")))
                simulate(config,inputs,id,original,80,radius,[](double,double){return 0.0;},
                    Json::Object{{"forcing_type","free"},{"controller_kind","none"},
                        {"data_role",split=="train"?"free_training":"free_validation"},{"initial_radius",radius}});
            const auto observations=readMeasurements(inputs/(id+".measurements.npz"));
            (split=="train"?train:validation).push_back(projectFree(observations,original));
            freeRecords.emplace_back(Json::Object{{"id",id},{"initial_radius",radius},{"split",split}});
        }
    }
    const auto oldCheck=inputFolder.empty()?validate(original.polynomial,config,validation)
        :Json("not applicable: a new basis needs its own polynomial");
    std::cout<<"original psi free validation "<<core::json::dump(oldCheck)<<std::endl;
    for(std::size_t rank:{2,4}) {
        double best=1e300; core::KoopmanSINDyResult bestPolynomial; Json chosen;
        for(double start:{5.0,15.0,35.0,50.0}) for(std::size_t degree:{1,3,5}) {
            if(rank==4 && degree==5) continue;
            std::vector<std::vector<double>> runs;
            for(const auto& run:train) runs.push_back(subset(run,rank,start).coordinates);
            core::PooledPolynomialEDMDOptions options;
            options.degree=degree; options.pairStride=1; options.sampleTime=original.controlStep;
            options.relativeRidge=1e-10;
            const auto edmd=core::fitPooledPolynomialEDMD(runs,rank,options);
            core::EDMDModeSearchOptions search;
            search.minimumContinuousReal=-0.5;search.maximumContinuousReal=-0.001;
            search.minimumContinuousImaginary=4.0;search.maximumContinuousImaginary=7.0;
            const auto modes=core::selectPolynomialEDMDModes(edmd,runs,search);
            for(const auto& mode:modes) {
                auto p=core::selectedPolynomialEigenfunction(edmd,mode);
                p.dynamicalStateDimension=2;
                const auto quality=validate(p,config,validation);
                const double score=quality.at("one_tau").number()+quality.at("four_tau").number();
                Json row=Json::Object{{"rank",double(rank)},{"degree",double(degree)},{"fit_start",start},
                    {"lambda_re",p.continuousEigenvalues[0].real()},{"lambda_im",p.continuousEigenvalues[0].imag()},
                    {"heldout_free",quality},{"selection_score",score}};
                candidates.push_back(row);
                if(score<best) {best=score;bestPolynomial=p;chosen=row;}
            }
        }
        if(bestPolynomial.coefficients.empty()) throw std::runtime_error("no eligible free-only eigenfunction");
        selected.push_back(chosen);
        std::cout<<"selected "<<core::json::dump(chosen)<<std::endl;
        Frozen f=original; f.polynomial=bestPolynomial;
        f.multiplier=std::exp(f.polynomial.continuousEigenvalues[0]*f.controlStep);
        const auto name="refit-xi"+std::to_string(rank);
        std::map<std::string,core::NpyArray> extra;
        for(const auto& key:{"history_delay_samples","history_sample_time","singular_values",
            "left_singular_vectors","observation_scale"}) extra[key]=archive.at(key);
        const auto output=folder/(name+".krom.npz");
        if(!std::filesystem::exists(output)) core::savePolynomialKoopmanModel(output,
            core::KoopmanModel(bestPolynomial,config),core::json::dump(Json::Object{
                {"data_role","reduced_order_model"},{"method","echo-history-koopman"},
                {"representation_training","unforced echo-only; frozen prior U; see eigenfunction-sweep.json"},
                {"selection",chosen}}),extra);
        retest(name,config,folder,inputs,f,fits,controls);
    }
    saveJson(folder/"eigenfunction-sweep.json",Json::Object{{"free_records",freeRecords},
        {"original_free_validation",oldCheck},{"candidates",candidates},{"selected",selected},
        {"maps",fits},{"controls",controls},{"input_records_directory",inputs.generic_string()}});
}
void shortHistoryStudy(const core::SystemConfiguration& config,const std::filesystem::path& folder,
    const Frozen& original) {
    const auto output=folder/"short-history-0p05";
    std::filesystem::create_directories(output);
    std::vector<std::vector<double>> freeEchoes;
    for(int k=1;k<=8;++k) freeEchoes.push_back(readMeasurements(
        folder/("free-train-"+std::to_string(k)+".measurements.npz")).echo);
    core::PooledHankelOptions options;
    options.delays=101; options.rank=4; options.subspaceIterations=8;
    Frozen f=original; f.basis=core::fitPooledHankelBasis(freeEchoes,options);
    core::NpzRecord basis;
    basis.arrays["history_delay_samples"]=core::NpyArray::fromDoubles({1},{double(f.basis.delays)});
    basis.arrays["history_sample_time"]=core::NpyArray::fromDoubles({1},{f.sensorStep});
    basis.arrays["singular_values"]=core::NpyArray::fromDoubles({4},f.basis.singularValues);
    basis.arrays["cumulative_energy"]=core::NpyArray::fromDoubles({4},f.basis.cumulativeEnergy);
    basis.arrays["left_singular_vectors"]=core::NpyArray::fromDoubles({101,4},f.basis.leftSingularVectors);
    basis.arrays["observation_scale"]=core::NpyArray::fromDoubles({1},{f.basis.observationScale});
    basis.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(Json::Object{
        {"schema","dynsyslab.echo-history-basis/v1"},{"window_tau",0.05},
        {"training","eight unforced free-train measurement records only"},
        {"input_records_directory",folder.generic_string()}}));
    if(!std::filesystem::exists(output/"basis.npz")) core::writeNpzAtomicNew(output/"basis.npz",basis);
    eigenfunctionSweep(config,output,f,basis,folder);
}
} // namespace echo_experiment
