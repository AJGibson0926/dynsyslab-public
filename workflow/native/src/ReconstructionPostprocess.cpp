#include "PostprocessInternal.hpp"
#include "SampledDomain.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/ReconstructionQuality.hpp>
#include <dynsyslab/core/Text.hpp>
#include <stdexcept>
#include <tuple>
namespace dynsyslab::workflow::postprocess_detail {
std::filesystem::path runObservableReconstruction(const std::filesystem::path& sourcePath,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    const auto krom=settings.at("krom").string();
    if(krom.empty())throw std::invalid_argument("observable reconstruction requires a KROM");
    const auto archive=dynsyslab::core::readNpz(sourcePath);
    SourceData source;source.path=sourcePath;source.metadataJson=archive.metadataJson();
    const auto array=settings.at("observable_array").string();
    const auto name=settings.at("observable").string(array);
    if(name.empty())throw std::invalid_argument("observable reconstruction requires a name");
    auto model=dynsyslab::core::loadKoopmanModel(dynsyslab::core::wide(krom));
    dynsyslab::core::ReconstructionResult fit;
    // settings.quality: the held-out gate (ReconstructionQuality.hpp). Keys validation_fraction,
    // split interleaved|tail, validation_source (an archive scored instead of a split),
    // minimum_r2, maximum_normalized_rmse, maximum_condition_number, allow_quality_failure.
    // The readout is then fitted on the training rows only; a failing fit is refused unless
    // allow_quality_failure, which publishes it marked unvalidated.
    const auto& qualitySettings=settings.at("quality");
    json::Value qualityReport;
    const auto gate=[&](const std::vector<std::complex<double>>& psiValues,std::size_t rows,
                        const std::vector<std::complex<double>>& yValues,std::size_t outputs,
                        const std::vector<std::complex<double>>& validationPsi,
                        const std::vector<std::complex<double>>& validationY,std::size_t validationRows) {
        dynsyslab::core::ReconstructionQualityPolicy policy;
        policy.validationFraction=qualitySettings.at("validation_fraction").number(policy.validationFraction);
        policy.split=qualitySettings.at("split").string(policy.split);
        policy.minimumR2=qualitySettings.at("minimum_r2").number(policy.minimumR2);
        policy.maximumNormalizedRmse=qualitySettings.at("maximum_normalized_rmse").number(policy.maximumNormalizedRmse);
        policy.maximumConditionNumber=qualitySettings.at("maximum_condition_number").number(policy.maximumConditionNumber);
        const auto quality=dynsyslab::core::assessReconstruction(psiValues,yValues,rows,model.liftedDimension(),
            outputs,settings.at("ridge").number(),policy,validationPsi,validationY,validationRows);
        const bool allow=qualitySettings.at("allow_quality_failure").boolean(false);
        json::Value::Array r2,nrmse,failures;
        for(double v:quality.r2)r2.emplace_back(v);
        for(double v:quality.normalizedRmse)nrmse.emplace_back(v);
        for(const auto& f:quality.failures)failures.emplace_back(f);
        qualityReport=json::Value::Object{{"passed",quality.passed},{"validated",quality.passed},
            {"r2",json::Value(std::move(r2))},{"normalized_rmse",json::Value(std::move(nrmse))},
            {"condition_number",quality.conditionNumber},{"training_samples",double(quality.trainingSamples)},
            {"validation_samples",double(quality.validationSamples)},{"failures",json::Value(std::move(failures))},
            {"policy",qualitySettings}};
        if(!quality.passed&&!allow)
            throw std::runtime_error("observable reconstruction failed its held-out quality gate: "
                +quality.failures.front()+" (set quality.allow_quality_failure to publish it unvalidated)");
        return quality.fit;
    };
    const auto validationArrays=[&](std::size_t outputs,bool stateRoute,const dynsyslab::core::Observable* observable)
        -> std::tuple<std::vector<std::complex<double>>,std::vector<std::complex<double>>,std::size_t> {
        const auto path=qualitySettings.at("validation_source").string();
        if(path.empty())return {};
        auto other=loadSource(std::filesystem::u8path(path));
        if(!stateRoute)throw std::invalid_argument("quality.validation_source needs the state route (no observable_array)");
        normalizeSourceCoordinates(other);
        std::vector<std::complex<double>> psiValues,yValues;
        for(std::size_t k=0;k<other.samples();++k) {
            const auto* x=other.values().data()+k*other.dimension();
            const auto p=model.evaluate(x);psiValues.insert(psiValues.end(),p.begin(),p.end());
            const auto y=observable->evaluate(x);yValues.insert(yValues.end(),y.begin(),y.end());
        }
        (void)outputs;
        return {psiValues,yValues,other.samples()};
    };
    if(array.empty()) {
        source=loadSource(sourcePath);normalizeSourceCoordinates(source);
        // Loading with source physics checks system/parameter/coordinate compatibility.
        model=dynsyslab::core::loadKoopmanModel(dynsyslab::core::wide(krom),systemConfiguration(source));
        auto system=dynsyslab::core::createSystem(systemConfiguration(source));
        const auto observable=system->observableCatalog().at(name);
        if(qualitySettings.isNull())
            fit=dynsyslab::core::fitAndSetObservableReconstruction(model,name,source.values(),
                source.samples(),observable,settings.at("ridge").number());
        else {
            std::vector<std::complex<double>> psiValues,yValues;
            for(std::size_t k=0;k<source.samples();++k) {
                const auto* x=source.values().data()+k*source.dimension();
                const auto p=model.evaluate(x);psiValues.insert(psiValues.end(),p.begin(),p.end());
                const auto y=observable.evaluate(x);yValues.insert(yValues.end(),y.begin(),y.end());
            }
            const auto outputs=observable.outputDimension();
            const auto [vp,vy,vr]=validationArrays(outputs,true,&observable);
            const auto trained=gate(psiValues,source.samples(),yValues,outputs,vp,vy,vr);
            fit=dynsyslab::core::fitAndSetObservableReconstruction(model,name,psiValues,source.samples(),yValues,
                outputs,settings.at("ridge").number());
            if(vr==0) {  // a split: attach the training-rows fit that was scored
                fit=trained;
            }
            model.setReconstruction(name,fit);
        }
    } else {
        // A common source archive supplies time-aligned measurements and Psi.
        // This route is usable with purely measured data and history coordinates.
        // A krom source supplies its own samples (eigenfunction_values) and, since it keeps no state
        // copy, its "training_states" through its origin.
        std::optional<dynsyslab::core::KoopmanModel> sourceKrom;
        if(dynsyslab::core::isKrom(archive))sourceKrom.emplace(dynsyslab::core::loadKoopmanModel(sourcePath));
        const auto psiName=settings.at("eigenfunction_array").string("eigenfunction_values");
        const auto psi=!archive.contains(psiName)&&psiName=="eigenfunction_values"&&sourceKrom
            ? dynsyslab::core::NpyArray::fromComplex({sourceKrom->samples().samples(),sourceKrom->liftedDimension()},
                                                     sourceKrom->samples().values)
            : archive.at(psiName);
        const auto values=!archive.contains(array)&&array=="training_states"&&sourceKrom
            ? [&]{
                const auto states=trainingStates(*sourceKrom);
                return dynsyslab::core::NpyArray::fromDoubles({sourceKrom->samples().samples(),states.dimension},states.values);
            }()
            : archive.at(array);
        if(psi.shape.size()!=2||psi.shape[1]!=model.liftedDimension()
            ||values.shape.empty()||values.shape.size()>2||values.shape[0]!=psi.shape[0])
            throw std::invalid_argument("measured reconstruction arrays are not sample aligned");
        const std::size_t rows=psi.shape[0];
        auto complexValues=[](const dynsyslab::core::NpyArray& array) {
            if(array.type==dynsyslab::core::ArrayType::complex128)return array.complexDoubles();
            const auto real=array.doubles();return std::vector<std::complex<double>>(real.begin(),real.end());
        };
        const std::size_t outputs=values.shape.size()==2?values.shape[1]:1;
        fit=dynsyslab::core::fitAndSetObservableReconstruction(model,name,complexValues(psi),
            rows,complexValues(values),outputs,settings.at("ridge").number());
        if(!qualitySettings.isNull()) {
            fit=gate(complexValues(psi),rows,complexValues(values),outputs,{},{},0);
            model.setReconstruction(name,fit);
        }
    }
    const auto details="\"source_krom\":\""+escaped(krom)+"\",\"observable\":\""+escaped(name)
        +"\",\"reconstruction_fit\":\"least-squares\",\"ridge\":"+std::to_string(settings.at("ridge").number());
    // A child of the krom: its modes and the readout C fitted here.
    std::map<std::string,dynsyslab::core::NpyArray> extra;
    extra["reconstruction_fit_errors"]=dynsyslab::core::NpyArray::fromDoubles(
        {3},{fit.residualRms,fit.residualMaximum,fit.relativeResidual});
    auto meta=metadata(source,"observable-reconstruction",label,details);
    if(!qualityReport.isNull()) {
        auto parsed=json::parse(meta);parsed.object()["reconstruction_quality"]=qualityReport;
        parsed.object()["validated"]=qualityReport.at("passed");
        meta=json::dump(parsed);
    }
    return publishKrom(root,source,model,"koopman-reconstruction",meta,std::move(extra));
}
} // namespace dynsyslab::workflow::postprocess_detail
