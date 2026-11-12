#include <dynsyslab/workflow/KoopmanSpectralQuery.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include "KoopmanSpectralSettings.hpp"
namespace dynsyslab::workflow {
dynsyslab::core::json::Value evaluateKoopmanSpectralQuery(const dynsyslab::core::json::Value& request) {
    using namespace spectral_settings;
    if(request.at("schema").string("dynsyslab.koopman-spectral-query/v1")!="dynsyslab.koopman-spectral-query/v1")
        throw std::invalid_argument("Unknown Koopman spectral query schema");
    const auto action=request.at("action").string("fit");
    if(action=="specrkhs") {
        const auto data=pairs(request.at("data"));const auto& settings=request.at("settings");
        c::SpectralData result;result.rkhsAdjoints.push_back(c::rkhsAdjointSpectrum(data,settings.at("bandwidth").number(1),
            count(settings.at("rank")),settings.at("svd_tolerance").number(1e-10),complexes(settings.at("residual_grid"))));
        return c::spectralDataJson(result);
    }
    if(action=="correlation_measure") {
        c::SpectralData data;data.measures.push_back(c::correlationSpectralMeasure(complexes(request.at("signal")),
            count(request.at("lags")),numbers(request.at("angles")),request.at("assume_stationary").boolean()));
        return c::spectralDataJson(data);
    }
    if(action=="pagerank")return J::Object{{"scores",array(c::koopmanPageRank(numbers(request.at("transition")),
        count(request.at("size")),numbers(request.at("personalization")),request.at("damping").number(.85)))}};
    if(action=="observable") {
        const auto& v=request.at("reference");c::ObservableDecompositionReference ref;
        ref.artifact=v.at("artifact").string();ref.modelIdentity=v.at("model_identity").string();ref.component=v.at("component").string();
        ref.kind=v.at("kind").string("reconstruction");ref.index=count(v.at("index"));
        auto observable=c::existingKromObservable(std::filesystem::u8path(ref.artifact),ref,request.at("name").string("krom-observable"));
        auto points=numbers(request.at("points"));
        if(points.empty()||points.size()%observable.stateDimension())throw std::invalid_argument("Referenced observable points shape mismatch");
        std::vector<std::complex<double>> values;
        for(std::size_t i=0;i<points.size();i+=observable.stateDimension()) {
            auto sample=observable.evaluate(points.data()+i);values.insert(values.end(),sample.begin(),sample.end());
        }
        return J::Object{{"values",array(values)},{"outputs",double(observable.outputDimension())},{"model_identity",ref.modelIdentity}};
    }
    c::SpectralKROM model;
    if(action=="fit") {
        const auto data=pairs(request.at("data"));const auto& settings=request.at("settings");
        model=c::fitKoopmanSpectral(data,dictionary(settings,data),options(settings));
        model.coordinateContract=settings.at("coordinate_contract").string("measured");
        c::fitSpectralReconstruction(model,"coordinates",data.x,data.x,data.dimension,data.weights);
        for(const auto& r:request.at("readouts").object())c::fitSpectralReconstruction(model,r.first,data.x,
            numbers(r.second.at("values")),count(r.second.at("outputs")),data.weights);
        measures(model,settings);
    } else {
        model=request.at("model").isNull()?c::loadSpectralKrom(std::filesystem::u8path(request.at("archive").string())):
            c::spectralKromFromJson(request.at("model"));
        if(action=="evaluate")return J::Object{{"values",array(model.evaluateEigenfunctions(numbers(request.at("points"))))},{"modes",double(model.rank)}};
        if(action=="reconstruct")return J::Object{{"values",array(model.reconstruct(request.at("reconstruction").string("coordinates"),numbers(request.at("points"))))}};
        if(action=="rollout")return J::Object{{"values",array(c::rolloutSpectralKrom(model,request.at("reconstruction").string("coordinates"),
            numbers(request.at("initial")),count(request.at("steps"))))},{"sample_interval",model.sampleInterval}};
        if(action=="residuals") {
            c::SpectralData d;d.residuals=c::koopmanResiduals(model,pairs(request.at("data")),
                complexes(request.at("residual_grid")));return c::spectralDataJson(d);
        }
        if(action=="rigged_measure")measures(model,request);
        else if(action!="inspect")throw std::invalid_argument("Unknown Koopman spectral query action");
    }
    const auto output=request.at("output_archive").string();
    if(!output.empty())c::writeNpzAtomicNew(std::filesystem::u8path(output),c::makeSpectralKromArchive(model));
    J::Object result{{"schema","dynsyslab.koopman-spectral-result/v1"},{"method",model.method},
        {"model_identity",c::spectralKromIdentity(model)},{"rank",double(model.rank)},
        {"assumptions",model.assumptions},{"spectral",c::spectralDataJson(model.spectral)}};
    if(!output.empty())result["archive"]=output;
    if(request.at("include_model").boolean(true))result["model"]=c::spectralKromJson(model);
    return result;
}
}
