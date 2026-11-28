// The koopman-spectral query (dynsyslab.koopman-spectral-query/v1): dictionaries, spectral fits and the
// analyses of a fitted krom, as JSON requests. A fit is a krom (KoopmanSpectralAnalysis.hpp); actions on
// an existing fit name its archive ("krom"). Results report the spectrum as JSON: eigenvalues and mode
// ids, residuals and pseudospectrum, the spectral bands, and the observables' measures (each with
// its generalized eigenfunctions, "packets", grid x features, and its atoms, derived from the
// observable's Koopman modes).
#include <dynsyslab/workflow/KoopmanSpectralQuery.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KromSpectralObservables.hpp>
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/Text.hpp>
#include "KoopmanSpectralSettings.hpp"
namespace dynsyslab::workflow {
namespace {
using namespace spectral_settings;
// The spectral parts of a krom as a JSON view.
J spectralJson(const c::KoopmanModel& model) {
    auto out=c::spectralDiagnosticsJson(model.spectralDiagnostics());
    const auto& spectrum=model.spectrum();
    J::Array bands,measures;
    for(const auto& band:spectrum.bands())bands.push_back(c::spectralBandJson(band));
    for(const auto& measure:spectrum.measures()) {
        auto item=c::spectralMeasureJson(measure);
        if(measure.packets)item.object()["packets"]=array(measure.packets->coefficients());
        J::Array atoms;
        for(const auto& atom:spectrum.atoms(measure.observable,measure.output))
            atoms.emplace_back(J::Object{{"mode",double(atom.mode)},{"mass",atom.mass}});
        item.object()["atoms"]=J(std::move(atoms));
        measures.push_back(std::move(item));
    }
    std::vector<double> residuals;std::string definition,scope;
    for(std::size_t j=0;j<spectrum.size();++j)if(const auto& r=spectrum.residual(j)){residuals.push_back(r->value);definition=r->definition;scope=r->scope;}
    J residualBlock=J::Object{{"definition",definition},{"scope",scope},{"eigenpair",array(residuals)}};
    if(const auto& p=spectrum.pseudospectrum()) {
        residualBlock.object()["grid"]=array(p->grid);residualBlock.object()["minimum"]=array(p->minimum);
        if(p->minimizers)residualBlock.object()["minimizing_coefficients"]=array(p->minimizers->coefficients());
    }
    out.object()["schema"]="dynsyslab.krom-spectral/v1";
    out.object()["eigenvalues"]=array(spectrum.discreteMultipliers());
    out.object()["mode_ids"]=[&]{J::Array ids;for(auto id:spectrum.modeIds())ids.emplace_back(double(id));return J(ids);}();
    out.object()["bands"]=bands;out.object()["measures"]=measures;out.object()["residuals"]=residualBlock;
    return out;
}
c::KoopmanModel kromOf(const J& request) {
    if(!request.at("model").isNull())
        throw std::invalid_argument("Spectral fits are kroms: pass the krom archive path as \"krom\" (write one with output_archive)");
    const auto path=request.at("krom").string(request.at("archive").string());
    if(path.empty())throw std::invalid_argument("This action needs a krom archive (\"krom\")");
    return c::loadKoopmanModel(std::filesystem::u8path(path));
}
}
dynsyslab::core::json::Value evaluateKoopmanSpectralQuery(const dynsyslab::core::json::Value& request) {
    if(request.at("schema").string("dynsyslab.koopman-spectral-query/v1")!="dynsyslab.koopman-spectral-query/v1")
        throw std::invalid_argument("Unknown Koopman spectral query schema");
    const auto action=request.at("action").string("fit");
    if(action=="dictionary") {
        c::SpectralDictionary d;std::vector<double> points;
        if(!request.at("dictionary").isNull()) {
            d=c::spectralDictionaryFromJson(request.at("dictionary"));points=numbers(request.at("points"));
        } else {
            const auto data=pairs(request.at("data"));d=dictionary(request.at("settings"),data);
            points=request.at("points").isNull()?data.x:numbers(request.at("points"));
        }
        const auto values=d.evaluate(points);
        return J::Object{{"schema","dynsyslab.spectral-dictionary-result/v1"},
            {"dictionary",c::spectralDictionaryJson(d)},{"values",array(values)},
            {"samples",double(points.size()/d.dimension)},{"features",double(d.size())}};
    }
    if(action=="specrkhs") {
        const auto data=pairs(request.at("data"));const auto& settings=request.at("settings");
        c::SpectralDiagnostics result;result.rkhsAdjoints.push_back(c::rkhsAdjointSpectrum(data,settings.at("bandwidth").number(1),
            count(settings.at("rank")),settings.at("svd_tolerance").number(1e-10),complexes(settings.at("residual_grid"))));
        return c::spectralDiagnosticsJson(result);
    }
    if(action=="correlation_measure") {
        const auto estimate=c::correlationSpectralMeasure(complexes(request.at("signal")),
            count(request.at("lags")),numbers(request.at("angles")),request.at("assume_stationary").boolean());
        return J::Object{{"schema","dynsyslab.krom-spectral/v1"},{"bands",J::Array{c::spectralBandJson(estimate.band)}},
            {"measures",J::Array{c::spectralMeasureJson(estimate.measure)}}};
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
    std::optional<c::KoopmanModel> fitted;
    if(action=="fit") {
        const auto data=pairs(request.at("data"));const auto& settings=request.at("settings");
        fitted=c::fitKoopmanSpectral(data,dictionary(settings,data),options(settings));
        auto compression=*fitted->compression();compression.coordinateContract=settings.at("coordinate_contract").string("measured");
        fitted->setCompression(std::move(compression));
        c::fitFeatureReadout(*fitted,"coordinates",data.x,data.x,data.dimension,data.weights);
        for(const auto& r:request.at("readouts").object())c::fitFeatureReadout(*fitted,r.first,data.x,
            numbers(r.second.at("values")),count(r.second.at("outputs")),data.weights);
        riggedMeasures(*fitted,settings.at("rigged_measures"));
    } else {
        fitted=kromOf(request);
        const auto& model=*fitted;
        if(!model.compression())throw std::invalid_argument("This krom has no compressed operator (a spectral fit's)");
        const auto& compression=*model.compression();
        if(action=="evaluate") {  // the eigenfunctions on the frozen dictionary's coordinates
            const auto points=numbers(request.at("points"));const auto n=compression.dictionary.size(),modes=model.liftedDimension();
            const auto phi=compression.dictionary.evaluate(points);const auto& coefficients=c::frozenDictionaryExpansion(model).coefficients();
            std::vector<std::complex<double>> values(points.size()/compression.dictionary.dimension*modes);
            for(std::size_t i=0;i<values.size()/modes;++i)for(std::size_t j=0;j<modes;++j)for(std::size_t k=0;k<n;++k)
                values[i*modes+j]+=phi[i*n+k]*coefficients[j*n+k];
            return J::Object{{"values",array(values)},{"modes",double(modes)}};
        }
        if(action=="reconstruct")return J::Object{{"values",array(compression.reconstruct(request.at("reconstruction").string("coordinates"),numbers(request.at("points"))))}};
        if(action=="rollout")return J::Object{{"values",array(c::rolloutCompression(compression,request.at("reconstruction").string("coordinates"),
            numbers(request.at("initial")),count(request.at("steps"))))},{"sample_interval",compression.sampleInterval}};
        if(action=="rollout_batch"){
            const auto points=numbers(request.at("initial_points"));const auto steps=count(request.at("steps"));
            const auto name=request.at("reconstruction").string("coordinates");
            const auto outputs=compression.readouts.at(name).outputs,dimension=compression.dictionary.dimension;
            if(points.empty()||points.size()%dimension||steps>=10000000||
               points.size()/dimension>10000000/outputs/(steps+1))
                throw std::invalid_argument("Batch rollout requires complete points and at most ten million output values");
            return J::Object{{"values",array(c::rolloutCompressionBatch(compression,name,points,steps))},
                {"sample_interval",compression.sampleInterval},{"trajectories",double(points.size()/dimension)},
                {"samples",double(steps+1)},{"outputs",double(outputs)},
                {"layout","trajectory/time/output"}};
        }
        if(action=="residuals") {  // held-out evidence: the residual block of a copy, the krom unchanged
            auto view=model;std::vector<std::complex<double>> minimizers;
            view.setSpectralEvidence(c::koopmanResiduals(model,pairs(request.at("data")),complexes(request.at("residual_grid")),
                                                         "empirical-held-out",&minimizers));
            if(!minimizers.empty())view.setSpectralFunction(c::pseudospectrumFunctionKey(),c::FunctionExpansion(
                std::make_shared<const c::ExpansionDictionary>(c::frozenDictionaryExpansion(model).dictionary()),std::move(minimizers),
                c::frozenDictionaryExpansion(model).coordinates()));
            return spectralJson(view);
        }
        if(action=="rigged_measure")riggedMeasures(*fitted,request.at("rigged_measures"));
        else if(action!="inspect")throw std::invalid_argument("Unknown Koopman spectral query action");
    }
    const auto& model=*fitted;const auto& compression=*model.compression();
    const auto output=request.at("output_archive").string();
    if(!output.empty())c::saveKoopmanModel(std::filesystem::u8path(output),model,"{}");
    J::Object result{{"schema","dynsyslab.koopman-spectral-result/v1"},{"method",compression.method},
        {"model_identity",c::koopmanCompressionIdentity(compression)},{"rank",double(compression.rank)},
        {"sample_interval",compression.sampleInterval},{"assumptions",compression.assumptions},{"spectral",spectralJson(model)}};
    if(!output.empty())result["archive"]=output;
    if(request.at("include_model").boolean(true))result["compression"]=c::koopmanCompressionJson(compression);
    return result;
}
}
