#include "HankelAnalysis.hpp"
#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/ObservableDerivative.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using J = json::Value;

std::size_t integer(const J& v,const char* key,std::size_t fallback,std::size_t minimum=0) {
    if (v.isNull()) return fallback;
    double x=v.number(-1.);
    if (!std::isfinite(x) || x<double(minimum) || x>9007199254740991. || std::floor(x)!=x)
        throw std::invalid_argument(std::string("invalid Hankel integer: ")+key);
    return static_cast<std::size_t>(x);
}

struct Signal {
    std::string name, base, derivative, units;
    std::vector<double> values;
};

class Signals {
public:
    explicit Signals(const SourceData& s):source(s) {}
    Signal resolve(const J& selector,std::size_t order) {
        if (std::holds_alternative<double>(selector.storage())) {
            const auto k=integer(selector,"observable",0);
            if (k>=source.dimension || order) throw std::invalid_argument("use a named observable for its derivative; check component index");
            Signal result{"state_"+std::to_string(k),"state_"+std::to_string(k),"none",{},{}};
            result.values.reserve(source.samples);
            for (std::size_t i=0;i<source.samples;++i) result.values.push_back(source.values[i*source.dimension+k]);
            return result;
        }
        if(std::holds_alternative<J::Object>(selector.storage())) {
            if(order)throw std::invalid_argument("A KROM decomposition reference has no implicit time derivative");
            c::ObservableDecompositionReference ref;ref.artifact=selector.at("artifact").string();
            ref.modelIdentity=selector.at("model_identity").string();ref.component=selector.at("component").string();
            ref.kind=selector.at("kind").string("reconstruction");ref.index=integer(selector.at("index"),"index",0);
            auto observable=c::existingKromObservable(std::filesystem::u8path(ref.artifact),ref,selector.at("name").string("krom-observable"));
            if(observable.stateDimension()!=source.dimension)throw std::invalid_argument("KROM observable source coordinate dimension differs");
            if(selector.at("coordinate_contract").string()!=observable.descriptor().stateCoordinateSystem)
                throw std::invalid_argument("KROM reference requires explicit matching coordinate_contract; dimensions alone do not identify coordinates");
            const auto output=integer(selector.at("output"),"output",0);
            if(output>=observable.outputDimension())throw std::invalid_argument("KROM observable output index out of range");
            const auto part=selector.at("part").string("real");
            if(part!="real"&&part!="imag")throw std::invalid_argument("KROM channel part must be real or imag");
            Signal result{observable.descriptor().name,ref.component,"none","",{}};
            for(std::size_t i=0;i<source.samples;++i) {
                const auto value=observable.evaluate(source.values.data()+i*source.dimension)[output];
                result.values.push_back(part=="real"?value.real():value.imag());
            }
            return result;
        }
        if (!system) {
            system=c::createSystem(systemConfiguration(source));
            if (system->descriptor().stateDimension!=source.dimension)
                throw std::invalid_argument("Hankel observable and source state dimensions differ");
            catalog=system->observableCatalog();
        }
        const auto& base=catalog->at(selector.string());
        if (base.outputDimension()!=1 || base.descriptor().complexValued)
            throw std::invalid_argument("Hankel channels require real scalar observables");
        const bool missingInput=source.inputs.empty() && datasetDynamics(source)=="forced"
            && system->configuration().driver.amplitude==0.;
        std::optional<double> constantInput;
        if (!source.inputs.empty()) {
            if (source.inputDimension==1 && std::all_of(source.inputs.begin(),source.inputs.end(),
                [&](double u){return u==source.inputs.front();})) constantInput=source.inputs.front();
        } else if (!missingInput && system->configuration().driver.amplitude==0.) constantInput=0.;
        auto plan=c::resolveObservableDerivative(*system,*catalog,base.descriptor().name,order,constantInput);
        Signal result{plan.name,base.descriptor().name,plan.method,plan.pushForward?"":plan.map.descriptor().units,{}};
        if (plan.pushForward) {
            if (missingInput) throw std::invalid_argument("Hankel derivative needs the recorded forcing or a prescribed driver");
            if (!velocities) velocities=exactVelocities(source,true);
            result.values=plan.map.pushForwardRealBatch(source.values,*velocities,source.samples,c::ObservableJacobianBackend::analytic);
        } else result.values=plan.map.evaluateRealBatch(source.values,source.samples);
        return result;
    }
private:
    const SourceData& source;
    std::unique_ptr<c::ODESystem> system;
    std::optional<c::ObservableCatalog> catalog;
    std::optional<std::vector<double>> velocities;
};
}

HankelAnalysis prepareHankel(const SourceData& source,const J& settings,bool paired) {
    const auto& saved=settings.at("hankel_archive");
    if (!saved.isNull()) {
        if (saved.string().empty()) throw std::invalid_argument("expected a Hankel matrix archive path");
        for (auto key:{"observables","observable","window","window_alignment","matrix_size","sampling_window",
                       "sample_stride","delay_stride","snapshot_stride","channel_scales"})
            if (!settings.at(key).isNull()) throw std::invalid_argument(std::string("stored Hankel layout already determines ")+key);
        return readHankel(source,std::filesystem::u8path(saved.string()),paired);
    }
    const auto& w=settings.at("window");
    const bool modern=!w.isNull();
    if (modern && !std::holds_alternative<J::Object>(w.storage()))
        throw std::invalid_argument("Hankel window must be an object");
    c::HankelWindow spec;
    if (modern) {
        for (const auto& field:w.object())
            if (!std::set<std::string>{"delays","columns","sample_stride","delay_stride","snapshot_stride"}.count(field.first))
                throw std::invalid_argument("unknown Hankel window setting: "+field.first);
        for (auto key:{"matrix_size","sampling_window","delay_stride","sample_stride","snapshot_stride"})
            if (!settings.at(key).isNull()) throw std::invalid_argument(std::string("conflicting legacy Hankel setting: ")+key);
        spec.delays=integer(w.at("delays"),"delays",0,1);
        if (!spec.delays) throw std::invalid_argument("window.delays is required");
        spec.columns=integer(w.at("columns"),"columns",0);
        spec.sampleStride=integer(w.at("sample_stride"),"sample_stride",1,1);
        spec.delayStride=integer(w.at("delay_stride"),"delay_stride",1,1);
        spec.snapshotStride=integer(w.at("snapshot_stride"),"snapshot_stride",1,1);
    } else if (paired) {
        spec.delays=integer(settings.at("sampling_window"),"sampling_window",20000,2);
        spec.columns=integer(settings.at("matrix_size"),"matrix_size",2000,2);
        spec.sampleStride=integer(settings.at("sample_stride"),"sample_stride",1,1);
        spec.snapshotStride=integer(settings.at("snapshot_stride"),"snapshot_stride",1,1);
    } else {
        spec.delays=integer(settings.at("matrix_size"),"matrix_size",10,2);
        spec.delayStride=integer(settings.at("delay_stride"),"delay_stride",1,1);
    }
    HankelAnalysis out;
    const auto alignment=settings.at("window_alignment").string(modern?"end":"start");
    if (alignment!="end" && alignment!="start") throw std::invalid_argument("invalid Hankel window_alignment");
    out.alignEnd=alignment=="end";
    J::Array requested;
    const auto& list=settings.at("observables");
    const auto& single=settings.at("observable");
    if (!list.isNull() && !std::holds_alternative<J::Array>(list.storage()))
        throw std::invalid_argument("observables must be an array");
    if (!list.array().empty()) {
        requested=list.array();
        if (modern && !single.isNull() && !(std::holds_alternative<std::string>(single.storage()) && single.string().empty()))
            throw std::invalid_argument("specify a nonempty observables list or one observable");
    } else if (!single.isNull() && !(std::holds_alternative<std::string>(single.storage()) && single.string().empty()))
        requested.push_back(single);
    else for (std::size_t k=0;k<(paired?1:source.dimension);++k) requested.emplace_back(double(k));
    Signals resolver(source);
    std::map<std::pair<std::string,std::size_t>,std::size_t> cache;
    std::vector<Signal> signals;
    std::vector<c::HankelChannel> channels;
    std::set<std::string> unique;
    J::Array metadata;
    const auto& scales=settings.at("channel_scales");
    if (!scales.isNull() && !std::holds_alternative<J::Array>(scales.storage()))
        throw std::invalid_argument("channel_scales must be an array");
    std::size_t maxLag=0;
    for (const auto& item:requested) {
        bool object=std::holds_alternative<J::Object>(item.storage());
        if (object) for (const auto& entry:item.object())
            if (entry.first!="observable" && entry.first!="derivative_order" && entry.first!="taps")
                throw std::invalid_argument("unknown Hankel observable field: "+entry.first);
        const auto& selector=object?item.at("observable"):item;
        if (!std::holds_alternative<double>(selector.storage()) && selector.string().empty())
            throw std::invalid_argument("expected a Hankel observable name or index");
        const auto order=integer(item.at("derivative_order"),"derivative_order",0);
        if (order>8) throw std::invalid_argument("Hankel derivative_order exceeds 8");
        auto key=std::make_pair(json::dump(selector),order);
        auto entry=cache.emplace(key,signals.size());
        if (entry.second) signals.push_back(resolver.resolve(selector,order));
        const auto index=entry.first->second;
        const auto& s=signals[index];
        const auto& tv=item.at("taps");
        const J::Array taps=tv.isNull()?J::Array{0.}:tv.array();
        if (taps.empty()) throw std::invalid_argument("Hankel taps must be a nonempty array");
        for (const auto& tap:taps) {
            if (tap.isNull()) throw std::invalid_argument("null Hankel lag");
            const auto lag=integer(tap,"tap",0);
            const auto name=s.name+(lag?"_lag"+std::to_string(lag):"");
            if (!unique.insert(name).second) throw std::invalid_argument("duplicate Hankel channel: "+name);
            const auto j=channels.size();
            const double scale=scales.array().empty()?1.:
                (j<scales.array().size()?scales.array()[j].number(-1.):-1.);
            channels.push_back({index,lag,scale}); out.names.push_back(name); out.scales.push_back(scale);
            metadata.emplace_back(J::Object{{"name",name},{"base",s.base},{"observable",s.name},
                {"derivative_order",double(order)},{"derivative_source",s.derivative},
                {"units",s.units},{"selector",selector},{"measurement_column",double(index)},{"lag_samples",double(lag)},{"scale",scale}});
            maxLag=std::max(maxLag,lag);
        }
    }
    if (!scales.array().empty() && scales.array().size()!=channels.size())
        throw std::invalid_argument("one scale is required per expanded Hankel channel");
    // Legacy H-DMD fixes the number of pairs and shortens the history when needed.
    if (!modern && paired)
        spec=c::fixedPairHankelWindow(source.samples,maxLag,spec.delays,spec.columns,spec.sampleStride,spec.snapshotStride);
    std::vector<double> values;
    if (signals.size()==1) values=std::move(signals.front().values);
    else {
        values.resize(source.samples*signals.size());
        for (std::size_t k=0;k<signals.size();++k)
            for (std::size_t i=0;i<source.samples;++i) values[i*signals.size()+k]=signals[k].values[i];
    }
    out.matrix=c::buildHankelMatrix(values,source.samples,signals.size(),channels,spec,paired);
    // Preserve the legacy POD by-index operation on irregular clocks. A
    // physical one-step DMD map and new time-based windows require uniform time.
    J analysisDt;
    if (!modern && !paired) {
        try { analysisDt=uniformDeltaTime(source)*spec.sampleStride; }
        catch (const std::invalid_argument&) { analysisDt=J{}; }
    } else analysisDt=uniformDeltaTime(source)*spec.sampleStride;
    out.provenance=J::Object{{"schema","dynsyslab.hankel-window/v1"},{"channels",metadata},
        {"delays",double(spec.delays)},{"columns",double(out.matrix.columns)},{"rows",double(out.matrix.rows)},
        {"sample_stride",double(spec.sampleStride)},{"delay_stride",double(spec.delayStride)},
        {"snapshot_stride",double(spec.snapshotStride)},{"analysis_dt",analysisDt},
        {"successor_stride",paired?J(1.):J{}},{"alignment",alignment},{"row_order","oldest-delay-then-channel"},
        {"consecutive_blocks",spec.delayStride==spec.snapshotStride}};
    return out;
}
} // namespace dynsyslab::workflow::postprocess_detail
