#include "HankelAnalysis.hpp"
#include "CoordinateChannels.hpp"
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

// The Hankel module never evaluates physics: a channel is a column of the
// time series it receives, chosen by index or by the archive's column name.
// Observables, their derivatives and KROM observables are coordinate changes
// (coordinate_transform) that must run before the series reaches this module.
class Columns {
public:
    explicit Columns(const SourceData& s):source(s),names(timeSeriesColumnNames(s)),
        named(!(s.metadataJson.empty()?J{}:json::parse(s.metadataJson)).at("coordinate_names").isNull()) {}
    Signal resolve(const J& item,const J& selector,std::size_t order) {
        if (order || selectorNeedsCoordinateChange(item,names))
            throw std::invalid_argument("Hankel channel "+json::dump(item)+" is not a column of the input time series ("
                +join()+"); evaluate it with a coordinate_transform stage and pass that output to this method");
        std::size_t k=0;
        if (std::holds_alternative<double>(selector.storage())) {
            k=integer(selector,"observable",0);
            if (k>=source.dimension) throw std::invalid_argument("Hankel column index exceeds the time-series width");
        } else {
            const auto name=selector.string();
            const auto found=std::find(names.begin(),names.end(),name);
            k=found!=names.end()?std::size_t(found-names.begin()):std::size_t(std::stoul(name.substr(6)));  // state_k alias
        }
        // Index selectors on a physical record keep the historical state_k label;
        // otherwise a channel carries its column name (state_0 -> position).
        const bool numeric=std::holds_alternative<double>(selector.storage());
        Signal result{numeric&&!named?"state_"+std::to_string(k):names[k],names[k],"none",{},{}};
        result.values.reserve(source.samples);
        for (std::size_t i=0;i<source.samples;++i) result.values.push_back(source.values[i*source.dimension+k]);
        return result;
    }
private:
    std::string join() const { std::string s; for (const auto& n:names) s+=(s.empty()?"":", ")+n; return s; }
    const SourceData& source;
    std::vector<std::string> names;
    bool named{};
};
}

SourceData hankelTimeWindow(const SourceData& full,const J& settings,std::size_t& first) {
    first=0;
    const auto& window=settings.at("time_window");
    if (window.isNull()) return full;
    if (!std::holds_alternative<J::Array>(window.storage()) || window.array().size()!=2)
        throw std::invalid_argument("Hankel time_window must be [start, stop]");
    const double start=window.array()[0].number(), stop=window.array()[1].number();
    if (!(stop>start)) throw std::invalid_argument("Hankel time_window needs stop > start");
    std::size_t last=0; bool any=false;
    for (std::size_t i=0;i<full.samples;++i) if (full.time[i]>=start && full.time[i]<=stop) {
        if (!any) first=i; last=i; any=true; }
    if (!any || last-first<2) throw std::invalid_argument("Hankel time_window leaves fewer than three samples");
    SourceData rows=full;
    const auto n=last-first+1;
    rows.time.assign(full.time.begin()+std::ptrdiff_t(first),full.time.begin()+std::ptrdiff_t(last+1));
    rows.values.assign(full.values.begin()+std::ptrdiff_t(first*full.dimension),
                       full.values.begin()+std::ptrdiff_t((last+1)*full.dimension));
    if (!full.inputs.empty())
        rows.inputs.assign(full.inputs.begin()+std::ptrdiff_t(first*full.inputDimension),
                           full.inputs.begin()+std::ptrdiff_t((last+1)*full.inputDimension));
    rows.samples=n;
    return rows;
}

HankelAnalysis prepareHankel(const SourceData& source,const J& settings,bool paired,bool windowHandled) {
    if (!windowHandled && !settings.at("time_window").isNull())
        throw std::invalid_argument("time_window is supported only by hdmd among the Hankel methods");
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
    Columns resolver(source);
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
        if (entry.second) signals.push_back(resolver.resolve(item,selector,order));
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
