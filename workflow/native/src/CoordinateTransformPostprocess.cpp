#include "PostprocessInternal.hpp"
#include "ObservationSource.hpp"
#include "KromDatasetSupport.hpp"
#include "CoordinateChangeKinds.hpp"
#include "CoordinateDiagnostics.hpp"
#include "CoordinateInverse.hpp"
#include "DerivativeSettings.hpp"
#include <dynsyslab/core/AcousticReceiver.hpp>
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
core::Observable receiverMap(const SourceData& source, const json::Value& settings, json::Value& meta) {
    const auto& spec=settings.at("receiver");
    const auto& pos=spec.at("position_m").array();
    if(pos.size()!=3)throw std::invalid_argument("receiver.position_m needs three offsets in metres");
    const double missing=std::numeric_limits<double>::quiet_NaN();
    const core::AcousticReceiver receiver{{pos[0].number(missing),pos[1].number(missing),pos[2].number(missing)},
        spec.at("sound_speed_m_s").number(missing)};
    const auto system=core::createSystem(systemConfiguration(source));
    const auto catalog=core::rpeReceiverPressures(*system,receiver);
    const auto name=settings.at("observable").string();
    const auto& map=catalog.at(name);
    const double range=receiver.rangeMetres();
    const double radius0=system->configuration().parameters.at("equilibrium_radius_m");
    if(source.dimension!=2)throw std::invalid_argument("Receiver pressures require RPE radius and velocity");
    for(std::size_t row=0;row<source.samples;++row) {
        const double radius=radius0*source.values[row*2];
        if(!(radius>0 && radius<range))
            throw std::invalid_argument("Receiver must be outside the bubble at every recorded state");
    }
    auto saved=spec;
    saved.object()["term"]=name;
    saved.object()["range_m"]=range;
    saved.object()["travel_time_s"]=receiver.travelTimeSeconds();
    saved.object()["clock"]="emission-time";
    meta.object()["acoustic_receiver"]=std::move(saved);
    return map;
}
}
namespace {
// kind "measured_derivatives": measured channels, optional measurement noise, and their time
// derivatives from a noise-robust differentiator. Coordinates are [y_1, y_1', y_2, y_2', ...] and
// coordinate_derivatives [y_1', y_1'', ...] (the same differentiator applied twice), so coordinate
// SINDy can use gamma_method "measured". Channels are catalog observables evaluated on a physical
// source (the simulated sensor) or columns of a coordinate archive.
core::CoordinateSeries measuredDerivativeCoordinates(const SourceData& source, const json::Value& settings,
                                                     json::Value& meta, bool physicalSource) {
    namespace c = dynsyslab::core;
    const std::size_t n = source.samples;
    if (n < 5) throw std::invalid_argument("measured_derivatives needs at least 5 samples");
    const double dt = source.time[1] - source.time[0];
    for (std::size_t i = 2; i < n; ++i)
        if (std::abs((source.time[i] - source.time[i - 1]) - dt) > 1e-9 * std::max(1.0, std::abs(dt)))
            throw std::invalid_argument("measured_derivatives needs uniformly sampled data");
    std::vector<std::vector<double>> channels, clean; std::vector<std::string> names;
    std::vector<std::vector<double>> reference;
    std::vector<std::string> observables;
    if (!settings.at("observable").isNull()) observables.push_back(settings.at("observable").string());
    for (const auto& v : settings.at("observables").array()) observables.push_back(v.string());
    if (!observables.empty()) {
        if (!physicalSource) throw std::invalid_argument("measured_derivatives observables need a physical state source; use columns");
        const auto catalog = c::createSystem(systemConfiguration(source))->observableCatalog();
        auto evaluate = [&](const std::string& name) {
            const auto& map = catalog.at(name);
            if (map.outputDimension() != 1) throw std::invalid_argument("measured_derivatives channels must be scalar: " + name);
            std::vector<double> y(n);
            for (std::size_t i = 0; i < n; ++i) y[i] = map.evaluateReal(&source.values[i * source.dimension]).front();
            return y;
        };
        for (const auto& name : observables) { channels.push_back(evaluate(name)); names.push_back(name); }
        for (const auto& v : settings.at("reference_derivatives").array()) reference.push_back(evaluate(v.string()));
    } else {
        for (const auto& v : settings.at("columns").array()) {
            const double k = v.number(-1);
            if (!(k >= 0) || k >= double(source.dimension) || k != std::floor(k)) throw std::invalid_argument("measured_derivatives column out of range");
            std::vector<double> y(n); for (std::size_t i = 0; i < n; ++i) y[i] = source.values[i * source.dimension + std::size_t(k)];
            channels.push_back(y); names.push_back("q_" + std::to_string(std::size_t(k) + 1));
        }
    }
    if (channels.empty()) throw std::invalid_argument("measured_derivatives needs observable(s) or columns");
    if (!reference.empty() && reference.size() != channels.size())
        throw std::invalid_argument("reference_derivatives needs one entry per channel");
    const auto derivative = parseDerivativeSpec(settings.at("derivative"));
    const auto noiseBase = parseNoiseSpec(settings.at("noise"));
    const auto valuesMode = settings.at("values").string("smoothed");
    if (valuesMode != "smoothed" && valuesMode != "measured") throw std::invalid_argument("measured_derivatives values must be smoothed or measured");
    const std::size_t m = channels.size();
    std::vector<double> values(n * 2 * m), rates(n * 2 * m);
    json::Value::Array noiseStd, derivativeError, valueError;
    for (std::size_t k = 0; k < m; ++k) {
        auto noise = noiseBase; noise.seed += 1000003u * k;
        const auto noisy = c::addMeasurementNoise(channels[k], dt, noise, settings.at("noise_scale").number(0));
        const auto first = c::differentiate(noisy, dt, derivative);
        const auto second = c::differentiate(first.derivative, dt, derivative);
        const auto& y = valuesMode == "measured" ? noisy : first.values;
        for (std::size_t i = 0; i < n; ++i) {
            values[i * 2 * m + 2 * k] = y[i]; values[i * 2 * m + 2 * k + 1] = first.derivative[i];
            rates[i * 2 * m + 2 * k] = first.derivative[i]; rates[i * 2 * m + 2 * k + 1] = second.derivative[i];
        }
        noiseStd.emplace_back(c::estimateWhiteNoiseStd(noisy));
        auto relative = [&](const std::vector<double>& a, const std::vector<double>& b) {
            const std::size_t skip = n / 50; double e = 0, r = 0;
            for (std::size_t i = skip; i + skip < n; ++i) { e += (a[i] - b[i]) * (a[i] - b[i]); r += b[i] * b[i]; }
            return r > 0 ? std::sqrt(e / r) : 0.0;
        };
        valueError.emplace_back(relative(y, channels[k]));
        if (!reference.empty()) derivativeError.emplace_back(relative(first.derivative, reference[k]));
    }
    json::Value record = json::Value::Object{};
    record.object()["noise"] = settings.at("noise"); record.object()["derivative"] = settings.at("derivative");
    record.object()["values"] = valuesMode; record.object()["noise_std_estimate"] = std::move(noiseStd);
    record.object()["relative_value_error"] = std::move(valueError);
    if (!reference.empty()) record.object()["relative_derivative_error"] = std::move(derivativeError);
    meta.object()["measured_derivative"] = std::move(record);
    meta.object()["derivative_evidence"] = "measured/" + derivative.method;
    if (settings.at("names").isNull()) {
        json::Value::Array out;
        for (const auto& name : names) { out.emplace_back(name); out.emplace_back(name + "_dot"); }
        meta.object()["coordinate_names"] = std::move(out);
    }
    return {c::TimeSeries(source.time, values, 2 * m), c::TimeSeries(source.time, rates, 2 * m)};
}

std::filesystem::path runCoordinateStep(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    namespace c=dynsyslab::core;
    const auto kind=settings.at("kind").string("identity");
    const bool modelRates=settings.at("system_derivatives").boolean(false);
    const bool integral=kind=="window_integral";
    const bool mapped=!settings.at("observable").isNull()||!settings.at("definition").isNull();
    if(integral&&modelRates)
        throw std::invalid_argument("Window integration produces coordinate values only");
    if(integral&&mapped&&!settings.at("measurement_array").isNull())
        throw std::invalid_argument("Choose either a measured array or an analytic observable for window integration");
    const auto original=source.metadataJson.empty()?json::Value{}:json::parse(source.metadataJson);
    const bool receiver=!settings.at("receiver").isNull();
    if(receiver&&((kind!="observable"&&!integral)||settings.at("observable").string().empty()||
       !settings.at("definition").isNull()||!settings.at("measurement_array").isNull()||
       !original.at("coordinate_provider").isNull()))
        throw std::invalid_argument("A receiver requires one named acoustic term and a physical state source");
    if(integral&&mapped&&!original.at("coordinate_provider").isNull())
        throw std::invalid_argument("Use the scalar measurement array on a transformed trajectory");
    if(kind=="observable"&&!receiver&&settings.at("definition").isNull()&&!original.at("coordinate_provider").isNull())
        throw std::invalid_argument("A catalog observable maps the physical state; put this step first in the chain or use a column of the input");
    if(modelRates&&datasetDynamics(source)!="unforced")
        throw std::invalid_argument("Autonomous derivative reference requires confirmed unforced data");
    const bool physicalSource=original.at("coordinate_provider").isNull();
    // History coordinates of an earlier coordinate archive take its stored exact rates.
    if(modelRates&&!physicalSource&&kind!="history")
        throw std::invalid_argument("System derivatives require physical-state source data");
    const auto q=c::TimeSeries::view(source.time,source.values.data(),source.samples,source.dimension);
    std::optional<c::TimeSeries> rates;
    if(modelRates&&physicalSource)rates=c::TimeSeries(source.time,exactVelocities(source,false),source.dimension);
    auto result=c::CoordinateSeries{q,rates};
    auto meta=json::parse(metadata(source,"coordinate-transform",label,""));
    meta.object()["coordinate_provider"]="measured";
    meta.object()["value_space"]="observation-coordinates";
    {   // keep forcing provenance through coordinate changes so later stages can check data roles
        const auto dynamics=datasetDynamics(source);
        meta.object()["forcing_type"]=dynamics=="unforced"?"free":dynamics=="forced"?"forced":"unknown";
    }
    meta.object()["derivative_evidence"]=modelRates?"system-rhs/chain-rule":"none";
    if(!settings.at("names").isNull())meta.object()["coordinate_names"]=settings.at("names");
    if(kind=="channels")result=channelCoordinates(source,settings,meta);
    if(kind=="columns")result=selectedColumns(source,settings,original,meta);
    if(kind=="measured_derivatives")result=measuredDerivativeCoordinates(source,settings,meta,physicalSource);
    std::optional<c::Observable> analyticMap;  // for the closure diagnostic
    if(kind=="observable"||(integral&&mapped)) {
        auto map=[&] {
            if(receiver)return receiverMap(source,settings,meta);
            if(!settings.at("definition").isNull())return c::Observable::fromExpressionDefinition(settings.at("definition"));
            const auto catalog=c::createSystem(systemConfiguration(source))->observableCatalog();
            const auto name=settings.at("observable").string();
            return catalog.coordinatePresets().count(name)?catalog.coordinatePreset(name):catalog.at(name);
        }();
        result=c::transformCoordinates(q,map,rates?&*rates:nullptr);
        meta.object()["coordinate_definition"]=map.expressionDefinition();
        if(kind=="observable"&&!receiver&&original.at("coordinate_provider").isNull())analyticMap=map;
        // Name the analytic map so model fits can re-attach its exact Jacobian.
        if(!receiver&&settings.at("definition").isNull()&&kind=="observable")
            meta.object()["coordinate_observable"]=settings.at("observable").string();
        const auto& outputs=map.descriptor().outputNames;
        if(settings.at("names").isNull()&&outputs.size()==result.coordinates.componentCount()&&kind=="observable") {
            json::Value::Array names,labels;
            const auto catalog=c::createSystem(systemConfiguration(source))->observableCatalog();
            for(const auto& n:outputs) {
                names.emplace_back(n);
                std::string label;  // display label (TeX) beside the identifier
                try{label=catalog.at(n).descriptor().label;}catch(const std::exception&){}
                labels.emplace_back(label.empty()?n:label);
            }
            meta.object()["coordinate_names"]=std::move(names);
            meta.object()["coordinate_labels"]=std::move(labels);
        }
    }
    if(integral) {
        const double duration=settings.at("window_duration").number();
        result.coordinates=c::windowIntegralCoordinates(result.coordinates,duration);
        json::Value::Object definition{{"kind","window_integral"},{"window_duration",duration},
            {"quadrature","piecewise-linear"},{"alignment","trailing-complete"}};
        if(mapped) definition["observable"]=meta.at("coordinate_definition");
        else {
            definition["measurement_array"]=settings.at("measurement_array").string("values");
            for(const auto key:{"coordinate_definition","coordinate_transform_identity","coordinate_names"})
                if(!original.at(key).isNull())definition[key]=original.at(key);
        }
        for(const auto key:{"measurement_units","measurement_time_units","model_time_units"})
            if(!original.at(key).isNull())definition[key]=original.at(key);
        meta.object()["coordinate_provider"]="causal-history";
        meta.object()["value_space"]="history-coordinates";
        meta.object()["coordinate_definition"]=definition;
        meta.object()["coordinate_transform_identity"]=json::dump(definition);
        meta.object()["history_span"]=duration;
        if(settings.at("names").isNull())
            meta.object()["coordinate_names"]=json::Value::Array{"signal","window_integral"};
    } else if(kind=="history"&&settings.at("transform_archive").string().empty()) {
        // The history basis is part of this coordinate change, never a separate stage.
        auto fitted=fitHistoryCoordinates(source,settings,root,label);
        meta.object()["derivative_evidence"]=modelRates?"frozen-projection-of-exact-measurement-rate":"none";
        result=std::move(fitted.series);
        meta.object()["coordinate_provider"]="causal-history";
        meta.object()["value_space"]="history-coordinates";
        meta.object()["coordinate_transform_identity"]=fitted.transform.identity();
        meta.object()["coordinate_transform_archive"]=std::filesystem::absolute(fitted.transformArchive).generic_u8string();
        meta.object()["history_span"]=double(fitted.transform.basis().delays-1)*fitted.transform.sampleInterval();
        meta.object()["warmup_samples"]=double(fitted.transform.basis().delays-1);
        meta.object()["controller_kind"]="none";
        meta.object()["input_history_source"]=std::filesystem::absolute(source.path).generic_u8string();
        json::Value::Array names,labels;
        for(std::size_t j=1;j<=result.coordinates.componentCount();++j) {
            names.emplace_back("xi_"+std::to_string(j));
            labels.emplace_back("\\xi_{"+std::to_string(j)+"}");
        }
        if(settings.at("names").isNull())meta.object()["coordinate_names"]=std::move(names);
        meta.object()["coordinate_labels"]=std::move(labels);
    } else if(kind=="history") {
        const auto path=historyTransformArchive(std::filesystem::u8path(settings.at("transform_archive").string()));
        const auto transform=c::loadFrozenHistoryTransform(path);
        const auto observable=settings.at("observable").string(settings.at("measurement").string(transform.measurement()));
        if(physicalSource&&observable!=transform.measurement()) {
            const auto system=c::createSystem(systemConfiguration(source));
            if(system->observableCatalog().at(observable).descriptor().name!=
                system->observableCatalog().at(transform.measurement()).descriptor().name)
                throw std::invalid_argument("Measurement differs from the frozen history transform");
        }
        for(const auto& pair:std::vector<std::pair<std::string,std::string>>{
            {"measurement_time_units",transform.timeUnits()},{"measurement_units",transform.observationUnits()}})
            if(!original.at(pair.first).isNull()&&original.at(pair.first).string()!=pair.second)
                throw std::invalid_argument("History measurement units differ from the frozen transform");
        std::optional<HistoryMeasurement> column;  // a column of an earlier coordinate change
        if(!physicalSource) {
            auto pick=settings; if(pick.at("measurement").isNull())pick.at("measurement")=transform.measurement();
            column=historyMeasurement(source.path,pick,modelRates);
            if(column->name!=transform.measurement())
                throw std::invalid_argument("Measurement column differs from the frozen history transform");
        }
        const auto measured=column?column->series
            :loadObservationSource(source.path,settings.at("measurement_array").string("echo_values"),observable);
        auto measurement=c::TimeSeries::view(measured.time,measured.values.data(),measured.samples,measured.dimension);
        std::optional<c::TimeSeries> measurementRates;
        if(column&&column->rates)measurementRates=c::TimeSeries(measured.time,*column->rates,1);
        else if(modelRates) {
            const auto map=c::createSystem(systemConfiguration(source))->observableCatalog().at(observable);
            // Exact reference is meaningful only when samples are this analytic observable.
            const auto analytic=c::transformCoordinates(q,map,&*rates);
            if(analytic.coordinates.componentCount()!=1||measurement.componentCount()!=1)
                throw std::invalid_argument("History transform requires one scalar observable");
            for(std::size_t k=0;k<measurement.sampleCount();++k)
                if(std::abs(measurement.value(k,0)-analytic.coordinates.value(k,0))>
                    1e-10*std::max(1.,std::abs(analytic.coordinates.value(k,0))))
                    throw std::invalid_argument("Measured signal differs from requested analytic derivative reference");
            measurementRates=analytic.derivatives;
        }
        result=c::transformHistoryCoordinates(measurement,transform,measurementRates?&*measurementRates:nullptr);
        meta.object()["coordinate_provider"]="causal-history";
        meta.object()["value_space"]="history-coordinates";
        meta.object()["coordinate_transform_identity"]=transform.identity();
        meta.object()["coordinate_transform_archive"]=std::filesystem::absolute(path).generic_u8string();
        meta.object()["history_span"]=double(transform.basis().delays-1)*transform.sampleInterval();
        meta.object()["derivative_evidence"]=modelRates?"frozen-projection-of-analytic-measurement-rate":"none";
        meta.object()["warmup_samples"]=double(transform.basis().delays-1);
        json::Value::Array names,labels;
        for(std::size_t j=1;j<=result.coordinates.componentCount();++j) {
            names.emplace_back("xi_"+std::to_string(j));labels.emplace_back("\\xi_{"+std::to_string(j)+"}");}
        if(settings.at("names").isNull())meta.object()["coordinate_names"]=std::move(names);
        meta.object()["coordinate_labels"]=std::move(labels);
    } else if(kind!="identity"&&kind!="observable"&&kind!="channels"&&kind!="columns"&&kind!="measured_derivatives")throw std::invalid_argument("Unknown coordinate transform kind");
    const auto& z=result.coordinates;const auto n=z.sampleCount(),d=z.componentCount();
    c::NpzRecord archive;
    archive.arrays["time"]=c::NpyArray::fromDoubles({n},{z.timeData(),z.timeData()+n});
    archive.arrays["values"]=c::NpyArray::fromDoubles({n,d},{z.valuesData(),z.valuesData()+n*d});
    const auto states=alignedSourceStates(source,z);
    if(states)
        archive.arrays["source_states"]=c::NpyArray::fromDoubles({n,states->size()/std::max<std::size_t>(n,1)},*states);
    appendCoordinateDiagnostics(source,settings,result,states,analyticMap?&*analyticMap:nullptr,meta,archive);
    // Retain the inverse when one exists (see CoordinateInverse.hpp).
    if(kind=="identity")recordCoordinateInverse(source,"exact","identity",z,meta,archive);
    else if(kind=="history")recordCoordinateInverse(source,"readout","",z,meta,archive);
    else {
        const auto& found=meta.at("coordinate_diagnostics").at("crossings").at("noninjective_witnesses");
        recordCoordinateInverse(source,"none",kind=="columns"?"a column subset discards the other coordinates"
            :!found.isNull()&&found.number()>0?"the diagnostics found distinct states with equal coordinates (not one-to-one)"
            :"no inverse is implemented for this map",z,meta,archive);
    }
    // Measured coordinates of a forced record keep its inputs, so input_map can regress on them.
    if(kind=="measured_derivatives"&&source.inputDimension&&source.inputs.size()==n*source.inputDimension)
        archive.arrays["control_values"]=c::NpyArray::fromDoubles({n,source.inputDimension},source.inputs);
    if(result.derivatives)archive.arrays["coordinate_derivatives"]=c::NpyArray::fromDoubles({n,d},
        {result.derivatives->valuesData(),result.derivatives->valuesData()+n*d});
    // Optional display labels (TeX), one per coordinate, e.g. ["e","e'"].
    if(const auto& labels=settings.at("labels");!labels.isNull()) {
        if(labels.array().size()!=result.coordinates.componentCount())
            throw std::invalid_argument("Coordinate labels need one entry per coordinate");
        meta.object()["coordinate_labels"]=labels;
    }
    meta.object()["coordinate_transform_kind"]=kind;
    if(!settings.at("chain_inputs").isNull())meta.object()["coordinate_chain"]=settings.at("chain_inputs");
    meta.object()["coordinate_dimension"]=double(d);
    meta.object()["delta_t"]=n>1&&z.uniformTimeSteps()?z.deltaTime():0;
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return publish(root,source,"coordinate-transform",label,std::move(archive));
}
}  // namespace

// One coordinate change, or a chain of them: each step reads the previous
// step's output, so (r,v) -> (e,e') -> xi(e) is one stage. Intermediate
// coordinates are archived for provenance; the stage's output is the last.
std::filesystem::path runCoordinateTransform(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    if(settings.at("kind").string()!="chain")return runCoordinateStep(source,settings,root,label);
    const auto& steps=settings.at("steps").array();
    if(steps.empty())throw std::invalid_argument("A coordinate chain needs a nonempty steps array");
    std::optional<SourceData> current;
    json::Value::Array produced;
    std::filesystem::path output;
    for(std::size_t k=0;k<steps.size();++k) {
        if(steps[k].at("kind").string()=="chain")throw std::invalid_argument("Nest coordinate chains by listing their steps");
        auto step=steps[k];
        if(!produced.empty())step.at("chain_inputs")=produced;
        // Chain-level diagnostics apply to every step that sets none of its own.
        if(step.at("diagnostics").isNull()&&!settings.at("diagnostics").isNull())step.at("diagnostics")=settings.at("diagnostics");
        output=runCoordinateStep(current?*current:source,step,root,label+" / step "+std::to_string(k+1));
        produced.emplace_back(std::filesystem::absolute(output).generic_u8string());
        if(k+1<steps.size())current=loadSource(output);
    }
    return output;
}
}
