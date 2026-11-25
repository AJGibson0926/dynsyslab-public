#include "HankelAnalysis.hpp"
#include "HdmdObjects.hpp"
#include "SampledDomain.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <algorithm>
#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <set>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::workflow::postprocess_detail {
namespace {
std::size_t countSetting(const json::Value& settings,const char* name,double fallback,std::size_t minimum=0) {
    const double value=settings.at(name).number(fallback);
    if(!std::isfinite(value)||value<double(minimum)||std::floor(value)!=value
        ||value>=double((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument(std::string(name)+" must be a valid integer count");
    return static_cast<std::size_t>(value);
}
}
std::filesystem::path runHdmd(const SourceData& full, const json::Value& settings,
                              const std::filesystem::path& root, const std::string& label) {
    namespace c = dynsyslab::core;
    // Input: one or more time series on a common clock, and nothing else.
    std::vector<c::TimeSeriesRecord> parts{timeSeriesRecord(full)};
    for (const auto& path : settings.at("additional_time_series").array())
        parts.push_back(timeSeriesRecord(loadSource(std::filesystem::u8path(path.string()))));
    const auto record = c::joinTimeSeriesRecords(parts);
    const auto view = sourceView(record, full);

    c::HankelDMDOptions options;
    // "time_window" [start, stop] restricts the Hankel matrix to that stretch of the record.
    std::size_t first=0;
    const SourceData source=hankelTimeWindow(view,settings,first);
    auto input=prepareHankel(source,settings,true,true);
    if(!settings.at("time_window").isNull())input.provenance.object()["time_window"]=settings.at("time_window");
    options.matrixSize=input.matrix.columns;
    options.samplingWindow=countSetting(settings,"sampling_window",input.matrix.window.delays,1);
    options.snapshotStride=input.matrix.window.snapshotStride;
    options.sampleStride=input.matrix.window.sampleStride;
    options.alignToWindowEnd=input.alignEnd;
    options.channelScales=input.scales;
    options.rank=countSetting(settings,"rank",0);
    options.singularValueTolerance=settings.at("svd_tolerance").number(1.0e-10);
    {   // "rank_rule": "optimal_hard_threshold" (Gavish-Donoho) truncates the noise floor of noisy data
        const auto rule=settings.at("rank_rule").string("tolerance");
        if(rule!="tolerance"&&rule!="optimal_hard_threshold")
            throw std::invalid_argument("hdmd rank_rule must be tolerance or optimal_hard_threshold");
        options.optimalHardThreshold=rule=="optimal_hard_threshold";
    }
    options.terminalAmplitudeEnergySort=settings.at("mode_sort").string("terminal-amplitude-energy")!="none";
    // H-DMD publishes data objects only: the eigenfunctions (with their spectrum and properties),
    // the reconstructions written over them, and on request the Hankel matrix it read.
    c::NpzRecord archive, reconstructions, matrix;
    if(settings.at("save_matrix").boolean()) archiveHankelLayout(matrix,input,true);
    const std::size_t labelOffset=input.alignEnd?input.matrix.ends.at(0)-input.matrix.starts.at(0):0;
    // The fit consumes the matrix data; input.matrix.window and input.provenance keep the layout.
    const auto result=c::fitHankelDMD(std::move(input.matrix),input.names,uniformDeltaTime(source),options);

    const HdmdFitContext context{&record,first,labelOffset,&input,&settings};
    const auto set=hdmdEigenfunctionSet(result,context);
    auto objects=c::writeEigenfunctionSet(archive,set,record.channels);
    archive.arrays["sample_interval"]=c::NpyArray::fromDoubles({1},{result.deltaTime});

    // Errors of the default representation on its default domain (first in each list).
    const auto errors=set.errors({}, {}, &record);
    archive.arrays["eigenfunction_value_error"]=c::NpyArray::fromDoubles({result.rank},errors.valueRelative);
    archive.arrays["eigenfunction_koopman_residual"]=c::NpyArray::fromDoubles({result.rank},errors.koopmanResidual);
    objects.object()["errors"]=json::Value::Object{{"representation",errors.representation},
        {"domain",errors.domain},{"samples",double(errors.samples)},
        {"arrays",json::Value::Array{json::Value("eigenfunction_value_error"),json::Value("eigenfunction_koopman_residual")}}};

    const auto channels=result.observableNames.size();
    reconstructions.arrays["reconstruction_matrix"]=c::NpyArray::fromComplex({channels,result.rank},result.reconstructionMatrix);
    std::string reconstructionNames;
    std::set<std::string> keys;
    for(std::size_t channel=0;channel<channels;++channel) {
        const auto& name=result.observableNames[channel];
        const auto key=safeArrayKey(name);
        if(!keys.insert(key).second)throw std::invalid_argument("observable archive keys collide");
        if(!reconstructionNames.empty())reconstructionNames+="\n";
        reconstructionNames+=name;
        reconstructions.arrays["reconstruction__"+key]=c::NpyArray::fromComplex({1,result.rank},
            std::vector<std::complex<double>>(result.reconstructionMatrix.begin()+channel*result.rank,
                result.reconstructionMatrix.begin()+(channel+1)*result.rank));
    }
    if(!settings.at("reconstruction_observables").array().empty()) {
        // Catalog observables are functions of the state, found through the input record's origin,
        // never a copy carried by H-DMD.
        const auto found=statesAtSamples(record,set.samples().time);
        if(!found)throw std::invalid_argument("reconstruction_observables needs the producing system's states; this input has none");
        const auto& states=found->values;
        auto system=c::createSystem(found->system);
        const auto catalog=system->observableCatalog();
        for(const auto& item:settings.at("reconstruction_observables").array()) {
            const auto name=item.string();const auto key=safeArrayKey(name);
            if(!keys.insert(key).second)throw std::invalid_argument("duplicate reconstruction name");
            const auto& observable=catalog.at(name);
            const auto real=observable.evaluateRealBatch(states,result.samples);
            const std::vector<std::complex<double>> values(real.begin(),real.end());
            const auto fitted=c::fitObservableReconstruction(result.coordinates,
                result.samples,result.rank,values,observable.outputDimension(),
                settings.at("reconstruction_ridge").number());
            reconstructions.arrays["reconstruction__"+key]=c::NpyArray::fromComplex({fitted.outputs,result.rank},fitted.matrix);
            reconstructions.arrays["reconstruction_error__"+key]=c::NpyArray::fromDoubles(
                {3},{fitted.residualRms,fitted.residualMaximum,fitted.relativeResidual});
            reconstructionNames+="\n"+name;
        }
    }
    reconstructions.arrays["reconstruction_names"]=c::NpyArray::fromString(reconstructionNames);
    archive.arrays["hankel_singular_values"]=c::NpyArray::fromDoubles({result.singularValues.size()},result.singularValues);
    archive.arrays["terminal_modal_amplitudes"]=c::NpyArray::fromComplex({result.rank},result.modalAmplitudes);
    std::vector<double> indices; for (auto value : result.sourceModeIndices) indices.push_back(static_cast<double>(value));
    archive.arrays["energy_sort_source_indices"]=c::NpyArray::fromDoubles({result.rank},indices);
    auto meta=json::parse(metadata(source, "hankel-dmd", label,
        "\"eigenvalue_domain\":\"discrete-time\",\"full_spectrum_retained\":"
        + std::string(result.rank==result.singularValues.size()?"true":"false")
        + ",\"mode_sort\":\"" + std::string(options.terminalAmplitudeEnergySort ? "terminal-amplitude-energy" : "none")
        + "\",\"matrix_size\":" + std::to_string(options.matrixSize)
        + ",\"observable\":\"" + escaped(result.observableName) + "\""
        + ",\"eigenfunction_representation\":\"scattered-interpolation\""
        + ",\"gradient_backend\":\"finite-difference\",\"interpolation_neighbors\":"
        + std::to_string(std::size_t(settings.at("interpolation_neighbors").number(8)))
        + ",\"gradient_relative_step\":0.000001"
        + ",\"eigenfunction_provenance\":\"H-DMD eigenfunction objects (dynsyslab.eigenfunction-objects/v1)\""
        + ",\"hankel_window\":"+json::dump(input.provenance)
        + ",\"delta_t\":" + json::dump(json::Value(result.deltaTime))
        + ",\"sampling_window\":" + std::to_string(options.samplingWindow)
        + ",\"matrix_rows\":" + std::to_string(result.matrixRows)
        + ",\"matrix_columns\":" + std::to_string(result.matrixColumns)
        + ",\"hankel_contract\":\"joint-hankel/v1\",\"row_layout\":\"delay-major\""
        + ",\"observables\":"+jsonStringArray(result.observableNames)
        + ",\"channel_scales\":"+jsonArray(options.channelScales)
        + ",\"snapshot_stride\":"+std::to_string(options.snapshotStride)
        + ",\"reference_delay\":"+std::to_string(result.referenceDelay)
        + ",\"delay_count\":"+std::to_string(result.delayCount)
        + ",\"sample_stride\":"+std::to_string(options.sampleStride)
        + ",\"rank\":" + std::to_string(result.rank)));
    meta.object()["eigenfunction_objects"]=std::move(objects);
    // Properties of the eigenfunctions as H-DMD produced them (per eigenfunction, except the input's
    // Hankel singular values, which describe the data the method read).
    meta.object()["object_properties"]=json::Value::Object{
        {"per_eigenfunction",json::Value::Array{"eigenfunction_value_error","eigenfunction_koopman_residual",
                                                 "terminal_modal_amplitudes","energy_sort_source_indices"}},
        {"of_input",json::Value::Array{"hankel_singular_values"}}};
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    const auto eigenfunctions=publish(root, source, "koopman-hankel-dmd", label, std::move(archive));

    auto reconstructionMeta=json::parse(metadata(source,"hankel-dmd-reconstructions",label,""));
    reconstructionMeta.object()["reconstructions"]=c::reconstructionsBlock(eigenfunctions,result.rank);
    reconstructions.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(reconstructionMeta));
    registerAdditionalObject(publish(root,source,"reconstructions",label,std::move(reconstructions)));
    if(!matrix.arrays.empty()){
        auto matrixMeta=json::parse(metadata(source,"hankel-matrix",label,""));
        matrixMeta.object()["hankel_window"]=input.provenance;
        matrix.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(matrixMeta));
        registerAdditionalObject(publish(root,source,"hankel_matrix",label,std::move(matrix)));
    }
    return eigenfunctions;
}


}
