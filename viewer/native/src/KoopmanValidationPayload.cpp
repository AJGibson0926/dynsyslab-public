#include "KoopmanValidationPayload.hpp"

#include <stdexcept>

namespace dsl {
namespace {

Trace metricTrace(const dynsyslab::core::NpyArray& matrix,std::size_t column,
                  std::string label,std::string color,const std::vector<double>& available = {}) {
    if(matrix.shape.size()!=2||column>=matrix.shape[1])
        throw std::invalid_argument("Koopman validation metric matrix is invalid");
    const auto values=matrix.doubles();Trace trace;trace.label=std::move(label);
    trace.color=std::move(color);trace.marker=true;
    if(!available.empty()&&available.size()!=matrix.shape[0])
        throw std::invalid_argument("Koopman metric availability must be mode aligned");
    for(std::size_t row=0;row<matrix.shape[0];++row){
        if(!available.empty()&&available[row]==0)continue;
        trace.x.push_back(static_cast<double>(row));
        trace.y.push_back(values[row*matrix.shape[1]+column]);
    }
    return trace;
}

Figure metricFigure(const dynsyslab::core::NpzRecord& record,std::size_t column,
                    const std::string& title,const std::string& yLabel) {
    Figure figure{"koopman_validation",title,"retained mode",yLabel,false};
    figure.traces.push_back(metricTrace(record.at("one_step_metrics"),column,
        "one step","#2dd4bf"));
    figure.traces.push_back(metricTrace(record.at("anchored_long_term_metrics"),column,
        "anchored long term","#f59e0b"));
    figure.traces.push_back(metricTrace(record.at("generator_metrics"),column,
        "generator","#a78bfa",record.contains("generator_available")?
        record.at("generator_available").doubles():std::vector<double>{}));
    return figure;
}

Figure horizonFigure(const dynsyslab::core::NpzRecord& record) {
    const auto& values=record.at("horizon_relative_l2");
    const auto& horizons=record.at("horizons");
    if(values.shape.size()!=2||horizons.shape.size()!=1
        ||values.shape[1]!=horizons.shape[0])
        throw std::invalid_argument("Koopman validation horizon grid is invalid");
    Trace trace;trace.label="semigroup relative L2";trace.heatmap=true;
    trace.x=horizons.doubles();trace.z=values.doubles();
    trace.zShape={static_cast<int>(values.shape[0]),static_cast<int>(values.shape[1])};
    for(std::size_t mode=0;mode<values.shape[0];++mode)
        trace.y.push_back(static_cast<double>(mode));
    Figure figure{"koopman_validation_horizons","Finite-horizon semigroup error",
        "horizon","retained mode",false};
    figure.traces.push_back(std::move(trace));return figure;
}

}  // namespace

Payload buildKoopmanValidationPayload(
    const std::filesystem::path& source,const std::string& view,
    const dynsyslab::core::NpzRecord& record,
    const dynsyslab::core::json::Value& metadata) {
    for(const auto* name:{"one_step_metrics","anchored_long_term_metrics",
            "generator_metrics","horizons","horizon_relative_l2"})
        if(!record.contains(name))
            throw std::runtime_error("archive lacks native Koopman validation arrays");
    Payload payload{source.u8string(),view,
        metadata.at("label").string(source.stem().u8string()),2};
    payload.figures.push_back(metricFigure(record,3,"Relative semigroup defects","relative L2"));
    payload.figures.push_back(metricFigure(record,5,"Tail pointwise defects","95th percentile relative"));
    payload.figures.push_back(metricFigure(record,1,"Absolute RMS defects","absolute RMS"));
    payload.figures.push_back(horizonFigure(record));
    if(record.contains("generator_available"))
        for(double available:record.at("generator_available").doubles())if(available==0) {
            payload.warnings.push_back("Generator defects are unavailable for modes without finite continuous rates; they are omitted, not zero-error results.");
            break;
        }
    const double outside=metadata.at("outside_training_fraction").number();
    if(outside>0.0)payload.warnings.push_back(
        "Validation includes samples outside the KROM training domain (fraction "
        +std::to_string(outside)+").");
    return payload;
}

}  // namespace dsl
