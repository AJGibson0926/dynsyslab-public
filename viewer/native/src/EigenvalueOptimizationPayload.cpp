#include "EigenvalueOptimizationPayload.hpp"

#include <cmath>
#include <stdexcept>

namespace dsl {
namespace {

constexpr double kPi=3.14159265358979323846;

Trace spectrumTrace(const dynsyslab::core::NpyArray& array,
                    std::string label,std::string color){
    Trace trace;trace.label=std::move(label);trace.color=std::move(color);
    trace.marker=true;
    for(const auto value:array.complexDoubles()){
        trace.x.push_back(value.real());trace.y.push_back(value.imag());
    }
    return trace;
}

Figure spectrumFigure(const dynsyslab::core::NpzRecord& record){
    Figure figure{"eigenvalue_optimization_spectrum",
        "Initial and optimized discrete spectra (unitless)",
        "Re(mu) (unitless)","Im(mu) (unitless)",true};
    Trace circle;circle.label="unit circle";circle.color="#64748b";
    for(int point=0;point<=360;++point){const double angle=2.0*kPi*point/360.0;
        circle.x.push_back(std::cos(angle));circle.y.push_back(std::sin(angle));}
    figure.traces.push_back(std::move(circle));
    figure.traces.push_back(spectrumTrace(record.at("initial_discrete_eigenvalues"),
        "initial spectrum","#94a3b8"));
    figure.traces.push_back(spectrumTrace(record.at("discrete_eigenvalues"),
        "optimized spectrum","#2dd4bf"));
    return figure;
}

Figure residualFigure(const dynsyslab::core::NpyArray& metrics){
    const auto values=metrics.doubles();Trace residual;residual.label="relative residual";
    residual.color="#f59e0b";residual.marker=true;
    for(std::size_t output=0;output<metrics.shape[0];++output){
        residual.x.push_back(static_cast<double>(output));
        residual.y.push_back(values[output*metrics.shape[1]+1]);
    }
    Figure figure{"eigenvalue_optimization_residual",
        "Variable-projection fit residual by observable output",
        "observable output","relative L2",false};
    figure.traces.push_back(std::move(residual));return figure;
}

Figure effortFigure(const dynsyslab::core::NpyArray& metrics){
    const auto values=metrics.doubles();Trace iterations,evaluations;
    iterations.label="iterations";iterations.color="#a78bfa";iterations.marker=true;
    evaluations.label="objective evaluations";evaluations.color="#38bdf8";
    evaluations.marker=true;
    for(std::size_t output=0;output<metrics.shape[0];++output){
        const auto base=output*metrics.shape[1];
        iterations.x.push_back(static_cast<double>(output));
        iterations.y.push_back(values[base+2]);
        evaluations.x.push_back(static_cast<double>(output));
        evaluations.y.push_back(values[base+3]);
    }
    Figure figure{"eigenvalue_optimization_effort","Optimization effort",
        "observable output","count",false};
    figure.traces.push_back(std::move(iterations));
    figure.traces.push_back(std::move(evaluations));return figure;
}

}  // namespace

Payload buildEigenvalueOptimizationPayload(
    const std::filesystem::path& source,const std::string& view,
    const dynsyslab::core::NpzRecord& record,
    const dynsyslab::core::json::Value& metadata){
    for(const auto* name:{"initial_discrete_eigenvalues","discrete_eigenvalues",
            "optimization_metrics"})if(!record.contains(name))
        throw std::runtime_error("archive lacks native eigenvalue optimization diagnostics");
    const auto& metrics=record.at("optimization_metrics");
    if(metrics.shape.size()!=2||metrics.shape[1]!=5)
        throw std::runtime_error("optimization metric matrix is invalid");
    Payload payload{source.u8string(),view,
        metadata.at("label").string(source.stem().u8string()),2};
    payload.figures.push_back(spectrumFigure(record));
    payload.figures.push_back(residualFigure(metrics));
    payload.figures.push_back(effortFigure(metrics));
    const auto values=metrics.doubles();
    for(std::size_t output=0;output<metrics.shape[0];++output)
        if(values[output*metrics.shape[1]+4]<0.5)payload.warnings.push_back(
            "Observable output "+std::to_string(output)+" did not meet the optimizer convergence test.");
    payload.warnings.push_back(
        "This archive contains eigenvalues only; construct eigenfunctions with SINDy or fitted boundary functions.");
    return payload;
}

}  // namespace dsl
