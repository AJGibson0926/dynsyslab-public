#include "TrajectoryPresentation.hpp"
#include "ArchiveSystem.hpp"
#include <algorithm>
#include <cmath>
namespace dsl {
namespace {
double positive(const json::Value& value,const char* key) {
    const auto number=value.at(key).number();
    return std::isfinite(number)&&number>0?number:0;
}
}
void decorateTrajectory(Payload& payload,const json::Value& metadata,
    const DisplayCoordinateTransform& current,const std::vector<double>& time) {
    const auto dimension=current.stateLabels.size();
    const auto other=displayCoordinateTransform(metadata,dimension,
        current.coordinates=="dimensionless"?"dimensional":"dimensionless");
    const auto system=metadata.at("system").string(metadata.at("system_configuration").at("identifier").string());
    const bool bubble=system=="rpe"||system=="rayleigh-plesset"||system=="rayleigh_plesset";
    double inputScale=positive(metadata,"control_input_scale_to_dimensional");
    if(!inputScale) inputScale=positive(metadata.at("control_history"),"scale_to_dimensional");
    if(!inputScale) inputScale=positive(metadata.at("nondimensionalization"),"control_pressure_scale_pa");
    if(!inputScale) inputScale=positive(metadata,"control_scale_to_dimensional");
    const auto controlCoordinates=metadata.at("control_input_coordinate_system").string(
        metadata.at("control_history").at("coordinate_system").string(
        metadata.at("control_coordinate_system").string(current.coordinates)));
    std::size_t spectrumComponent=0;
    for(auto& figure:payload.figures) {
        const bool spectrum=figure.kind=="fourier"||figure.kind=="control_fourier";
        const bool control=figure.title=="Control history"||figure.kind=="control_fourier";
        const bool phase=figure.kind=="phase";
        const bool spatial=phase||figure.kind=="three_d";
        const auto found=std::find(current.stateLabels.begin(),current.stateLabels.end(),figure.yLabel);
        const auto component=spectrum?spectrumComponent++:found==current.stateLabels.end()?std::size_t(0):
            static_cast<std::size_t>(found-current.stateLabels.begin());
        figure.coordinates=current.coordinates;
        if(control) {
            if(controlCoordinates!=current.coordinates&&inputScale>0) {
                const auto scale=current.coordinates=="dimensional"?inputScale:1/inputScale;
                for(auto& trace:figure.traces) for(auto& value:trace.y) value*=scale;
            } else if(controlCoordinates!=current.coordinates) figure.coordinates=controlCoordinates;
            figure.yLabel=spectrum?"control amplitude":bubble&&figure.coordinates=="dimensional"?"u\\; (Pa)":"u";
        }
        if(spatial&&!figure.traces.empty()) {
            figure.square=false;
            figure.traces.front().encoding="time";
            figure.traces.front().colorValues=time;
            figure.colorLabel=current.timeLabel;
        }
        if(spatial&&!figure.traces.empty()&&!figure.traces.front().x.empty()&&!figure.traces.front().y.empty()) {
            const auto& trace=figure.traces.front();
            figure.references.push_back({trace.x.front(),trace.y.front(),"IC",false});
            if(phase) appendVerifiedEquilibria(figure,metadata);
        }
        if(other.coordinates==current.coordinates||(control&&!inputScale)) continue;
        Figure::CoordinateAxes alternate;
        alternate.name=other.coordinates;
        alternate.timeScale=other.timeMultiplier/current.timeMultiplier;
        alternate.xScale=spatial?other.stateMultipliers[0]/current.stateMultipliers[0]:alternate.timeScale;
        if(spectrum) alternate.xScale=1/alternate.timeScale;
        alternate.yScale=control?(other.coordinates=="dimensional"?inputScale:1/inputScale):
            other.stateMultipliers[component]/current.stateMultipliers[component];
        if(figure.kind=="three_d"&&dimension>=3)
            alternate.zScale=other.stateMultipliers[2]/current.stateMultipliers[2];
        alternate.xLabel=spatial?other.stateLabels[0]:other.timeLabel;
        alternate.yLabel=control?(bubble&&other.coordinates=="dimensional"?"u\\; (Pa)":"u"):other.stateLabels[component];
        if(figure.kind=="three_d"&&dimension>=3) alternate.zLabel=other.stateLabels[2];
        if(spectrum) {
            alternate.xLabel=other.coordinates=="dimensionless"?"frequency (cycles / tau)":
                positive(metadata,"time_scale_to_seconds")>0?"frequency (Hz)":"frequency (cycles / time)";
            alternate.yLabel=figure.yLabel;
        }
        for(const auto& trace:figure.traces)
            alternate.traceLabels.push_back(!spatial&&!control&&trace.role=="data"
                ?other.stateLabels[component]:trace.label);
        figure.alternateCoordinates=std::move(alternate);
    }
}
}
