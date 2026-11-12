#include "KromControlCampaignPayload.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dsl {
namespace {

std::vector<double> column(const dynsyslab::core::NpyArray& array,std::size_t index){
    if(array.shape.size()!=2||index>=array.shape[1])throw std::runtime_error(
        "campaign metric matrix is malformed");const auto source=array.doubles();
    std::vector<double> result(array.shape[0]);for(std::size_t row=0;row<array.shape[0];++row)
        result[row]=source[row*array.shape[1]+index];return result;}
Trace trace(std::string name,const std::vector<double>& x,std::vector<double> y,
    std::string color,std::string style="solid"){Trace result;result.label=std::move(name);
    result.x=x;result.y=std::move(y);result.color=std::move(color);result.lineStyle=std::move(style);
    result.marker=true;return result;}
Figure metricFigure(const std::string& title,const std::string& yLabel,
    const std::vector<double>& attempts,const dynsyslab::core::NpyArray& matrix){
    Figure figure{"krom_control_campaign",title,"attempt",yLabel,false};
    const char* names[]{"exact B_psi","constant B","bilinear B"};
    const char* colors[]{"#38bdf8","#f59e0b","#22c55e"};
    for(std::size_t index=0;index<3;++index)
        figure.traces.push_back(trace(names[index],attempts,column(matrix,index),colors[index]));
    return figure;}
void threshold(Figure& figure,const std::vector<double>& attempts,double value,
    const std::string& name){if(attempts.empty()||!std::isfinite(value))return;
    figure.traces.push_back(trace(name,{attempts.front(),attempts.back()},{value,value},
        "#ef4444","dash"));}
std::vector<double> data(const dynsyslab::core::NpzRecord& record,const char* name){
    if(!record.contains(name))throw std::runtime_error(std::string("campaign archive is missing ")+name);
    return record.at(name).doubles();}

}  // namespace

Payload buildKromControlCampaignPayload(const std::filesystem::path& source,
    const std::string& view,const dynsyslab::core::NpzRecord& record,const json::Value& metadata){
    const auto attempts=data(record,"attempt_number");Payload payload{source.u8string(),view,
        metadata.at("label").string(source.stem().u8string()),2};
    auto release=metricFigure("Target error at control release","state error",attempts,
        record.at("release_error"));threshold(release,attempts,metadata.at("success_policy")
        .at("release_error_tolerance").number(std::numeric_limits<double>::quiet_NaN()),"success tolerance");
    payload.figures.push_back(std::move(release));auto finalError=metricFigure(
        "Final error after release","state error",attempts,record.at("final_error"));
    threshold(finalError,attempts,metadata.at("success_policy").at("final_error_tolerance")
        .number(std::numeric_limits<double>::quiet_NaN()),"success tolerance");
    payload.figures.push_back(std::move(finalError));payload.figures.push_back(metricFigure(
        "Peak control demand","|u| maximum",attempts,record.at("control_peak")));
    Figure maps{"krom_control_campaign","Mean B_psi approximation error","attempt",
        "mean relative matrix error",false};maps.traces.push_back(trace("constant B",attempts,
        data(record,"constant_error_mean"),"#f59e0b"));maps.traces.push_back(trace("bilinear B",
        attempts,data(record,"bilinear_error_mean"),"#22c55e"));payload.figures.push_back(std::move(maps));
    Figure screening{"krom_control_campaign","Candidate pre-screening","attempt","selection score",false};
    screening.traces.push_back(trace("selection score",attempts,data(record,"candidate_selection_score"),
        "#a78bfa"));payload.figures.push_back(std::move(screening));
    Figure modes{"krom_control_campaign","Extracted KROM size","attempt","retained modes",false};
    modes.traces.push_back(trace("mode count",attempts,data(record,"candidate_mode_count"),"#94a3b8"));
    payload.figures.push_back(std::move(modes));const auto success=data(record,"case_success");
    bool bilinear=false;for(std::size_t row=0;row<attempts.size();++row)
        if(success[row*3+2]>0.5)bilinear=true;if(!bilinear)payload.warnings.push_back(
            "No bilinear approximation met both retained success tolerances.");
    return payload;}

}  // namespace dsl
