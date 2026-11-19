#include "ArchiveMatrixPayload.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dsl {
namespace {
using dynsyslab::core::NpzRecord;
namespace json = dynsyslab::core::json;

std::string archiveLabel(const json::Value& metadata, const std::filesystem::path& source) {
    return metadata.at("label").string(source.stem().u8string());
}

Figure errorFigure(const NpzRecord& record, const char* arrayName,
                   const std::string& title) {
    const auto& values = record.at(arrayName);
    if (values.shape.size() != 2) throw std::invalid_argument("control error field must be a matrix");
    Trace trace; trace.label = title; trace.heatmap = true;
    trace.x = record.at("control_error_axis0").doubles();
    trace.y = record.at("control_error_axis1").doubles();
    trace.z = values.doubles();
    trace.zShape = {static_cast<int>(values.shape[0]), static_cast<int>(values.shape[1])};
    Figure figure{"control_map_error", title, "state 0", "state 1", true};
    figure.traces.push_back(std::move(trace)); return figure;
}

}  // namespace

Payload buildControlErrorPayload(const std::filesystem::path& source, const std::string& view,
                            const NpzRecord& record, const json::Value& metadata) {
    if (!record.contains("control_error_axis0") || !record.contains("control_error_axis1"))
        throw std::runtime_error("archive has no native control-map error grid");
    Payload payload{source.u8string(), view, archiveLabel(metadata, source), 2};
    const std::pair<const char*, const char*> fields[]{
        {"constant_absolute_error", "Constant B absolute error"},
        {"constant_relative_error", "Constant B relative error"},
        {"bilinear_absolute_error", "Bilinear B absolute error"},
        {"bilinear_relative_error", "Bilinear B relative error"}};
    for (const auto& field : fields) if (record.contains(field.first))
        payload.figures.push_back(errorFigure(record, field.first, field.second));
    if (payload.figures.empty()) throw std::runtime_error("archive has no control-map error fields");
    return payload;
}

Payload buildObservabilityPayload(const std::filesystem::path& source,const std::string& view,
                             const NpzRecord& record,const json::Value& metadata) {
    const auto& singular=record.at("observability_singular_values");
    const auto& directions=record.contains("observable_directions")?record.at("observable_directions")
        :record.at("observability_directions");
    if (singular.shape.size()!=1 || directions.shape.size()!=2)
        throw std::runtime_error("archive has no native observability analysis");
    Trace spectrum; spectrum.label="observability singular values"; spectrum.marker=true;
    const auto values=singular.doubles();
    for (std::size_t index=0;index<values.size();++index) {
        spectrum.x.push_back(static_cast<double>(index));
        spectrum.y.push_back(std::log10(std::max(values[index],1.0e-300)));
    }
    Figure singularFigure{"observability_spectrum","Observability singular spectrum",
                          "direction","log10(singular value)",false};
    singularFigure.traces.push_back(std::move(spectrum));
    Trace map; map.label="right singular directions"; map.heatmap=true;
    map.z=directions.doubles();
    map.zShape={static_cast<int>(directions.shape[0]),static_cast<int>(directions.shape[1])};
    for (std::size_t i=0;i<directions.shape[1];++i) map.x.push_back(static_cast<double>(i));
    for (std::size_t i=0;i<directions.shape[0];++i) map.y.push_back(static_cast<double>(i));
    Figure directionFigure{"observability_directions","Observable lifted-state directions",
                           "singular direction","lifted component",true};
    directionFigure.traces.push_back(std::move(map));
    Payload payload{source.u8string(),view,archiveLabel(metadata,source),2};
    payload.figures.push_back(std::move(singularFigure));
    payload.figures.push_back(std::move(directionFigure));
    return payload;
}

}  // namespace dsl
