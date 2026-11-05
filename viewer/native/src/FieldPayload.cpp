#include "FieldPayload.hpp"
#include <algorithm>
#include <stdexcept>
namespace dsl {
namespace {
using dynsyslab::core::NpyArray;
const NpyArray *select(const dynsyslab::core::NpzRecord &record,
                       std::initializer_list<const char *> names, std::size_t dimensions) {
    for (auto name : names)
        if (record.contains(name) && record.at(name).shape.size() == dimensions)
            return &record.at(name);
    return nullptr;
}
Trace line(std::string label, const std::vector<double> &x, std::vector<double> y) {
    Trace trace;
    trace.label = std::move(label);
    trace.x = x;
    trace.y = std::move(y);
    trace.role = "data";
    return trace;
}
} // namespace
using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;
namespace json = dynsyslab::core::json;
Payload buildFieldPayload(const std::filesystem::path &source, const std::string &view,
                          const NpzRecord &record, const json::Value &metadata) {
    const NpyArray *fields = select(record, {"values", "fields"}, 2);
    const NpyArray *space = select(record, {"space", "x"}, 1);
    const NpyArray *times = select(record, {"time"}, 1);
    if (!fields || !space || !times || fields->shape[0] == 0 || fields->shape[1] == 0 ||
        fields->shape[0] != times->elementCount() || fields->shape[1] != space->elementCount())
        throw std::runtime_error("PDE archive has incompatible field coordinates");
    const auto rowIndices = [&] {
        const std::size_t count = std::min<std::size_t>(fields->shape[0], 480);
        std::vector<std::size_t> result(count);
        for (std::size_t i = 0; i < count; ++i)
            result[i] = count == 1 ? 0 : i * (fields->shape[0] - 1) / (count - 1);
        return result;
    }();
    const auto columnIndices = [&] {
        const std::size_t count = std::min<std::size_t>(fields->shape[1], 320);
        std::vector<std::size_t> result(count);
        for (std::size_t i = 0; i < count; ++i)
            result[i] = count == 1 ? 0 : i * (fields->shape[1] - 1) / (count - 1);
        return result;
    }();
    const auto fieldValues = fields->doubles();
    const auto xValues = space->doubles();
    const auto timeValues = times->doubles();
    Trace map;
    map.label = "field";
    map.heatmap = true;
    map.zShape = {static_cast<int>(rowIndices.size()), static_cast<int>(columnIndices.size())};
    for (const auto column : columnIndices)
        map.x.push_back(xValues[column]);
    for (const auto row : rowIndices)
        map.y.push_back(timeValues[row]);
    map.z.reserve(rowIndices.size() * columnIndices.size());
    for (const auto row : rowIndices)
        for (const auto column : columnIndices)
            map.z.push_back(fieldValues[row * fields->shape[1] + column]);
    Figure field{"field", "Space-time field", "space", "time", false};
    field.traces.push_back(std::move(map));
    Figure profile{"field_profile", "Final spatial profile", "space", "field", false};
    std::vector<double> finalValues;
    finalValues.reserve(columnIndices.size());
    const std::size_t finalRow = fields->shape[0] - 1;
    for (const auto column : columnIndices)
        finalValues.push_back(fieldValues[finalRow * fields->shape[1] + column]);
    profile.traces.push_back(line("final field", field.traces[0].x, std::move(finalValues)));
    Payload payload{source.u8string(), view, metadata.at("label").string(source.stem().u8string()),
                    2};
    if (record.contains("analytic_values")) {
        const auto &analytic = record.at("analytic_values");
        bool valid = analytic.shape == fields->shape;
        if (record.contains("analytic_valid")) {
            const auto mask = record.at("analytic_valid").doubles();
            valid = valid && mask.size() == fields->shape[0] && mask.back() == 1;
        }
        if (valid) {
            const auto values = analytic.doubles();
            std::vector<double> y;
            for (const auto column : columnIndices)
                y.push_back(values[finalRow * fields->shape[1] + column]);
            auto overlay = line("Analytic", profile.traces.front().x, std::move(y));
            overlay.role = "analytic_solution";
            overlay.color = "#7b6d9c";
            overlay.lineStyle = "dash";
            profile.traces.push_back(std::move(overlay));
        } else
            payload.warnings.push_back("Analytic reference is unavailable at the final field time "
                                       "or has an incompatible shape.");
    }
    payload.figures.push_back(std::move(field));
    payload.figures.push_back(std::move(profile));
    return payload;
}

} // namespace dsl
