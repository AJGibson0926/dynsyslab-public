#include "AnalyticOverlay.hpp"
#include <cmath>
#include <dynsyslab/core/Json.hpp>
#include <limits>
namespace dsl {
void appendAnalyticOverlay(Payload &payload, const dynsyslab::core::NpzRecord &record,
                           const DisplayCoordinateTransform &coordinates,
                           const std::vector<double> &time, const std::vector<std::size_t> &indices,
                           std::size_t samples, std::size_t dimension) {
    const auto metadataText = record.metadataJson();
    const auto metadata = metadataText.empty() ? dynsyslab::core::json::Value{}
                                               : dynsyslab::core::json::parse(metadataText);
    if (metadata.at("analytic_solution").at("status").string() == "partial")
        payload.warnings.push_back("Analytic reference: " +
                                   metadata.at("analytic_solution").at("reason").string());
    if (!record.contains("analytic_values"))
        return;
    const auto &raw = record.at("analytic_values");
    if (raw.shape != std::vector<std::size_t>{samples, dimension}) {
        payload.warnings.push_back("Analytic comparison shape does not match this trajectory.");
        return;
    }
    auto values = transformDisplayStates(raw, coordinates).doubles();
    if (record.contains("analytic_valid")) {
        const auto valid = record.at("analytic_valid").doubles();
        if (valid.size() != samples) {
            payload.warnings.push_back("Analytic validity mask has an incompatible shape.");
            return;
        }
        for (std::size_t row = 0; row < samples; ++row)
            if (valid[row] != 1)
                for (std::size_t col = 0; col < dimension; ++col)
                    values[row * dimension + col] = std::numeric_limits<double>::quiet_NaN();
    }
    auto column = [&](std::size_t c) {
        std::vector<double> out;
        out.reserve(indices.size());
        for (auto row : indices)
            out.push_back(values[row * dimension + c]);
        return out;
    };
    std::size_t component = 0;
    for (auto &figure : payload.figures) {
        if (figure.traces.empty() || figure.traces.front().role != "data")
            continue;
        Trace trace;
        trace.label = "Analytic";
        trace.role = "analytic_solution";
        trace.color = "#7b6d9c";
        trace.lineStyle = "dash";
        if (figure.kind == "phase" && dimension >= 2) {
            trace.x = column(0);
            trace.y = column(1);
        } else if (figure.kind == "three_d" && dimension >= 3) {
            trace.x = column(0);
            trace.y = column(1);
            trace.z = column(2);
        } else if (figure.kind == "time" && component < dimension) {
            trace.x = time;
            trace.y = column(component++);
        } else
            continue;
        figure.traces.push_back(std::move(trace));
    }
}
} // namespace dsl
