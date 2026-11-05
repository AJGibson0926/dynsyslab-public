#include "LatticePayload.hpp"

#include <stdexcept>

namespace dsl {
namespace {

Figure latticeFigure(const dynsyslab::core::NpyArray& values,
                      const std::string& title) {
    if (values.shape.size() != 2 && values.shape.size() != 3)
        throw std::invalid_argument("lattice view requires a 2-D or 3-D site array");
    const std::size_t rows = values.shape[values.shape.size() - 2];
    const std::size_t columns = values.shape.back();
    Trace trace; trace.label = title; trace.heatmap = true;
    for (std::size_t column = 0; column < columns; ++column)
        trace.x.push_back(static_cast<double>(column));
    for (std::size_t row = 0; row < rows; ++row)
        trace.y.push_back(static_cast<double>(row));
    const auto all = values.doubles();
    const std::size_t offset = values.shape.size() == 3
        ? (values.shape.front() / 2) * rows * columns : 0;
    trace.z.assign(all.begin() + static_cast<std::ptrdiff_t>(offset),
                   all.begin() + static_cast<std::ptrdiff_t>(offset + rows * columns));
    trace.zShape = {static_cast<int>(rows), static_cast<int>(columns)};
    const std::string shown = values.shape.size() == 3 ? title + " (central slice)" : title;
    Figure figure{"lattice", shown, "site column", "site row", true};
    figure.traces.push_back(std::move(trace)); return figure;
}

}  // namespace

Payload buildLatticePayload(const std::filesystem::path& source, const std::string& view,
    const dynsyslab::core::NpzRecord& record, const json::Value& metadata) {
    if (!record.contains("lattice_initial") || !record.contains("lattice_final"))
        throw std::runtime_error("archive has no retained lattice configurations");
    Payload payload{source.u8string(), view,
        metadata.at("label").string(source.stem().u8string()), 2};
    payload.figures.push_back(latticeFigure(record.at("lattice_initial"), "Initial lattice"));
    payload.figures.push_back(latticeFigure(record.at("lattice_final"), "Final lattice"));
    return payload;
}

}  // namespace dsl
