#include <dynsyslab/core/GridFieldArchive.hpp>

#include <stdexcept>

namespace dynsyslab::core {
namespace {
constexpr const char* kFormat = "dynsyslab.grid-field/v1";
}

bool isGridFieldArchive(const json::Value& metadata) { return metadata.at("format").string() == kFormat; }

NpzRecord gridFieldArchive(const GridFieldSeries& series) {
    if (series.frames.empty() || series.frameTime.size() != series.frames.size())
        throw std::invalid_argument("A grid-field archive needs frames, each with its time");
    const auto& first = series.frames.front();
    const auto& grid = first.grid;
    grid.validate();
    const std::size_t points = grid.size(), outputs = first.outputs, frames = series.frames.size();
    std::vector<std::complex<double>> values;
    std::vector<double> status, spread, counts;
    const bool hasSpread = !first.spread.empty(), hasCounts = !first.counts.empty();
    for (const auto& f : series.frames) {
        if (f.outputs != outputs || f.values.size() != points * outputs || f.status.size() != points)
            throw std::invalid_argument("Grid-field frames must share the grid and outputs");
        values.insert(values.end(), f.values.begin(), f.values.end());
        for (const auto s : f.status) status.push_back(double(s));
        if (hasSpread) spread.insert(spread.end(), f.spread.begin(), f.spread.end());
        if (hasCounts) for (const auto c : f.counts) counts.push_back(double(c));
    }
    NpzRecord out;
    json::Value::Array axes, names, outputNames;
    for (std::size_t j = 0; j < grid.axes.size(); ++j) {
        const auto& a = grid.axes[j];
        out.arrays["axis_" + std::to_string(j)] = NpyArray::fromDoubles({a.points}, grid.axisValues(j));
        axes.push_back(json::Value::Object{{"coordinate", double(a.coordinate)}, {"lower", a.lower}, {"upper", a.upper},
                                           {"points", double(a.points)}});
    }
    for (const auto& n : grid.domain.names) names.emplace_back(n);
    for (const auto& n : first.outputNames) outputNames.emplace_back(n);
    out.arrays["anchor"] = NpyArray::fromDoubles({grid.anchor.size()}, grid.anchor);
    out.arrays["field_values"] = NpyArray::fromComplex({frames, points, outputs}, values);
    out.arrays["status"] = NpyArray::fromDoubles({frames, points}, status);
    if (hasSpread) out.arrays["spread"] = NpyArray::fromDoubles({frames, points, outputs}, spread);
    if (hasCounts) out.arrays["counts"] = NpyArray::fromDoubles({frames, points}, counts);
    out.arrays["frame_time"] = NpyArray::fromDoubles({frames}, series.frameTime);
    auto meta = series.metadata.isNull() ? json::Value(json::Value::Object{}) : series.metadata;
    meta.object()["format"] = kFormat;
    meta.object()["grid"] = json::Value::Object{{"domain", functionDomainJson(grid.domain)}, {"axes", axes},
                                                {"coordinate_names", names}};
    meta.object()["outputs"] = double(outputs);
    meta.object()["output_names"] = outputNames;
    meta.object()["complex"] = first.complexValued;
    meta.object()["frames"] = double(frames);
    out.arrays["metadata_json"] = NpyArray::fromString(json::dump(meta));
    return out;
}

GridFieldSeries readGridFieldArchive(const NpzRecord& archive) {
    GridFieldSeries out;
    out.metadata = json::parse(archive.metadataJson());
    if (!isGridFieldArchive(out.metadata)) throw std::invalid_argument("not a grid-field archive");
    const auto& gridMeta = out.metadata.at("grid");
    DomainGrid grid;
    grid.domain = functionDomainFromJson(gridMeta.at("domain"));
    grid.anchor = archive.at("anchor").doubles();
    for (const auto& a : gridMeta.at("axes").array())
        grid.axes.push_back({std::size_t(a.at("coordinate").number()), a.at("lower").number(), a.at("upper").number(),
                             std::size_t(a.at("points").number())});
    grid.validate();
    const auto& valuesArray = archive.at("field_values");
    if (valuesArray.shape.size() != 3 || valuesArray.shape[1] != grid.size())
        throw std::invalid_argument("grid-field values do not match the grid");
    const std::size_t frames = valuesArray.shape[0], points = valuesArray.shape[1], outputs = valuesArray.shape[2];
    const auto values = valuesArray.complexDoubles();
    const auto status = archive.at("status").doubles();
    const auto spread = archive.contains("spread") ? archive.at("spread").doubles() : std::vector<double>{};
    const auto counts = archive.contains("counts") ? archive.at("counts").doubles() : std::vector<double>{};
    out.frameTime = archive.at("frame_time").doubles();
    std::vector<std::string> names;
    for (const auto& n : out.metadata.at("output_names").array()) names.push_back(n.string());
    for (std::size_t f = 0; f < frames; ++f) {
        GridField field;
        field.grid = grid;
        field.outputs = outputs;
        field.outputNames = names;
        field.complexValued = out.metadata.at("complex").boolean();
        field.values.assign(values.begin() + std::ptrdiff_t(f * points * outputs), values.begin() + std::ptrdiff_t((f + 1) * points * outputs));
        for (std::size_t p = 0; p < points; ++p) field.status.push_back(GridStatus(int(status[f * points + p])));
        if (!spread.empty())
            field.spread.assign(spread.begin() + std::ptrdiff_t(f * points * outputs), spread.begin() + std::ptrdiff_t((f + 1) * points * outputs));
        if (!counts.empty())
            for (std::size_t p = 0; p < points; ++p) field.counts.push_back(std::size_t(counts[f * points + p]));
        out.frames.push_back(std::move(field));
    }
    return out;
}

}  // namespace dynsyslab::core
