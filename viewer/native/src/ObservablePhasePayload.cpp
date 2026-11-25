#include "ObservablePhasePayload.hpp"
#include "DisplayCoordinates.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace dsl {
namespace c = dynsyslab::core;

Payload buildObservablePhasePayload(const std::filesystem::path& source, const c::NpzRecord& record,
                                    const json::Value& metadata) {
    const auto sourcePath = std::filesystem::u8path(metadata.at("source_archive").string());
    const auto trajectory = c::readNpz(sourcePath);
    const auto trajectoryMeta = trajectory.metadataJson().empty() ? json::Value{} : json::parse(trajectory.metadataJson());
    const auto& rawStates = trajectory.contains("values") ? trajectory.at("values") : trajectory.at("states");
    if (rawStates.shape.size() != 2 || rawStates.shape[1] < 2)
        throw std::runtime_error("the source trajectory has no two-dimensional phase plane");
    const auto transform = displayCoordinateTransform(trajectoryMeta, rawStates.shape[1]);
    const auto states = transformDisplayStates(rawStates, transform).doubles();
    const auto trajectoryTime = trajectory.at("time").doubles();
    const std::size_t n = rawStates.shape[0], width = rawStates.shape[1];

    const auto& valueArray = record.at("values");
    const auto time = record.at("time").doubles();
    const auto values = valueArray.doubles();
    const std::size_t columns = valueArray.shape.size() == 2 ? valueArray.shape[1] : 1;
    const bool complex = columns == 2 && !metadata.at("eigenfunction").isNull();
    const auto& names = metadata.at("coordinate_names");
    std::string name = !metadata.at("eigenfunction").isNull()
        ? "\\psi_{" + std::to_string(std::size_t(metadata.at("eigenfunction").at("mode").number())) + "}"
        : (!names.isNull() && !names.array().empty() ? names.array()[0].string() : std::string("g"));

    Figure figure;
    figure.kind = "phase";
    figure.square = false;
    const auto label = [&](std::size_t axis) {
        return (axis < transform.stateLabels.size() ? transform.stateLabels[axis] : "x_" + std::to_string(axis)) + transform.suffix;
    };
    figure.xLabel = label(0);
    figure.yLabel = label(1);
    const auto representation = metadata.at("eigenfunction").at("representation").string();
    figure.title = name + (representation.empty() ? "" : " (" + representation + ")") + " on the phase plane";
    figure.colorLabel = complex ? "arg " + name : name;

    Trace path;
    path.label = "trajectory";
    path.role = "reference";
    path.color = "#94a3b8";
    path.opacity = 0.45;
    for (std::size_t k = 0; k < n; ++k) { path.x.push_back(states[k * width]); path.y.push_back(states[k * width + 1]); }
    figure.traces.push_back(std::move(path));

    Trace coloured;
    coloured.label = name;
    coloured.role = "state";
    coloured.encoding = complex ? "domain-color" : "time";
    coloured.lineWidthPoints = 2.2;
    double largest = 0;
    std::size_t row = 0, unmatched = 0;
    for (std::size_t k = 0; k < time.size(); ++k) {  // the source row at each evaluated time
        const double tolerance = 1e-9 * std::max(1.0, std::abs(time[k]));
        while (row < n && trajectoryTime[row] < time[k] - tolerance) ++row;
        if (row >= n || std::abs(trajectoryTime[row] - time[k]) > tolerance) { ++unmatched; continue; }
        coloured.x.push_back(states[row * width]);
        coloured.y.push_back(states[row * width + 1]);
        if (complex) {
            const std::complex<double> g(values[k * columns], values[k * columns + 1]);
            coloured.colorValues.push_back(std::arg(g));
            coloured.magnitudes.push_back(std::abs(g));
            largest = std::max(largest, std::abs(g));
        } else {
            coloured.colorValues.push_back(values[k * columns]);
        }
    }
    if (coloured.x.empty()) throw std::runtime_error("the observable's times do not match its source trajectory");
    if (complex && largest > 0) for (auto& m : coloured.magnitudes) m /= largest;
    figure.references.push_back({coloured.x.front(), coloured.y.front(), "start", false});
    figure.traces.push_back(std::move(coloured));

    Payload payload{source.u8string(), "observable_phase", metadata.at("label").string(source.stem().u8string()), 1};
    if (unmatched) payload.warnings.push_back(std::to_string(unmatched) + " evaluated times are not samples of the source trajectory.");
    if (complex) payload.warnings.push_back("Hue is the argument and brightness the modulus (relative to its largest value here).");
    payload.figures.push_back(std::move(figure));
    return payload;
}

}  // namespace dsl
