#include "ObservablePhasePayload.hpp"
#include "ArchiveSystem.hpp"
#include "DisplayCoordinates.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace dsl {
namespace c = dynsyslab::core;

namespace {
// Display name of an eigenfunction's representation. "expansion" names both a dictionary
// expansion and a window map; the window map is shown as "time-delay". The choice follows
// KoopmanModel::representation: the first representation of the requested kind in the krom.
std::string representationLabel(const json::Value& eigenfunction) {
    const auto requested = eigenfunction.at("representation").string();
    if (requested != "expansion" || eigenfunction.at("object").isNull()) return requested;
    try {
        const auto index = c::readNpzIndex(std::filesystem::u8path(eigenfunction.at("object").string()));
        const auto meta = json::parse(index.metadataJson);
        for (const auto& rep : meta.at("krom").at("representations").array()) {
            const auto kind = rep.at("kind").string();
            if (kind == "window_map") return "time-delay";
            if (kind != "expansion") continue;
            return rep.at("dictionary").isNull() && rep.at("family").string() == "linear-window" ? "time-delay" : requested;
        }
    } catch (const std::exception&) {}  // object moved or unreadable: keep the stored name
    return requested;
}
}  // namespace

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

    // Training runs are drawn bolder than test runs (the source's recorded data role).
    const auto role = trajectoryMeta.at("data_role").string();
    const bool training = role.find("training") != std::string::npos;
    const bool test = role.find("test") != std::string::npos;

    Figure figure;
    figure.kind = "phase";
    figure.square = false;
    figure.coordinates = transform.coordinates;
    const auto label = [&](std::size_t axis) {
        return (axis < transform.stateLabels.size() ? transform.stateLabels[axis] : "x_" + std::to_string(axis)) + transform.suffix;
    };
    figure.xLabel = label(0);
    figure.yLabel = label(1);
    const auto representation = representationLabel(metadata.at("eigenfunction"));
    figure.title = name + (representation.empty() ? "" : " (" + representation + ")") + " on the phase plane";
    figure.colorLabel = complex ? "arg " + name : name;

    Trace path;
    path.label = training ? "training trajectory" : test ? "test trajectory" : "trajectory";
    path.role = "reference";
    path.color = "#94a3b8";
    path.opacity = training ? 0.7 : 0.45;
    path.lineWidthPoints = training ? 3.4 : 1.6;
    for (std::size_t k = 0; k < n; ++k) { path.x.push_back(states[k * width]); path.y.push_back(states[k * width + 1]); }
    figure.traces.push_back(std::move(path));

    Trace coloured;
    coloured.label.clear();  // the colour bar names the observable
    coloured.role = "state";
    coloured.encoding = complex ? "domain-color" : "time";
    coloured.lineWidthPoints = training ? 3.4 : 1.8;
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
        } else {
            coloured.colorValues.push_back(values[k * columns]);
        }
    }
    if (coloured.x.empty()) throw std::runtime_error("the observable's times do not match its source trajectory");
    // Magnitudes stay raw |g|: the plot scales brightness by the largest modulus drawn
    // (shared across runs when several are overlaid).
    figure.references.push_back({states[0], states[1], "IC", false});  // drawn with its coordinates
    figure.traces.push_back(std::move(coloured));
    appendVerifiedEquilibria(figure, trajectoryMeta);

    Payload payload{source.u8string(), "observable_phase", metadata.at("label").string(source.stem().u8string()), 1};
    if (unmatched) payload.warnings.push_back(std::to_string(unmatched) + " evaluated times are not samples of the source trajectory.");
    if (complex) payload.warnings.push_back("Hue is the argument and brightness the modulus (relative to the largest modulus drawn).");
    payload.figures.push_back(std::move(figure));
    return payload;
}

}  // namespace dsl
