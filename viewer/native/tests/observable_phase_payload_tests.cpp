// View "observable_phase": an eigenfunction's values along a trajectory drawn on the trajectory's own
// phase plane, coloured by argument (hue) and modulus (brightness), over the faint whole trajectory.
#include "dynsyslab_viewer/ComparisonPayload.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "../src/ObservablePhasePayload.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-observable-phase-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        constexpr std::size_t n = 40;
        std::vector<double> time, states, psiTime, psi;
        for (std::size_t k = 0; k < n; ++k) {
            const double t = 0.1 * double(k);
            time.push_back(t);
            states.push_back(std::cos(t));
            states.push_back(-std::sin(t));
            if (k % 2 == 0 && k < 30) {  // evaluated on part of the trajectory only (e.g. complete windows)
                const auto value = std::exp(std::complex<double>(-0.1, 1.0) * t);
                psiTime.push_back(t);
                psi.push_back(value.real());
                psi.push_back(value.imag());
            }
        }
        c::NpzRecord trajectory;
        trajectory.arrays["time"] = c::NpyArray::fromDoubles({n}, time);
        trajectory.arrays["values"] = c::NpyArray::fromDoubles({n, 2}, states);
        trajectory.arrays["metadata_json"] = c::NpyArray::fromString(R"({"system":"duffing","integration_form":"dimensionless","data_role":"free_training"})");
        const auto trajectoryPath = root / "trajectory.npz";
        c::writeNpzAtomic(trajectoryPath, trajectory);
        c::NpzRecord series;
        series.arrays["time"] = c::NpyArray::fromDoubles({psiTime.size()}, psiTime);
        series.arrays["values"] = c::NpyArray::fromDoubles({psiTime.size(), 2}, psi);
        dynsyslab::core::json::Value::Object meta{{"value_space", "observable-values"},
            {"source_archive", trajectoryPath.generic_u8string()}, {"label", "psi along a trajectory"},
            {"eigenfunction", dynsyslab::core::json::Value::Object{{"mode", 2.0}, {"representation", "expansion"}}}};
        series.arrays["metadata_json"] = c::NpyArray::fromString(dynsyslab::core::json::dump(dynsyslab::core::json::Value(meta)));
        const auto seriesPath = root / "series.npz";
        c::writeNpzAtomic(seriesPath, series);

        require(dsl::isObservableSeries(dynsyslab::core::json::Value(meta)), "an observable series is recognized");
        const auto payload = dsl::buildNativePayload(seriesPath, "observable_phase");
        require(payload.figures.size() == 1 && payload.figures.front().kind == "phase", "one phase-plane figure");
        const auto& figure = payload.figures.front();
        require(figure.traces.size() == 2 && figure.traces[0].x.size() == n, "the whole trajectory underneath");
        require(figure.traces[0].label == "training trajectory" && figure.traces[0].lineWidthPoints > 3
                && figure.traces[1].lineWidthPoints > 3, "a training run is drawn bolder");
        require(figure.traces[1].label.empty(), "the colour bar, not the legend, names the observable");
        require(!figure.references.empty() && figure.references.front().label == "IC"
                && !figure.references.front().equilibrium, "the initial condition is an open marker");
        const auto& coloured = figure.traces[1];
        require(coloured.encoding == "domain-color" && coloured.x.size() == psiTime.size(), "one coloured state per value");
        for (std::size_t k = 0; k < psiTime.size(); ++k) {
            require(std::abs(coloured.x[k] - std::cos(psiTime[k])) < 1e-12, "states matched by time");
            require(std::abs(coloured.colorValues[k] - std::arg(std::complex<double>(psi[2 * k], psi[2 * k + 1]))) < 1e-12,
                    "hue is the argument");
        }
        require(std::abs(coloured.magnitudes.front() - 1.0) < 1e-12, "brightness relative to the largest modulus");
        require(figure.title.find("(expansion)") != std::string::npos, "a dictionary expansion keeps its name");

        // A window map is stored as an "expansion" but shown as "time-delay".
        c::NpzRecord object;
        object.arrays["metadata_json"] = c::NpyArray::fromString(
            R"({"krom":{"representations":[{"kind":"window_map","causal":true}]}})");
        const auto objectPath = root / "eigenfunctions.npz";
        c::writeNpzAtomic(objectPath, object);
        meta["eigenfunction"] = dynsyslab::core::json::Value::Object{{"mode", 2.0}, {"representation", "expansion"},
            {"object", objectPath.generic_u8string()}};
        series.arrays["metadata_json"] = c::NpyArray::fromString(dynsyslab::core::json::dump(dynsyslab::core::json::Value(meta)));
        const auto windowSeriesPath = root / "window-series.npz";
        c::writeNpzAtomic(windowSeriesPath, series);
        const auto windowPayload = dsl::buildNativePayload(windowSeriesPath, "observable_phase");
        require(windowPayload.figures.front().title.find("(time-delay)") != std::string::npos,
                "a window map is labelled time-delay");

        // Overlaid runs (training and test trajectories) share one phase plane: each keeps its domain
        // colouring and grey path, and brightness uses the largest modulus among all of them.
        auto doubled = payload;
        for (auto& m : doubled.figures.front().traces[1].magnitudes) m *= 2.0;
        const auto overlay = dsl::comparePayloads({{"train", payload, "#1f77b4"}, {"test", doubled, "#d62728"}});
        require(overlay.figures.size() == 1 && overlay.figures.front().traces.size() == 4, "one overlaid figure");
        double largest = 0.0;
        for (const double m : doubled.figures.front().traces[1].magnitudes) largest = std::max(largest, m);
        for (const auto& trace : overlay.figures.front().traces) {
            if (trace.role == "reference") require(trace.color == "#94a3b8" && trace.lineStyle.empty(), "grey path kept");
            else require(trace.encoding == "domain-color" && std::abs(trace.magnitudeScale - largest) < 1e-12,
                         "domain colouring kept with one shared brightness scale");
        }
        // One legend entry per kind of line; each run is named at its start marker.
        std::size_t legend = 0;
        for (const auto& trace : overlay.figures.front().traces) legend += !trace.label.empty();
        const auto& markers = overlay.figures.front().references;
        require(legend == 1, "one legend entry per kind of line, no run names");
        require(overlay.figures.front().title.find("train") == std::string::npos, "no run names in the title");
        std::size_t initial = 0, equilibria = 0;
        for (const auto& marker : markers) {
            if (marker.equilibrium) { ++equilibria; continue; }
            ++initial;
            require(marker.label == "IC" && std::abs(marker.x - 1.0) < 1e-12, "the initial condition is an IC marker");
        }
        require(initial == 1, "a shared initial condition is drawn once");
        require(equilibria == payload.figures.front().references.size() - 1, "each equilibrium is drawn once");
        std::filesystem::remove_all(root);
        std::cout << "observable phase payload tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
