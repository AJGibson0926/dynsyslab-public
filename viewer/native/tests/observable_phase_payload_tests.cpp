// View "observable_phase": an eigenfunction's values along a trajectory drawn on the trajectory's own
// phase plane, coloured by argument (hue) and modulus (brightness), over the faint whole trajectory.
#include "dynsyslab_viewer/NativePayload.hpp"
#include "../src/ObservablePhasePayload.hpp"

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
        trajectory.arrays["metadata_json"] = c::NpyArray::fromString(R"({"system":"duffing","integration_form":"dimensionless"})");
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
        const auto& coloured = figure.traces[1];
        require(coloured.encoding == "domain-color" && coloured.x.size() == psiTime.size(), "one coloured state per value");
        for (std::size_t k = 0; k < psiTime.size(); ++k) {
            require(std::abs(coloured.x[k] - std::cos(psiTime[k])) < 1e-12, "states matched by time");
            require(std::abs(coloured.colorValues[k] - std::arg(std::complex<double>(psi[2 * k], psi[2 * k + 1]))) < 1e-12,
                    "hue is the argument");
        }
        require(std::abs(coloured.magnitudes.front() - 1.0) < 1e-12, "brightness relative to the largest modulus");
        std::filesystem::remove_all(root);
        std::cout << "observable phase payload tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
