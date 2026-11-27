// Receiver arrays through jobs (receiver_signals, state_readout, ring_errors) on free RPE ringdowns:
// the far receiver holds the near one's range-corrected monopole delayed by (D_2 - D_1)/(c T); the
// neighbour readout recovers a held-out ringdown from two receivers; the linearized delay inverse
// recovers a small-amplitude one; ring errors pool to the whole-record score.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
using O = json::Value::Object;
using A = json::Value::Array;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path ringdown(const std::filesystem::path& root, const std::string& name, double r0) {
    core::SystemConfiguration spec;
    spec.identifier = "rpe";
    auto system = core::createSystem(spec);
    core::IntegrationOptions o;
    o.finalTime = 12; o.outputStep = 1e-3; o.method = core::IntegrationMethod::bdf; o.relativeTolerance = o.absoluteTolerance = 1e-10;
    const auto s = core::integrate(*system, {r0, 0}, 1, o);
    std::vector<double> values(s.values.begin(), s.values.begin() + std::ptrdiff_t(s.sampleCount * 2));
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({s.sampleCount}, s.time);
    record.arrays["values"] = core::NpyArray::fromDoubles({s.sampleCount, 2}, values);
    record.arrays["metadata_json"] = core::NpyArray::fromString(
        R"({"schema":"dynsyslab.trajectory/v1","label":")" + name + R"(","system":"rpe","archive_coordinate_system":"dimensionless","integration_form":"dimensionless","forcing_type":"free"})");
    const auto path = root / (name + ".npz");
    core::writeNpzAtomic(path, record);
    return path;
}

std::filesystem::path run(const std::filesystem::path& source, const std::filesystem::path& root, const std::string& method, O settings) {
    const json::Value job = O{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", method},
        {"sources", A{source.generic_u8string()}}, {"methods", A{O{{"name", method}, {"settings", std::move(settings)}}}}};
    return dynsyslab::workflow::executeJob(job, root).archives.front();
}

double score(const std::filesystem::path& rings, std::size_t component) {
    return json::parse(core::readNpz(rings).metadataJson()).at("validation").at("relative_rms_percent").array().at(component).number();
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-receivers-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const O geometry{{"receivers", A{A{0.005, 0.0, 0.0}, A{0.0, 0.006, 0.0}}}, {"sound_speed", 1480.0}, {"stride", 2.0},
                         {"terms", A{"monopole_pressure"}}};
        const auto training = run(ringdown(root, "training", 1.3), root / "a", "receiver_signals", geometry);
        const auto heldOut = run(ringdown(root, "held-out", 1.22), root / "b", "receiver_signals", geometry);

        // Retardation: the range-corrected monopole is the same echo at both receivers, delayed.
        const auto features = core::readNpz(heldOut);
        const auto t = features.at("time").doubles(), y = features.at("values").doubles(), lag = features.at("lag_tau").doubles();
        require(lag[0] == 0 && lag[1] > 0, "the far receiver lags the near one");
        double worst = 0, size = 0;
        for (std::size_t k = 0; k < t.size(); ++k) {
            const double back = t[k] - lag[1];
            if (back < t.front()) continue;
            const auto j = std::size_t(std::upper_bound(t.begin(), t.end(), back) - t.begin()) - 1;
            const double f = (back - t[j]) / (t[j + 1] - t[j]);
            worst = std::max(worst, std::abs(y[2 * k + 1] - (y[2 * j] + f * (y[2 * j + 2] - y[2 * j]))));
            size = std::max(size, std::abs(y[2 * k]));
        }
        require(worst < 2e-3 * size, "the far receiver holds the retarded near-receiver echo");
        require(features.at("ring_bounds").doubles().size() > 8, "radius maxima of the ringdown");

        // Neighbour readout of the held-out ringdown from the training library.
        const auto knn = run(heldOut, root / "c", "state_readout", {{"training", A{training.generic_u8string()}}, {"neighbors", 4.0}});
        const auto readout = core::readNpz(knn);
        require(readout.contains("neighbor_rows") && readout.at("neighbor_weights").shape[1] == 4, "neighbour provenance");
        const auto rings = run(knn, root / "d", "ring_errors", {{"center", A{1.0, 0.0}}, {"windows", O{{"early", A{0.0, 3.0}}}}});
        require(score(rings, 0) < 15 && score(rings, 1) < 15, "two receivers recover the held-out ringdown");
        const auto ringArchive = core::readNpz(rings);
        require(ringArchive.at("relative_rms_percent").shape[1] == 2 && ringArchive.at("count").shape[0] > 6, "per-ring errors");

        // Affine least-squares readout from the same features.
        const auto affine = run(heldOut, root / "h", "state_readout", {{"training", A{training.generic_u8string()}}, {"mode", "affine"}});
        require(core::readNpz(affine).at("readout_weights").shape[0] == 2, "affine readout weights");
        // Linear delay inverse about rest for a small ringdown.
        const auto small = run(run(ringdown(root, "small", 1.01), root / "e", "receiver_signals", geometry), root / "f", "state_readout",
                               {{"training", A{training.generic_u8string()}}, {"mode", "linear"}, {"equilibrium", A{1.0, 0.0}}});
        const auto linearRings = run(small, root / "g", "ring_errors", {{"center", A{1.0, 0.0}}});
        require(score(linearRings, 0) < 5 && score(linearRings, 1) < 5, "the linear delay inverse recovers a small ringdown");
        std::filesystem::remove_all(root);
        std::cout << "receiver pipeline tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
