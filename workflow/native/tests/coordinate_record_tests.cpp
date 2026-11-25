// Coordinate changes are recorded, not copied (core/Archive.hpp recorded archives): the stage writes
// only its metadata and the record; reading it replays the change from its source; a chain replays
// through every step; a changed source is refused; "store_values" keeps the arrays.
#include <dynsyslab/core/Archive.hpp>
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
std::filesystem::path transform(const std::filesystem::path& source, const std::filesystem::path& root, O settings) {
    const json::Value job = O{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", "coords"},
        {"sources", A{source.generic_u8string()}},
        {"methods", A{O{{"name", "coordinate_transform"}, {"enabled", true}, {"settings", std::move(settings)}}}},
        {"archive", O{{"enabled", true}}}};
    return dynsyslab::workflow::executeJob(job, root).archives.front();
}
}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-coordinate-record-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto sourcePath = root / "trajectory.npz";
        constexpr std::size_t samples = 50;
        std::vector<double> time(samples), states;
        for (std::size_t k = 0; k < samples; ++k) {
            time[k] = 0.05 * double(k);
            states.push_back(std::cos(time[k]));
            states.push_back(-std::sin(time[k]));
        }
        core::NpzRecord source;
        source.arrays["time"] = core::NpyArray::fromDoubles({samples}, time);
        source.arrays["values"] = core::NpyArray::fromDoubles({samples, 2}, states);
        source.arrays["metadata_json"] = core::NpyArray::fromString(
            R"({"format":"dynsyslab.trajectory/v1","label":"source","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{"delta":0},"integration_form":"dimensionless"})");
        core::writeNpzAtomic(sourcePath, source);

        // An observable map: stored as metadata and a record only; reading replays it.
        const auto forcePath = transform(sourcePath, root / "force", {{"kind", "observable"}, {"observable", "force"}});
        const auto stored = core::readNpzStored(forcePath);
        require(stored.arrays.size() == 1 && stored.contains("metadata_json"), "the coordinate change stores no arrays");
        require(core::isRecordedArchive(stored.metadataJson()), "the stored metadata holds the record");
        const auto replayed = core::readNpz(forcePath);
        const auto values = replayed.at("values").doubles();
        const auto force = core::createSystem([] {
            core::SystemConfiguration c; c.identifier = "duffing"; c.preset = "double_well_reference"; c.parameters = {{"delta", 0.0}};
            return c;
        }())->observableCatalog().at("force");
        for (std::size_t k = 0; k < samples; ++k)
            require(std::abs(values[k] - force.evaluateReal(states.data() + 2 * k)[0]) < 1e-12, "replayed values are the map of the source");
        require(core::readNpzIndex(forcePath).contains("values"), "the index lists the replayed arrays");

        // A chain: identity, then a column of the identity; every step is replayed from the root.
        const auto chained = transform(sourcePath, root / "chain", {{"kind", "chain"}, {"steps", A{
            O{{"kind", "identity"}}, O{{"kind", "columns"}, {"columns", A{1.0}}}}}});
        require(core::readNpzStored(chained).arrays.size() == 1, "chained steps store no arrays");
        const auto column = core::readNpz(chained).at("values").doubles();
        for (std::size_t k = 0; k < samples; ++k) require(column[k] == states[2 * k + 1], "the chain replays through both steps");

        // An explicit copy keeps the arrays.
        const auto copied = transform(sourcePath, root / "copy", {{"kind", "identity"}, {"store_values", true}});
        require(core::readNpzStored(copied).contains("values"), "store_values keeps the arrays");

        // A changed source is refused rather than silently replayed.
        source.arrays["values"] = core::NpyArray::fromDoubles({samples, 2}, std::vector<double>(2 * samples, 0.5));
        core::writeNpzAtomic(sourcePath, source);
        bool refused = false;
        try { (void)core::readNpz(forcePath); } catch (const std::runtime_error&) { refused = true; }
        require(refused, "a changed source is refused");
        std::filesystem::remove_all(root);
        std::cout << "coordinate record tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
