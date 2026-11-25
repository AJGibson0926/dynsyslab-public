// The domain_field job: a function on a grid of its domain, the Koopman action of the flow on it,
// and a push-forward through a coordinate map, archived as dynsyslab.grid-field/v1.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/GridFieldArchive.hpp>
#include <dynsyslab/core/Json.hpp>
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

core::GridFieldSeries run(const std::filesystem::path& source, const std::filesystem::path& root, O settings) {
    const json::Value job = O{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", "domain"},
        {"sources", A{source.generic_u8string()}},
        {"methods", A{O{{"name", "domain_field"}, {"enabled", true}, {"settings", std::move(settings)}}}},
        {"archive", O{{"enabled", true}}}};
    const auto result = dynsyslab::workflow::executeJob(job, root);
    return core::readGridFieldArchive(core::readNpz(result.archives.front()));
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-domain-field-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto sourcePath = root / "trajectory.npz";
        core::NpzRecord source;
        source.arrays["time"] = core::NpyArray::fromDoubles({3}, std::vector<double>{0, 0.1, 0.2});
        source.arrays["values"] = core::NpyArray::fromDoubles({3, 2}, std::vector<double>{1, 0, 0.99, -0.1, 0.98, -0.2});
        source.arrays["metadata_json"] = core::NpyArray::fromString(
            R"({"format":"dynsyslab.trajectory/v1","label":"source","system":"harmonic_oscillator","system_preset":"undamped","resolved_system_parameters":{"mass_kg":1,"stiffness_n_m":4,"damping_kg_s":0},"integration_form":"dimensional"})");
        core::writeNpzAtomic(sourcePath, source);
        const json::Value box = O{{"lower", A{-1.0, -1.0}}, {"upper", A{1.0, 1.0}}, {"points", 11.0}};

        // Evaluate: g = x^2 + v on the grid.
        const auto plain = run(sourcePath, root / "evaluate", {{"function", O{{"expression", "position^2 + velocity"}}}, {"domain", box}});
        const auto& field = plain.frames.front();
        require(field.grid.size() == 121 && std::abs(field.value(0).real() - 0.0) < 1e-12, "evaluate on the grid");
        require(plain.metadata.at("format").string() == "dynsyslab.grid-field/v1", "grid-field archive");

        // Koopman action: the energy k x^2 + m v^2 is invariant, and x o Phi^t is the exact rotation.
        const auto energy = run(sourcePath, root / "action", {{"function", O{{"expression", "stiffness_n_m*position^2 + mass_kg*velocity^2"}}},
            {"domain", box}, {"operation", "koopman_action"}, {"time_step", 0.25}, {"steps", 3.0}, {"maximum_step", 1e-3}});
        require(energy.frames.size() == 4 && std::abs(energy.frameTime[3] - 0.75) < 1e-12, "one frame per step");
        for (std::size_t k = 0; k < 121; ++k)
            require(std::abs(energy.frames[3].value(k) - energy.frames[0].value(k)) < 1e-8, "invariant under the Koopman action");

        // Push-forward of x through h = (x^2, v): two-valued fibres are conflicts; u < 0 has no preimage.
        const auto pushed = run(sourcePath, root / "push", {{"function", O{{"expression", "position"}}},
            {"map", O{{"expression", A{"position^2", "velocity"}}}}, {"operation", "push_forward"}, {"domain", box},
            {"target", O{{"lower", A{-1.0, -1.0}}, {"upper", A{1.0, 1.0}}, {"points", 11.0}}}, {"tolerance", 0.05}});
        const auto& p = pushed.frames.front();
        require(p.count(core::GridStatus::conflict) > 0 && p.count(core::GridStatus::noPreimage) > 0, "push-forward status");
        require(!p.spread.empty() && !p.counts.empty(), "push-forward keeps the fibre spread and counts");
        std::filesystem::remove_all(root);
        std::cout << "domain field tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
