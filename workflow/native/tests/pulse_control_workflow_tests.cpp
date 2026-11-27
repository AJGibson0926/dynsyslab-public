// Pulse-branch acquisition and input-map audits through executeJob: pulse_branches solves restarted
// from archived states (archive_states, repeat), their paired causal sensitivities (pulse_response),
// the symmetric_pulse input-map fit on the paired first responses, and control_map_audit (path
// representation error against the exact map, design rank and condition).
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
using J = json::Value;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

J solve(double stop, J controller, J initial) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"pulses",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":1,"output_step":0.01,"coordinate_system":"automatic"},
      "solver":{"method":"RK4"},"archive":{"enabled":true,"data_role":"forced_training"}})");
    job.at("time").at("stop") = stop;
    job.at("controller") = controller;
    job.at("initial_conditions") = initial;
    return job;
}

std::vector<std::filesystem::path> run(const J& job, const std::filesystem::path& root) {
    const auto result = dynsyslab::workflow::executeJob(job, root);
    require(!result.archives.empty(), "job published no archive");
    return result.archives;
}

std::filesystem::path postprocess(const std::string& method, const J& settings, const std::vector<std::filesystem::path>& sources,
                                  const std::filesystem::path& root) {
    J::Array list;
    for (const auto& path : sources) list.emplace_back(path.generic_u8string());
    const J job = J::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", method},
        {"sources", J(std::move(list))}, {"methods", J::Array{J::Object{{"name", method}, {"enabled", true}, {"settings", settings}}}},
        {"archive", J::Object{{"enabled", true}}}};
    return run(job, root).front();
}

// psi = x on the Duffing state: the exact input map grad psi B is [0, 1].
std::filesystem::path kromFixture(const std::filesystem::path& root) {
    auto fit = core::polynomialEigenfunctions(2, {1, 0, 0, 1}, {{1, 0}, {0, 0}, {0, 0}, {1, 0}},
        core::Spectrum::fromContinuousRates({{-1.0, 0.0}, {-2.0, 0.0}}), {0.0, 0.0}, {1.0, 1.0});
    core::SystemConfiguration system;
    system.identifier = "duffing"; system.preset = "double_well_reference";
    core::KoopmanModel model(std::move(fit), std::move(system));
    model.setControlMap("constant", core::ConstantControlMap{{{0.0, 0.0}, {1.0, 0.0}}, 2, 1, "fixture"});
    model.setControlMap("biased", core::ConstantControlMap{{{0.0, 0.0}, {1.1, 0.0}}, 2, 1, "fixture"});
    const auto path = root / "fixture-krom.npz";
    core::saveKoopmanModel(path, model, R"({"schema":"dynsyslab.derived/v1","label":"fixture","delta_t":0.01})");
    return path;
}

J metadataOf(const std::filesystem::path& path) { return json::parse(core::readNpz(path).metadataJson()); }

void pulseTest(const std::filesystem::path& root) {
    const auto krom = kromFixture(root);
    const auto base = run(solve(1.0, J::Object{{"kind", "none"}},
        json::parse(R"({"mode":"explicit","count":1,"values":[[-0.8,0.1]]})")), root).front();
    const auto anchors = J::Object{{"mode", "archive_states"}, {"archive", base.generic_u8string()}, {"times", J::Array{0.2, 0.5}},
                                   {"count", 10.0}, {"repeat", 5.0}};
    auto pulseJob = solve(0.3, J::Object{{"kind", "pulse_branches"}, {"settings",
        J::Object{{"amplitudes", J::Array{0.0, 0.1, -0.1, 0.2, -0.2}}, {"widths", 0.01}}}}, anchors);
    pulseJob.at("solver").object()["control_evaluation"] = "sampled";  // held pulses
    const auto pulses = run(pulseJob, root);
    require(pulses.size() == 10, "five branches from each of two anchors");
    const auto first = core::readNpz(pulses.front()), baseRecord = core::readNpz(base);
    require(std::abs(first.at("values").doubles()[0] - baseRecord.at("values").doubles()[2 * 20]) < 1e-12, "the anchor is the archived state at t = 0.2");

    const auto response = metadataOf(postprocess("pulse_response", J::Object{{"krom", krom.generic_u8string()},
        {"bands", J::Array{J::Array{0.01, 0.1}, J::Array{0.1, 0.31}}}}, pulses, root));
    const auto& anchorRows = response.at("anchors").array();
    require(anchorRows.size() == 2, "two anchors analysed");
    const auto& velocity = anchorRows[0].at("modes").array().at(1);
    require(std::abs(velocity.at("immediate_sensitivity_magnitude").number() - 0.01) < 5e-4, "G(t_1) of the rate is dt B: " + json::dump(anchorRows[0]));
    require(velocity.at("amplitude_checks").array().at(0).at("sensitivity_change").number() < 0.05, "doubling a small pulse barely changes G");
    require(!velocity.at("even_nonlinear_over_odd_response").isNull(), "the zero branch gives the even response");

    (void)postprocess("input_map", J::Object{{"krom", krom.generic_u8string()}, {"representation", "constant"},
        {"fit", "symmetric_pulse"}, {"map_name", "pulse"}}, pulses, root);

    const auto audit = metadataOf(postprocess("control_map_audit", J::Object{{"krom", krom.generic_u8string()},
        {"maps", J::Array{"constant", "biased"}}, {"training_sources", J::Array{pulses[1].generic_u8string()}}, {"audit_stride", 1.0}},
        {base}, root));
    const auto& path = audit.at("path_error").at(base.generic_u8string());
    require(path.at("constant").at("relative_frobenius").number() < 1e-12, "the exact constant map has no path error");
    require(std::abs(path.at("biased").at("maximum_pointwise_relative").number() - 0.1) < 1e-9, "a 10 percent bias is a 10 percent error");
    const auto& design = audit.at("design").at(pulses[1].generic_u8string());
    require(design.at("constant_design_rank").number() == 1, "one input column");
    require(design.at("bilinear_design_rank").number() >= 1, "bilinear design rank reported");
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-pulses-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        pulseTest(root);
        std::filesystem::remove_all(root);
        std::cout << "pulse and control-map workflow tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
