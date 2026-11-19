// `empirical_observability` workflow stage: a forced linear Duffing solve, then local
// observability through `position` along it (applied input replayed), compared with the
// exact finite-horizon gramian, which does not depend on the input or the state.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path() / ("dynsyslab-empirical-observability-" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto solve = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"linear duffing",
            "system":{"name":"duffing","parameters":{"alpha":-1,"beta":0,"delta":0.5}},
            "time":{"start":0,"stop":8,"output_step":0.01},
            "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[0.4,0]]},
            "controller":{"kind":"waveform","settings":{"frequencies":[0.3],"model_amplitudes":[0.5],"phases_rad":[0]}},
            "solver":{"method":"RK45","rtol":1e-9,"atol":1e-11},"archive":{"enabled":true}})");
        const auto run = dynsyslab::workflow::executeJob(solve, root);
        const auto source = run.archives.front();
        const json::Value job = json::Value::Object{
            {"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", "observability"},
            {"sources", json::Value::Array{source.generic_u8string()}},
            {"methods", json::Value::Array{json::Value::Object{{"name", "empirical_observability"}, {"enabled", true},
                {"settings", json::Value::Object{{"observable", "position"}, {"start", 1.0}, {"stop", 5.0},
                    {"points", 5.0}, {"horizon", 2.0}, {"samples", 401.0}, {"substeps", 10.0},
                    {"noise_intensity", 0.01}}}}}}};
        const auto result = dynsyslab::workflow::executeJob(job, root);
        const auto record = core::readNpz(result.archives.front());
        for (const char* name : {"sample_times", "sample_states", "gramians", "eigenvalues", "eigenvectors",
                                 "condition_number", "scaled_eigenvalues", "estimation_bound_std",
                                 "relative_estimation_bound", "output_jacobian", "valid", "state_scale"})
            require(record.contains(name), name);
        require(record.at("sample_times").shape[0] == 5 && record.at("gramians").shape[1] == 2, "five 2x2 gramians");
        const auto exact = core::finiteObservabilityGramian(
            core::makeStateSpace({0, 1, -1, -0.5}, {0, 1}, {1, 0}, {}, 2, 1, 1), 2.0);
        const auto g = record.at("gramians").doubles();
        for (std::size_t m = 0; m < 5; ++m)
            for (std::size_t k = 0; k < 4; ++k)
                require(std::abs(g[m * 4 + k] - exact[k]) < 1e-5 * std::abs(exact[0]), "linear plant: input-independent exact gramian");
        const auto meta = json::parse(record.metadataJson());
        require(meta.at("schema").string() == "dynsyslab.empirical-observability/v1", "schema");
        require(meta.at("input").string() == "recorded applied input", "forced run replays its input");
        require(meta.at("invalid_points").number() == 0, "all points valid");
        bool rejected = false;  // the window plus horizon must lie inside the run
        try {
            auto bad = job; bad.at("methods").array()[0].at("settings").at("stop") = 7.5;
            (void)dynsyslab::workflow::executeJob(bad, root);
        } catch (const std::exception&) { rejected = true; }
        require(rejected, "window past the end of the run rejected");
        std::filesystem::remove_all(root);
        std::cout << "dynsyslab_empirical_observability_postprocess_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        if (!root.empty()) { std::error_code ignored; std::filesystem::remove_all(root, ignored); }
        std::cerr << "dynsyslab_empirical_observability_postprocess_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
