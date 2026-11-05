#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TimeSeriesArchive.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-solve-set-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    try {
        const auto job = dynsyslab::core::json::parse(R"({
          "schema":"dynsyslab.compute-job/v1","kind":"solve_set","label":"tree-job",
          "batches":[{"label":"Batch 1","systems":[
            {"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"Duffing",
             "ui_system_label":"Duffing oscillator",
             "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
             "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
             "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,"values":[],"ranges":[]},
             "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0,0]}},
             "campaign":{"enabled":false},"preflight":{"enabled":false},
             "solver":{"method":"RK4"},"archive":{"enabled":true,"data_role":"free"}},
            {"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"Van der Pol",
             "ui_system_label":"Van der Pol oscillator",
             "system":{"name":"van_der_pol","preset":null,"parameters":{}},
             "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
             "initial_conditions":{"mode":"default","count":2,"clusters":1,"seed":0,"values":[],"ranges":[]},
             "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0,0]}},
             "campaign":{"enabled":false},"preflight":{"enabled":false},
             "solver":{"method":"RK4"},"archive":{"enabled":true,"data_role":"free"}}
          ]}]})");
        const auto plan = dynsyslab::workflow::planJob(job);
        require(!std::filesystem::exists(root),
                "planning a proposed solve tree created archive state");
        require(plan.at("kind").string() == "solve_set", "solve-set plan kind");
        require(plan.at("batch_count").number() == 1.0, "solve-set batch count");
        require(plan.at("configured_system_count").number() == 2.0,
                "solve-set configured system count");
        require(plan.at("trajectory_count").number() == 3.0,
                "solve-set trajectory count");
        const auto run = dynsyslab::workflow::executeJob(job, root);
        require(run.kind == "solve_set" && run.archives.size() == 3,
                "solve-set execution result");
        for (const auto& path : run.archives) {
            const auto archive = dynsyslab::core::readNpz(path);
            const auto metadata = dynsyslab::core::json::parse(
                archive.at("metadata_json").string());
            const auto series = dynsyslab::core::loadTimeSeries(archive);
            require(metadata.at("time_series_format").string() ==
                    "dynsyslab.TimeSeries/v2" && series.sampleCount() == 3,
                    "native solve did not publish a canonical TimeSeries");
            require(!archive.contains("frequencies") &&
                    !archive.contains("fourier_amplitudes"),
                    "native solve persisted derived Fourier data");
        }
        const auto manifest = dynsyslab::core::json::load(run.manifest.wstring());
        require(manifest.at("format").string() == "dynsyslab.solve-set/v1",
                "solve-set manifest format");
        require(manifest.at("hierarchy_manifests").array().size() == 1,
                "one hierarchy per batch");
        const auto hierarchy = dynsyslab::core::json::load(std::filesystem::u8path(
            manifest.at("hierarchy_manifests").array()[0].string()).wstring());
        require(hierarchy.at("systems").array().size() == 2,
                "batch hierarchy did not merge its systems");
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_solve_set_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_solve_set_tests: FAIL: " << error.what() << '\n';
        std::filesystem::remove_all(root, ignored);
        return 1;
    }
}
