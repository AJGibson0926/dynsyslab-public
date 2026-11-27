// Experiment sweep and select stages with the orbit_geometry method: three orbits that cross in
// the plane of coordinates (0, 1) and are separated in (0, 2); a sweep over the plane and a
// validation-only selection must choose (0, 2), and its record keeps every candidate.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Experiment.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
using O = json::Value::Object;
using A = json::Value::Array;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path orbit(const std::filesystem::path& root, int k) {
    std::vector<double> time, values, states;
    for (std::size_t s = 0; s < 4000; ++s) {
        const double t = double(s) * 0.01;
        time.push_back(t);
        values.insert(values.end(), {std::cos(t), std::sin(t) + 0.8 * k, 3.0 * k + 0.1 * std::sin(t)});
        states.insert(states.end(), {double(k), std::cos(t)});
    }
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({time.size()}, time);
    record.arrays["values"] = core::NpyArray::fromDoubles({time.size(), 3}, values);
    record.arrays["source_states"] = core::NpyArray::fromDoubles({time.size(), 2}, states);
    record.arrays["metadata_json"] = core::NpyArray::fromString(R"({"schema":"dynsyslab.trajectory/v1","label":"orbit"})");
    const auto path = root / ("orbit-" + std::to_string(k) + ".npz");
    core::writeNpzAtomic(path, record);
    return path;
}

json::Value readJson(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    std::stringstream text; text << stream.rdbuf();
    return json::parse(text.str());
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-sweep-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        A orbits;
        for (int k = 0; k < 3; ++k) orbits.push_back(orbit(root, k).generic_u8string());
        const json::Value geometry = O{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"},
            {"sources", "${orbits.archives}"},
            {"methods", A{O{{"name", "orbit_geometry"}, {"settings", O{{"role", "validation"}, {"axes", "${sweep.axes}"},
                {"groups", A{"inside", "inside", "outside"}}, {"cycle_samples", "${sweep.samples}"}}}}}}};
        const json::Value experiment = O{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "experiment"}, {"label", "sweep-test"},
            {"stages", A{
                O{{"id", "orbits"}, {"job", O{{"kind", "archive"}, {"archives", orbits}}}},
                O{{"id", "planes"}, {"sweep", O{{"variables", O{{"axes", A{A{0.0, 1.0}, A{0.0, 2.0}}}, {"samples", A{200.0, 300.0}}}},
                    {"mode", "zip"}, {"stages", A{O{{"id", "score"}, {"job", geometry}}}}}}},
                O{{"id", "best"}, {"select", O{{"sweep", "planes"}, {"stage", "score"}, {"metric", "validation.crossings"},
                    {"goal", "minimize"}}}},
                O{{"id", "after"}, {"job", O{{"kind", "archive"}, {"archives", A{"${best.archive}"}}}}}}}};

        const auto plan = dynsyslab::workflow::planExperiment(experiment);
        require(plan.at("stage_count").number() == 5, "a two-point sweep expands to two stages");
        const auto result = dynsyslab::workflow::executeJob(experiment, root);
        const auto manifest = readJson(result.manifest);
        const auto& stages = manifest.at("stages").array();
        require(stages.size() == 5 && stages[1].at("id").string() == "planes_1_score", "sweep stage ids");
        const auto crossing = json::parse(core::readNpz(std::filesystem::u8path(stages[1].at("archives").array()[0].string())).metadataJson());
        require(crossing.at("validation").at("crossings").number() >= 4, "orbits cross in the (0, 1) plane");
        require(crossing.at("validation").at("crossings_by_group").at("inside_inside").number() >= 2, "crossings by group");
        require(crossing.at("validation").at("state_ambiguity_maximum").number() > 0.1, "crossings join different states");
        const auto& select = stages[3];
        require(select.at("kind").string() == "select", "select stage");
        const auto record = readJson(std::filesystem::u8path(select.at("manifest").string()));
        require(record.at("winner").at("stage").string() == "planes_2_score", "the separated plane wins");
        require(record.at("winner").at("score").number() == 0, "no crossings in the (0, 2) plane");
        require(record.at("candidates").array()[0].at("rejected").string() == "outscored", "rejects recorded");
        require(stages[4].at("archives").array()[0].string() == select.at("archives").array()[0].string(),
                "${best.archive} is the winner's archive");

        {   // stands_for: the winner's archive of another stage of its point
            auto other = experiment;
            other.at("stages").array()[1].at("sweep").at("stages").array().push_back(
                O{{"id", "copy"}, {"job", O{{"kind", "archive"}, {"archives", A{"${orbits.archives[0]}"}}}}});
            other.at("stages").array()[2].at("select").at("stands_for") = json::Value("copy");
            const auto again = dynsyslab::workflow::executeJob(other, root);
            const auto stagesAgain = readJson(again.manifest).at("stages").array();
            const auto chosen = std::find_if(stagesAgain.begin(), stagesAgain.end(), [](const json::Value& s) { return s.at("id").string() == "best"; });
            require(chosen != stagesAgain.end() && chosen->at("archives").array()[0].string() == orbits[0].string(), "stands_for result");
        }
        bool refused = false;  // training scores need an explicit opt-in
        auto training = experiment;
        training.at("stages").array()[2].at("select").at("metric") = json::Value("training.crossings");
        try { (void)dynsyslab::workflow::planExperiment(training); } catch (const std::invalid_argument&) { refused = true; }
        require(refused, "select refuses a training score");
        std::filesystem::remove_all(root);
        std::cout << "experiment sweep tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
