#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TrainingDomain.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Solve.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace json = dynsyslab::core::json;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value plannedJob() {
    auto job = json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve",
      "label":"duffing-domain-expansion",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,
                            "values":[],"ranges":[]},
      "controller":{"kind":"none","settings":{},
                    "target":{"kind":"constant","state":[0,0]}},
      "preflight":{"enabled":true,"kind":"krom_training_domain",
        "source_krom":"parent-krom.npz","coordinate_system":"dimensionless",
        "current_domain":{"minimum":[-1,-0.25],"maximum":[1,0.25]},
        "required_domain":{"minimum":[-1.2,-0.25],"maximum":[1,0.5]},
        "points_per_axis":[3,3],"anchors":[[0,0]],
        "margin_fraction":0.1,"outside_current_only":true,
        "screening":{"time_end":0.01,"output_step":0.005,"maximum_step":0.002,
          "method":"RK4","rtol":1e-9,"atol":1e-12,
          "state_minimum":[],"state_maximum":[],
          "maximum_absolute_state":[10,10],"require_all":true}},
      "solver":{"method":"RK4","step_size":0.001},
      "archive":{"enabled":true,"data_role":"free_training"}
    })");
    dynsyslab::core::StateDomainExpansionOptions options;
    options.pointsPerAxis = {3,3}; options.marginFraction = 0.1;
    options.outsideCurrentOnly = true; options.anchors = {0,0}; options.anchorCount = 1;
    const auto plan = dynsyslab::core::planStateDomainExpansion(
        {{-1,-0.25},{1,0.25},"dimensionless"},
        {{-1.2,-0.25},{1,0.5},"dimensionless"}, options);
    json::Value::Array rows;
    for (std::size_t row = 0; row < plan.initialConditionCount; ++row) {
        json::Value::Array values;
        values.emplace_back(plan.initialConditions[row * 2]);
        values.emplace_back(plan.initialConditions[row * 2 + 1]);
        rows.emplace_back(std::move(values));
    }
    job.at("initial_conditions").at("count") =
        static_cast<double>(plan.initialConditionCount);
    job.at("initial_conditions").at("values") = json::Value(std::move(rows));
    return job;
}
}

int main() {
    try {
        auto job = plannedJob();
        const auto plan = dynsyslab::workflow::planJob(job);
        require(plan.at("trajectory_count").number() > 1.0,
                "preflight plan validates generated initial conditions");
        const auto root = std::filesystem::temp_directory_path()
            / "dynsyslab-training-preflight-test";
        std::error_code ignored;
        std::filesystem::remove_all(root, ignored);
        double lastFraction = 0.0;
        bool monotonic = true;
        dynsyslab::core::ExecutionContext execution;
        execution.progress = [&](const dynsyslab::core::ExecutionProgress& event) {
            monotonic = monotonic && event.fraction + 1.0e-12 >= lastFraction;
            lastFraction = event.fraction;
        };
        const auto solved = dynsyslab::workflow::executeSolve(
            job, root, root / "progress.json", &execution);
        require(!solved.archives.empty(), "screened solve publishes archives");
        require(monotonic && lastFraction == 1.0, "screening and solve progress is monotonic");
        const auto metadata = json::parse(
            dynsyslab::core::readNpz(solved.archives.front()).metadataJson());
        require(metadata.at("training_preflight_result").at("all_accepted").boolean(false),
                "solution retains passing preflight evidence");
        require(metadata.at("training_preflight_result").at("trajectory")
                    .at("accepted").boolean(false),
                "solution retains its own screening result");
        const auto hierarchy = json::load(solved.hierarchyManifest.wstring());
        require(hierarchy.at("properties").at("training_preflight").string() == "passed",
                "hierarchy retains training preflight provenance");

        auto rejected = plannedJob();
        rejected.at("preflight").at("screening").at("maximum_absolute_state") =
            json::parse("[0.1,10]");
        bool blocked = false;
        try {
            (void)dynsyslab::workflow::executeSolve(
                rejected, root / "rejected", root / "rejected-progress.json");
        } catch (const std::runtime_error&) { blocked = true; }
        require(blocked, "unsafe training-domain solve is blocked");
        require(!std::filesystem::exists(root / "rejected" / "solutions"),
                "blocked preflight publishes no solution archives");

        auto stale = plannedJob();
        stale.at("initial_conditions").at("values").array()[0].array()[0] = 99.0;
        bool staleRejected = false;
        try { (void)dynsyslab::workflow::executeSolve(stale, root, {}); }
        catch (const std::invalid_argument&) { staleRejected = true; }
        require(staleRejected, "stale edited ICs cannot bypass the generated plan");
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_training_preflight_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_training_preflight_tests: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
