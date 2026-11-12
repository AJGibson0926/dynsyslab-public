#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include <chrono>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/DelaySolve.hpp>
#include <dynsyslab/workflow/PDESolve.hpp>
#include <iostream>
#include <stdexcept>
namespace {
namespace j = dynsyslab::core::json;
void check(bool value, const char *reason) {
    if (!value)
        throw std::runtime_error(reason);
}
void ownsData(const dynsyslab::core::NpzRecord &record) {
    const auto m = j::parse(record.metadataJson());
    check(record.contains("analytic_values") &&
              record.at("analytic_values").shape == record.at("values").shape,
          "Analytic arrays missing or incompatible");
    check(m.at("dynamical_knowledge").isNull() &&
              m.at("analytic_solution").at("solution").at("references").isNull(),
          "Papers leaked into a solution archive");
}
} // namespace
int main() {
    namespace core = dynsyslab::core;
    using namespace dynsyslab::workflow;
    const auto root = std::filesystem::temp_directory_path() /
                      ("dsl-nonode-" +
                       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        auto job =
            j::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"analytic-heat",
            "system":{"name":"heat_equation_1d","parameters":{"points":17}},
            "time":{"start":0,"stop":0.04,"output_step":0.02},
            "initial_conditions":{"mode":"default","count":1,"clusters":1},
            "controller":{"kind":"none"},"solver":{"method":"RK4","step_size":0.002}})");
        const auto solved = executePDESolve(job, root);
        const auto record = core::readNpz(solved.archives.front());
        ownsData(record);
        const auto metadata = j::parse(record.metadataJson());
        check(metadata.at("analytic_solution").at("status").string() == "complete",
              "Heat reference is incomplete");
        auto payload = dsl::buildNativePayload(solved.archives.front(), "field");
        check(payload.figures.size() == 2 && payload.figures[1].traces.size() == 2,
              "Field profile lacks analytic overlay");
        const auto &overlay = payload.figures[1].traces.back();
        check(overlay.role == "analytic_solution" && overlay.lineStyle == "dash" &&
                  overlay.label == "Analytic",
              "Field analytic style missing");
        check(payload.figures[1].references.empty(), "Field profile has state markers");
        (void)dsl::renderPlotPng(payload, root / "heat-analytic.png", {});
        auto query = j::parse(
            R"({"schema":"dynsyslab.analytic-query/v1","system":{"name":"heat_equation_1d","parameters":{"points":17}},"time":[0,0.02,0.04]})");
        j::Value::Array initial;
        for (double x : record.at("values").doubles()) {
            initial.emplace_back(x);
            if (initial.size() == 17)
                break;
        }
        query.object()["initial_state"] = initial;
        auto result = evaluateAnalyticQuery(query);
        check(result.at("status").string() == "complete" &&
                  result.at("system_family").string() == "pde",
              "PDE analytic CLI query failed");
        job.object()["label"] = "analytic-delay";
        job.object()["system"] = j::parse(
            R"({"name":"mackey_glass","parameters":{"exponent":0,"beta":0.3,"delay":0.02}})");
        job.at("time").object()["stop"] = .1;
        job.at("solver").object()["method"] = "Euler";
        auto delayed = executeDelaySolve(job, root);
        auto delayRecord = core::readNpz(delayed.archives.front());
        ownsData(delayRecord);
        check(j::parse(delayRecord.metadataJson()).at("analytic_solution").at("status").string() ==
                  "complete",
              "Linear delay reference is incomplete");
        auto delayPayload = dsl::buildNativePayload(delayed.archives.front(), "time");
        check(delayPayload.figures.front().traces.size() == 2 &&
                  delayPayload.figures.front().references.empty(),
              "Delay curves have missing reference or unwanted IC dots");
        job.object()["label"] = "analytic-delay-partial";
        job.at("system").at("parameters").object()["exponent"] = 10.0;
        delayed = executeDelaySolve(job, root);
        delayRecord = core::readNpz(delayed.archives.front());
        ownsData(delayRecord);
        check(j::parse(delayRecord.metadataJson()).at("analytic_solution").at("status").string() ==
                  "partial",
              "Nonlinear first interval was extrapolated");
        delayPayload = dsl::buildNativePayload(delayed.archives.front(), "time");
        check(!delayPayload.warnings.empty(), "Viewer omitted the first-delay validity boundary");
        job.object()["label"] = "analytic-delay-off";
        job.at("solver").object()["analytic_solution"] = "off";
        delayed = executeDelaySolve(job, root);
        delayRecord = core::readNpz(delayed.archives.front());
        check(!delayRecord.contains("analytic_values") && j::parse(delayRecord.metadataJson())
                                                                  .at("analytic_solution")
                                                                  .at("status")
                                                                  .string() == "disabled",
              "Delay analytic off setting was ignored");
        std::cout << "PDE/DDE query, raw archiving, bibliography ownership and Viewer overlays "
                     "passed\nEvidence: "
                  << root.u8string() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\nEvidence: " << root.u8string() << '\n';
        return 1;
    }
}
