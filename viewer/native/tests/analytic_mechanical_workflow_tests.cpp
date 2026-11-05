#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include <chrono>
#include <cmath>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <iostream>
#include <stdexcept>
namespace core = dynsyslab::core;
namespace j = core::json;
using namespace dynsyslab::workflow;
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
int main() {
    const auto root = std::filesystem::temp_directory_path() /
                      ("dsl-analytic-mechanics-" +
                       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        auto job = j::parse(
            R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"gravity-free-double-pendulum",
            "system":{"name":"double_pendulum","parameters":{"gravity_m_s2":0,"damping1_n_m_s":0,"damping2_n_m_s":0}},
            "time":{"start":0.37,"stop":8.37,"output_step":0.02},
            "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[0.6,0.5,0,0.5]]},
            "controller":{"kind":"none"},"solver":{"method":"RK45","rtol":1e-10,"atol":1e-12}})");
        const auto solved = executeSolve(job, root);
        const auto archive = core::readNpz(solved.archives.front());
        const auto metadata = j::parse(archive.metadataJson());
        check(archive.at("analytic_values").shape == archive.at("values").shape,
              "Raw analytic shape mismatch");
        check(metadata.at("analytic_solution").at("status").string() == "complete",
              "Incomplete analytic reference");
        check(metadata.at("dynamical_knowledge").isNull() &&
                  metadata.at("analytic_solution").at("solution").at("references").isNull(),
              "Paper metadata leaked into solution");
        const auto numeric = archive.at("values").doubles(),
                   analytic = archive.at("analytic_values").doubles();
        for (std::size_t i = 0; i < numeric.size(); ++i)
            check(std::abs(numeric[i] - analytic[i]) < 2e-7,
                  "Archived reference disagrees with numerical solve");
        const auto time = dsl::buildNativePayload(solved.archives.front(), "time");
        check(time.figures.size() == 4, "Expected four state time series");
        check(time.figures[0].yLabel == R"(\theta_1)" &&
                  time.figures[1].yLabel == R"(\frac{d\theta_1}{dt})",
              "System-owned LaTeX labels missing");
        for (const auto &figure : time.figures) {
            check(figure.traces.size() == 2 && figure.references.empty(),
                  "Time series must show two curves and no IC/equilibrium dots");
            const auto &trace = figure.traces.back();
            check(trace.role == "analytic_solution" && trace.lineStyle == "dash" &&
                      trace.label == "Analytic" && trace.color == "#7b6d9c",
                  "Analytic Viewer style missing");
        }
        const auto phase = dsl::buildNativePayload(solved.archives.front(), "phase");
        check(!phase.figures.empty(), "Phase plot missing");
        for (const auto &figure : phase.figures) {
            bool found = false;
            for (const auto &trace : figure.traces)
                found = found || trace.role == "analytic_solution";
            check(found, "Phase plot missing analytic reference");
        }
        (void)dsl::renderPlotPng(time, root / "double-pendulum-time.png", {1400, 900});
        (void)dsl::renderPlotPng(phase, root / "double-pendulum-phase.png", {1400, 900});
        const auto research = evaluateResearchCatalogueQuery(
            j::parse(
                R"({"schema":"dynsyslab.research-query/v1","action":"list","system":{"name":"double_pendulum"}})"),
            root);
        check(research.at("papers").array().size() == 2,
              "System bibliography missing double-pendulum papers");
        const auto special = evaluateSpecialFunctionQuery(j::parse(
            R"({"schema":"dynsyslab.special-function-query/v1","function":"lauricella_fd","a":1,"b":[1],"c":2,"x":[0.5]})"));
        check(std::abs(special.at("value").number() - 2 * std::log(2.0)) < 1e-12,
              "Lauricella public query mismatch");
        auto query = j::parse(
            R"({"schema":"dynsyslab.phase-space-query/v1","mode":"compute","system":{"name":"double_pendulum"},"options":{"lower":[-4,-1,-4,-1],"upper":[4,1,4,1],"manifolds":false,"periodic_orbits":false,"vector_field":false}})");
        const auto context = evaluatePhaseSpaceQuery(query, root);
        check(context.at("cache").at("status").string() == "computed", "Phase metadata not cached");
        std::size_t verified = 0;
        for (const auto &guide : context.at("knowledge_checks").array())
            if (guide.at("status").string() == "verified")
                ++verified;
        check(verified == 9,
              "Nine vertical equilibrium candidates must be verified against native equations");
        query.object()["mode"] = "lookup";
        check(evaluatePhaseSpaceQuery(query, root).at("cache").at("status").string() == "hit",
              "Phase metadata cache was not reused");
        job.object()["label"] = "gravity-free-disabled";
        job.at("solver").object()["analytic_solution"] = "off";
        const auto disabled = executeSolve(job, root);
        check(!core::readNpz(disabled.archives.front()).contains("analytic_values"),
              "Analytic opt-out ignored");
        std::cout << "Analytic mechanics, raw archives, Viewer, bibliography ownership and checked "
                     "phase guides passed\nEvidence: "
                  << root.u8string() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\nEvidence: " << root.u8string() << '\n';
        return 1;
    }
}
