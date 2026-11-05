#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include <chrono>
#include <cmath>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <iostream>
#include <stdexcept>
namespace core = dynsyslab::core;
namespace j = core::json;
using namespace dynsyslab::workflow;
void check(bool b, const char *text) {
    if (!b)
        throw std::runtime_error(text);
}
int main() {
    const auto root =
        std::filesystem::temp_directory_path() /
        ("dsl-cat-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        const auto cases = j::parse(R"([
            {"name":"kuramoto","parameters":{"oscillators":2,"coupling":0.7,"frequency_spread":0.5},"initial":[-0.3,0.4],"duration":20},
            {"name":"hodgkin_huxley","parameters":{"g_na":0,"g_k":0,"g_l":0,"bias_current":0},"initial":[-65,0.2,0.3,0.4],"duration":0.003},
            {"name":"spring_coupled_pendulums","parameters":{"damping1_n_m_s":0,"damping2_n_m_s":0,"spring_constant_n_m":40},"initial":[0.7,0,0.7,0],"duration":3},
            {"name":"marmottant","parameters":{},"initial":[1,0],"duration":2},
            {"name":"reddy_szeri","parameters":{"vapor_pressure_pa":0},"initial":[1,0,0,0,0.1,0.02,0.05,0,-0.03,0.01,0.02,-0.01],"duration":3}
        ])");
        for (const auto &test : cases.array()) {
            const auto name = test.at("name").string();
            const double duration = test.at("duration").number();
            auto job = j::parse(
                R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","controller":{"kind":"none"},"solver":{"method":"RK45","rtol":1e-10,"atol":1e-12}})");
            job.object()["label"] = "catalogue-limit-" + name;
            job.object()["system"] =
                j::Value::Object{{"name", name}, {"parameters", test.at("parameters")}};
            job.object()["time"] = j::Value::Object{
                {"start", 0.}, {"stop", duration}, {"output_step", duration / 100}};
            job.object()["initial_conditions"] =
                j::Value::Object{{"mode", "explicit"},
                                 {"count", 1.},
                                 {"clusters", 1.},
                                 {"values", j::Value::Array{test.at("initial")}}};
            const auto solve = executeSolve(job, root);
            const auto archive = core::readNpz(solve.archives.front());
            const auto text = archive.metadataJson();
            const auto meta = j::parse(text);
            check(meta.at("analytic_solution").at("status").string() == "complete",
                  "Missing complete raw analytic archive");
            check(text.find("paper_catalogue") == std::string::npos &&
                      text.find("doi.org") == std::string::npos &&
                      text.find("references") == std::string::npos,
                  "Bibliography leaked into a solve archive");
            check(archive.at("analytic_values").shape == archive.at("values").shape,
                  "Raw analytic state shape differs");
            const auto a = archive.at("analytic_values").doubles(),
                       b = archive.at("values").doubles();
            for (std::size_t i = 0; i < a.size(); ++i)
                check(std::abs(a[i] - b[i]) < 2e-7, "Raw analytic and numerical arrays disagree");
            const auto plots = dsl::buildNativePayload(solve.archives.front(), "time");
            check(plots.figures.size() == test.at("initial").array().size(), "State plot missing");
            if (name == "kuramoto")
                check(plots.figures[0].yLabel == R"(\theta_{1})" &&
                          plots.figures[1].yLabel == R"(\theta_{2})",
                      "Native Kuramoto component symbols missing");
            if (name == "spring_coupled_pendulums")
                check(plots.figures[1].yLabel == R"(\frac{d\theta_1}{dt})",
                      "Pendulum derivative label missing");
            for (const auto &figure : plots.figures) {
                check(figure.references.empty(), "Time series contains IC or equilibrium dots");
                check(figure.traces.size() == 2, "Missing numerical or analytic curve");
                const auto &trace = figure.traces.back();
                check(trace.role == "analytic_solution" && trace.label == "Analytic" &&
                          trace.lineStyle == "dash" && trace.color == "#7b6d9c",
                      "Viewer reference style differs");
            }
            const auto research = evaluateResearchCatalogueQuery(
                j::Value::Object{{"schema", "dynsyslab.research-query/v1"},
                                 {"action", "list"},
                                 {"system", job.at("system")}},
                root);
            check(!research.at("papers").array().empty(), "System bibliography missing");
            if (name == "kuramoto")
                (void)dsl::renderPlotPng(plots, root / "kuramoto-time.png", {1400, 900});
            std::cout << name << ": raw arrays, system-only papers and Viewer reference passed\n";
        }
        std::cout << "Evidence: " << root.u8string() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\nEvidence: " << root.u8string() << '\n';
        return 1;
    }
}
