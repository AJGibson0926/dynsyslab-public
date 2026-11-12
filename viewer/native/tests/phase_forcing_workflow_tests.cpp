#include "../src/PhaseContext.hpp"
#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include <chrono>
#include <cmath>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <iostream>
#include <stdexcept>
namespace core = dynsyslab::core;
namespace j = core::json;
using namespace dynsyslab::workflow;
void check(bool ok, const char *text) {
    if (!ok)
        throw std::runtime_error(text);
}
int main() {
    const auto root = std::filesystem::temp_directory_path() /
                      ("dsl-force-" +
                       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        auto job = j::parse(
            R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"Periodic forcing fixture",
          "system":{"name":"harmonic_oscillator","parameters":{"damping_kg_s":0.4}},
          "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[0,0]]},
          "controller":{"kind":"waveform","settings":{"amplitude":0.3,"frequency":1,"phase_radians":0.4}},
          "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12},"time":{"start":0,"stop":6,"output_step":0.01}})");
        const auto solved = executeSolve(job, root);
        const auto archived = core::readNpz(solved.archives.front());
        const auto metadata = j::parse(archived.metadataJson());
        const auto originalMetadata = archived.metadataJson();
        check(originalMetadata.find("paper_catalogue") == std::string::npos &&
                  originalMetadata.find("references") == std::string::npos &&
                  originalMetadata.find("doi.org") == std::string::npos,
              "Papers leaked into raw solution metadata");
        dsl::Backend backend({root, root, {}, root / "viewer-cache", root / "no-mask.json"});
        dsl::CatalogItem item;
        item.source = solved.archives.front().u8string();
        item.kind = "trajectory";
        auto phase = backend.loadPayload(item, "phase");
        check(!std::filesystem::exists(root / "system_metadata"),
              "Opening a trajectory computed the sidecar");
        auto options = j::parse(
            R"({"lower":[-0.2,-0.3],"upper":[0.2,0.3],"periodic_seeds":[{"state":[0,0]}],"forcing_section_time":0.137,"manifolds":false,"vector_grid_per_axis":3})");
        auto missing = backend.queryPhaseSpace(item, "lookup", std::nullopt, options);
        check(missing.at("results").array()[0].at("result").at("cache").at("status").string() ==
                  "miss",
              "Expected uncomputed forcing context");
        auto computed = backend.queryPhaseSpace(item, "compute", std::nullopt, options);
        const auto &entry = computed.at("results").array()[0];
        const auto &result = entry.at("result");
        check(result.at("cache").at("status").string() == "computed",
              "Periodic context not computed");
        check(result.at("options").at("dynamics").string() == "periodic-forcing",
              "Viewer lost archived periodic forcing");
        check(result.at("fixed_points").array().empty(),
              "Periodic response mislabeled equilibrium");
        check(result.at("periodic_orbits").array().size() == 1, "Verified response missing");
        const auto &orbit = result.at("periodic_orbits").array()[0];
        constexpr double w = 6.283185307179586476925286766559;
        const double t = .137, k = 1 - w * w, c = .4 * w, denominator = k * k + c * c;
        const double exact =
            .3 * (k * std::sin(w * t + .4) - c * std::cos(w * t + .4)) / denominator;
        const double velocity =
            .3 * w * (k * std::cos(w * t + .4) + c * std::sin(w * t + .4)) / denominator;
        check(std::abs(orbit.at("section_state").array()[0].number() - exact) < 1e-7 &&
                  std::abs(orbit.at("section_state").array()[1].number() - velocity) < 1e-7,
              "Archived sine controller phase or units changed during reconstruction");
        check(result.at("forcing_clock").at("section_time").number() == t &&
                  std::abs(result.at("forcing_clock").at("period").number() - 1) < 1e-14,
              "Forcing clock metadata absent");
        auto again = backend.queryPhaseSpace(item, "lookup", std::nullopt, options);
        check(again.at("results").array()[0].at("result").at("cache").at("status").string() ==
                  "hit",
              "Forced context cache not reused");
        options.object()["forcing_section_time"] = .2;
        auto other = backend.queryPhaseSpace(item, "lookup", std::nullopt, options);
        check(other.at("results").array()[0].at("result").at("cache").at("status").string() ==
                  "miss",
              "Forcing section missing from cache identity");
        dsl::appendPhaseContext(phase.figures[0], metadata, result);
        const auto plotBounds = dsl::automaticAxisLimits(phase.figures[0]);
        for (const auto &trace : phase.figures[0].traces)
            if (trace.role == "phase_vector_field")
                for (std::size_t i = 0; i + 1 < trace.x.size(); i += 3)
                    check(std::abs(trace.x[i + 1] - trace.x[i]) <=
                                  .151 * (plotBounds.xMax - plotBounds.xMin) &&
                              std::abs(trace.y[i + 1] - trace.y[i]) <=
                                  .151 * (plotBounds.yMax - plotBounds.yMin),
                          "Arrows obscure a plot smaller than the analysis domain");
        bool response = false;
        for (const auto &trace : phase.figures[0].traces)
            if (trace.role == "phase_periodic_orbit")
                response = trace.label == "periodic response";
        check(response, "Viewer lacks periodic-response overlay");
        for (const auto &point : phase.figures[0].references)
            check(!point.equilibrium, "Viewer called a forced response an equilibrium");
        auto series = backend.loadPayload(item, "time");
        for (auto &figure : series.figures) {
            const auto count = figure.traces.size();
            dsl::appendPhaseContext(figure, metadata, result);
            check(figure.traces.size() == count && figure.references.empty(),
                  "Phase context leaked into time plots");
        }
        dsl::Figure heatmap;
        heatmap.kind = "observable_heatmap";
        dsl::appendPhaseContext(heatmap, metadata, result);
        check(heatmap.traces.empty(), "Phase context leaked into observable heat map");
        auto request = entry.at("request");
        request.object()["mode"] = "lookup";
        request.object()["controller"].object()["settings"].object()["active_time_window"] =
            j::Value::Array{0., 10.};
        bool rejected = false;
        try {
            (void)evaluatePhaseSpaceQuery(request, root);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        check(rejected, "Windowed control silently became globally periodic");
        check(core::readNpz(solved.archives.front()).metadataJson() == originalMetadata,
              "Sidecar changed an existing solve archive");
        job.object()["label"] = "zero-waveform";
        job.object()["controller"].object()["settings"].object()["amplitude"] = 0.;
        auto zeroSolve = executeSolve(job, root);
        auto zeroItem = item; zeroItem.source = zeroSolve.archives.front().u8string();
        auto zeroOptions = j::parse(R"({"lower":[-1,-1],"upper":[1,1],"periodic_orbits":false,"manifolds":false,"vector_field":false,"grid_per_axis":3})");
        auto zeroContext = backend.queryPhaseSpace(zeroItem, "compute", std::nullopt, zeroOptions);
        const auto &zero = zeroContext.at("results").array()[0].at("result");
        check(zero.at("options").at("dynamics").string() == "autonomous-zero-input" &&
              zero.at("fixed_points").array().size() == 1, "Zero-amplitude control lost autonomous context");
        (void)dsl::renderPlotPng(phase, root / "forced-phase.png", {1200, 900});
        std::cout << "Archived waveform, fixed-phase cache, raw metadata and Viewer overlays: "
                     "PASS\nEvidence: "
                  << root.u8string() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\nEvidence: " << root.u8string() << '\n';
        return 1;
    }
}
