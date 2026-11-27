#include "../src/PhaseContext.hpp"
#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include <chrono>
#include <cmath>
#include <dynsyslab/core/Archive.hpp>
#include <iostream>
#include <stdexcept>
namespace {
void require(bool ok, const char *text) {
    if (!ok)
        throw std::runtime_error(text);
}
} // namespace
int main() {
    namespace core = dynsyslab::core;
    namespace j = core::json;
    const auto root =
        std::filesystem::temp_directory_path() /
        ("dynsyslab-phase-viewer-" +
         std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    try {
        const auto source = root / "oscillator.npz";
        std::vector<double> time, states;
        for (size_t i = 0; i <= 1000; ++i) {
            const double t = i * .01;
            time.push_back(t);
            states.push_back(std::cos(t));
            states.push_back(-std::sin(t));
        }
        auto metadata = j::parse(
            R"({"schema":"dynsyslab.solution/v1","system":"harmonic_oscillator","integration_form":"dimensional","label":"Analytic and phase context fixture"})");
        core::NpzRecord record;
        record.arrays["values"] = core::NpyArray::fromDoubles({time.size(), 2}, states);
        record.arrays["time"] = core::NpyArray::fromDoubles({time.size()}, time);
        record.arrays["analytic_values"] = record.at("values");
        record.arrays["metadata_json"] = core::NpyArray::fromString(j::dump(metadata));
        core::writeNpzAtomicNew(source, record);
        dsl::Backend backend({root, root, {}, root / "cache", root / "no-mask.json"});
        dsl::CatalogItem item;
        item.source = source.u8string();
        item.kind = "trajectory";
        auto phase = backend.loadPayload(item, "phase");
        require(phase.figures.size() == 1, "Phase figure missing");
        require(phase.figures[0].traces.size() == 2, "Analytic overlay missing");
        require(phase.figures[0].traces[1].role == "analytic_solution" &&
                    phase.figures[0].traces[1].lineStyle == "dash",
                "Analytic style missing");
        require(!std::filesystem::exists(root / "system_metadata"),
                "Viewing triggered phase computation");
        const auto options = j::parse(
            R"({"lower":[-1.1,-1.1],"upper":[1.1,1.1],"periodic_orbits":false,"manifolds":false,"grid_per_axis":3,"vector_grid_per_axis":5})");
        auto computed = backend.queryPhaseSpace(item, "compute", std::nullopt, options);
        const auto &entry = computed.at("results").array().at(0);
        const auto &result = entry.at("result");
        require(result.at("cache").at("status").string() == "computed",
                "Explicit context did not compute");
        require(result.at("fixed_points").array().size() == 1, "Oscillator equilibrium missing");
        auto again = backend.queryPhaseSpace(item, "lookup", std::nullopt, options);
        require(again.at("results").array()[0].at("result").at("cache").at("status").string() ==
                    "hit",
                "Phase cache was not reused");
        const auto bounds = dsl::automaticAxisLimits(phase.figures[0]);
        dsl::appendPhaseContext(phase.figures[0], metadata, result);
        require(phase.figures[0].traces.back().role == "phase_vector_field",
                "Sparse arrows missing");
        require(std::any_of(phase.figures[0].references.begin(),
            phase.figures[0].references.end(), [](const auto& reference) {
                return reference.equilibrium && reference.label == "EC";
            }), "phase-context equilibria are not labeled EC");
        const auto after = dsl::automaticAxisLimits(phase.figures[0]);
        require(bounds.xMin == after.xMin && bounds.xMax == after.xMax &&
                    bounds.yMin == after.yMin && bounds.yMax == after.yMax,
                "Context changed trajectory bounds");
        auto series = backend.loadPayload(item, "time");
        for (auto &figure : series.figures) {
            const auto count = figure.traces.size();
            dsl::appendPhaseContext(figure, metadata, result);
            require(figure.traces.size() == count && figure.references.empty(),
                    "Context leaked into time series");
        }
        dsl::Figure heatmap;
        heatmap.kind = "observable_heatmap";
        dsl::appendPhaseContext(heatmap, metadata, result);
        require(heatmap.traces.empty() && heatmap.references.empty(),
                "Context leaked into observable heat map");
        auto meshResult=j::parse(R"({"manifold_surfaces":[{"kind":"stable_manifold_surface","status":"resolved","vertices":[0,0,0,1,0,0,0,1,1,-1,0,0,2,0,0,0,2,4,-2,0,0,3,0,0,0,3,9,-3,0,0],"front_offsets":[1,4,7]}]})");
        auto surfaceFigure=phase.figures[0];const auto beforeSurface=dsl::automaticAxisLimits(surfaceFigure);
        const auto previousCount=surfaceFigure.traces.size();
        dsl::appendPhaseSurfaces(surfaceFigure,meshResult,3,1,1);
        require(surfaceFigure.traces.size()==previousCount+1 && surfaceFigure.traces.back().label=="stable surface (projection)","Sparse surface projection missing");
        require(dsl::automaticAxisLimits(surfaceFigure).xMax==beforeSurface.xMax,"Surface projection changed trajectory bounds");
        meshResult.at("manifold_surfaces").array()[0].object()["status"]="unresolved";
        dsl::appendPhaseSurfaces(surfaceFigure,meshResult,3,1,1);
        require(surfaceFigure.traces.size()==previousCount+1,"Unresolved surface was drawn");
        dsl::appendPhaseSurfaces(heatmap,meshResult,3,1,1);
        require(heatmap.traces.empty(),"Surface leaked into observable heat map");
        const auto rendered = dsl::renderPlotPng(phase, root / "phase-context.png", {});
        require(std::filesystem::file_size(root / "phase-context.png") > 1000,
                "Background render absent");
        std::cout << "Phase context cache, analytic overlays, state-only scope and background "
                     "rendering passed\nEvidence: "
                  << root.u8string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << "\nEvidence: " << root.u8string() << '\n';
        return 1;
    }
}
