// Document payloads (periodic comparison, phase covers, research figures), SVG and interactive HTML
// export, multi-page reports (HTML and PDF), and the bifurcation, observability and system-default
// presentation.
#include "../src/ArchiveMatrixPayload.hpp"
#include "../src/DocumentPayload.hpp"
#include "../src/PlotSvg.hpp"
#include "../src/ReportDocument.hpp"
#include "../src/SystemViewDefaults.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
namespace fs = std::filesystem;
namespace json = dynsyslab::core::json;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
constexpr double kPi = 3.14159265358979323846;

json::Value circle(std::size_t samples) {
    json::Value::Array values;
    for (std::size_t i = 0; i < samples; ++i) {
        values.emplace_back(10 + std::cos(2 * kPi * i / samples));
        values.emplace_back(20 + std::sin(2 * kPi * i / samples));
    }
    return json::Value(values);
}

void documents() {
    const auto comparison = json::parse(R"({"method":"periodic-spectral-comparison","angles":[-1,0,1],
        "validation_frequencies":[1.0,2.0],"validation_families":["circle","circle"],"packet_centers":[0.2,0.4],
        "references":[{"name":"position","channel":0,"density":[0,1,0],"packet_rms":[1,2,3,4]}],
        "comparisons":[{"channel":0,"model":"mpEDMD","rank":4,"density":[0,0.9,0],"packet_rms":[1,2,3,4]}]})");
    require(dsl::isDocumentPayload(comparison), "comparison document recognized");
    const auto pages = dsl::buildDocumentPayload(comparison, "frequency");
    require(pages.pageLabels.size() == 2 && pages.figures.size() == 4, "one page per channel, family and centre");
    require(pages.figures[0].valueScale == "symlog" && pages.figures[1].traces.size() == 2, "density and packet panels");

    json::Value::Object geometry{{"schema", json::Value("dynsyslab.phase-cover-geometry/v1")}, {"prime", json::Value(2.0)},
                                 {"samples", json::Value(32.0)}};
    json::Value::Array orbits;
    for (int j = 0; j < 3; ++j)
        orbits.emplace_back(json::Value::Object{{"state_dimension", json::Value(2.0)}, {"values", circle(32)}});
    geometry["orbits"] = json::Value(orbits);
    const auto cover = dsl::buildDocumentPayload(json::Value(geometry), {});
    require(cover.figures.size() == 3 && cover.figures[2].traces.size() == 1 + 4, "level 2 has four phase points");
    geometry["display"] = json::Value(json::Value::Object{{"digits", json::parse("[1,1]")},
                                                          {"elapsed_base_periods", json::parse("[0, 1, 2]")}});
    const auto flow = dsl::buildDocumentPayload(json::Value(geometry), {});
    require(flow.animation && flow.pageLabels.size() == 3 && flow.figures.size() == 9, "the clock animates one frame per time");
    require(flow.pageLabels[1].find("residues 0, 0") != std::string::npos, "one base turn carries digits (1, 1) to (0, 0)");

    const char* kinds[] = {
        R"({"figure":"cycle_spectrum_sequence","records":[{"label":"c","atoms":[0,3.14],"eigenvalues":{"resdmd":[[1,0],[-1,0]]}}]})",
        R"({"figure":"return_sections","records":[{"label":"r","states":[[0,1],[1,0],[0,-1]]}]})",
        R"({"figure":"prediction_comparison","record":{"label":"p","truth":[0,1,2],"predictions":{"mpedmd":[0,1,2.1]}}})",
        R"({"figure":"regularized_measures","records":[{"label":"m","measures":{"mpedmd":{"grid":[0,1],"density":[1,-1]}}}]})",
        R"({"figure":"multiplier_branch","branch":{"parameter":[0,1],"multiplier":[-0.5,-1.2],"bracket":[0.4,0.6]},"local":{"parameter":[0.4,0.6],"multiplier":[-0.9,-1.1]}})",
        R"({"figure":"perturbation","records":[{"label":"q","time":[0,1],"reference_deviation":[[0,0],[1,1]],"models":[{"label":"m","deviation":[[0,0],[1,1.1]],"error_norm":[0,0.1]}]}]})",
        R"({"figure":"truncation_amplitude","records":[{"label":"t","quotients":[2,4],"spectral_radii":[0.9,1.1],"perturbations":[0.01,0.1],"rms_errors":[0.001,0.1]}]})",
        R"({"figure":"cyclic_partition","depth":1,"record":{"label":"k","states":[[0,0],[1,1],[2,0]],"geometry":{"levels":[{"depth":0,"cells":1,"maximum_diameter":2,"minimum_sample_separation":null,"singleton_terminal_level":false},{"depth":1,"cells":2,"maximum_diameter":1,"minimum_sample_separation":0.5,"singleton_terminal_level":false}]}}})",
        R"({"figure":"coordinate_reconstruction","records":[{"label":"x","time":[0,1],"reference":[[0,1],[1,0]],"decoded":[[0,1.01],[1,0]]}]})",
        R"({"figure":"character_mass","records":[{"label":"a","conductor_masses":[{"denominator":1,"mass":1},{"denominator":2,"mass":0.1}]}]})",
        R"({"figure":"error_comparison","x":[8,16],"curves":{"L2":[0.1,0.01]}})",
        R"({"figure":"solenoid_suspension","base":2,"depth":3,"phase_turns":0.4,"address":3})",
        R"({"figure":"cycle_clock","periods":[2,4,8],"initial_step":6,"steps":8})",
        R"({"figure":"clock_character","evolution":{"samples":[{"time":0,"bare_characters":[[1,0]],"corrected_characters":[[1,0]]},{"time":1,"bare_characters":[[0,1]],"corrected_characters":[[0.5,0.5]]}]},"linear_reference":[[1,0],[0.5,0.5]]})"};
    for (const char* text : kinds) {
        auto document = json::parse(text);
        document.at("schema") = json::Value("dynsyslab.research-figure/v1");
        require(!dsl::buildDocumentPayload(document, {}).figures.empty(), "every research figure kind draws");
    }
}

dsl::Payload samplePayload() {
    dsl::Payload payload{"memory", "time", "sample", 2};
    dsl::Figure line; line.kind = "time"; line.title = "x versus time"; line.xLabel = "t"; line.yLabel = "x";
    dsl::Trace trace; trace.label = "x"; for (int i = 0; i < 50; ++i) { trace.x.push_back(i); trace.y.push_back(std::sin(0.2 * i)); }
    line.traces.push_back(trace);
    dsl::Figure heat; heat.kind = "heatmap"; heat.title = "map";
    dsl::Trace map; map.heatmap = true; map.z = {0, 1, 2, 3}; map.zShape = {2, 2}; map.x = {0, 1}; map.y = {0, 1};
    heat.traces.push_back(map);
    payload.figures = {line, heat};
    return payload;
}

void exports(const fs::path& root) {
    const auto payload = samplePayload();
    const auto svg = dsl::payloadSvg(payload, 0, 900, 400);
    require(svg.rfind("<svg", 0) == 0 && svg.find("<path") != std::string::npos && svg.find("<rect x=") != std::string::npos,
            "SVG has a curve and heat-map cells");
    const auto html = dsl::payloadHtml(payload, "sample </script> title");
    require(html.find("id=\"dsl-data\"") != std::string::npos && html.find("sample &lt;/script>") != std::string::npos,
            "HTML embeds the data and escapes the title");

    const auto spec = json::parse(R"({"schema":"dynsyslab.report/v1","title":"Report","footer":"test","pages":[
        {"title":"First","blocks":[{"kind":"paragraph","text":"Hello report."},
          {"kind":"table","columns":["name","value"],"rows":[["a","1"],["b","2"]]},
          {"kind":"figure","source":"memory","caption":"A figure"}]},
        {"title":"Second","blocks":[{"kind":"svg","markup":"<svg></svg>"}]}]})");
    const dsl::ReportFigureSource figures = [&](const json::Value&) { return payload; };
    const auto report = dsl::reportHtml(spec, figures);
    require(report.find("<table>") != std::string::npos && report.find("class=\"plot\"") != std::string::npos
            && report.find("<section>") != std::string::npos, "HTML report with table and inline SVG figure");
    const auto pdf = root / "report.pdf";
    dsl::writeReportPdf(spec, figures, pdf);
    std::ifstream file(pdf, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    // Two report pages; the figure does not fit under the table and continues on a third PDF page.
    require(bytes.rfind("%PDF-1.4", 0) == 0 && bytes.find("/Count 3") != std::string::npos
            && bytes.find("continued") != std::string::npos && bytes.find("/DCTDecode") != std::string::npos
            && bytes.find("%%EOF") != std::string::npos, "PDF pages with an embedded figure");
}

void presentation() {
    dynsyslab::core::NpzRecord record;
    record.arrays["parameter"] = dynsyslab::core::NpyArray::fromDoubles({3}, {1, 1, 2});
    record.arrays["maximum"] = dynsyslab::core::NpyArray::fromDoubles({3}, {0.5, 0.7, 1.2});
    const auto bifurcation = dsl::buildBifurcationPayload("b.npz", "bifurcation", record,
        json::parse(R"({"parameter_label":"A","value_label":"max r"})"));
    require(bifurcation.figures.size() == 1 && bifurcation.figures[0].traces[0].marker
            && bifurcation.figures[0].xLabel == "A", "bifurcation scatter");

    dynsyslab::core::NpzRecord observability;
    observability.arrays["observability_singular_values"] = dynsyslab::core::NpyArray::fromDoubles({2}, {2, 1});
    observability.arrays["observable_directions"] = dynsyslab::core::NpyArray::fromDoubles({2, 2}, {1, 0, 0, 1});
    observability.arrays["mode_participation"] = dynsyslab::core::NpyArray::fromDoubles({1, 2}, {1, 1});
    observability.arrays["observer_poles"] = dynsyslab::core::NpyArray::fromComplex({1}, {{-1.0, 2.0}});
    observability.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({2}, {0, 1});
    observability.arrays["state_error_norm"] = dynsyslab::core::NpyArray::fromDoubles({2}, {1, 0.5});
    const auto analysis = dsl::buildObservabilityPayload("o.npz", "observability", observability, json::Value{});
    require(analysis.figures.size() == 5, "spectrum, directions, participation, poles and error norm");

    std::vector<std::string> views{"time", "phase", "three_d", "fourier"};
    dsl::orderViewsForSystem("lorenz", views);
    require(views.front() == "three_d" && views.size() == 4, "Lorenz opens on its three-state view");
    require(dsl::systemCamera("lorenz").has_value() && !dsl::systemCamera("duffing"), "oblique camera for attractors only");
}

}  // namespace

int main() {
    const auto root = fs::temp_directory_path() / ("dynsyslab-document-export-"
        + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        fs::create_directories(root);
        documents();
        exports(root);
        presentation();
        fs::remove_all(root);
        std::cout << "dynsyslab_document_export_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored; fs::remove_all(root, ignored);
        std::cerr << "dynsyslab_document_export_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
