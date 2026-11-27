// dynsyslab-cli experiment report <manifest.experiment.json> [--output report.html|report.pdf]
//
// Builds a dynsyslab.report/v1 spec from an experiment manifest (stage dependency graph, then one
// page per stage: the graph with that stage highlighted, its method table, measured results and its
// Viewer plot) and renders it with DynSysLabViewerCli report into one self-contained file. Writes the
// spec beside the output (<output>.report.json). Prints dynsyslab.experiment-report/v1.
#include "CliExperimentReport.hpp"
#include "CliProcess.hpp"
#include "CliRuntimePaths.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace dynsyslab::cli {
namespace {
namespace json = dynsyslab::core::json;

std::string compact(const json::Value& value) {
    auto text = json::dump(value);
    return text.size() <= 160 ? text : text.substr(0, 157) + "...";
}

std::string escaped(const std::string& text) {
    std::string out;
    for (const char c : text) { if (c == '&') out += "&amp;"; else if (c == '<') out += "&lt;"; else if (c == '>') out += "&gt;"; else out += c; }
    return out;
}

std::vector<std::string> methodNames(const json::Value& job) {
    std::vector<std::string> names;
    for (const auto& method : job.at("methods").array()) names.push_back(method.at("name").string());
    return names;
}

bool uses(const std::vector<std::string>& methods, std::initializer_list<const char*> names) {
    return std::any_of(methods.begin(), methods.end(), [&](const std::string& m) {
        return std::any_of(names.begin(), names.end(), [&](const char* n) { return m == n; }); });
}

struct Role { std::string layer, text; std::vector<std::string> views; };

Role roleOf(const json::Value& job) {
    const auto kind = job.at("kind").string();
    if (kind == "solve") {
        const auto controller = job.at("controller").at("kind").string("none");
        if (controller == "none" || controller.empty()) return {"solver", "training data", {"time", "phase"}};
        return {"controller", controller + " control run", {"control", "time"}};
    }
    const auto methods = methodNames(job);
    std::string text;
    for (const auto& m : methods) text += (text.empty() ? "" : " + ") + m;
    if (text.empty()) text = kind;
    if (uses(methods, {"krom_extract"})) return {"extracted data", text, {"krom_dashboard", "koopman_spectrum"}};
    if (uses(methods, {"sindy"})) return {"post-processing", text, {"koopman_spectrum", "krom_dashboard"}};
    if (uses(methods, {"control_regression", "history_control_regression", "input_map"}))
        return {"identification", text, {"control_map_error", "koopman_spectrum"}};
    if (uses(methods, {"coordinate_transform", "history_basis", "history_projection"}))
        return {"coordinate change", text, {"time", "observable_time"}};
    if (uses(methods, {"klqe", "observability"})) return {"estimation", text, {"observability", "time"}};
    if (uses(methods, {"hdmd", "dmd", "edmd", "resdmd", "mpedmd", "koopman_spectral"}))
        return {"post-processing", text, {"koopman_spectrum"}};
    return {kind == "postprocess" ? "post-processing" : kind, text, {}};
}

// Edges from ${stage.archive|manifest|archives[k]|output:x} references inside each stage's job.
std::vector<std::pair<std::string, std::string>> edges(const json::Value::Array& stages) {
    std::set<std::string> ids;
    for (const auto& stage : stages) ids.insert(stage.at("id").string());
    std::vector<std::pair<std::string, std::string>> links;
    for (const auto& stage : stages) {
        const auto text = json::dump(stage.at("job"));
        std::set<std::string> found;
        for (std::size_t at = text.find("${"); at != std::string::npos; at = text.find("${", at + 2)) {
            const auto dot = text.find('.', at);
            if (dot == std::string::npos) continue;
            const auto id = text.substr(at + 2, dot - at - 2);
            if (ids.count(id) && id != stage.at("id").string()) found.insert(id);
        }
        for (const auto& id : found) links.emplace_back(id, stage.at("id").string());
    }
    return links;
}

std::string graphSvg(const json::Value::Array& stages, const std::vector<std::pair<std::string, std::string>>& links,
                     const std::string& active, const std::map<std::string, std::string>& status) {
    std::map<std::string, int> depth;
    for (const auto& stage : stages) depth[stage.at("id").string()] = 0;
    for (std::size_t pass = 0; pass < stages.size(); ++pass)  // longest-path layering (stages are ordered)
        for (const auto& [a, b] : links) depth[b] = std::max(depth[b], depth[a] + 1);
    std::map<int, std::vector<std::string>> layers;
    for (const auto& stage : stages) layers[depth[stage.at("id").string()]].push_back(stage.at("id").string());
    constexpr int w = 190, h = 58, gx = 46, gy = 22;
    std::map<std::string, std::pair<int, int>> position;
    std::size_t tallest = 1;
    for (const auto& [d, members] : layers) {
        tallest = std::max(tallest, members.size());
        for (std::size_t k = 0; k < members.size(); ++k) position[members[k]] = {20 + d * (w + gx), 20 + static_cast<int>(k) * (h + gy)};
    }
    const int width = 40 + (layers.rbegin()->first + 1) * (w + gx) - gx;
    const int height = 40 + static_cast<int>(tallest) * (h + gy) - gy;
    std::ostringstream out;
    out << "<svg viewBox=\"0 0 " << width << ' ' << height << "\" width=\"100%\" style=\"max-width:" << width
        << "px\" role=\"img\"><defs><marker id=\"a\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"7\" "
           "markerHeight=\"7\" orient=\"auto\"><path d=\"M0,0L10,5L0,10z\" fill=\"#5b6876\"/></marker></defs>";
    for (const auto& [a, b] : links) {
        const auto [x1, y1] = position[a]; const auto [x2, y2] = position[b];
        out << "<path d=\"M" << x1 + w << ',' << y1 + h / 2 << " C" << x1 + w + gx / 2 << ',' << y1 + h / 2 << ' '
            << x2 - gx / 2 << ',' << y2 + h / 2 << ' ' << x2 - 2 << ',' << y2 + h / 2
            << "\" fill=\"none\" stroke=\"#5b6876\" stroke-width=\"1.4\" marker-end=\"url(#a)\"/>";
    }
    for (const auto& stage : stages) {
        const auto id = stage.at("id").string();
        const auto [x, y] = position[id];
        const auto role = roleOf(stage.at("job"));
        const auto state = status.count(id) ? status.at(id) : "";
        const bool on = id == active;
        auto label = role.text.size() <= 30 ? role.text : role.text.substr(0, 29) + "...";
        out << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << h << "\" rx=\"8\" fill=\"#ffffff\" stroke=\""
            << (state == "failed" ? "#b42318" : on ? "#0f766e" : "#c9d1da") << "\" stroke-width=\"" << (on ? 3.2 : 1.2) << "\"/>"
            << "<text x=\"" << x + 10 << "\" y=\"" << y + 19 << "\" font-size=\"11.5\" fill=\"#5b6876\" font-weight=\""
            << (on ? 700 : 400) << "\">(" << escaped(role.layer) << ") " << escaped(id)
            << (state == "failed" || state == "reused" || state == "not run" ? " - " + state : "") << "</text>"
            << "<text x=\"" << x + 10 << "\" y=\"" << y + 40 << "\" font-size=\"13\" fill=\"#16202b\" font-weight=\""
            << (on ? 700 : 400) << "\">" << escaped(label) << "</text>";
    }
    out << "</svg>";
    return out.str();
}

json::Value row(const std::string& name, const std::string& value) {
    return json::Value(json::Value::Array{json::Value(name), json::Value(value)});
}

json::Value::Array methodRows(const json::Value& job) {
    json::Value::Array rows;
    if (job.at("kind").string() == "solve") {
        rows.push_back(row("system", compact(job.at("system").at("name"))));
        rows.push_back(row("time", compact(job.at("time"))));
        rows.push_back(row("initial condition", compact(job.at("initial_conditions").at("values"))));
        rows.push_back(row("solver", compact(job.at("solver"))));
        if (!job.at("controller").isNull()) {
            auto settings = job.at("controller").at("settings");
            settings.object().erase("krom");
            rows.push_back(row("controller", compact(job.at("controller").at("kind"))));
            rows.push_back(row("controller settings", compact(settings)));
            rows.push_back(row("target", compact(job.at("controller").at("target"))));
        }
        return rows;
    }
    for (const auto& method : job.at("methods").array()) rows.push_back(row(method.at("name").string(), compact(method.at("settings"))));
    return rows;
}

// Measured results of a stage's first archive: arrays and their shapes, the leading eigenvalues,
// and the scalar metadata entries.
json::Value::Array resultRows(const std::filesystem::path& archive) {
    json::Value::Array rows;
    const auto record = core::readNpz(archive);
    std::string arrays;
    for (const auto& [name, array] : record.arrays) {
        if (name == "metadata_json") continue;
        std::string shape;
        for (const auto extent : array.shape) shape += (shape.empty() ? "" : " x ") + std::to_string(extent);
        arrays += (arrays.empty() ? "" : ", ") + name + " [" + shape + "]";
    }
    rows.push_back(row("arrays", arrays));
    if (record.contains("continuous_eigenvalues")) {
        const auto values = record.at("continuous_eigenvalues").complexDoubles();
        std::ostringstream text;
        for (std::size_t k = 0; k < values.size() && k < 6; ++k)
            text << (k ? ", " : "") << values[k].real() << (values[k].imag() < 0 ? " - " : " + ") << std::abs(values[k].imag()) << "i";
        rows.push_back(row("modes", std::to_string(values.size())));
        rows.push_back(row("leading eigenvalues", text.str()));
    }
    const auto text = record.metadataJson();
    if (!text.empty()) {
        std::size_t shown = 0;
        const auto metadata = json::parse(text);  // kept alive: the loop iterates into it
        for (const auto& [key, value] : metadata.object()) {
            const auto& storage = value.storage();
            const bool scalar = std::holds_alternative<double>(storage) || std::holds_alternative<bool>(storage)
                || (std::holds_alternative<std::string>(storage) && value.string().size() <= 80);
            if (!scalar || key == "schema" || key == "created_utc" || ++shown > 24) continue;
            rows.push_back(row(key, value.string(json::dump(value))));
        }
    }
    return rows;
}

json::Value table(const std::string& title, json::Value::Array rows) {
    return json::Value(json::Value::Object{{"kind", json::Value("table")}, {"title", json::Value(title)},
                                           {"columns", json::Value(json::Value::Array{})}, {"rows", json::Value(std::move(rows))}});
}

json::Value paragraph(const std::string& text, const std::string& style = {}) {
    return json::Value(json::Value::Object{{"kind", json::Value("paragraph")}, {"text", json::Value(text)}, {"style", json::Value(style)}});
}

}  // namespace

json::Value experimentReportSpec(const json::Value& manifest, const std::filesystem::path& manifestPath) {
    const auto& stages = manifest.at("expanded_job").at("stages").array();
    std::map<std::string, const json::Value*> records;
    for (const auto& record : manifest.at("stages").array()) records[record.at("id").string()] = &record;
    std::string failed;
    if (manifest.at("status").string() == "failed") {
        const auto error = manifest.at("error").string();
        const auto open = error.find("stage '");
        if (open != std::string::npos) failed = error.substr(open + 7, error.find('\'', open + 7) - open - 7);
    }
    std::map<std::string, std::string> status;
    for (const auto& stage : stages) {
        const auto id = stage.at("id").string();
        status[id] = records.count(id) ? (records[id]->at("reused").boolean() ? "reused" : "ran") : id == failed ? "failed" : "not run";
    }
    const auto links = edges(stages);
    const auto svg = [&](const std::string& active) {
        return json::Value(json::Value::Object{{"kind", json::Value("svg")}, {"markup", json::Value(graphSvg(stages, links, active, status))}});
    };
    json::Value::Array pages;
    json::Value::Array overview{paragraph("Experiment manifest " + manifestPath.filename().u8string() + " - status "
                                          + manifest.at("status").string("complete"), "muted")};
    if (!manifest.at("error").isNull()) overview.push_back(paragraph(manifest.at("error").string(), "bad"));
    overview.push_back(svg({}));
    pages.emplace_back(json::Value::Object{{"title", json::Value("Stages")}, {"blocks", json::Value(std::move(overview))}});
    for (const auto& stage : stages) {
        const auto id = stage.at("id").string();
        const auto& job = records.count(id) ? records[id]->at("job") : stage.at("job");
        const auto role = roleOf(job);
        json::Value::Array blocks{svg(id)};
        if (status[id] == "failed") blocks.push_back(paragraph("Failed: " + manifest.at("error").string(), "bad"));
        if (status[id] == "not run") blocks.push_back(paragraph("Not run (an earlier stage failed or the run stopped before it).", "muted"));
        blocks.push_back(table("Method", methodRows(job)));
        const auto& archives = records.count(id) ? records[id]->at("archives").array() : json::Value::Array{};
        if (!archives.empty()) {
            const auto archive = std::filesystem::u8path(archives.front().string());
            try {
                auto rows = resultRows(archive);
                rows.push_back(row("archive", archive.filename().u8string()));
                blocks.push_back(table(status[id] == "reused" ? "Results (reused from an identical earlier run)" : "Results", std::move(rows)));
            } catch (const std::exception& error) {
                blocks.push_back(paragraph(std::string("Results unavailable: ") + error.what(), "bad"));
            }
            json::Value::Array views;
            for (const auto& view : role.views) views.emplace_back(view);
            blocks.emplace_back(json::Value::Object{{"kind", json::Value("figure")}, {"source", json::Value(archive.generic_u8string())},
                {"views", json::Value(std::move(views))}, {"width", json::Value(1500.0)}, {"height", json::Value(1000.0)},
                {"caption", json::Value("Viewer plot of " + archive.filename().u8string())}});
        }
        pages.emplace_back(json::Value::Object{{"title", json::Value("(" + role.layer + ") " + id + ": " + role.text)},
                                               {"blocks", json::Value(std::move(blocks))}});
    }
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.report/v1")},
        {"title", json::Value(manifest.at("label").string(manifestPath.stem().u8string()))},
        {"footer", json::Value("DynSysLab experiment report - " + manifestPath.filename().u8string())},
        {"pages", json::Value(std::move(pages))}});
}

int experimentCommand(const Arguments& args) {
    if (args.positional.size() < 3 || args.positional[1] != L"report")
        throw std::invalid_argument("experiment requires report <manifest.experiment.json> [--output report.html|report.pdf]");
    const std::filesystem::path manifestPath(args.positional[2]);
    const auto manifest = json::load(manifestPath.wstring());
    if (manifest.at("schema").string() != "dynsyslab.experiment/v1")
        throw std::invalid_argument("not an experiment manifest (dynsyslab.experiment/v1)");
    auto stem = manifestPath.stem().wstring();
    if (const auto dot = stem.find(L".experiment"); dot != std::wstring::npos) stem.resize(dot);
    const auto output = core::prepareOutputFile(args.output.empty()
        ? std::filesystem::current_path() / L"experiment-reports" / (stem + L".html") : args.output, "--output");
    const auto specPath = output.parent_path() / (output.stem().wstring() + L".report.json");
    {
        std::ofstream file(specPath, std::ios::binary);
        file << json::dump(experimentReportSpec(manifest, manifestPath));
        if (!file) throw std::runtime_error("unable to write " + specPath.u8string());
    }
    const auto viewer = cliExecutablePath().parent_path() / L"DynSysLabViewerCli.exe";
    if (!std::filesystem::is_regular_file(viewer)) throw std::runtime_error("DynSysLabViewerCli.exe is not next to dynsyslab-cli");
    const auto run = runProcess(viewer, {L"report", specPath.wstring(), L"--output", output.wstring(),
                                         L"--archives", args.archiveRoot.wstring()}, 900000);
    if (!run.started || run.exitCode != 0)
        throw std::runtime_error("the Viewer CLI could not render the report: " + run.err);
    std::cout << json::dump(json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.experiment-report/v1")},
        {"manifest", json::Value(std::filesystem::absolute(manifestPath).generic_u8string())},
        {"output", json::Value(output.generic_u8string())}, {"spec", json::Value(specPath.generic_u8string())},
        {"stages", json::Value(static_cast<double>(manifest.at("expanded_job").at("stages").array().size()))}})) << '\n';
    return 0;
}

}  // namespace dynsyslab::cli
