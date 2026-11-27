// DynSysLabViewerCli report SPEC.json --output report.html|report.pdf [--project P --archives A --cache C]
// Renders a dynsyslab.report/v1 spec (src/ReportDocument.hpp). Figure blocks name an archive
// ("source"), a route or panel label ("view", "panel"), optional "coordinates", "start", "stop",
// "scale" and "title"; they are drawn by the same payload services as `render`.
#include "../src/ReportDocument.hpp"
#include "../src/DocumentPayload.hpp"

#include "dynsyslab_viewer/RuntimePaths.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"

#include <dynsyslab/core/OutputPath.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>

namespace {

std::map<std::wstring, std::wstring> options(int argc, wchar_t** argv, std::filesystem::path& spec) {
    std::map<std::wstring, std::wstring> values;
    for (int i = 2; i < argc; ++i) {
        const std::wstring key = argv[i];
        if (key.rfind(L"--", 0) == 0) {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + dsl::narrow(key));
            values[key] = argv[++i];
        } else if (spec.empty()) spec = key;
        else throw std::runtime_error("Unexpected argument: " + dsl::narrow(key));
    }
    return values;
}

}  // namespace

int viewerReportCommand(int argc, wchar_t** argv) {
    try {
        std::filesystem::path specPath;
        auto values = options(argc, argv, specPath);
        if (specPath.empty() || values[L"--output"].empty())
            throw std::runtime_error("report requires SPEC.json and --output report.html|report.pdf");
        const auto spec = dsl::json::load(specPath.wstring());
        const auto paths = dsl::discoverRuntimePaths(L"Viewer", values[L"--project"], values[L"--archives"], values[L"--cache"]);
        dsl::BackendConfig config;
        config.projectRoot = paths.projectRoot; config.archiveRoot = paths.archiveRoot; config.cache = paths.cacheRoot;
        dsl::ViewerSession session(config);
        const auto base = std::filesystem::absolute(specPath).parent_path();
        const dsl::ReportFigureSource figures = [&](const dsl::json::Value& block) {
            auto source = std::filesystem::u8path(block.at("source").string());
            if (source.is_relative()) source = base / source;
            if (dsl::isDocumentPayloadFile(source)) return dsl::buildDocumentPayload(source, block.at("view").string());
            session.selectSource(source);
            dsl::ViewRequest request;
            request.view = block.at("view").string(); request.panel = block.at("panel").string();
            if (request.view.empty() && !block.at("views").array().empty()) {  // first route the archive offers
                const auto routes = session.routes();
                for (const auto& view : block.at("views").array())
                    if (std::find(routes.begin(), routes.end(), view.string()) != routes.end()) { request.view = view.string(); break; }
            }
            request.coordinates = block.at("coordinates").string();
            if (!block.at("start").isNull()) request.window.start = block.at("start").number();
            if (!block.at("stop").isNull()) request.window.stop = block.at("stop").number();
            auto payload = session.workspace(request);
            if (!block.at("scale").isNull()) dsl::setValueScale(payload, block.at("scale").string());
            if (!block.at("title").isNull())
                for (auto& figure : payload.figures) figure.title = block.at("title").string();
            return payload;
        };
        const auto output = dynsyslab::core::prepareOutputFile(values[L"--output"], "--output");
        const bool pdf = output.extension() == L".pdf";
        if (pdf) dsl::writeReportPdf(spec, figures, output);
        else {
            std::ofstream file(output, std::ios::binary);
            file << dsl::reportHtml(spec, figures);
            if (!file) throw std::runtime_error("Unable to write " + output.u8string());
        }
        std::cout << dsl::json::dump(dsl::json::Value::Object{{"schema", "dynsyslab.viewer.report/v1"},
            {"output", output.u8string()}, {"format", pdf ? "pdf" : "html"},
            {"pages", static_cast<double>(spec.at("pages").array().size())}}) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << dsl::json::dump(dsl::json::Value::Object{{"schema", "dynsyslab.viewer.error/v1"},
            {"error", error.what()}}) << '\n';
        return 1;
    }
}
