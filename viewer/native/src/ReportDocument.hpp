#pragma once
// Multi-page reports from a JSON spec (DynSysLabViewerCli report SPEC --output out.html|out.pdf).
//
// {"schema":"dynsyslab.report/v1","title":"...","subtitle":"...","footer":"...",
//  "pages":[{"title":"...","blocks":[
//     {"kind":"paragraph","text":"..."},
//     {"kind":"table","columns":["name","value"],"rows":[["a","1"],...]},
//     {"kind":"figure","source":"run.npz","view":"time","panel":"label","page":0,
//      "width":1400,"height":800,"caption":"..."},
//     {"kind":"image","path":"plot.png","caption":"..."},
//     {"kind":"svg","markup":"<svg ...>...</svg>"}]}]}       (svg: HTML only)
// HTML output is one self-contained file (figures as inline SVG, images embedded), one section per
// page with print page breaks. PDF output is landscape US Letter (792 x 612 pt), one PDF page per
// report page (figures rasterized through the Viewer canvas, continued onto further pages when a
// page overflows), text in the standard Helvetica fonts (ISO 32000-1:2008, PDF 1.7).
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <functional>
#include <string>

namespace dsl {

inline constexpr const char* kReportSchema = "dynsyslab.report/v1";

// Builds the Viewer payload a figure block names (the CLI resolves source, view and panel).
using ReportFigureSource = std::function<Payload(const dynsyslab::core::json::Value& block)>;

[[nodiscard]] std::string reportHtml(const dynsyslab::core::json::Value& spec, const ReportFigureSource& figures);
// Writes the PDF; figure blocks are rendered to PNG beside `output` and removed afterwards.
void writeReportPdf(const dynsyslab::core::json::Value& spec, const ReportFigureSource& figures,
                    const std::filesystem::path& output);

}  // namespace dsl
