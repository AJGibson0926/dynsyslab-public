// HTML reports (ReportDocument.hpp).
#include "ReportDocument.hpp"

#include "PlotSvg.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace dsl {
namespace {
namespace json = dynsyslab::core::json;

std::string escaped(const std::string& text) {
    std::string out;
    for (const char c : text) {
        if (c == '&') out += "&amp;"; else if (c == '<') out += "&lt;"; else if (c == '>') out += "&gt;";
        else if (c == '"') out += "&quot;"; else out += c;
    }
    return out;
}

std::string base64(const std::string& bytes) {
    static constexpr char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    std::size_t i = 0;
    for (; i + 2 < bytes.size(); i += 3) {
        const unsigned v = (static_cast<unsigned char>(bytes[i]) << 16) | (static_cast<unsigned char>(bytes[i + 1]) << 8)
                           | static_cast<unsigned char>(bytes[i + 2]);
        out += {table[v >> 18], table[(v >> 12) & 63], table[(v >> 6) & 63], table[v & 63]};
    }
    if (i < bytes.size()) {
        const unsigned v = (static_cast<unsigned char>(bytes[i]) << 16)
                           | (i + 1 < bytes.size() ? static_cast<unsigned char>(bytes[i + 1]) << 8 : 0);
        out += {table[v >> 18], table[(v >> 12) & 63], i + 1 < bytes.size() ? table[(v >> 6) & 63] : '=', '='};
    }
    return out;
}

std::string mimeType(const std::filesystem::path& path) {
    const auto extension = path.extension().u8string();
    if (extension == ".png") return "image/png";
    if (extension == ".jpg" || extension == ".jpeg") return "image/jpeg";
    if (extension == ".svg") return "image/svg+xml";
    throw std::invalid_argument("report images must be PNG, JPEG or SVG: " + path.u8string());
}

constexpr const char* kCss =
    ":root{--bg:#f6f7f9;--card:#fff;--ink:#16202b;--muted:#5b6876;--line:#c9d1da;--accent:#0f766e;--bad:#b42318}"
    "@media (prefers-color-scheme:dark){:root{--bg:#11161c;--card:#19212a;--ink:#e6ebf0;--muted:#9aa7b4;"
    "--line:#33404d;--accent:#2dd4bf;--bad:#f97066}}"
    "body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,Segoe UI,sans-serif}"
    "main{max-width:1180px;margin:0 auto;padding:24px 16px}h1{font-size:24px;margin:0 0 4px}"
    "h2{margin:36px 0 8px;font-size:19px}.card{background:var(--card);border:1px solid var(--line);"
    "border-radius:10px;padding:16px;margin:12px 0;overflow-x:auto}table{border-collapse:collapse;width:100%}"
    "th,td{border-top:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}"
    "th{color:var(--muted);font-weight:600}.muted{color:var(--muted)}.bad{color:var(--bad);font-weight:600}"
    "figure{margin:0}figcaption{color:var(--muted);font-size:13px;margin-top:6px}"
    "svg.plot,img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff}"
    "section{break-before:page}section:first-of-type{break-before:auto}footer{color:var(--muted);font-size:12px;margin:24px 0}";

std::string block(const json::Value& item, const ReportFigureSource& figures) {
    const auto kind = item.at("kind").string();
    std::ostringstream out;
    if (kind == "paragraph") {
        out << "<p" << (item.at("style").string() == "bad" ? " class=\"bad\"" : item.at("style").string() == "muted"
            ? " class=\"muted\"" : "") << '>' << escaped(item.at("text").string()) << "</p>";
    } else if (kind == "table") {
        out << "<div class=\"card\">";
        if (!item.at("title").isNull()) out << "<b>" << escaped(item.at("title").string()) << "</b>";
        out << "<table>";
        if (!item.at("columns").array().empty()) {
            out << "<tr>";
            for (const auto& column : item.at("columns").array()) out << "<th>" << escaped(column.string()) << "</th>";
            out << "</tr>";
        }
        for (const auto& row : item.at("rows").array()) {
            out << "<tr>";
            for (const auto& cell : row.array()) out << "<td>" << escaped(cell.string(json::dump(cell))) << "</td>";
            out << "</tr>";
        }
        out << "</table></div>";
    } else if (kind == "figure" || kind == "image" || kind == "svg") {
        out << "<div class=\"card\"><figure>";
        if (kind == "figure") {
            try {
                const auto payload = figures(item);
                auto svg = payloadSvg(payload, static_cast<int>(item.at("page").number(0)),
                                      static_cast<int>(item.at("width").number(1400)),
                                      static_cast<int>(item.at("height").number(800)));
                svg.insert(4, " class=\"plot\"");
                out << svg;
            } catch (const std::exception& error) {
                out << "<p class=\"bad\">Figure unavailable: " << escaped(error.what()) << "</p>";
            }
        } else if (kind == "image") {
            const auto path = std::filesystem::u8path(item.at("path").string());
            std::ifstream file(path, std::ios::binary);
            if (!file) throw std::invalid_argument("report image is missing: " + path.u8string());
            const std::string bytes((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
            out << "<img alt=\"" << escaped(item.at("caption").string()) << "\" src=\"data:" << mimeType(path)
                << ";base64," << base64(bytes) << "\">";
        } else {
            out << item.at("markup").string();  // trusted markup composed by the report writer
        }
        if (!item.at("caption").isNull()) out << "<figcaption>" << escaped(item.at("caption").string()) << "</figcaption>";
        out << "</figure></div>";
    } else {
        throw std::invalid_argument("unknown report block kind: " + kind);
    }
    return out.str();
}

}  // namespace

std::string reportHtml(const json::Value& spec, const ReportFigureSource& figures) {
    if (spec.at("schema").string() != kReportSchema)
        throw std::invalid_argument(std::string("report spec must have schema ") + kReportSchema);
    const auto title = spec.at("title").string("DynSysLab report");
    std::ostringstream out;
    out << "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" "
           "content=\"width=device-width,initial-scale=1\"><title>" << escaped(title) << "</title><style>" << kCss
        << "</style></head><body><main><h1>" << escaped(title) << "</h1>";
    if (!spec.at("subtitle").isNull()) out << "<p class=\"muted\">" << escaped(spec.at("subtitle").string()) << "</p>";
    for (const auto& page : spec.at("pages").array()) {
        out << "<section>";
        if (!page.at("title").isNull()) out << "<h2>" << escaped(page.at("title").string()) << "</h2>";
        for (const auto& item : page.at("blocks").array()) out << block(item, figures);
        out << "</section>";
    }
    if (!spec.at("footer").isNull()) out << "<footer>" << escaped(spec.at("footer").string()) << "</footer>";
    out << "</main></body></html>";
    return out.str();
}

}  // namespace dsl
