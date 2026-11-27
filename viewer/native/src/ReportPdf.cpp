// PDF reports (ReportDocument.hpp writeReportPdf). A small PDF 1.4 writer (ISO 32000-1:2008):
// landscape Letter pages, Helvetica / Helvetica-Bold text (WinAnsi), table rules, and figures as
// JPEG image XObjects rendered by the Viewer canvas (renderPlotPng), so a report figure is the
// on-screen plot.
#include "ReportDocument.hpp"

#include "dynsyslab_viewer/ViewerSession.hpp"

#include <windows.h>
#include <objidl.h>
#include <gdiplus.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace dsl {
namespace {
namespace json = dynsyslab::core::json;

constexpr double kWidth = 792, kHeight = 612, kMargin = 40;

struct Image { std::vector<unsigned char> jpeg; int width = 0, height = 0; };

std::string pdfText(const std::string& text) {  // WinAnsi subset: ASCII, others as '?'
    std::string out;
    for (const unsigned char c : text) {
        if (c == '(' || c == ')' || c == '\\') { out += '\\'; out += static_cast<char>(c); }
        else if (c >= 32 && c < 127) out += static_cast<char>(c);
        else if (c >= 0xC0) out += '?';  // first byte of a multi-byte UTF-8 character
    }
    return out;
}

bool jpegEncoder(CLSID& result) {
    UINT count = 0, bytes = 0;
    if (Gdiplus::GetImageEncodersSize(&count, &bytes) != Gdiplus::Ok || !bytes) return false;
    std::vector<BYTE> storage(bytes);
    auto* encoders = reinterpret_cast<Gdiplus::ImageCodecInfo*>(storage.data());
    if (Gdiplus::GetImageEncoders(count, bytes, encoders) != Gdiplus::Ok) return false;
    for (UINT i = 0; i < count; ++i)
        if (std::wstring(encoders[i].MimeType) == L"image/jpeg") { result = encoders[i].Clsid; return true; }
    return false;
}

Image jpegImage(const std::filesystem::path& path) {
    Gdiplus::GdiplusStartupInput input; ULONG_PTR token = 0;
    if (Gdiplus::GdiplusStartup(&token, &input, nullptr) != Gdiplus::Ok) throw std::runtime_error("GDI+ is unavailable");
    Image image;
    {
        Gdiplus::Bitmap bitmap(path.c_str());
        CLSID encoder{}; IStream* stream = nullptr;
        if (bitmap.GetLastStatus() == Gdiplus::Ok && jpegEncoder(encoder)
                && CreateStreamOnHGlobal(nullptr, TRUE, &stream) == S_OK) {
            if (bitmap.Save(stream, &encoder) == Gdiplus::Ok) {
                HGLOBAL memory = nullptr; STATSTG status{};
                if (stream->Stat(&status, STATFLAG_NONAME) == S_OK && GetHGlobalFromStream(stream, &memory) == S_OK) {
                    const auto* bytes = static_cast<const unsigned char*>(GlobalLock(memory));
                    image.jpeg.assign(bytes, bytes + static_cast<std::size_t>(status.cbSize.QuadPart));
                    GlobalUnlock(memory);
                    image.width = static_cast<int>(bitmap.GetWidth()); image.height = static_cast<int>(bitmap.GetHeight());
                }
            }
            stream->Release();
        }
    }
    Gdiplus::GdiplusShutdown(token);
    if (image.jpeg.empty()) throw std::runtime_error("unable to encode report image " + path.u8string());
    return image;
}

class Writer {
public:
    void newPage() { pages_.emplace_back(); y_ = kHeight - kMargin; }
    double& y() { return y_; }
    void text(double x, double y, double size, bool bold, const std::string& value) {
        page() << "BT /" << (bold ? "F2 " : "F1 ") << size << " Tf " << x << ' ' << y << " Td (" << pdfText(value) << ") Tj ET\n";
    }
    void line(double x0, double y0, double x1, double y1) { page() << x0 << ' ' << y0 << " m " << x1 << ' ' << y1 << " l S\n"; }
    void image(const Image& picture, double x, double y, double width, double height) {
        images_.push_back(picture); imagePages_.push_back(pages_.size() - 1);
        page() << "q " << width << " 0 0 " << height << ' ' << x << ' ' << y << " cm /Im" << images_.size() - 1 << " Do Q\n";
    }
    // Ensures `height` points remain above the bottom margin, else starts a continuation page.
    void need(double height, const std::string& header) {
        if (y_ - height >= kMargin + 20) return;
        newPage(); text(kMargin, y_ - 12, 11, true, header + " (continued)"); y_ -= 30;
    }
    void save(const std::filesystem::path& path, const std::string& title, const std::string& footer) {
        std::ofstream out(path, std::ios::binary);
        if (!out) throw std::runtime_error("unable to write " + path.u8string());
        out << "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n";
        std::vector<std::streamoff> offsets;
        const auto begin = [&](std::size_t id) { if (offsets.size() <= id) offsets.resize(id + 1); offsets[id] = out.tellp(); out << id << " 0 obj\n"; };
        const std::size_t pageCount = pages_.size(), firstPage = 5, firstImage = firstPage + 2 * pageCount;
        begin(1); out << "<< /Type /Catalog /Pages 2 0 R >>\nendobj\n";
        begin(2); out << "<< /Type /Pages /Count " << pageCount << " /Kids [";
        for (std::size_t p = 0; p < pageCount; ++p) out << firstPage + 2 * p << " 0 R ";
        out << "] >>\nendobj\n";
        begin(3); out << "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>\nendobj\n";
        begin(4); out << "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>\nendobj\n";
        for (std::size_t p = 0; p < pageCount; ++p) {
            std::ostringstream content; content << std::fixed << std::setprecision(2) << "0.35 G 0.6 w\n" << pages_[p].str()
                << "BT /F1 8 Tf " << kMargin << ' ' << 18 << " Td (" << pdfText(footer + "  -  page " + std::to_string(p + 1)
                + " of " + std::to_string(pageCount)) << ") Tj ET\n";
            const auto stream = content.str();
            begin(firstPage + 2 * p);
            out << "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 " << kWidth << ' ' << kHeight << "] /Resources << /Font << /F1 3 0 R /F2 4 0 R >> /XObject << ";
            for (std::size_t i = 0; i < images_.size(); ++i)
                if (imagePages_[i] == p) out << "/Im" << i << ' ' << firstImage + i << " 0 R ";
            out << ">> >> /Contents " << firstPage + 2 * p + 1 << " 0 R >>\nendobj\n";
            begin(firstPage + 2 * p + 1);
            out << "<< /Length " << stream.size() << " >>\nstream\n" << stream << "endstream\nendobj\n";
        }
        for (std::size_t i = 0; i < images_.size(); ++i) {
            begin(firstImage + i);
            out << "<< /Type /XObject /Subtype /Image /Width " << images_[i].width << " /Height " << images_[i].height
                << " /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length " << images_[i].jpeg.size() << " >>\nstream\n";
            out.write(reinterpret_cast<const char*>(images_[i].jpeg.data()), static_cast<std::streamsize>(images_[i].jpeg.size()));
            out << "\nendstream\nendobj\n";
        }
        const std::size_t infoId = offsets.size();
        begin(infoId); out << "<< /Title (" << pdfText(title) << ") /Producer (DynSysLab) >>\nendobj\n";
        const auto xref = out.tellp();
        out << "xref\n0 " << offsets.size() << "\n0000000000 65535 f \n";
        for (std::size_t id = 1; id < offsets.size(); ++id) out << std::setw(10) << std::setfill('0') << offsets[id] << " 00000 n \n";
        out << "trailer\n<< /Size " << offsets.size() << " /Root 1 0 R /Info " << infoId << " 0 R >>\nstartxref\n" << xref << "\n%%EOF\n";
        if (!out) throw std::runtime_error("unable to write " + path.u8string());
    }

private:
    std::ostringstream& page() { if (pages_.empty()) newPage(); return pages_.back(); }
    std::vector<std::ostringstream> pages_;
    std::vector<Image> images_;
    std::vector<std::size_t> imagePages_;
    double y_ = kHeight - kMargin;
};

std::vector<std::string> wrapped(const std::string& text, double size, double width) {
    const auto perLine = std::max<std::size_t>(10, static_cast<std::size_t>(width / (0.5 * size)));
    std::vector<std::string> lines; std::istringstream words(text); std::string word, line;
    while (words >> word) {
        if (!line.empty() && line.size() + 1 + word.size() > perLine) { lines.push_back(line); line.clear(); }
        line += (line.empty() ? "" : " ") + word;
    }
    if (!line.empty()) lines.push_back(line);
    return lines;
}

void placeImage(Writer& pdf, const Image& image, const std::string& header, const std::string& caption) {
    const double width = kWidth - 2 * kMargin;
    double height = width * image.height / std::max(1, image.width);
    const double maximum = kHeight - 2 * kMargin - 60;
    double drawnWidth = width;
    if (height > maximum) { drawnWidth *= maximum / height; height = maximum; }
    pdf.need(height + 20, header);
    pdf.y() -= height;
    pdf.image(image, kMargin + (width - drawnWidth) / 2, pdf.y(), drawnWidth, height);
    if (!caption.empty()) { pdf.y() -= 13; pdf.text(kMargin, pdf.y(), 9, false, caption); }
    pdf.y() -= 12;
}

}  // namespace

void writeReportPdf(const json::Value& spec, const ReportFigureSource& figures, const std::filesystem::path& output) {
    if (spec.at("schema").string() != kReportSchema)
        throw std::invalid_argument(std::string("report spec must have schema ") + kReportSchema);
    const auto title = spec.at("title").string("DynSysLab report");
    Writer pdf;
    std::size_t index = 0;
    for (const auto& page : spec.at("pages").array()) {
        pdf.newPage();
        const auto header = page.at("title").string(title);
        if (index++ == 0) { pdf.text(kMargin, pdf.y() - 18, 18, true, title); pdf.y() -= 30;
            if (!spec.at("subtitle").isNull()) { pdf.text(kMargin, pdf.y() - 10, 10, false, spec.at("subtitle").string()); pdf.y() -= 18; } }
        pdf.text(kMargin, pdf.y() - 14, 13, true, header); pdf.y() -= 26;
        for (const auto& block : page.at("blocks").array()) {
            const auto kind = block.at("kind").string();
            if (kind == "paragraph") {
                for (const auto& line : wrapped(block.at("text").string(), 10, kWidth - 2 * kMargin)) {
                    pdf.need(14, header); pdf.text(kMargin, pdf.y() - 10, 10, false, line); pdf.y() -= 14;
                }
                pdf.y() -= 6;
            } else if (kind == "table") {
                const auto& columns = block.at("columns").array();
                std::size_t count = columns.size();
                for (const auto& row : block.at("rows").array()) count = std::max(count, row.array().size());
                if (!count) continue;
                const double cell = (kWidth - 2 * kMargin) / static_cast<double>(count);
                const auto row = [&](const json::Value::Array& cells, bool bold) {
                    pdf.need(15, header);
                    for (std::size_t c = 0; c < cells.size(); ++c) {
                        auto value = cells[c].string(json::dump(cells[c]));
                        const auto limit = static_cast<std::size_t>(cell / 4.6);
                        if (value.size() > limit && limit > 3) value = value.substr(0, limit - 3) + "...";
                        pdf.text(kMargin + c * cell + 3, pdf.y() - 10, 9, bold, value);
                    }
                    pdf.y() -= 14; pdf.line(kMargin, pdf.y(), kWidth - kMargin, pdf.y());
                };
                if (!block.at("title").isNull()) { pdf.need(16, header); pdf.text(kMargin, pdf.y() - 10, 10, true, block.at("title").string()); pdf.y() -= 16; }
                if (!columns.empty()) row(columns, true);
                for (const auto& cells : block.at("rows").array()) row(cells.array(), false);
                pdf.y() -= 8;
            } else if (kind == "figure" || kind == "image") {
                std::filesystem::path path;
                const auto scratch = output.parent_path() / (output.stem().wstring() + L".figure-" + std::to_wstring(index) + L".png");
                try {
                    if (kind == "image") path = std::filesystem::u8path(block.at("path").string());
                    else {
                        PlotRenderOptions options;
                        options.width = static_cast<int>(block.at("width").number(1400));
                        options.height = static_cast<int>(block.at("height").number(800));
                        options.page = static_cast<int>(block.at("page").number(0)); options.toolbar = false;
                        (void)renderPlotPng(figures(block), scratch, options);
                        path = scratch;
                    }
                    placeImage(pdf, jpegImage(path), header, block.at("caption").string());
                } catch (const std::exception& error) {
                    pdf.need(14, header); pdf.text(kMargin, pdf.y() - 10, 10, false, std::string("Figure unavailable: ") + error.what());
                    pdf.y() -= 16;
                }
                std::error_code ignored; std::filesystem::remove(scratch, ignored);
            } else if (kind == "svg") {
                pdf.need(14, header); pdf.text(kMargin, pdf.y() - 10, 9, false, "(diagram in the HTML report)"); pdf.y() -= 14;
            } else throw std::invalid_argument("unknown report block kind: " + kind);
        }
    }
    pdf.save(output, title, spec.at("footer").string(title));
}

}  // namespace dsl
