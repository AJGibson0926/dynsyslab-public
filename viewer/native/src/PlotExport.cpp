#include "dynsyslab_viewer/PlotExport.hpp"

#include <commdlg.h>
#include <objidl.h>
#include <gdiplus.h>

#include <algorithm>
#include <cwchar>
#include <cwctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace dsl {
namespace {

class CapturedBitmap {
public:
    explicit CapturedBitmap(HWND window) {
        RECT area{};
        GetClientRect(window, &area);
        width_ = std::max(1L, area.right - area.left);
        height_ = std::max(1L, area.bottom - area.top);
        HDC source = GetDC(window);
        HDC memory = CreateCompatibleDC(source);
        bitmap_ = CreateCompatibleBitmap(source, width_, height_);
        HGDIOBJ previous = SelectObject(memory, bitmap_);
        SendMessageW(window, WM_PRINTCLIENT, reinterpret_cast<WPARAM>(memory),
                     PRF_CLIENT | PRF_ERASEBKGND);
        SelectObject(memory, previous);
        DeleteDC(memory);
        ReleaseDC(window, source);
    }

    ~CapturedBitmap() { if (bitmap_) DeleteObject(bitmap_); }
    HBITMAP get() const { return bitmap_; }
    int width() const { return width_; }
    int height() const { return height_; }
    HBITMAP release() { HBITMAP result = bitmap_; bitmap_ = nullptr; return result; }

private:
    HBITMAP bitmap_ = nullptr;
    int width_ = 1;
    int height_ = 1;
};

bool imageEncoder(const wchar_t* mimeType,CLSID& result) {
    UINT count = 0, bytes = 0;
    if (Gdiplus::GetImageEncodersSize(&count, &bytes) != Gdiplus::Ok || !bytes) return false;
    std::vector<BYTE> storage(bytes);
    auto* encoders = reinterpret_cast<Gdiplus::ImageCodecInfo*>(storage.data());
    if (Gdiplus::GetImageEncoders(count, bytes, encoders) != Gdiplus::Ok) return false;
    for (UINT index = 0; index < count; ++index) {
        if (encoders[index].MimeType && std::wstring(encoders[index].MimeType) == mimeType) {
            result = encoders[index].Clsid;
            return true;
        }
    }
    return false;
}

bool savePng(const std::filesystem::path& path,const CapturedBitmap& capture){
    CLSID encoder{};if(!imageEncoder(L"image/png",encoder))return false;
    Gdiplus::Bitmap image(capture.get(),nullptr);
    return image.Save(path.c_str(),&encoder)==Gdiplus::Ok;
}

std::vector<unsigned char> jpegBytes(HBITMAP bitmap){
    CLSID encoder{};if(!imageEncoder(L"image/jpeg",encoder))return {};
    IStream* stream=nullptr;
    if(CreateStreamOnHGlobal(nullptr,TRUE,&stream)!=S_OK)return {};
    Gdiplus::Bitmap image(bitmap,nullptr);
    if(image.Save(stream,&encoder)!=Gdiplus::Ok){stream->Release();return {};}
    STATSTG status{};HGLOBAL memory=nullptr;
    if(stream->Stat(&status,STATFLAG_NONAME)!=S_OK
        ||GetHGlobalFromStream(stream,&memory)!=S_OK){stream->Release();return {};}
    const auto size=static_cast<std::size_t>(status.cbSize.QuadPart);
    const auto* bytes=static_cast<const unsigned char*>(GlobalLock(memory));
    if(!bytes){stream->Release();return {};}
    std::vector<unsigned char> result(bytes,bytes+size);GlobalUnlock(memory);
    stream->Release();return result;
}

bool savePdf(const std::filesystem::path& path,HBITMAP bitmap,
             int bitmapWidth,int bitmapHeight){
    const auto jpeg=jpegBytes(bitmap);if(jpeg.empty())return false;
    constexpr double pageWidth=842.0,pageHeight=595.0,margin=24.0;
    const double scale=std::min((pageWidth-2.0*margin)/bitmapWidth,
        (pageHeight-2.0*margin)/bitmapHeight);
    const double width=bitmapWidth*scale,height=bitmapHeight*scale;
    const double x=(pageWidth-width)/2.0,y=(pageHeight-height)/2.0;
    std::ostringstream commands;commands<<std::fixed<<std::setprecision(3)
        <<"q "<<width<<" 0 0 "<<height<<' '<<x<<' '<<y<<" cm /Im0 Do Q\n";
    const auto content=commands.str();std::ofstream output(path,std::ios::binary);
    if(!output)return false;output<<"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n";
    std::vector<std::streamoff> offsets(6);
    auto object=[&](int id,const std::string& value){offsets[id]=output.tellp();
        output<<id<<" 0 obj\n"<<value<<"\nendobj\n";};
    object(1,"<< /Type /Catalog /Pages 2 0 R >>");
    object(2,"<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    object(3,"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 842 595] "
        "/Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>");
    offsets[4]=output.tellp();output<<"4 0 obj\n<< /Type /XObject /Subtype /Image /Width "
        <<bitmapWidth<<" /Height "<<bitmapHeight
        <<" /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length "
        <<jpeg.size()<<" >>\nstream\n";
    output.write(reinterpret_cast<const char*>(jpeg.data()),
        static_cast<std::streamsize>(jpeg.size()));output<<"\nendstream\nendobj\n";
    offsets[5]=output.tellp();output<<"5 0 obj\n<< /Length "<<content.size()
        <<" >>\nstream\n"<<content<<"endstream\nendobj\n";
    const auto xref=output.tellp();output<<"xref\n0 6\n0000000000 65535 f \n";
    for(int id=1;id<=5;++id)output<<std::setw(10)<<std::setfill('0')
        <<offsets[id]<<" 00000 n \n";
    output<<"trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n"
        <<xref<<"\n%%EOF\n";return static_cast<bool>(output);
}

std::wstring safeStem(std::wstring value) {
    if (value.empty()) value = L"dynsyslab-plot";
    constexpr wchar_t forbidden[] = L"<>:\"/\\|?*";
    for (wchar_t& character : value) {
        if (wcschr(forbidden, character) || character < 32) character = L'-';
    }
    while (!value.empty() && (value.back() == L'.' || value.back() == L' ')) value.pop_back();
    return value.empty() ? L"dynsyslab-plot" : value;
}

}  // namespace

bool saveBitmapAsPdf(const std::filesystem::path& path,HBITMAP bitmap,
                     int width,int height){
    return bitmap&&width>0&&height>0&&savePdf(path,bitmap,width,height);
}

PlotExportResult copyPlotToClipboard(HWND owner, HWND plotWindow) {
    CapturedBitmap image(plotWindow);
    if (!image.get()) return {false, false, L"Unable to capture the visible plot page"};
    if (!OpenClipboard(owner)) return {false, false, L"The Windows clipboard is currently unavailable"};
    EmptyClipboard();
    const bool copied = SetClipboardData(CF_BITMAP, image.get()) != nullptr;
    if (copied) image.release();
    CloseClipboard();
    return copied ? PlotExportResult{true, false, L"Visible plot page copied to the clipboard"} :
                    PlotExportResult{false, false, L"Windows could not copy the plot image"};
}

PlotExportResult savePlotAsPng(HWND owner, HWND plotWindow,
                               const std::wstring& suggestedStem) {
    std::wstring initial = safeStem(suggestedStem) + L".png";
    std::vector<wchar_t> path(32768, L'\0');
    std::copy(initial.begin(), initial.end(), path.begin());
    OPENFILENAMEW dialog{};
    dialog.lStructSize = sizeof(dialog);
    dialog.hwndOwner = owner;
    dialog.lpstrFilter = L"PNG image (*.png)\0*.png\0All files (*.*)\0*.*\0\0";
    dialog.lpstrFile = path.data();
    dialog.nMaxFile = static_cast<DWORD>(path.size());
    dialog.lpstrDefExt = L"png";
    dialog.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetSaveFileNameW(&dialog)) {
        if (CommDlgExtendedError() == 0) return {false, true, L"Save cancelled"};
        return {false, false, L"Windows could not open the plot save dialog"};
    }

    CapturedBitmap capture(plotWindow);
    if (!capture.get()) return {false, false, L"Unable to capture the visible plot page"};
    if (!savePng(path.data(),capture))
        return {false, false, L"The plot image could not be written"};
    return {true, false, L"Saved plot page to " + std::filesystem::path(path.data()).filename().wstring()};
}

PlotExportResult savePlotToDisk(HWND owner,HWND plotWindow,
                                const std::wstring& suggestedStem){
    std::wstring initial=safeStem(suggestedStem)+L".png";
    std::vector<wchar_t> path(32768,L'\0');std::copy(initial.begin(),initial.end(),path.begin());
    OPENFILENAMEW dialog{};dialog.lStructSize=sizeof(dialog);dialog.hwndOwner=owner;
    dialog.lpstrFilter=L"PNG image (*.png)\0*.png\0PDF report (*.pdf)\0*.pdf\0All files (*.*)\0*.*\0\0";
    dialog.lpstrFile=path.data();dialog.nMaxFile=static_cast<DWORD>(path.size());
    dialog.lpstrDefExt=L"png";dialog.Flags=OFN_OVERWRITEPROMPT|OFN_PATHMUSTEXIST|OFN_NOCHANGEDIR;
    if(!GetSaveFileNameW(&dialog)){if(CommDlgExtendedError()==0)return {false,true,L"Save cancelled"};
        return {false,false,L"Windows could not open the plot save dialog"};}
    std::filesystem::path destination(path.data());auto extension=destination.extension().wstring();
    std::transform(extension.begin(),extension.end(),extension.begin(),::towlower);
    CapturedBitmap capture(plotWindow);if(!capture.get())
        return {false,false,L"Unable to capture the visible plot page"};
    const bool pdf=extension==L".pdf"||dialog.nFilterIndex==2;
    if(pdf&&extension!=L".pdf")destination.replace_extension(L".pdf");
    if(!pdf&&extension!=L".png")destination.replace_extension(L".png");
    const bool saved=pdf?saveBitmapAsPdf(destination,capture.get(),
        capture.width(),capture.height()):savePng(destination,capture);
    if(!saved)return {false,false,pdf?L"The PDF report could not be written":
        L"The PNG image could not be written"};
    return {true,false,L"Saved visible plot page to "+destination.filename().wstring()};
}

}  // namespace dsl
