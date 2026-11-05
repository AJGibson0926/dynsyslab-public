#include "dynsyslab_viewer/Workbench.hpp"
#include <gdiplus.h>
#include <stdexcept>
namespace dsl {
namespace {
struct Snapshot { HWND parent; HDC dc; };
BOOL CALLBACK paintChild(HWND child, LPARAM data) {
    if (!(GetWindowLongPtrW(child,GWL_STYLE)&WS_VISIBLE)) return TRUE;
    auto& snapshot=*reinterpret_cast<Snapshot*>(data);
    RECT area{};GetClientRect(child,&area);POINT origin{};MapWindowPoints(child,snapshot.parent,&origin,1);
    const int saved=SaveDC(snapshot.dc);SetViewportOrgEx(snapshot.dc,origin.x,origin.y,nullptr);
    IntersectClipRect(snapshot.dc,0,0,area.right,area.bottom);
    SendMessageW(child,WM_PRINTCLIENT,reinterpret_cast<WPARAM>(snapshot.dc),PRF_CLIENT);
    RestoreDC(snapshot.dc,saved);return TRUE;
}
}
void Workbench::saveBackgroundWorkspaceSnapshot(const wchar_t* filename) {
    RECT area{};GetClientRect(window_,&area);
    Gdiplus::Bitmap bitmap(area.right,area.bottom,PixelFormat32bppARGB);
    Gdiplus::Graphics graphics(&bitmap);graphics.Clear(theme::color(theme::background));
    HDC dc=graphics.GetHDC();paint(dc);Snapshot snapshot{window_,dc};
    EnumChildWindows(window_,paintChild,reinterpret_cast<LPARAM>(&snapshot));graphics.ReleaseHDC(dc);
    UINT count=0,bytes=0;Gdiplus::GetImageEncodersSize(&count,&bytes);
    std::vector<unsigned char> storage(bytes);auto* codecs=reinterpret_cast<Gdiplus::ImageCodecInfo*>(storage.data());
    Gdiplus::GetImageEncoders(count,bytes,codecs);
    std::filesystem::create_directories(arguments_.backend.cache);
    const auto output=arguments_.backend.cache/filename;
    for(UINT i=0;i<count;++i) if(wcscmp(codecs[i].MimeType,L"image/png")==0) {
        if(bitmap.Save(output.c_str(),&codecs[i].Clsid)!=Gdiplus::Ok)
            throw std::runtime_error("Unable to save hidden workspace snapshot");
        return;
    }
    throw std::runtime_error("PNG encoder unavailable");
}
}
