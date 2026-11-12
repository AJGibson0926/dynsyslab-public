#include "dynsyslab_compute/Workbench.hpp"

#include "dynsyslab_viewer/Json.hpp"

#include <commctrl.h>
#include <gdiplus.h>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace dsl::compute {
namespace {

struct Snapshot { HWND parent; HDC dc; };

BOOL CALLBACK paintChild(HWND child,LPARAM data) {
    auto& snapshot=*reinterpret_cast<Snapshot*>(data);
    // The root is deliberately hidden, so IsWindowVisible cannot be used.
    // Descendants of a hidden panel must still be excluded from this capture.
    for(HWND item=child;item&&item!=snapshot.parent;item=GetParent(item))
        if(!(GetWindowLongPtrW(item,GWL_STYLE)&WS_VISIBLE))return TRUE;
    RECT area{};GetClientRect(child,&area);POINT origin{};
    MapWindowPoints(child,snapshot.parent,&origin,1);
    const int saved=SaveDC(snapshot.dc);SetViewportOrgEx(snapshot.dc,origin.x,origin.y,nullptr);
    IntersectClipRect(snapshot.dc,0,0,area.right,area.bottom);
    SendMessageW(child,WM_PRINTCLIENT,reinterpret_cast<WPARAM>(snapshot.dc),PRF_CLIENT);
    RestoreDC(snapshot.dc,saved);return TRUE;
}

bool tabStop(HWND control) {
    return control && (GetWindowLongPtrW(control,GWL_STYLE)&WS_TABSTOP)!=0;
}

}  // namespace

void Workbench::saveBackgroundSnapshot(const wchar_t* filename) {
    RECT area{};GetClientRect(window_,&area);
    Gdiplus::Bitmap bitmap(std::max(1L,area.right),std::max(1L,area.bottom),PixelFormat32bppARGB);
    Gdiplus::Graphics graphics(&bitmap);graphics.Clear(dsl::theme::color(dsl::theme::background));
    HDC dc=graphics.GetHDC();paint(dc);Snapshot snapshot{window_,dc};
    EnumChildWindows(window_,paintChild,reinterpret_cast<LPARAM>(&snapshot));graphics.ReleaseHDC(dc);
    UINT count=0,bytes=0;Gdiplus::GetImageEncodersSize(&count,&bytes);
    std::vector<unsigned char> storage(bytes);
    auto* codecs=reinterpret_cast<Gdiplus::ImageCodecInfo*>(storage.data());
    Gdiplus::GetImageEncoders(count,bytes,codecs);
    std::filesystem::create_directories(arguments_.backend.cache);
    const auto output=arguments_.backend.cache/filename;
    for(UINT index=0;index<count;++index)if(wcscmp(codecs[index].MimeType,L"image/png")==0){
        if(bitmap.Save(output.c_str(),&codecs[index].Clsid)!=Gdiplus::Ok)
            throw std::runtime_error("Unable to save hidden Solver snapshot");
        return;
    }
    throw std::runtime_error("PNG encoder unavailable");
}

void Workbench::finishBackgroundCheck(int status,const std::wstring& error) {
    KillTimer(window_,5051);
    std::wstring finalError=error;
    RECT system{},conditions{},control{},integrator{},review{},rail{};
    GetWindowRect(stageSystem_,&system);GetWindowRect(stageInitialConditions_,&conditions);
    GetWindowRect(stageControl_,&control);GetWindowRect(stageIntegrator_,&integrator);
    GetWindowRect(stageReview_,&review);GetWindowRect(reviewRailTitle_,&rail);
    const bool ordered=system.right<=conditions.left&&conditions.right<=control.left&&
        control.right<=integrator.left&&integrator.right<=review.left;
    const bool stagesFocusable=tabStop(stageSystem_)&&tabStop(stageInitialConditions_)&&
        tabStop(stageControl_)&&tabStop(stageIntegrator_)&&tabStop(stageReview_);
    const bool panesSeparated=review.right<=rail.left;
    const bool tooltips=headerTooltip_&&SendMessageW(headerTooltip_,TTM_GETTOOLCOUNT,0,0)>=7;
    if(status==0&&(!ordered||!stagesFocusable||!panesSeparated||!tooltips)){
        status=1;finalError=L"Solver workflow accessibility or pane layout contract failed";
    }
    std::filesystem::create_directories(arguments_.backend.cache);
    try{saveBackgroundSnapshot(L"solver-main.png");}
    catch(const std::exception& exception){status=1;finalError=dsl::wide(exception.what());}
    dsl::json::Value::Object result{{"schema","dynsyslab.compute.workbench-check/v1"},
        {"process_id",static_cast<double>(GetCurrentProcessId())},
        {"status",static_cast<double>(status)},{"error",dsl::narrow(finalError)},
        {"catalog_items",static_cast<double>(catalog_.items.size())},
        {"visible",IsWindowVisible(window_)!=FALSE},{"solver_active",solverActive_},
        {"workflow_ordered",ordered},{"workflow_tabstops",stagesFocusable},
        {"panes_separated",panesSeparated},{"header_tooltip_count",static_cast<double>(
            headerTooltip_?SendMessageW(headerTooltip_,TTM_GETTOOLCOUNT,0,0):0)}};
    std::ofstream output(arguments_.backend.cache/L"workbench-check.json",std::ios::binary);
    output<<dsl::json::dump(result);PostQuitMessage(status);
}

}  // namespace dsl::compute
