#include "dynsyslab_viewer/ViewerWindow.hpp"
#include "dynsyslab_viewer/Model.hpp"
#include <windows.h>
#include <commctrl.h>
namespace dsl {
namespace {
struct Reload { unsigned long pid; bool requested = false; };
BOOL CALLBACK reload(HWND window,LPARAM context) {
    auto& request = *reinterpret_cast<Reload*>(context);
    DWORD pid = 0; GetWindowThreadProcessId(window,&pid);
    wchar_t type[128]{}; GetClassNameW(window,type,128);
    if (pid == request.pid && std::wstring(type) == L"DynSysLabViewerWorkbench") {
        request.requested = PostMessageW(window,WM_COMMAND,102,0) != FALSE;
        return FALSE;
    }
    return TRUE;
}
BOOL CALLBACK inspect(HWND window,LPARAM context) {
    wchar_t type[128]{};
    GetClassNameW(window,type,128);
    if (std::wstring(type) != L"DynSysLabViewerWorkbench") return TRUE;
    DWORD pid = 0; GetWindowThreadProcessId(window,&pid);
    wchar_t title[1024]{}; GetWindowTextW(window,title,1024);
    wchar_t path[32768]{}; DWORD size = 32768;
    const auto process = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION,FALSE,pid);
    if (process) { QueryFullProcessImageNameW(process,0,path,&size); CloseHandle(process); }
    DWORD_PTR count = 0;
    const auto tree = GetDlgItem(window,100);
    const bool responds = tree && SendMessageTimeoutW(tree,TVM_GETCOUNT,0,0,
        SMTO_ABORTIFHUNG | SMTO_BLOCK,1000,&count);
    auto& rows = *reinterpret_cast<json::Value::Array*>(context);
    rows.emplace_back(json::Value::Object{{"process_id",static_cast<double>(pid)},
        {"executable",narrow(path)},{"title",narrow(title)},
        {"visible",IsWindowVisible(window)!=FALSE},{"responding",responds},
        {"navigation_nodes",responds?json::Value(static_cast<double>(count)):json::Value()}});
    return TRUE;
}
}
json::Value runningViewerWindows() {
    json::Value::Array rows;
    EnumWindows(inspect,reinterpret_cast<LPARAM>(&rows));
    return json::Value::Object{{"schema","dynsyslab.viewer.live-windows/v1"},{"windows",std::move(rows)}};
}
bool reloadRunningViewer(unsigned long processId) {
    Reload request{processId};
    EnumWindows(reload,reinterpret_cast<LPARAM>(&request));
    return request.requested;
}
}
