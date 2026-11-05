#include "dynsyslab_viewer/Workbench.hpp"
#include "dynsyslab_viewer/StartupDiagnostics.hpp"

#include <commctrl.h>
#include <gdiplus.h>
#include <shobjidl.h>
#include <windows.h>

#include <exception>

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int showCommand) {
    if (wcsstr(GetCommandLineW(), L"--smoke-test") != nullptr) {
        try {
            const auto arguments = dsl::parseArguments();
            const dsl::Backend backend(arguments.backend);
            const auto catalog = backend.loadCatalog();
            dsl::writeStartupDiagnostics("Viewer smoke test", arguments.backend.cache,
                arguments.backend.projectRoot, arguments.backend.archiveRoot,
                arguments.backend.catalog, catalog.items.size(), 0,
                arguments.backend.archiveFocus, catalog.focusLabel, catalog.hiddenItemCount);
            return 0;
        }
        catch (...) { return 1; }
    }
    SetCurrentProcessExplicitAppUserModelID(L"DynSysLab.Viewer");
    SetProcessDPIAware();
    INITCOMMONCONTROLSEX controls{sizeof(controls), ICC_TREEVIEW_CLASSES |
                                  ICC_STANDARD_CLASSES | ICC_BAR_CLASSES};
    InitCommonControlsEx(&controls);
    Gdiplus::GdiplusStartupInput graphicsInput;
    ULONG_PTR graphicsToken = 0;
    Gdiplus::GdiplusStartup(&graphicsToken, &graphicsInput, nullptr);
    int result = 1;
    try {
        if (!dsl::Workbench::registerClass(instance)) throw std::runtime_error("Unable to register Viewer classes");
        const auto arguments = dsl::parseArguments();
        dsl::Workbench workbench(arguments);
        result = workbench.run(instance, showCommand);
    } catch (const std::exception& error) {
        if (wcsstr(GetCommandLineW(),L"--background-check") == nullptr)
            MessageBoxA(nullptr, error.what(), "DynSysLab Viewer", MB_OK | MB_ICONERROR);
    }
    Gdiplus::GdiplusShutdown(graphicsToken);
    return result;
}
