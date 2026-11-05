#include "dynsyslab_compute/Workbench.hpp"

#include <commctrl.h>
#include <gdiplus.h>
#include <shobjidl.h>
#include <exception>

int WINAPI wWinMain(HINSTANCE instance,HINSTANCE,PWSTR,int showCommand) {
    if(wcsstr(GetCommandLineW(),L"--smoke-test")!=nullptr){
        try{dsl::compute::Workbench workbench(dsl::compute::parseArguments());return 0;}
        catch(...){return 1;}
    }
    SetCurrentProcessExplicitAppUserModelID(L"DynSysLab.App");
    SetProcessDPIAware(); INITCOMMONCONTROLSEX controls{sizeof(controls),ICC_TREEVIEW_CLASSES|ICC_STANDARD_CLASSES|ICC_BAR_CLASSES};
    InitCommonControlsEx(&controls); Gdiplus::GdiplusStartupInput input; ULONG_PTR token=0;Gdiplus::GdiplusStartup(&token,&input,nullptr);
    int result=1;try{if(!dsl::compute::Workbench::registerClass(instance))throw std::runtime_error("Unable to register compute UI classes");
        const auto arguments=dsl::compute::parseArguments();dsl::compute::Workbench workbench(arguments);
        result=workbench.run(instance,showCommand);
    }catch(const std::exception& error){MessageBoxA(nullptr,error.what(),"DynSysLab Solver and Post-processor",MB_OK|MB_ICONERROR);}
    Gdiplus::GdiplusShutdown(token);return result;
}
