#include "dynsyslab_viewer/PlotExport.hpp"

#include <objidl.h>
#include <gdiplus.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(){ULONG_PTR token=0;Gdiplus::GdiplusStartupInput input;
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-native-plot.pdf";
    try{if(Gdiplus::GdiplusStartup(&token,&input,nullptr)!=Gdiplus::Ok)
            throw std::runtime_error("GDI+ startup failed");
        HDC screen=GetDC(nullptr);HBITMAP bitmap=CreateCompatibleBitmap(screen,64,32);
        HDC memory=CreateCompatibleDC(screen);const auto previous=SelectObject(memory,bitmap);
        RECT area{0,0,64,32};FillRect(memory,&area,static_cast<HBRUSH>(GetStockObject(WHITE_BRUSH)));
        SelectObject(memory,previous);DeleteDC(memory);ReleaseDC(nullptr,screen);
        if(!dsl::saveBitmapAsPdf(path,bitmap,64,32))
            throw std::runtime_error("native PDF export failed");
        DeleteObject(bitmap);std::ifstream stream(path,std::ios::binary);
        const std::string contents((std::istreambuf_iterator<char>(stream)),{});
        if(contents.rfind("%PDF-1.4",0)!=0||contents.find("/DCTDecode")==std::string::npos
            ||contents.size()<6||contents.substr(contents.size()-6)!="%%EOF\n")
            throw std::runtime_error("native PDF structure is incomplete");
        Gdiplus::GdiplusShutdown(token);token=0;std::error_code ignored;
        std::filesystem::remove(path,ignored);
        std::cout<<"dynsyslab_plot_export_tests: PASS\n";return 0;
    }catch(const std::exception& error){if(token)Gdiplus::GdiplusShutdown(token);
        std::error_code ignored;std::filesystem::remove(path,ignored);
        std::cerr<<"dynsyslab_plot_export_tests: FAIL: "<<error.what()<<'\n';return 1;}
}
